import os
import re
import pysam
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
from collections import defaultdict

# =====================================================================
# 1. 基础绘图与安全解析工具
# =====================================================================
def draw_bezier_ribbon(ax, x_top, w_top, y_top, x_bot, w_bot, y_bot, color, alpha=0.4, invert=False):
    if not invert:
        verts = [(x_top, y_top), (x_top, y_bot + (y_top - y_bot)/2), (x_bot, y_top - (y_top - y_bot)/2),
                 (x_bot, y_bot), (x_bot + w_bot, y_bot), (x_bot + w_bot, y_top - (y_top - y_bot)/2),
                 (x_top + w_top, y_bot + (y_top - y_bot)/2), (x_top + w_top, y_top), (x_top, y_top)]
    else:
        verts = [(x_top, y_top), (x_top, y_bot + (y_top - y_bot)/2), (x_bot + w_bot, y_top - (y_top - y_bot)/2),
                 (x_bot + w_bot, y_bot), (x_bot, y_bot), (x_bot, y_top - (y_top - y_bot)/2),
                 (x_top + w_top, y_bot + (y_top - y_bot)/2), (x_top + w_top, y_top), (x_top, y_top)]
    codes = [Path.MOVETO] + [Path.CURVE4]*3 + [Path.LINETO] + [Path.CURVE4]*3 + [Path.CLOSEPOLY]
    ax.add_patch(patches.PathPatch(Path(verts, codes), facecolor=color, lw=0, alpha=alpha))


def safe_get(info, key, default_val):
    if key not in info: return default_val
    val = info.get(key)
    if val is None: return default_val
    if isinstance(val, (tuple, list)):
        return val[0] if len(val) > 0 and val[0] is not None else default_val
    return val


def _parse_srcreg(srcreg):
    """Parse 'chr:start-end' → (chr, start, end). Returns (None,0,0) on failure."""
    try:
        chrom, coords = srcreg.rsplit(":", 1)
        start_s, end_s = coords.split("-", 1)
        return chrom, int(start_s), int(end_s)
    except (ValueError, AttributeError):
        return None, 0, 0


def _coerce_bool(val):
    if isinstance(val, bool): return val
    if isinstance(val, int):  return bool(val)
    if isinstance(val, str):  return val.lower() not in ('false', '0', '')
    return bool(val)


# =====================================================================
# 2. VCF 解析（新格式）
# =====================================================================
def parse_asmsv_complex_event(vcf_path, target_event, hap_name):
    """
    解析 bam_info.py 生成的 COMPLEX 事件 VCF（新命名约定）。

    同一事件的所有记录在 INFO.EVENT 中均为 target_event，
    通过 VCF ID 字段（rec.id）区分记录类型：
      Summary:  <target_event>_SUMMARY   FILTER=COMPLEX, SVTYPE=CSV, SUBTYPE=COMPLEX
      Span:     <target_event>_SPAN_<N>  SVTYPE=CSV, CPLX_IDX=N, SRCREG, INS_LEN
      Nest:     <target_event>_NEST_<N>  SVTYPE=DUP, DUP:TANDEM（嵌套于 backbone 内）
      BND:      <target_event>_BND_<N>_ENTRY / _EXIT
      TailBND:  <target_event>_TAIL_BND1 / _TAIL_BND2

    DUP:TANDEM 块有两种来源：
      1. NEST 记录（SVTYPE=DUP）：直接给出 pos/end，按 BND_N_ENTRY 的 to_pos 匹配 CPLX_IDX
      2. 无 NEST 记录时：从 BND_N_ENTRY chain 推断 src_enter/src_leave
    """
    if not os.path.exists(vcf_path):
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    vcf = pysam.VariantFile(vcf_path)
    parent_event = None
    sub_blocks   = {}   # cplx_idx (int) → block dict
    bnd_entries  = {}   # cplx_idx (int) → {from_chr, from_pos, to_chr, to_pos, to_strand}
    # NEST 记录：以 (chrom, pos) 为键，供后续按 bnd_entry["to_pos"] 匹配
    nest_by_pos  = {}   # (chrom, pos) → block dict
    tail_bnd_chr = None
    tail_bnd_pos = None

    valid_gts = [(1, 1)]
    if hap_name == "Hap1":   valid_gts.append((1, 0))
    elif hap_name == "Hap2": valid_gts.append((0, 1))

    te_upper     = target_event.upper()
    tail_bnd1_up = te_upper + "_TAIL_BND1"
    bnd_entry_re = re.compile(rf'{re.escape(te_upper)}_BND_(\d+)_ENTRY$')
    nest_re      = re.compile(rf'{re.escape(te_upper)}_NEST_\d+$')

    for rec in vcf:
        sample_name = list(rec.samples.keys())[0]
        if rec.samples[sample_name]['GT'] not in valid_gts:
            continue

        info       = rec.info
        event_id   = str(safe_get(info, "EVENT", "")).upper()
        if event_id != te_upper:
            continue

        svtype_raw = str(safe_get(info, "SVTYPE", ""))
        subtype    = str(safe_get(info, "SUBTYPE", ""))
        rec_id     = (rec.id or "").upper()

        # ── Summary record ────────────────────────────────────────────────
        if svtype_raw == "CSV" and subtype == "COMPLEX":
            true_end = rec.stop if rec.stop is not None else int(safe_get(info, "END", rec.pos + 1))
            parent_event = {
                "chrom":     rec.chrom,
                "pos":       int(rec.pos),
                "end":       int(true_end),
                "query_len": int(safe_get(info, "INS_LEN", 0)),
                "strand":    str(safe_get(info, "BACKBONE_STRAND", "++")),
                "up_down":   _coerce_bool(safe_get(info, "BACKBONE_FLIP", False)),
            }

        # ── SPAN records (非 TANDEM 孤儿簇，有 CPLX_IDX) ─────────────────
        elif svtype_raw == "CSV" and "CPLX_IDX" in info:
            idx    = int(safe_get(info, "CPLX_IDX", 0))
            srcreg = str(safe_get(info, "SRCREG", ""))
            src_chr, src_pos, src_end = _parse_srcreg(srcreg)
            if src_chr is None:
                src_chr, src_pos, src_end = rec.chrom, int(rec.pos), int(rec.pos)

            q_len   = int(safe_get(info, "INS_LEN", abs(src_end - src_pos)))
            ref_len = abs(src_end - src_pos)

            sub_blocks[idx] = {
                "subtype":  subtype,
                "src_chr":  src_chr,
                "src_pos":  src_pos,
                "src_end":  src_end,
                "q_len":    max(1, q_len),
                "ref_len":  max(1, ref_len),
                "is_inv":   subtype.startswith("Inv"),
            }

        # ── NEST records（SVTYPE=DUP，嵌套 TANDEM，直接给出坐标）─────────
        elif svtype_raw == "DUP" and nest_re.match(rec_id):
            src_pos = int(rec.pos)
            src_end = rec.stop if rec.stop is not None else int(safe_get(info, "END", src_pos + 1))
            ref_len = max(1, abs(src_end - src_pos))
            q_len   = int(safe_get(info, "INS_LEN", ref_len))
            nest_by_pos[(rec.chrom, src_pos)] = {
                "subtype":  subtype or "DUP:TANDEM",
                "src_chr":  rec.chrom,
                "src_pos":  src_pos,
                "src_end":  src_end,
                "q_len":    max(1, q_len),
                "ref_len":  ref_len,
                "is_inv":   False,  # 由 anchor_strand 检查后覆盖
            }

        # ── BND records（BND_N_ENTRY 用于定位 TANDEM / NEST 的 CPLX_IDX）─
        elif svtype_raw == "BND":
            m = bnd_entry_re.match(rec_id)
            if m:
                idx     = int(m.group(1))
                strands = str(safe_get(info, "STRANDS", "++"))
                bnd_entries[idx] = {
                    "from_chr":  rec.chrom,
                    "from_pos":  int(rec.pos),
                    "to_chr":    str(safe_get(info, "CHR2", "")),
                    "to_pos":    int(safe_get(info, "POS2", 0)),
                    "to_strand": strands[1] if len(strands) > 1 else "+",
                }
            elif rec_id == tail_bnd1_up:
                tail_bnd_chr = rec.chrom
                tail_bnd_pos = int(rec.pos)

    if not parent_event:
        raise ValueError(f"COMPLEX event '{target_event}' not found in {vcf_path}.")

    # 确定主干方向（用于判断 TANDEM 是否倒置）
    bb_strands   = parent_event["strand"]
    if len(bb_strands) < 2: bb_strands += bb_strands
    anchor_strand = bb_strands[1] if parent_event["up_down"] else bb_strands[0]

    # ── 重建缺失的 DUP:TANDEM 块 ─────────────────────────────────────────
    for idx in sorted(set(bnd_entries.keys()) - set(sub_blocks.keys())):
        entry     = bnd_entries[idx]
        to_chr    = entry["to_chr"]
        to_pos    = entry["to_pos"]
        to_strand = entry["to_strand"]

        # 优先：若有 NEST 记录在 bnd_entry 的落点（to_chr, to_pos），直接使用
        nest = nest_by_pos.get((to_chr, to_pos))
        if nest:
            is_inv = (to_strand != anchor_strand)
            sub_blocks[idx] = {**nest, "is_inv": is_inv,
                                "subtype": ("Inv" if is_inv else "") + "DUP:TANDEM"}
            continue

        # 退出坐标 = idx+1 的 BND_ENTRY 在该染色体上的 backbone 出发点
        src_enter  = to_pos
        src_leave  = src_enter   # zero-length 兜底
        next_entry = bnd_entries.get(idx + 1)
        if next_entry and next_entry["from_chr"] == to_chr:
            src_leave = next_entry["from_pos"]
        elif tail_bnd_chr == to_chr and tail_bnd_pos is not None:
            src_leave = tail_bnd_pos

        ref_len = abs(src_leave - src_enter)
        is_inv  = (to_strand != anchor_strand)

        sub_blocks[idx] = {
            "subtype":  ("Inv" if is_inv else "") + "DUP:TANDEM",
            "src_chr":  to_chr,
            "src_pos":  min(src_enter, src_leave),
            "src_end":  max(src_enter, src_leave),
            "q_len":    max(1, ref_len),
            "ref_len":  max(1, ref_len),
            "is_inv":   is_inv,
        }

    ordered_blocks = [sub_blocks[i] for i in sorted(sub_blocks.keys())]
    return parent_event, ordered_blocks


# =====================================================================
# 3. 绘图引擎
# =====================================================================
def plot_parsed_complex_sv(parent, blocks, event_id, hap_name, output_path="auto_complex_plot.png"):
    fig, ax = plt.subplots(figsize=(16, 8), dpi=300)
    main_chr    = parent["chrom"]
    flank_len   = 50000
    color_palette = ["#4A90E2", "#50E3C2", "#F5A623", "#D0021B", "#9013FE"]

    up_down    = parent.get("up_down", False)
    bb_strands = parent.get("strand", "++")
    if len(bb_strands) < 2: bb_strands += bb_strands

    # BACKBONE_STRAND 编码（bam_info.py 第1287行）：
    #   up_down=False → strands_str = anchor_strand   + bb_right.strand
    #   up_down=True  → strands_str = bb_right.strand + anchor_strand
    #
    # parent["pos"] 始终为较小坐标（genome-left），parent["end"] 为较大坐标（genome-right）：
    #   up_down=False: left_flank 在 bb_left(anchor) 侧 → strands_str[0]=anchor_strand
    #   up_down=True:  left_flank 在 bb_right 侧         → strands_str[0]=bb_right.strand
    # 两种情况下 strands_str[0] 恒对应 left_flank，strands_str[1] 恒对应 right_flank。
    # 修复：根据 assembly 遍历方向，将入口/出口的 strand 分配给对应的 Ref 侧
    # if up_down:
    #     # Assembly 反向遍历：入口(0)在 Ref 右侧，出口(1)在 Ref 左侧
    #     left_is_inv  = (bb_strands[1] == "-")
    #     right_is_inv = (bb_strands[0] == "-")
    # else:
    #     # Assembly 正向遍历：入口(0)在 Ref 左侧，出口(1)在 Ref 右侧
    #     left_is_inv  = (bb_strands[0] == "-")
    #     right_is_inv = (bb_strands[1] == "-")
    left_is_inv  = (bb_strands[0] == "-")
    right_is_inv = (bb_strands[1] == "-")

    # ── 构建所有块 ──────────────────────────────────────────────────────
    left_flank = {
        "name": "Left_Flank", "subtype": "Flank",
        "q_len": flank_len, "ref_len": flank_len,
        "src_chr": main_chr, "src_pos": parent["pos"] - flank_len,
        "is_inv": left_is_inv, "color": "#A9A9A9", "is_flank": True,
    }
    right_flank = {
        "name": "Right_Flank", "subtype": "Flank",
        "q_len": flank_len, "ref_len": flank_len,
        "src_chr": main_chr, "src_pos": parent["end"],
        "is_inv": right_is_inv, "color": "#A9A9A9", "is_flank": True,
    }

    core_blocks = []
    for idx, b in enumerate(blocks):
        core_blocks.append({
            **b,
            "name":     b["subtype"],
            "color":    color_palette[idx % len(color_palette)],
            "is_flank": False,
        })

    # Assembly 方向决定顺序[right_flank] + core_blocks + [left_flank]
                #   if up_down else
    all_blocks = (
                  [left_flank] + core_blocks + [right_flank])

    # ── Assembly 轴坐标（以 q_len 步进）──────────────────────────────────
    curr_asm_x = 0
    for b in all_blocks:
        b["asm_x"]  = curr_asm_x
        curr_asm_x += b["q_len"]

    # ── Reference 轴坐标（基因组坐标顺序）────────────────────────────────
    # main_chr 按真实基因组坐标排布（归一化到 0 起点），ribbon 在 BACKBONE_FLIP
    # 时自然产生交叉，交叉本身即表达了"assembly 反向遍历"的生物学含义。
    # 跨染色体 block 按各自染色体坐标排布，并居中对齐到对应 assembly 区间下方。
    ref_tracks = defaultdict(list)
    for b in all_blocks:
        ref_tracks[b["src_chr"]].append(b)
        b["len_on_ref"] = max(1, b["ref_len"])

    main_ref_min = min(b["src_pos"] for b in ref_tracks[main_chr])
    for b in ref_tracks[main_chr]:
        b["ref_x"] = b["src_pos"] - main_ref_min

    for chrom, blist in ref_tracks.items():
        if chrom == main_chr:
            continue
        chr_ref_min = min(b["src_pos"] for b in blist)
        for b in blist:
            b["ref_x"] = b["src_pos"] - chr_ref_min
        asm_lo = min(b["asm_x"]                for b in blist)
        asm_hi = max(b["asm_x"] + b["q_len"]   for b in blist)
        ref_lo = min(b["ref_x"]                for b in blist)
        ref_hi = max(b["ref_x"] + b["len_on_ref"] for b in blist)
        shift  = (asm_lo + asm_hi) / 2 - (ref_lo + ref_hi) / 2
        for b in blist:
            b["ref_x"] += shift

    # ref_x 可能出现负值（居中对齐时），统一上移
    min_ref_x = min(b["ref_x"] for b in all_blocks)
    if min_ref_x < 0:
        for b in all_blocks:
            b["ref_x"] -= min_ref_x

    max_plot_x = max(
        max(b["ref_x"]  + b["len_on_ref"] for b in all_blocks),
        max(b["asm_x"]  + b["q_len"]      for b in all_blocks),
    )
    x_span  = max_plot_x
    label_x = -(x_span * 0.05 + 5000)

    # ── 绘制轴线与标签 ────────────────────────────────────────────────────
    involved_chrs = [main_chr] + [c for c in ref_tracks if c != main_chr]
    Y_ASM, v_gap, bh = 0, 2.0, 0.35
    chr_y_map = {chrom: (len(involved_chrs) - i) * v_gap for i, chrom in enumerate(involved_chrs)}

    ax.axhline(Y_ASM, color='black', linewidth=1.2, zorder=1, alpha=0.6)
    ax.text(label_x, Y_ASM, f"Assembly\n({hap_name})", va='center', ha='right', fontsize=11, fontweight='bold')

    for chrom, y_pos in chr_y_map.items():
        ax.axhline(y_pos, color='black', linewidth=1.2, zorder=1, alpha=0.6)
        label = "Main Ref" if chrom == main_chr else "Ref"
        ax.text(label_x, y_pos, f"{label}\n{chrom}", va='center', ha='right', fontsize=11, fontweight='bold')

    # ── 绘制各块 ─────────────────────────────────────────────────────────
    for b in all_blocks:
        w_asm  = b["q_len"]
        w_ref  = b["len_on_ref"]
        color  = b["color"]
        is_inv = b["is_inv"]
        ref_y  = chr_y_map[b["src_chr"]]
        zo     = 2 if b.get("is_flank") else 3

        # Ref 块
        ax.add_patch(patches.Rectangle(
            (b["ref_x"], ref_y - bh/2), w_ref, bh,
            color=color, alpha=0.65, zorder=zo, ec='black', lw=0.8,
        ))
        mid_bp = b["src_pos"] + b["ref_len"] / 2
        ax.text(b["ref_x"] + w_ref / 2, ref_y + bh/2 + 0.1,
                f'{mid_bp / 1e6:.2f}M',
                ha='center', va='bottom', fontsize=8.5, color='#222222',
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=0.5), zorder=4)

        # Assembly 块
        ax.add_patch(patches.Rectangle(
            (b["asm_x"], Y_ASM - bh/2), w_asm, bh,
            color=color, alpha=0.65, zorder=3, ec='black', lw=0.8,
        ))
        if w_asm > x_span * 0.005:
            dx = -w_asm * 0.6 if is_inv else w_asm * 0.6
            sx =  b["asm_x"] + w_asm * 0.8 if is_inv else b["asm_x"] + w_asm * 0.2
            ax.arrow(sx, Y_ASM, dx, 0,
                     head_width=0.08, head_length=min(w_asm * 0.1, 5000),
                     fc='white', ec='white', zorder=4)

        if not b.get("is_flank"):
            ax.text(b["asm_x"] + w_asm / 2, Y_ASM - bh/2 - 0.2, b["name"],
                    ha='center', va='top', fontsize=9, fontweight='bold', color=color,
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=0.2))

        # Bezier 丝带（Ref 顶 → Asm 底）
        draw_bezier_ribbon(
            ax,
            b["ref_x"], w_ref, ref_y - bh/2,
            b["asm_x"], w_asm, Y_ASM + bh/2,
            color, alpha=0.35, invert=is_inv,
        )

    # ── 参考缺口红线 ──────────────────────────────────────────────────────
    # ref 轴按基因组坐标排布，缺口直接由 parent["pos"]/["end"] 定位，与 up_down 无关。
    ref_gap_len = parent["end"] - parent["pos"]
    if ref_gap_len > 0:
        # main_ref_min 可能因 min_ref_x 修正而偏移，需同步
        adjusted_ref_min = main_ref_min - (min_ref_x if min_ref_x < 0 else 0)
        lx = parent["pos"] - adjusted_ref_min
        rx = parent["end"] - adjusted_ref_min
        y_val = chr_y_map[main_chr]
        ax.plot([lx, rx], [y_val, y_val], color='#D0021B', linestyle='--', zorder=5, lw=2)

        l_ref = (f'Ref Gap: {ref_gap_len/1e6:.2f}Mb' if ref_gap_len >= 1e6
                 else f'Ref Gap: {ref_gap_len/1e3:.1f}kb')
        q_len = parent.get("query_len", 0)
        l_q   = (f'Asm Seq: {q_len/1e6:.2f}Mb' if q_len >= 1e6
                 else f'Asm Seq: {q_len/1e3:.1f}kb')
        ax.text((lx + rx) / 2, y_val, f"{l_ref}\n{l_q}",
                color='#D0021B', ha='center', va='center',
                fontweight='bold', fontsize=8.5,
                bbox=dict(facecolor='white', edgecolor='#D0021B', alpha=1.0,
                          boxstyle='round,pad=0.3', lw=1), zorder=6)

    # ── 收尾 ─────────────────────────────────────────────────────────────
    all_y = [Y_ASM] + list(chr_y_map.values())
    ax.set_xlim(label_x - x_span * 0.02, max_plot_x + x_span * 0.05)
    ax.set_ylim(min(all_y) - 1.2, max(all_y) + 1.5)
    ax.axis('off')
    ax.text(0.5, -0.05,
            f"* Flanking: {flank_len/1e3:.0f} kb"
            f" | Assembly width = query length, Ref width = source span"
            f" | Ref axis: genomic order (ribbon crossing = BACKBONE_FLIP, i.e. assembly traverses ref right→left)",
            transform=ax.transAxes, fontsize=9, color='#777777', style='italic', ha='center')

    plt.title(f"Structural Variant Architecture: {event_id}\nHaplotype: {hap_name}",
              fontsize=15, fontweight='bold', pad=30)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f">>> Plot saved: {output_path}")


if __name__ == "__main__":
    VCF_FILE    = "/io/wuxiaomia/sv_/my_tools/AsmSV/src/t2t.vcf"
    TARGET_EVENT = "haplotype2-0000141_CPLX_84"
    HAP_NAME    = "Hap2"

    try:
        parent_info, sub_blocks = parse_asmsv_complex_event(VCF_FILE, TARGET_EVENT, HAP_NAME)
        plot_parsed_complex_sv(parent_info, sub_blocks, TARGET_EVENT, HAP_NAME,
                               output_path=f"{TARGET_EVENT}_{HAP_NAME}_traced.png")
    except Exception as e:
        import traceback
        print(f"Error: {e}")
        traceback.print_exc()
