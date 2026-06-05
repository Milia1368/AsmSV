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
    解析新版 bam_info.py 生成的 COMPLEX 事件 VCF。

    Summary record (FILTER=COMPLEX):
      SVTYPE=CSV, EVENT=<target_event>, SUBTYPE=COMPLEX
      关键 INFO: END, INS_LEN (dist_q), DEL_LEN, BACKBONE_STRAND, BACKBONE_FLIP

    Sub-block records (每个非 TANDEM 孤儿簇):
      SVTYPE=CSV, EVENT=<target_event>_sub<N>
      关键 INFO: CPLX_IDX, SRCREG=chr:start-end, SUBTYPE, INS_LEN (true_q_len)

    DUP:TANDEM 块无 span 记录，从 BND entry 对重建坐标。
    """
    if not os.path.exists(vcf_path):
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    vcf = pysam.VariantFile(vcf_path)
    parent_event = None
    sub_blocks   = {}   # cplx_idx (int) → block dict
    bnd_entries  = {}   # cplx_idx (int) → {from_chr, from_pos, to_chr, to_pos, to_strand}
    tail_bnd_chr = None
    tail_bnd_pos = None

    valid_gts = [(1, 1)]
    if hap_name == "Hap1":   valid_gts.append((1, 0))
    elif hap_name == "Hap2": valid_gts.append((0, 1))

    sub_prefix   = target_event + "_sub"
    tail_bnd_key = target_event + "_tail_bnd1"

    for rec in vcf:
        sample_name = list(rec.samples.keys())[0]
        if rec.samples[sample_name]['GT'] not in valid_gts:
            continue

        info       = rec.info
        event_id   = str(safe_get(info, "EVENT", ""))
        svtype_raw = str(safe_get(info, "SVTYPE", ""))
        subtype    = str(safe_get(info, "SUBTYPE", ""))

        # ── Summary record ────────────────────────────────────────────────
        if event_id == target_event and subtype == "COMPLEX":
            true_end = rec.stop if rec.stop is not None else int(safe_get(info, "END", rec.pos + 1))
            parent_event = {
                "chrom":     rec.chrom,
                "pos":       int(rec.pos),
                "end":       int(true_end),
                "query_len": int(safe_get(info, "INS_LEN", 0)),
                "strand":    str(safe_get(info, "BACKBONE_STRAND", "++")),
                "up_down":   _coerce_bool(safe_get(info, "BACKBONE_FLIP", False)),
            }

        # ── Sub-block span records (non-TANDEM orphan clusters) ───────────
        elif event_id.startswith(sub_prefix) and "CPLX_IDX" in info:
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

        # ── BND entry records（用于重建 TANDEM 坐标）────────────────────────
        elif svtype_raw == "BND":
            m = re.match(rf'{re.escape(sub_prefix)}(\d+)_bnd_entry$', event_id)
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
            elif event_id == tail_bnd_key:
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
        entry    = bnd_entries[idx]
        src_chr  = entry["to_chr"]
        src_enter = entry["to_pos"]
        m_strand  = entry["to_strand"]

        # 退出坐标 = 下一个 entry BND 在该染色体上的出发点
        src_leave = src_enter
        next_entry = bnd_entries.get(idx + 1)
        if next_entry and next_entry["from_chr"] == src_chr:
            src_leave = next_entry["from_pos"]
        elif tail_bnd_chr == src_chr and tail_bnd_pos is not None:
            src_leave = tail_bnd_pos

        ref_len = abs(src_leave - src_enter)
        is_inv  = (m_strand != anchor_strand)

        sub_blocks[idx] = {
            "subtype":  ("Inv" if is_inv else "") + "DUP:TANDEM",
            "src_chr":  src_chr,
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
    MAX_VIS_GAP = 50000
    color_palette = ["#4A90E2", "#50E3C2", "#F5A623", "#D0021B", "#9013FE"]

    up_down    = parent.get("up_down", False)
    bb_strands = parent.get("strand", "++")
    if len(bb_strands) < 2: bb_strands += bb_strands

    # BACKBONE_STRAND 编码规则（来自 bam_info.py）:
    #   up_down=False → [0]=左断点strand, [1]=右断点strand
    #   up_down=True  → [0]=右断点strand, [1]=左断点strand
    if up_down:
        left_is_inv  = (bb_strands[1] == "-")
        right_is_inv = (bb_strands[0] == "-")
    else:
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

    # Assembly 方向决定顺序
    all_blocks = ([right_flank] + core_blocks + [left_flank]
                  if up_down else
                  [left_flank] + core_blocks + [right_flank])

    # ── Assembly 轴坐标（以 q_len 步进）──────────────────────────────────
    curr_asm_x = 0
    for b in all_blocks:
        b["asm_x"]  = curr_asm_x
        curr_asm_x += b["q_len"]

    # ── Reference 轴坐标（空白压缩）──────────────────────────────────────
    ref_tracks = defaultdict(list)
    for b in all_blocks:
        ref_tracks[b["src_chr"]].append(b)

    for chrom, t_blocks in ref_tracks.items():
        intervals = sorted(
            [[b["src_pos"], b["src_pos"] + b["ref_len"]] for b in t_blocks],
            key=lambda x: x[0],
        )
        merged_cov = []
        for iv in intervals:
            if not merged_cov or merged_cov[-1][1] < iv[0]:
                merged_cov.append(list(iv))
            else:
                merged_cov[-1][1] = max(merged_cov[-1][1], iv[1])

        shrink_zones = []
        for i in range(len(merged_cov) - 1):
            gap_len = merged_cov[i + 1][0] - merged_cov[i][1]
            if gap_len > MAX_VIS_GAP:
                shrink_zones.append((merged_cov[i][1], merged_cov[i + 1][0], gap_len - MAX_VIS_GAP))

        def map_pos(p, _sz=shrink_zones):
            shift = 0
            for sz_s, sz_e, shrink in _sz:
                if p > sz_e:
                    shift += shrink
                elif sz_s < p <= sz_e:
                    shift += shrink * (p - sz_s) / (sz_e - sz_s)
            return p - shift

        for b in t_blocks:
            b["visual_pos"] = map_pos(b["src_pos"])
            b["visual_len"] = map_pos(b["src_pos"] + b["ref_len"]) - b["visual_pos"]

        center_asm = sum(b["asm_x"] + b["q_len"] / 2   for b in t_blocks) / len(t_blocks)
        center_vis = sum(b["visual_pos"] + b["visual_len"] / 2 for b in t_blocks) / len(t_blocks)
        shift = center_asm - center_vis

        for b in t_blocks:
            b["ref_x"]      = b["visual_pos"] + shift
            b["len_on_ref"] = b["visual_len"]

    # 归一化最左起点到 0
    min_x = min(min(b["ref_x"] for b in all_blocks), min(b["asm_x"] for b in all_blocks))
    for b in all_blocks:
        b["ref_x"]  -= min_x
        b["asm_x"]  -= min_x

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
    ref_gap_len = parent["end"] - parent["pos"]
    if ref_gap_len > 0:
        lf = next((b for b in all_blocks if b["name"] == "Left_Flank"), None)
        rf = next((b for b in all_blocks if b["name"] == "Right_Flank"), None)
        if lf and rf:
            lx    = lf["ref_x"] + lf["len_on_ref"]
            rx    = rf["ref_x"]
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
            f"* Flanking: {flank_len/1e3:.0f} kb | Gaps > {MAX_VIS_GAP/1e3:.0f} kb compressed"
            f" | Assembly width = query length, Ref width = source span",
            transform=ax.transAxes, fontsize=9, color='#777777', style='italic', ha='center')

    plt.title(f"Structural Variant Architecture: {event_id}\nHaplotype: {hap_name}",
              fontsize=15, fontweight='bold', pad=30)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f">>> Plot saved: {output_path}")


if __name__ == "__main__":
    VCF_FILE    = "t2t_test.vcf"
    TARGET_EVENT = "haplotype1-0000010_CPLX_3"
    HAP_NAME    = "Hap1"

    try:
        parent_info, sub_blocks = parse_asmsv_complex_event(VCF_FILE, TARGET_EVENT, HAP_NAME)
        plot_parsed_complex_sv(parent_info, sub_blocks, TARGET_EVENT, HAP_NAME,
                               output_path=f"{TARGET_EVENT}_{HAP_NAME}_traced.png")
    except Exception as e:
        import traceback
        print(f"Error: {e}")
        traceback.print_exc()
