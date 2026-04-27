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

def parse_bnd_alt(rec, info):
    if "CHR2" in info and "POS2" in info:
        return str(safe_get(info, "CHR2", "")), int(safe_get(info, "POS2", 0))
    alt = rec.alts[0] if rec.alts else ""
    match = re.search(r'[\[\]](.+?):(\d+)[\[\]]', alt)
    if match: return match.group(1), int(match.group(2))
    return str(rec.chrom), int(rec.pos)

# =====================================================================
# 2. VCF 解析与【BND 图游走算法】
# =====================================================================
def parse_asmsv_complex_event(vcf_path, target_event, hap_name):
    if not os.path.exists(vcf_path): raise FileNotFoundError(f"VCF file not found: {vcf_path}")
    vcf = pysam.VariantFile(vcf_path)
    parent_event, blocks, bnds = None, [], []
    valid_gts = [(1, 1)] 
    if hap_name == "Hap1": valid_gts.append((1, 0))
    elif hap_name == "Hap2": valid_gts.append((0, 1))
    
    for rec in vcf:
        sample_name = list(rec.samples.keys())[0]
        if rec.samples[sample_name]['GT'] not in valid_gts: continue
            
        info = rec.info
        event_id = str(safe_get(info, "EVENT", ""))
        pevent_id = str(safe_get(info, "PEVENT", ""))
        svtype = str(safe_get(info, "SVTYPE", ""))
        true_end = rec.stop if rec.stop is not None else safe_get(info, "END", rec.pos + 1)
        
        if target_event in event_id.split(",") and svtype == "COMPLEX":
            parent_event = {
                "chrom": rec.chrom, "pos": int(rec.pos), 
                "end": int(true_end), 
                "query_len": int(safe_get(info, "REPLACED_QUERY_LEN", 0)),
                "strand": safe_get(info, "BACKBONE_STRAND", "++"),
                # 【修改点 1】：解析你新加的 UP_DOWN 字段
                "up_down": rec.info.get("INV_BACKBONE", False)
            }
        elif target_event in pevent_id.split(","):
            if svtype == "BND":
                c2, p2 = parse_bnd_alt(rec, info)
                bnds.append({"c1": str(rec.chrom), "p1": int(rec.pos), "c2": c2, "p2": p2})
            elif svtype in ["INV", "DUP", "INS"]:
                src_pos, src_end = safe_get(info, "SOURCE_POS", rec.pos), safe_get(info, "SOURCE_END", true_end)
                svlen = abs(int(safe_get(info, "SVLEN", src_end - src_pos)))
                blocks.append({
                    "type": svtype, "subtype": str(safe_get(info, "SUBTYPE", svtype)),
                    "len": svlen, "src_chr": str(safe_get(info, "SOURCE_CHR", rec.chrom)), 
                    "src_pos": int(src_pos), "src_end": int(src_end)
                })

    if not parent_event: raise ValueError(f"Event '{target_event}' not found.")
    return parent_event, blocks, bnds

def trace_alt_path_with_bnds(parent, blocks, bnds, tol=50):
    ordered_blocks = []
    unvisited = list(blocks)
    avail_bnds = list(bnds)
    curr_chr = parent["chrom"]
    curr_pos = parent["end"] if parent.get("up_down", False)  else parent["pos"]
    
    for _ in range(len(blocks) * 2): 
        if not unvisited: break
            
        bnd_idx, next_chr, next_pos = None, None, None
        for i, b in enumerate(avail_bnds):
            if b['c1'] == curr_chr and abs(b['p1'] - curr_pos) <= tol:
                bnd_idx, next_chr, next_pos = i, b['c2'], b['p2']; break
            elif b['c2'] == curr_chr and abs(b['p2'] - curr_pos) <= tol:
                bnd_idx, next_chr, next_pos = i, b['c1'], b['p1']; break
        
        if bnd_idx is not None:
            avail_bnds.pop(bnd_idx)
        else:
            next_chr, next_pos = curr_chr, curr_pos
            
        blk_idx, matched_blk, direction, exit_pos = None, None, None, None
        for i, blk in enumerate(unvisited):
            if blk['src_chr'] == next_chr:
                dist_fwd = abs(blk['src_pos'] - next_pos) 
                dist_rev = abs(blk['src_end'] - next_pos) 
                
                if dist_fwd <= tol and dist_fwd <= dist_rev:
                    blk_idx, matched_blk, direction, exit_pos = i, blk, 'forward', blk['src_end']
                    break
                elif dist_rev <= tol:
                    blk_idx, matched_blk, direction, exit_pos = i, blk, 'reverse', blk['src_pos']
                    break
                    
        if matched_blk is not None:
            unvisited.pop(blk_idx)
            matched_blk['is_inv_traversal'] = (direction == 'reverse')
            ordered_blocks.append(matched_blk)
            curr_chr, curr_pos = matched_blk['src_chr'], exit_pos
        else:
            if next_chr == parent["chrom"] and abs(next_pos - parent["end"]) <= tol: break
            break
            
    # 兜底方向逻辑
    primary_strand = parent.get("strand", "++")[1] if parent.get("up_down", False) else parent.get("strand", "++")[0]
    for blk in unvisited:
        name_indicates_inv = ('Inv' in blk.get('subtype', '') or blk.get('type') == 'INV')
        if primary_strand == '+':
            blk['is_inv_traversal'] = name_indicates_inv
        else:
            blk['is_inv_traversal'] = not name_indicates_inv
        ordered_blocks.append(blk)
        
    return ordered_blocks
# =====================================================================
# 3. 绘图引擎
# =====================================================================
def plot_parsed_complex_sv(parent, blocks, bnds, event_id, hap_name, output_path="auto_complex_plot.png"):
    fig, ax = plt.subplots(figsize=(16, 8), dpi=300)
    main_chr, flank_len = parent["chrom"], 50000
    color_palette = ["#4A90E2", "#50E3C2", "#F5A623", "#D0021B", "#9013FE"]

    ordered_core_blocks = trace_alt_path_with_bnds(parent, blocks, bnds)

    left_flank_pos = parent["pos"] - flank_len
    right_flank_pos = parent["end"]
    
    # 你的 VCF 修正保证了 bb_strands[0] 绝对是左物理断点，[1] 绝对是右物理断点
    bb_strands = parent.get("strand", "++")
    if len(bb_strands) == 1: bb_strands = bb_strands * 2
    left_is_inv = (bb_strands[0] == "-")
    right_is_inv = (bb_strands[1] == "-")
    
    up_down = parent.get("up_down", False)

    # 格式化所有核心块
    formatted_core = []
    for idx, b in enumerate(ordered_core_blocks):
        formatted_core.append({
            "name": b["subtype"], "len": b["len"], "chr": b["src_chr"], "pos": b["src_pos"],
            "is_inv": b.get("is_inv_traversal", False), 
            "color": color_palette[idx % len(color_palette)], "is_flank": False
        })

    left_flank = {"name": "Left_Flank", "len": flank_len, "chr": main_chr, "pos": left_flank_pos, "is_inv": left_is_inv, "color": "#A9A9A9", "is_flank": True}
    right_flank = {"name": "Right_Flank", "len": flank_len, "chr": main_chr, "pos": right_flank_pos, "is_inv": right_is_inv, "color": "#A9A9A9", "is_flank": True}

    all_blocks = []
    # 【修改点 3】：极度智能的 Assembly 轴模块排布
    if not up_down:
        # 正向组装：Assembly 从左向右走 (左侧翼 -> 核心 -> 右侧翼)
        all_blocks = [left_flank] + formatted_core + [right_flank]
    else:
        # 倒置组装：Assembly 从右向左走 (右侧翼 -> 核心 -> 左侧翼)
        all_blocks = [right_flank] + formatted_core + [left_flank]

    # --- 对齐与排版坐标 ---
    ref_tracks = defaultdict(list)
    for b in all_blocks: ref_tracks[b["chr"]].append(b)
        
    curr_asm_x = 0
    for b in all_blocks:
        b["asm_x"] = curr_asm_x
        curr_asm_x += b["len"]

    main_min_pos = min(b["pos"] for b in ref_tracks[main_chr])
    for chrom, t_blocks in ref_tracks.items():
        if chrom == main_chr:
            for b in t_blocks: b["ref_x"] = b["pos"] - main_min_pos
        else:
            first_b = min(t_blocks, key=lambda x: x["asm_x"])
            chr_offset = first_b["asm_x"] - first_b["pos"]
            for b in t_blocks: b["ref_x"] = b["pos"] + chr_offset

    alignment_shift = min(b["ref_x"] for b in all_blocks) - min(b["asm_x"] for b in all_blocks)
    for b in all_blocks: b["asm_x"] += alignment_shift

    # --- 绘图引擎开始 ---
    involved_chrs = [main_chr] + [c for c in ref_tracks.keys() if c != main_chr]
    Y_ASM, v_gap, bh = 0, 2.0, 0.35
    chr_y_map = {chrom: (len(involved_chrs) - i) * v_gap for i, chrom in enumerate(involved_chrs)}
    
    min_x = min(min((b["ref_x"] for b in all_blocks)), min((b["asm_x"] for b in all_blocks)))
    max_x = max(max((b["ref_x"] + b["len"] for b in all_blocks)), max((b["asm_x"] + b["len"] for b in all_blocks)))
    x_span = max_x - min_x
    label_x = min_x - x_span * 0.05 - 5000  
    
    ax.axhline(Y_ASM, color='black', linewidth=1.2, zorder=1, alpha=0.6)
    ax.text(label_x, Y_ASM, f"Assembly\n({hap_name})", va='center', ha='right', fontsize=11, fontweight='bold')
    
    for chrom, y_pos in chr_y_map.items():
        ax.axhline(y_pos, color='black', linewidth=1.2, zorder=1, alpha=0.6)
        label_prefix = "Main Ref" if chrom == main_chr else "Ref"
        ax.text(label_x, y_pos, f"{label_prefix}\n{chrom}", va='center', ha='right', fontsize=11, fontweight='bold')

    for b in all_blocks:
        w, color, is_inv, ref_y = b["len"], b["color"], b["is_inv"], chr_y_map[b["chr"]]
        zo = 2 if b.get("is_flank") else 3
        
        # Ref Block
        ax.add_patch(patches.Rectangle((b["ref_x"], ref_y - bh/2), w, bh, color=color, alpha=0.65, zorder=zo, ec='black', lw=0.8))
        ax.text(b["ref_x"] + w/2, ref_y + bh/2 + 0.1, f'{(b["pos"] + w/2) / 1e6:.2f}M', 
                ha='center', va='bottom', fontsize=8.5, color='#222222',
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=0.5), zorder=4)

        # Asm Block
        ax.add_patch(patches.Rectangle((b["asm_x"], Y_ASM - bh/2), w, bh, color=color, alpha=0.65, zorder=3, ec='black', lw=0.8))
        if w > x_span * 0.005: 
            dx = -w*0.6 if is_inv else w*0.6
            sx = b["asm_x"] + w*0.8 if is_inv else b["asm_x"] + w*0.2
            ax.arrow(sx, Y_ASM, dx, 0, head_width=0.08, head_length=min(w*0.1, 5000), fc='white', ec='white', zorder=4)
            
        if not b["is_flank"]:
            ax.text(b["asm_x"] + w/2, Y_ASM - bh/2 - 0.2, b["name"], 
                    ha='center', va='top', fontsize=9, fontweight='bold', color=color,
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=0.2))

        # Ribbon
        draw_bezier_ribbon(ax, b["ref_x"], w, ref_y - bh/2, b["asm_x"], w, Y_ASM + bh/2, color, alpha=0.35, invert=is_inv)

    # --- 绘制缺失红线 ---
    ref_missing_len = parent["end"] - parent["pos"]
    if ref_missing_len > 0:
        lx = parent["pos"] - main_min_pos
        rx = parent["end"] - main_min_pos
        y_val = chr_y_map[main_chr]
        
        ax.plot([lx, rx], [y_val, y_val], color='#D0021B', linestyle='--', zorder=5, lw=2)
        
        l_ref = f'Ref Gap: {ref_missing_len/1e6:.2f}Mb' if ref_missing_len >= 1e6 else f'Ref Gap: {ref_missing_len/1e3:.1f}kb'
        q_len = parent.get("query_len", 0)
        l_q = f'Asm Seq: {q_len/1e6:.2f}Mb' if q_len >= 1e6 else f'Asm Seq: {q_len/1e3:.1f}kb'
        
        ax.text((lx + rx)/2, y_val, f"{l_ref}\n{l_q}", 
                color='#D0021B', ha='center', va='center', fontweight='bold', fontsize=8.5,
                bbox=dict(facecolor='white', edgecolor='#D0021B', alpha=1.0, pad=1.0, boxstyle='round,pad=0.3', lw=1), zorder=6)

    # --- 收尾视觉美化 ---
    ax.set_xlim(min_x - x_span * 0.1 - 25000, max_x + x_span * 0.05)
    all_y = [Y_ASM] + list(chr_y_map.values())
    ax.set_ylim(min(all_y) - 1.2, max(all_y) + 1.5)
    ax.axis('off')
    
    ax.text(0.5, -0.05, f"* Default flanking sequence length: {flank_len/1e3:.0f} kb ", 
            transform=ax.transAxes, fontsize=9, color='#777777', style='italic', ha='center')

    plt.title(f"Structural Variant Architecture: {event_id}\nHaplotype: {hap_name}", fontsize=15, fontweight='bold', pad=30)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f">>> Plot saved: {output_path}")

if __name__ == "__main__":
    VCF_FILE = "t2t_test.vcf"  
    TARGET_EVENT = "haplotype1-0000010_CPLX_3"
    HAP_NAME = "Hap1"
    
    try:
        parent_info, sub_blocks, bnds = parse_asmsv_complex_event(VCF_FILE, TARGET_EVENT, HAP_NAME)
        plot_parsed_complex_sv(parent_info, sub_blocks, bnds, TARGET_EVENT, HAP_NAME, output_path=f"{TARGET_EVENT}_{HAP_NAME}_traced.png")
    except Exception as e:
        print(f"Error: {e}")