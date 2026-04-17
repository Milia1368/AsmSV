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
    """绘制连接 Assembly 轴和 Reference 轴的贝塞尔曲线带"""
    if not invert:
        verts = [(x_top, y_top), (x_top, y_bot + (y_top - y_bot)/2), (x_bot, y_top - (y_top - y_bot)/2),
                 (x_bot, y_bot), (x_bot + w_bot, y_bot), (x_bot + w_bot, y_top - (y_top - y_bot)/2),
                 (x_top + w_top, y_bot + (y_top - y_bot)/2), (x_top + w_top, y_top), (x_top, y_top)]
    else:
        # 倒位交叉 (呈现 X 型)
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
    if match:
        return match.group(1), int(match.group(2))
    return str(rec.chrom), int(rec.pos)

# =====================================================================
# 2. VCF 解析与【BND 图游走算法】
# =====================================================================
def parse_asmsv_complex_event(vcf_path, target_event, hap_name):
    if not os.path.exists(vcf_path): raise FileNotFoundError(f"VCF file not found: {vcf_path}")
    vcf = pysam.VariantFile(vcf_path)
    parent_event, blocks, bnds = None, [], []
    valid_gts = [(1, 1)] 
    if hap_name == "Hap1":
        valid_gts.append((1, 0))
    elif hap_name == "Hap2":
        valid_gts.append((0, 1))
    
    for rec in vcf:
        # 【核心过滤逻辑】：直接检查你的 SUPPORT 标签
        sample_name = list(rec.samples.keys())[0]
        gt_tuple = rec.samples[sample_name]['GT']
        
        # 过滤非当前单倍型的变异
        if gt_tuple not in valid_gts:
            continue
        info = rec.info
        event_id = str(safe_get(info, "EVENT", ""))
        pevent_id = str(safe_get(info, "PEVENT", ""))
        svtype = str(safe_get(info, "SVTYPE", ""))
        
       # 【终极防弹逻辑】：提到最前面，全局计算真实 END
        true_end = rec.stop if rec.stop is not None else safe_get(info, "END", rec.pos + 1)
        
        # 将原来的精确等于 `==` 修改为逗号分隔后的成员判断 `in`
        if target_event in event_id.split(",") and svtype == "COMPLEX":
            parent_event = {
                "chrom": rec.chrom, "pos": int(rec.pos), 
                "end": int(true_end),  # 直接使用
                "replaced_len": int(safe_get(info, "REPLACED_REF_LEN", 0)),
                "strand": safe_get(info, "BACKBONE_STRAND", "+") 
            }
        elif target_event in pevent_id.split(","):  # <--- 这里的改动最关键
            if svtype == "BND":
                c2, p2 = parse_bnd_alt(rec, info)
                bnds.append({"c1": str(rec.chrom), "p1": int(rec.pos), "c2": c2, "p2": p2})
            elif svtype in ["INV", "DUP", "INS"]:
                source_pos = safe_get(info, "SOURCE_POS", rec.pos)
                
                # 【关键修改】：兜底时直接使用上方的 true_end，不再去 info 里捞 "END"
                source_end = safe_get(info, "SOURCE_END", true_end)
                
                svlen = abs(int(safe_get(info, "SVLEN", source_end - source_pos)))
                blocks.append({
                    "type": svtype, "subtype": str(safe_get(info, "SUBTYPE", svtype)),
                    "len": svlen, "src_chr": str(safe_get(info, "SOURCE_CHR", rec.chrom)), 
                    "src_pos": int(source_pos), "src_end": int(source_end)
                })

    if not parent_event: raise ValueError(f"Event '{target_event}' not found.")
    return parent_event, blocks, bnds

def trace_alt_path_with_bnds(parent, blocks, bnds, tol=50):
    """通过 BND 连接完美还原 ALT 组装序列真实的组成顺序与方向"""
    ordered_blocks = []
    unvisited = list(blocks)
    avail_bnds = list(bnds)
    curr_chr = parent["chrom"]
    curr_pos = parent["end"] if parent.get("strand") == "-" else parent["pos"]
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
    backbone_strand = parent.get("strand", "+")
    for blk in unvisited:
        # 判断它的 SUBTYPE 名字里是否带有 INV
        name_indicates_inv = ('Inv' in blk.get('subtype', '') or blk.get('type') == 'INV')
        
        if backbone_strand == '+':
            # 主干是正链：带有 INV 才是反向游走 (<-)
            blk['is_inv_traversal'] = name_indicates_inv
        else:
            # 【负链魔法】：主干是负链！
            # 不发生 INV 的正常块，是顺着负链反向游走的 (<-)
            # 发生了 INV 的块，负负得正，反而变成了正链方向 (->)
            blk['is_inv_traversal'] = not name_indicates_inv
            
        ordered_blocks.append(blk)
        
    return ordered_blocks

# =====================================================================
# 3. 绘图引擎
# =====================================================================
def plot_parsed_complex_sv(parent, blocks, bnds, event_id, hap_name, output_path="auto_complex_plot.png"):
    fig, ax = plt.subplots(figsize=(16, 8), dpi=300)
    main_chr = parent["chrom"]
    flank_len = 50000
    color_palette = ["#4A90E2", "#50E3C2", "#F5A623", "#D0021B", "#9013FE"]

    ordered_core_blocks = trace_alt_path_with_bnds(parent, blocks, bnds)

    # --- 智能判定主干正负链的侧翼挂载逻辑 ---
    backbone_strand = parent.get("strand", "+")
    
    if backbone_strand == "+":
        left_flank_pos = parent["pos"] - flank_len
        right_flank_pos = parent["end"]
        flank_is_inv = False
    else:
        # 【负链魔法】：Assembly 的前端(Left)对应 REF 的右侧，末端(Right)对应 REF 的左侧！
        left_flank_pos = parent["end"]
        right_flank_pos = parent["pos"] - flank_len
        flank_is_inv = True # 负链侧翼必须呈现倒位交叉！

    all_blocks = []
    # 起始侧翼 (Left Flank)
    all_blocks.append({"name": "Left_Flank", "len": flank_len, "chr": main_chr, "pos": left_flank_pos, "is_inv": flank_is_inv, "color": "#A9A9A9", "is_flank": True})
    
    # 核心变异块 (顺序由 BND 游走算法决定)
    for idx, b in enumerate(ordered_core_blocks):
        all_blocks.append({
            "name": b["subtype"], "len": b["len"], "chr": b["src_chr"], "pos": b["src_pos"],
            "is_inv": b.get("is_inv_traversal", False), 
            "color": color_palette[idx % len(color_palette)], "is_flank": False
        })
        
    # 终止侧翼 (Right Flank)
    all_blocks.append({"name": "Right_Flank", "len": flank_len, "chr": main_chr, "pos": right_flank_pos, "is_inv": flank_is_inv, "color": "#A9A9A9", "is_flank": True})

    # ... 下方的坐标计算和绘图代码保持原样完全不变 ...
    # --- 阶段 A：计算 ALT (Assembly) 轴坐标 (连续排布) ---
    curr_asm_x = 0
    for b in all_blocks:
        b["asm_x"] = curr_asm_x
        curr_asm_x += b["len"]

    # --- 阶段 B：计算 REF (Reference) 轴坐标 (物理嵌套映射) ---
    ref_tracks = defaultdict(list)
    for b in all_blocks:
        ref_tracks[b["chr"]].append(b)
        
    for chrom, t_blocks in ref_tracks.items():
        # 1. 获取所有物理区间并排序
        intervals = [(b["pos"], b["pos"] + b["len"]) for b in t_blocks]
        intervals.sort(key=lambda x: x[0])
        
        # 2. 合并重叠或距离极近 (5kb以内) 的物理区间，构成真实连续的参考片段
        merged = []
        for iv in intervals:
            if not merged:
                merged.append(list(iv))
            else:
                last = merged[-1]
                if iv[0] <= last[1] + 5000:  
                    last[1] = max(last[1], iv[1])
                else:
                    merged.append(list(iv))
                    
        # 3. 为合并后的每个大片段分配视觉起始 X 坐标
        region_visual_start = {}
        curr_ref_x = 20000 
        for reg in merged:
            region_visual_start[reg[0]] = curr_ref_x
            curr_ref_x += (reg[1] - reg[0]) + 30000  # 不同大片段之间间隔 30kb
            
        # 4. 根据块的真实物理位置，反推它在视觉轴上的精确位置
        for b in t_blocks:
            for reg in merged:
                if reg[0] <= b["pos"] <= reg[1]:
                    # 按照 1:1 物理缩放嵌入到模块中
                    b["ref_x"] = region_visual_start[reg[0]] + (b["pos"] - reg[0])
                    break

    # --- 阶段 C：开始画图 ---
    involved_chrs = [main_chr] + [c for c in ref_tracks.keys() if c != main_chr]
    Y_ASM = len(involved_chrs) * 1.5
    chr_y_map = {chrom: Y_ASM - (i + 1) * 1.5 for i, chrom in enumerate(involved_chrs)}
    bh = 0.3
    
    ax.axhline(Y_ASM, color='black', linewidth=1.5, zorder=1)
    ax.text(-10000, Y_ASM, f"Assembly ({hap_name})", va='center', ha='right', fontsize=12, fontweight='bold')
    for chrom, y_pos in chr_y_map.items():
        ax.axhline(y_pos, color='black', linewidth=1.5, zorder=1)
        ax.text(-10000, y_pos, f"Ref {chrom}", va='center', ha='right', fontsize=12, fontweight='bold')
    # 用于记录已经在 REF 轴上标注过的坐标位置，防止相邻块断点文字重叠
    drawn_ref_labels = set()

   
    for b in all_blocks:
        w, color, is_inv = b["len"], b["color"], b["is_inv"]
        
        # 1. 绘制 Assembly 块
        ax.add_patch(patches.Rectangle((b["asm_x"], Y_ASM - bh/2), w, bh, color=color, zorder=3, ec='black'))
        if is_inv:
            ax.arrow(b["asm_x"] + w*0.8, Y_ASM, -w*0.6, 0, head_width=0.1, head_length=min(w*0.1, 5000), fc='white', ec='white', zorder=4)
        else:
            ax.arrow(b["asm_x"] + w*0.2, Y_ASM, w*0.6, 0, head_width=0.1, head_length=min(w*0.1, 5000), fc='white', ec='white', zorder=4)
        if not b["is_flank"]:
            ax.text(b["asm_x"] + w/2, Y_ASM + bh + 0.1, b["name"], ha='center', fontsize=9, fontweight='bold', color=color)

        # 2. 绘制 Reference 块 
        ref_y = chr_y_map[b["chr"]]
        # 巧思：如果是 Flank，设为底层(zorder=2)；如果是核心 SV 块，覆盖在上层(zorder=3)
        zo = 2 if b.get("is_flank") else 3
        ax.add_patch(patches.Rectangle((b["ref_x"], ref_y - bh/2), w, bh, color=color, zorder=zo, ec='black'))
        
        # 为了防止重叠，给文本加上了倾斜角度
        ref_label = f'{b["pos"]/1e6:.2f}M'
        ax.text(b["ref_x"] + w/2, ref_y - bh - 0.1, ref_label, ha='center', va='top', fontsize=9, rotation=20)

        # 3. 绘制相连 Ribbon
        draw_bezier_ribbon(ax, b["asm_x"], w, Y_ASM - bh/2, b["ref_x"], w, ref_y + bh/2, color, invert=is_inv)

    # 画巨型缺失红线
    if parent["replaced_len"] > 0:
        main_chr_blocks = sorted(ref_tracks[main_chr], key=lambda x: x["pos"])
        left_flank = next((b for b in main_chr_blocks if b["name"] == "Left_Flank"), None)
        right_flank = next((b for b in main_chr_blocks if b["name"] == "Right_Flank"), None)
        if left_flank and right_flank:
            del_start = left_flank["ref_x"] + left_flank["len"]
            del_end = right_flank["ref_x"]
            if del_end > del_start:
                ax.plot([del_start, del_end], [chr_y_map[main_chr], chr_y_map[main_chr]], color='red', linestyle='--', zorder=5, lw=2)
                ax.text((del_start + del_end)/2, chr_y_map[main_chr] + 0.1, f'{parent["replaced_len"]/1e6:.1f}Mb Replaced', color='red', ha='center', fontweight='bold', fontsize=9)

    max_x = max(curr_asm_x, max((b["ref_x"] + b["len"] for b in all_blocks)))
    ax.set_xlim(-20000, max_x + 20000)
    ax.set_ylim(min(chr_y_map.values()) - 1.5, Y_ASM + 1.5)
    ax.axis('off')
    ax.text(0.01, 0.02, f"* Note: Default display length for flanking reference sequences (Flanks): {flank_len} bp", 
            transform=ax.transAxes, fontsize=10, color='#666666', style='italic')

    plt.title(f"Automated Sequence Path Map: {event_id} ({hap_name})", fontsize=14, fontweight='bold', pad=20)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Plot correctly generated for {hap_name} and saved to {output_path}")
if __name__ == "__main__":
    VCF_FILE = "t2t_test.vcf"  
    TARGET_EVENT = "haplotype1-0000021_CPLX_24"
    HAP_NAME = "Hap1"
    
    try:
        parent_info, sub_blocks, bnds = parse_asmsv_complex_event(VCF_FILE, TARGET_EVENT, HAP_NAME)
        plot_parsed_complex_sv(parent_info, sub_blocks, bnds, TARGET_EVENT, HAP_NAME, output_path=f"{TARGET_EVENT}_{HAP_NAME}_traced.png")
    except Exception as e:
        print(f"Error: {e}")
