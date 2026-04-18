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
# --- 阶段 C：开始画图 (顶部 Ref, 底部 Asm) ---
    # 重新整理染色体顺序：Main Chr 在列表首位
    involved_chrs = [main_chr] + [c for c in ref_tracks.keys() if c != main_chr]
    num_refs = len(involved_chrs)
    
    # 1. 坐标轴高度定义
    Y_ASM = 0  # Assembly 固定在最底部
    
    # 计算 Ref 轨道坐标：
    # 列表首位 (Main Chr) 分配最高 Y 值，后续 Chr 依次减小，但都大于 Y_ASM
    # 间距设为 2.0 以留出足够的 Ribbon 摆动空间和标注空间
    v_gap = 2.0 
    chr_y_map = {}
    for i, chrom in enumerate(involved_chrs):
        # 比如有 2 条 Ref: i=0 (Main) Y=4.0; i=1 (Other) Y=2.0
        chr_y_map[chrom] = (num_refs - i) * v_gap

    bh = 0.3 # 色块厚度
    
    # 绘制横线和标签
    ax.axhline(Y_ASM, color='black', linewidth=1.2, zorder=1, alpha=0.6)
    ax.text(-15000, Y_ASM, f"Assembly\n({hap_name})", va='center', ha='right', fontsize=11, fontweight='bold')
    
    for chrom, y_pos in chr_y_map.items():
        ax.axhline(y_pos, color='black', linewidth=1.2, zorder=1, alpha=0.6)
        label_prefix = "Main Ref" if chrom == main_chr else "Ref"
        ax.text(-15000, y_pos, f"{label_prefix}\n{chrom}", va='center', ha='right', fontsize=11, fontweight='bold')

    # 遍历绘制所有块
    for b in all_blocks:
        w, color, is_inv = b["len"], b["color"], b["is_inv"]
        
        # -------------------------------------------------------
        # 2. 绘制 Reference 块 (上方)
        # -------------------------------------------------------
        ref_y = chr_y_map[b["chr"]]
        zo = 2 if b.get("is_flank") else 3
        # 绘制色块
        ax.add_patch(patches.Rectangle((b["ref_x"], ref_y - bh/2), w, bh, color=color, zorder=zo, ec='black', lw=0.8))
        
        # 标注物理坐标：在 Ref 块的正上方，纯水平排布
        ax_pos = b["pos"] + (b["len"] / 2)
        ref_label = f'{ax_pos / 1e6:.2f}M'
        
        # 【关键修改】：去掉了 rotation，文字变为纯水平
        # 增加了半透明白底 bbox，即使偶尔密集也能互相遮挡，保持数字清晰
        ax.text(b["ref_x"] + w/2, ref_y + bh/2 + 0.1, ref_label, 
                ha='center', va='bottom', fontsize=8.5, color='#222222',
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=0.5),
                zorder=4)

        # -------------------------------------------------------
        # 3. 绘制 Assembly 块 (最底部)
        # -------------------------------------------------------
        ax.add_patch(patches.Rectangle((b["asm_x"], Y_ASM - bh/2), w, bh, color=color, zorder=3, ec='black', lw=0.8))
        
        # 绘制方向箭头
        if is_inv:
            ax.arrow(b["asm_x"] + w*0.8, Y_ASM, -w*0.6, 0, head_width=0.08, head_length=min(w*0.1, 5000), fc='white', ec='white', zorder=4)
        else:
            ax.arrow(b["asm_x"] + w*0.2, Y_ASM, w*0.6, 0, head_width=0.08, head_length=min(w*0.1, 5000), fc='white', ec='white', zorder=4)
            
        # 标注变异名称：在 Assembly 块的正下方
        if not b["is_flank"]:
            ax.text(b["asm_x"] + w/2, Y_ASM - bh/2 - 0.2, b["name"], 
                    ha='center', va='top', fontsize=9, fontweight='bold', color=color,
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=0.2))

        # -------------------------------------------------------
        # 4. 绘制 Ribbon (连接上方的 Ref 和底部的 Asm)
        # -------------------------------------------------------
        draw_bezier_ribbon(ax, 
                           b["ref_x"], w, ref_y - bh/2,   # 上：Ref 的底部
                           b["asm_x"], w, Y_ASM + bh/2,   # 下：Asm 的顶部
                           color, invert=is_inv)

    # -------------------------------------------------------
    # 5. 绘制缺失红线 (在 Main Chr 轨道上)
    # -------------------------------------------------------
    if parent.get("replaced_len", 0) > 0:
        main_blocks = sorted(ref_tracks[main_chr], key=lambda x: x["pos"])
        left_f = next((b for b in main_blocks if b["name"] == "Left_Flank"), None)
        right_f = next((b for b in main_blocks if b["name"] == "Right_Flank"), None)
        
        if left_f and right_f:
            lx = left_f["ref_x"] + left_f["len"]
            rx = right_f["ref_x"]
            
            if rx > lx:
                y_val = chr_y_map[main_chr]
                # 绘制红色虚线
                ax.plot([lx, rx], [y_val, y_val], color='#D0021B', linestyle='--', zorder=5, lw=2)
                
                # 计算显示文本 (根据长度自动选择 Mb 或 kb)
                missing_val = parent["replaced_len"]
                label_text = f'Missing: {missing_val/1e6:.2f}Mb' if missing_val >= 1e6 else f'Missing: {missing_val/1e3:.1f}kb'
                
                # --- 关键修改点：va='center', Y坐标设为 y_val ---
                ax.text((lx + rx)/2, y_val, label_text, 
                        color='#D0021B', ha='center', va='center', 
                        fontweight='bold', fontsize=9,
                        # 使用白色背景 bbox 遮挡其下方的虚线，达到嵌入效果
                        bbox=dict(facecolor='white', edgecolor='#D0021B', 
                                  alpha=1.0, pad=1.5, boxstyle='round,pad=0.3', lw=1),
                        zorder=6)

    # --- 阶段 D：收尾与视觉美化 ---
    max_x = max(curr_asm_x, max((b["ref_x"] + b["len"] for b in all_blocks)))
    ax.set_xlim(-max_x*0.05, max_x * 1.05)
    
    # 动态计算 Y 轴范围，确保所有轨道和标注都能看到
    all_y = [Y_ASM] + list(chr_y_map.values())
    ax.set_ylim(min(all_y) - 1.2, max(all_y) + 1.5)
    
    ax.axis('off')
    
    # 添加图注说明
    ax.text(0.5, -0.05, f"* Note: Default display length for flanking reference sequences (Flanks): {flank_len/1e3:.0f} kb", 
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
