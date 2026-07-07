#!/usr/bin/env python3
"""
Plot AsmSV COMPLEX region architecture from VCF *_SUMMARY records (SUB_SVS field).
Updated with multi-track linear plotting, specific Ref Span indicator,
and strictly Reference-anchored Left-to-Right orientation mapping (Reverse projection for negative strand).

Usage:
  python plot.py -i results.vcf --list
  python plot.py -i results.vcf -e haplotype1-0000097_CPLX_2
  python plot.py -i results.vcf -e haplotype1-0000097_CPLX_2 --hap Hap1 -o out.png

Requires: pysam, matplotlib
"""
import argparse
import os
import sys

import pysam

VALID_GTS = {
    "Hap1": {(1, 1), (1, 0)},
    "Hap2": {(1, 1), (0, 1)},
}


def _require_matplotlib():
    try:
        import matplotlib.patches as patches
        import matplotlib.pyplot as plt
        from matplotlib.path import Path
    except ImportError as exc:
        sys.exit(
            "matplotlib is required for plotting. Install with: pip install matplotlib\n"
            f"Original error: {exc}"
        )
    return patches, plt, Path


def draw_bezier_ribbon(ax, x_top, w_top, y_top, x_bot, w_bot, y_bot, color, alpha=0.4, invert=False):
    """绘制连接 Top (Ref) 和 Bottom (Assembly) 的贝塞尔丝带"""
    patches, _, Path = _require_matplotlib()
    if not invert:
        verts = [(x_top, y_top), (x_top, y_bot + (y_top - y_bot) / 2), (x_bot, y_top - (y_top - y_bot) / 2),
                 (x_bot, y_bot), (x_bot + w_bot, y_bot), (x_bot + w_bot, y_top - (y_top - y_bot) / 2),
                 (x_top + w_top, y_bot + (y_top - y_bot) / 2), (x_top + w_top, y_top), (x_top, y_top)]
    else:  # 负链/倒位导致交叉反转
        verts = [(x_top, y_top), (x_top, y_bot + (y_top - y_bot) / 2), (x_bot + w_bot, y_top - (y_top - y_bot) / 2),
                 (x_bot + w_bot, y_bot), (x_bot, y_bot), (x_bot, y_top - (y_top - y_bot) / 2),
                 (x_top + w_top, y_bot + (y_top - y_bot) / 2), (x_top + w_top, y_top), (x_top, y_top)]

    codes = [Path.MOVETO] + [Path.CURVE4] * 3 + [Path.LINETO] + [Path.CURVE4] * 3 + [Path.CLOSEPOLY]
    ax.add_patch(patches.PathPatch(Path(verts, codes), facecolor=color, lw=0, alpha=alpha, zorder=2))


def safe_get(info, key, default_val):
    if key not in info:
        return default_val
    val = info.get(key)
    if val is None:
        return default_val
    if isinstance(val, (tuple, list)):
        return val[0] if len(val) > 0 and val[0] is not None else default_val
    return val


def _parse_sub_svs(raw, parent_event):
    if isinstance(raw, (tuple, list)):
        sub_svs_str = ",".join(map(str, raw))
    else:
        sub_svs_str = str(raw) if raw is not None else ""
    if not sub_svs_str.strip():
        return []

    blocks = []
    all_qs = []
    all_qe = []
    bb_q_start = None
    bb_q_end = None

    for item in sub_svs_str.split("|"):
        item = item.strip()
        if not item:
            continue
        parts = item.split(",")
        if len(parts) != 7:
            continue
        b_type, c_chr, rs, re, qs, qe, strand = parts
        try:
            rs_i, re_i, qs_i, qe_i = int(rs), int(re), int(qs), int(qe)
        except ValueError:
            continue

        src_pos, src_end = min(rs_i, re_i), max(rs_i, re_i)
        q_start, q_end = min(qs_i, qe_i), max(qs_i, qe_i)

        all_qs.append(q_start)
        all_qe.append(q_end)

        if b_type == "backbone":
            bb_q_start = q_start
            bb_q_end = q_end
            parent_event["strand"] = strand
            continue

        blocks.append({
            "subtype": b_type,
            "src_chr": c_chr,
            "src_pos": src_pos,
            "src_end": src_end,
            "qs": q_start,
            "qe": q_end,
            "strand": strand,
            "is_inv": strand == "-",
            "q_len": q_end - q_start,
            "ref_len": src_end - src_pos,
        })

    # 原生空间排序：始终按照 Contig 物理方向从低到高排列
    blocks.sort(key=lambda x: x["qs"])

    if all_qs and all_qe:
        parent_event["qs_left"] = min(all_qs)
        parent_event["qs_right"] = max(all_qe)
    else:
        parent_event["qs_left"] = bb_q_start if bb_q_start is not None else 0
        parent_event["qs_right"] = bb_q_end if bb_q_end is not None else 1000

    return blocks


def list_complex_events(vcf_path):
    if not os.path.isfile(vcf_path):
        raise FileNotFoundError(f"VCF not found: {vcf_path}")
    events = []
    vcf = pysam.VariantFile(vcf_path)
    for rec in vcf:
        if not rec.id or not rec.id.endswith("_SUMMARY"):
            continue
        if str(safe_get(rec.info, "SUBTYPE", "")).upper() != "COMPLEX":
            continue
        if not rec.info.get("SUB_SVS"):
            continue
        events.append(str(safe_get(rec.info, "EVENT", rec.id.replace("_SUMMARY", ""))))
    vcf.close()
    return events


def parse_asmsv_complex_event(vcf_path, target_event, hap_name):
    if hap_name not in VALID_GTS:
        raise ValueError(f"Invalid haplotype '{hap_name}'. Use Hap1 or Hap2.")
    if not os.path.isfile(vcf_path):
        raise FileNotFoundError(f"VCF not found: {vcf_path}")

    vcf = pysam.VariantFile(vcf_path)
    samples = list(vcf.header.samples)
    if not samples:
        vcf.close()
        raise ValueError("VCF has no sample columns.")
    sample_name = samples[0]

    parent_event = None
    blocks = []
    te_upper = target_event.upper()

    for rec in vcf:
        if rec.samples[sample_name]["GT"] not in VALID_GTS[hap_name]:
            continue
        if str(safe_get(rec.info, "EVENT", "")).upper() != te_upper:
            continue
        if not (rec.id and rec.id.endswith("_SUMMARY")):
            continue

        true_end = rec.stop if rec.stop is not None else int(safe_get(rec.info, "END", rec.pos + 1))
        parent_event = {
            "chrom": rec.chrom,
            "pos": int(rec.pos),
            "end": int(true_end),
            "query_len": int(safe_get(rec.info, "INS_LEN", 0)),
        }
        blocks = _parse_sub_svs(rec.info.get("SUB_SVS", ""), parent_event)
        break

    vcf.close()

    if not parent_event:
        raise ValueError(f"COMPLEX event '{target_event}' SUMMARY not found for {hap_name}.")
    return parent_event, blocks


def plot_parsed_complex_sv(parent, blocks, event_id, hap_name, output_path):
    patches, plt, _ = _require_matplotlib()
    flank_len = 50000
    main_chr = parent["chrom"]
    color_palette = ["#4A90E2", "#50E3C2", "#F5A623", "#D0021B", "#9013FE", "#8B572A"]

    qs_left = parent.get("qs_left", 0)
    qs_right = parent.get("qs_right", parent.get("query_len", 1000))
    bb_strand = parent.get("strand", "+")
    flip_asm = (bb_strand == "-")  # 核心变量：如果主骨架是负链，Assembly 将逆序绘制，贴合 Reference 坐标

    # 根据正反链严格定义 Flank 锚点：保证 Reference 从左到右恒为递增
    if not flip_asm:
        lf_ref = parent["pos"] - flank_len
        rf_ref = parent["end"]
        lf_strand, rf_strand = "+", "+"
    else:
        lf_ref = parent["end"]
        rf_ref = parent["pos"] - flank_len
        lf_strand, rf_strand = "-", "-"

    # flank_qs_low 对应 Assembly (Contig) 的极左端
    flank_qs_low = {
        "name": "Left_Flank" if not flip_asm else "Right_Flank", 
        "subtype": "Flank", "is_flank": True, "src_chr": main_chr,
        "qs": qs_left - flank_len, "qe": qs_left, "q_len": flank_len,
        "src_pos": lf_ref, "src_end": lf_ref + flank_len, "ref_len": flank_len,
        "strand": lf_strand, "is_inv": (lf_strand == "-"), "color": "#A9A9A9",
    }
    # flank_qs_high 对应 Assembly (Contig) 的极右端
    flank_qs_high = {
        "name": "Right_Flank" if not flip_asm else "Left_Flank",
        "subtype": "Flank", "is_flank": True, "src_chr": main_chr,
        "qs": qs_right, "qe": qs_right + flank_len, "q_len": flank_len,
        "src_pos": rf_ref, "src_end": rf_ref + flank_len, "ref_len": flank_len,
        "strand": rf_strand, "is_inv": (rf_strand == "-"), "color": "#A9A9A9",
    }

    core_blocks = []
    for idx, b in enumerate(blocks):
        core_blocks.append({
            **b,
            "name": b["subtype"],
            "color": color_palette[idx % len(color_palette)],
            "is_flank": False,
        })

    all_blocks = [flank_qs_low] + core_blocks + [flank_qs_high]

    min_asm_qs = min(b["qs"] for b in all_blocks)
    max_asm_qe = max(b["qe"] for b in all_blocks)

    # ==== 核心映射：负链进行完美逆序布局与渲染反转逻辑 ====
    for b in all_blocks:
        if not flip_asm:
            # 正链直接正序布局
            b["asm_x"] = b["qs"] - min_asm_qs
            b["visual_inv"] = b["is_inv"]
            b["visual_arrow_dir"] = "left" if b["is_inv"] else "right"
        else:
            # 负链时画布上严格从左往右画，变异顺序因 `qe` 递减而被自动逆序
            b["asm_x"] = max_asm_qe - b["qe"]
            # 视觉交叉与箭头均需与 Assembly 整体翻转状态进行 XOR 异或修正
            b["visual_inv"] = b["is_inv"] ^ True
            b["visual_arrow_dir"] = "right" if b["is_inv"] else "left"

    asm_total_len = max_asm_qe - min_asm_qs

    chroms_seen = []
    for b in all_blocks:
        if b["src_chr"] not in chroms_seen:
            chroms_seen.append(b["src_chr"])
    if main_chr in chroms_seen:
        chroms_seen.remove(main_chr)
    unique_chroms = [main_chr] + chroms_seen

    Y_ASM = 0.0
    track_spacing = 2.8
    chrom_y = {chrom: 2.5 + idx * track_spacing for idx, chrom in enumerate(unique_chroms)}

    # 取消任何形式的缩放，严格 1:1 Reference 对齐
    chrom_shifts = {}
    for chrom in unique_chroms:
        c_blks = [b for b in all_blocks if b["src_chr"] == chrom]
        if not c_blks:
            C_asm = asm_total_len / 2
            C_ref_coord = (parent["pos"] + parent["end"]) / 2
        else:
            C_asm = sum(b["asm_x"] + b["q_len"] / 2 for b in c_blks) / len(c_blks)
            C_ref_coord = sum(b["src_pos"] + b["ref_len"] / 2 for b in c_blks) / len(c_blks)
        
        chrom_shifts[chrom] = (C_ref_coord, C_asm)

        for b in all_blocks:
            if b["src_chr"] == chrom:
                b["ref_x"] = b["src_pos"] - C_ref_coord + C_asm

    max_plot_x = max(
        max(b["ref_x"] + b["ref_len"] for b in all_blocks),
        max(b["asm_x"] + b["q_len"] for b in all_blocks),
    )
    min_plot_x = min(
        min(b["ref_x"] for b in all_blocks),
        min(b["asm_x"] for b in all_blocks),
    )
    plot_width = max_plot_x - min_plot_x

    fig_height = 5 + len(unique_chroms) * 2.2
    fig, ax = plt.subplots(figsize=(16, fig_height), dpi=300)
    bh = 0.38

    ax.axhline(Y_ASM, color="black", linewidth=1.2, zorder=1, alpha=0.6)
    label_x = min_plot_x - (plot_width * 0.05 + 5000)
    
    # 动态标记当前视角的投影系方向
    asm_label = f"Assembly\n({hap_name})"
    if flip_asm:
        asm_label += "\n[RevComp mapped]"
    ax.text(label_x, Y_ASM, asm_label, va="center", ha="right", fontsize=11, fontweight="bold")

    for chrom in unique_chroms:
        y_pos = chrom_y[chrom]
        ax.axhline(y_pos, color="black", linestyle="--", linewidth=1.0, zorder=1, alpha=0.4)
        ax.text(label_x, y_pos, f"Reference\n({chrom})", va="center", ha="right", fontsize=11, fontweight="bold")

    for b in all_blocks:
        w_asm, w_ref = b["q_len"], b["ref_len"]
        x_asm, x_ref = b["asm_x"], b["ref_x"]
        color, visual_inv = b["color"], b["visual_inv"]
        strand_sym = b["strand"]
        zo = 2 if b.get("is_flank") else 3
        y_ref = chrom_y[b["src_chr"]]

        # 1. 绘制 Reference
        ax.add_patch(patches.Rectangle(
            (x_ref, y_ref - bh / 2), w_ref, bh, color=color, alpha=0.65, zorder=zo, ec="black", lw=0.8,
        ))
        mid_bp = b["src_pos"] + b["ref_len"] / 2
        ax.text(
            x_ref + w_ref / 2, y_ref + bh / 2 + 0.12, f"{mid_bp / 1e6:.3f}M ({strand_sym})",
            ha="center", va="bottom", fontsize=8.5, color="#222222",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=0.3), zorder=4,
        )

        # 2. 绘制 Assembly (已支持负链时按逆序排布)
        ax.add_patch(patches.Rectangle(
            (x_asm, Y_ASM - bh / 2), w_asm, bh, color=color, alpha=0.65, zorder=3, ec="black", lw=0.8,
        ))

        # 内部序列真实 5'->3' 箭头朝向
        if w_asm > plot_width * 0.005:
            if b["visual_arrow_dir"] == "left":
                dx = -w_asm * 0.4
                sx = x_asm + w_asm * 0.7
            else:
                dx = w_asm * 0.4
                sx = x_asm + w_asm * 0.3
            ax.arrow(sx, Y_ASM, dx, 0, head_width=0.1, head_length=min(w_asm * 0.15, 6000),
                     fc="white", ec="white", zorder=4)

        if not b.get("is_flank"):
            ax.text(
                x_asm + w_asm / 2, Y_ASM - bh / 2 - 0.22, b["name"],
                ha="center", va="top", fontsize=8.5, fontweight="bold", color=color,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=0.2),
            )

        # 3. 绘制重构后的连接带 (XOR 判定)
        draw_bezier_ribbon(
            ax, x_ref, w_ref, y_ref - bh / 2, x_asm, w_asm, Y_ASM + bh / 2, color, alpha=0.3, invert=visual_inv,
        )

    # Reference 全局覆盖跨度指示
    C_ref_coord_main, C_asm_main = chrom_shifts[main_chr]
    refspan_start_x = parent["pos"] - C_ref_coord_main + C_asm_main
    refspan_end_x = parent["end"] - C_ref_coord_main + C_asm_main
    Y_REF_MAIN = chrom_y[main_chr]

    ax.plot([refspan_start_x, refspan_end_x], [Y_REF_MAIN - 0.45, Y_REF_MAIN - 0.45],
            color="#D0021B", linewidth=2.0, linestyle="-", zorder=5)
    ax.plot([refspan_start_x, refspan_start_x], [Y_REF_MAIN - 0.55, Y_REF_MAIN - 0.35], color="#D0021B", linewidth=2.0, zorder=5)
    ax.plot([refspan_end_x, refspan_end_x], [Y_REF_MAIN - 0.55, Y_REF_MAIN - 0.35], color="#D0021B", linewidth=2.0, zorder=5)

    ref_span_len = parent["end"] - parent["pos"]
    ref_span_label = f"Ref Span: {ref_span_len / 1e3:.1f} kb" if ref_span_len < 1e6 else f"Ref Span: {ref_span_len / 1e6:.2f} Mb"
    ax.text((refspan_start_x + refspan_end_x) / 2, Y_REF_MAIN - 0.6, ref_span_label,
            color="#D0021B", ha="center", va="top", fontsize=9, fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="#D0021B", alpha=1.0, boxstyle="round,pad=0.2", lw=1), zorder=6)

    # 动态标记核心 Assembly 跨度区
    if core_blocks:
        core_start_x = min(b["asm_x"] for b in core_blocks)
        core_end_x = max(b["asm_x"] + b["q_len"] for b in core_blocks)
    else:
        core_start_x = flank_qs_high["asm_x"] + flank_qs_high["q_len"] if flip_asm else flank_qs_low["asm_x"] + flank_qs_low["q_len"]
        core_end_x = flank_qs_low["asm_x"] if flip_asm else flank_qs_high["asm_x"]

    ax.plot([core_start_x, core_end_x], [Y_ASM - bh / 8, Y_ASM - bh / 8],
            color="#4A90E2", linestyle="--", zorder=5, lw=2)
    
    q_len = parent["qs_right"] - parent["qs_left"]
    asm_label_text = (
        f"Total Query Span: {q_len / 1e6:.2f}Mb (Backbone: {bb_strand})" if q_len >= 1e6
        else f"Total Query Span: {q_len / 1e3:.1f}kb (Backbone: {bb_strand})"
    )
    ax.text(
        (core_start_x + core_end_x) / 2, Y_ASM - bh / 2 - 0.9, asm_label_text, color="#4A90E2", ha="center", va="center",
        fontweight="bold", fontsize=9,
        bbox=dict(facecolor="white", edgecolor="#4A90E2", alpha=1.0, boxstyle="round,pad=0.3", lw=1),
        zorder=6,
    )

    max_y = max(chrom_y.values())
    ax.set_xlim(label_x - plot_width * 0.02, max_plot_x + plot_width * 0.05)
    ax.set_ylim(Y_ASM - 1.8, max_y + 1.2)
    ax.axis("off")

    plt.title(
        f"Multi-Track Structural Variant Architecture: {event_id}\nHaplotype: {hap_name} | Main Backbone: {main_chr}",
        fontsize=15, fontweight="bold", pad=30,
    )
    plt.tight_layout()

    out_dir = os.path.dirname(os.path.abspath(output_path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f">>> Done. Unified-scale multi-track plot saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Plot AsmSV COMPLEX SV strictly mapped to Reference (5' -> 3') orientation.",
    )
    parser.add_argument("-i", "--vcf", required=True, help="Input VCF (.vcf or .vcf.gz)")
    parser.add_argument("-e", "--event", help="COMPLEX EVENT ID (exact match)")
    parser.add_argument("--hap", default="Hap1", choices=["Hap1", "Hap2"], help="Haplotype (default: Hap1)")
    parser.add_argument("-o", "--output", help="Output PNG path (default: ./plots/<event>_<hap>.png)")
    parser.add_argument("--list", action="store_true", help="List COMPLEX EVENT IDs and exit")
    args = parser.parse_args()

    try:
        if args.list:
            events = list_complex_events(args.vcf)
            if not events:
                print("No COMPLEX SUMMARY records with SUB_SVS found.")
            else:
                for ev in events:
                    print(ev)
                print(f">>> Done. Found {len(events)} event(s).")
            return

        if not args.event:
            parser.error("-e/--event is required unless --list is used")

        output = args.output or os.path.join("plots", f"{args.event}_{args.hap}.png")
        parent, blocks = parse_asmsv_complex_event(args.vcf, args.event, args.hap)
        plot_parsed_complex_sv(parent, blocks, args.event, args.hap, output)
    except (FileNotFoundError, ValueError, OSError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()