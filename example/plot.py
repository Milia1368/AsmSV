#!/usr/bin/env python3
"""
Plot AsmSV COMPLEX region architecture from VCF *_SUMMARY records (SUB_SVS field).

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
    else:
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

        if b_type == "backbone":
            parent_event["qs_left"] = qs_i
            parent_event["qs_right"] = qe_i
            parent_event["strand"] = strand
            continue

        bb_strand = parent_event.get("strand", "++")
        blocks.append({
            "subtype": b_type,
            "src_chr": c_chr,
            "src_pos": rs_i,
            "src_end": re_i,
            "qs": qs_i,
            "qe": qe_i,
            "strand": strand,
            "is_inv": strand != bb_strand,
            "q_len": qe_i - qs_i,
            "ref_len": re_i - rs_i,
        })
    return blocks


def list_complex_events(vcf_path):
    """列出 VCF 中所有 COMPLEX *_SUMMARY 事件的 EVENT ID。"""
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
    """从 VCF 解析指定 COMPLEX 事件的 SUMMARY / SUB_SVS 记录。"""
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
            "up_down": str(safe_get(rec.info, "BACKBONE_FLIP", "False")).lower() == "true",
            "strand": str(safe_get(rec.info, "BACKBONE_STRAND", "++")),
            "query_len": int(safe_get(rec.info, "INS_LEN", 0)),
        }
        blocks = _parse_sub_svs(rec.info.get("SUB_SVS", ""), parent_event)
        break

    vcf.close()

    if not parent_event:
        raise ValueError(f"COMPLEX event '{target_event}' SUMMARY not found for {hap_name}.")
    if not blocks:
        raise ValueError(f"Event '{target_event}' has no valid SUB_SVS blocks.")
    return parent_event, blocks


def plot_parsed_complex_sv(parent, blocks, event_id, hap_name, output_path):
    patches, plt, _ = _require_matplotlib()
    flank_len = 50000
    fig, ax = plt.subplots(figsize=(16, 8), dpi=300)
    main_chr = parent["chrom"]
    color_palette = ["#4A90E2", "#50E3C2", "#F5A623", "#D0021B", "#9013FE", "#8B572A"]

    up_down = False
    qs_left = parent.get("qs_left", 0)
    qs_right = parent.get("qs_right", parent.get("query_len", 1000))
    bb_strands = parent.get("strand", "++")
    if bb_strands != "++":
        up_down = True
    if len(bb_strands) < 2:
        bb_strands += bb_strands

    if not up_down:
        lf_ref = parent["pos"] - flank_len
        rf_ref = parent["end"]
        lf_inv = bb_strands[0] == "-"
        rf_inv = bb_strands[1] == "-"
    else:
        lf_ref = parent["end"]
        rf_ref = parent["pos"] - flank_len
        lf_inv = bb_strands[0] == "-"
        rf_inv = bb_strands[1] == "-"

    left_flank = {
        "name": "Left_Flank", "subtype": "Flank", "is_flank": True,
        "qs": qs_left - flank_len, "qe": qs_left, "q_len": flank_len,
        "src_pos": lf_ref, "src_end": lf_ref + flank_len, "ref_len": flank_len,
        "is_inv": lf_inv, "color": "#A9A9A9",
    }
    right_flank = {
        "name": "Right_Flank", "subtype": "Flank", "is_flank": True,
        "qs": qs_right, "qe": qs_right + flank_len, "q_len": flank_len,
        "src_pos": rf_ref, "src_end": rf_ref + flank_len, "ref_len": flank_len,
        "is_inv": rf_inv, "color": "#A9A9A9",
    }

    core_blocks = []
    for idx, b in enumerate(blocks):
        core_blocks.append({
            **b,
            "name": b["subtype"],
            "color": color_palette[idx % len(color_palette)],
            "is_flank": False,
        })

    all_blocks = [left_flank] + core_blocks + [right_flank]

    min_asm_x = min(b["qs"] for b in all_blocks)
    min_ref_x = min(b["src_pos"] for b in all_blocks)
    for b in all_blocks:
        b["asm_x"] = b["qs"] - min_asm_x
        b["ref_x"] = b["src_pos"] - min_ref_x

    max_plot_x = max(
        max(b["ref_x"] + b["ref_len"] for b in all_blocks),
        max(b["asm_x"] + b["q_len"] for b in all_blocks),
    )

    Y_REF = 2.0
    Y_ASM = 0.0
    bh = 0.35

    ax.axhline(Y_REF, color="black", linewidth=1.2, zorder=1, alpha=0.6)
    ax.axhline(Y_ASM, color="black", linewidth=1.2, zorder=1, alpha=0.6)

    label_x = -(max_plot_x * 0.05 + 5000)
    ax.text(label_x, Y_REF, f"Main Ref\n{main_chr}", va="center", ha="right", fontsize=11, fontweight="bold")
    ax.text(label_x, Y_ASM, f"Assembly\n({hap_name})", va="center", ha="right", fontsize=11, fontweight="bold")

    for b in all_blocks:
        w_asm, w_ref = b["q_len"], b["ref_len"]
        x_asm, x_ref = b["asm_x"], b["ref_x"]
        color, is_inv = b["color"], b["is_inv"]
        zo = 2 if b.get("is_flank") else 3

        ax.add_patch(patches.Rectangle(
            (x_ref, Y_REF - bh / 2), w_ref, bh, color=color, alpha=0.65, zorder=zo, ec="black", lw=0.8,
        ))
        mid_bp = b["src_pos"] + b["ref_len"] / 2
        ax.text(
            x_ref + w_ref / 2, Y_REF + bh / 2 + 0.1, f"{mid_bp / 1e6:.2f}M",
            ha="center", va="bottom", fontsize=8.5, color="#222222",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=0.5), zorder=4,
        )

        ax.add_patch(patches.Rectangle(
            (x_asm, Y_ASM - bh / 2), w_asm, bh, color=color, alpha=0.65, zorder=3, ec="black", lw=0.8,
        ))

        if w_asm > max_plot_x * 0.005:
            dx = -w_asm * 0.6 if is_inv else w_asm * 0.6
            sx = x_asm + w_asm * 0.8 if is_inv else x_asm + w_asm * 0.2
            ax.arrow(sx, Y_ASM, dx, 0, head_width=0.08, head_length=min(w_asm * 0.1, 5000),
                     fc="white", ec="white", zorder=4)

        if not b.get("is_flank"):
            ax.text(
                x_asm + w_asm / 2, Y_ASM - bh / 2 - 0.2, b["name"],
                ha="center", va="top", fontsize=9, fontweight="bold", color=color,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=0.2),
            )

        draw_bezier_ribbon(
            ax, x_ref, w_ref, Y_REF - bh / 2, x_asm, w_asm, Y_ASM + bh / 2, color, alpha=0.35, invert=is_inv,
        )

    ref_gap_len = parent["end"] - parent["pos"]
    q_len = parent["qs_right"] - parent["qs_left"]

    if ref_gap_len > 0:
        ref_start_x = parent["pos"] - min_ref_x
        ref_end_x = parent["end"] - min_ref_x
        gap_center_x = (ref_start_x + ref_end_x) / 2

        ax.plot([ref_start_x, ref_end_x], [Y_REF, Y_REF], color="#D0021B", linestyle="--", zorder=5, lw=2)
        ref_label = (
            f"Ref Gap: {ref_gap_len / 1e6:.2f}Mb" if ref_gap_len >= 1e6
            else f"Ref Gap: {ref_gap_len / 1e3:.1f}kb"
        )
        ax.text(
            gap_center_x, Y_REF, ref_label, color="#D0021B", ha="center", va="center",
            fontweight="bold", fontsize=8.5,
            bbox=dict(facecolor="white", edgecolor="#D0021B", alpha=1.0, boxstyle="round,pad=0.3", lw=1),
            zorder=6,
        )

        asm_start_x = left_flank["asm_x"] + left_flank["q_len"]
        asm_end_x = right_flank["asm_x"]
        asm_center_x = (asm_start_x + asm_end_x) / 2
        ax.plot([asm_start_x, asm_end_x], [Y_ASM - bh / 8, Y_ASM - bh / 8],
                color="#4A90E2", linestyle="--", zorder=5, lw=2)
        asm_label = (
            f"Query Span: {q_len / 1e6:.2f}Mb" if q_len >= 1e6
            else f"Query Span: {q_len / 1e3:.1f}kb"
        )
        ax.text(
            asm_center_x, Y_ASM - bh / 2, asm_label, color="#4A90E2", ha="center", va="center",
            fontweight="bold", fontsize=8.5,
            bbox=dict(facecolor="white", edgecolor="#4A90E2", alpha=1.0, boxstyle="round,pad=0.3", lw=1),
            zorder=6,
        )

    ax.set_xlim(label_x - max_plot_x * 0.02, max_plot_x * 1.05)
    ax.set_ylim(Y_ASM - 1.2, Y_REF + 1.5)
    ax.axis("off")
    ax.text(
        0.5, -0.05, f"*Default Backbone Flanking: {flank_len / 1e3:.0f} kb",
        transform=ax.transAxes, fontsize=9, color="#777777", style="italic", ha="center",
    )

    plt.title(
        f"Structural Variant Architecture: {event_id}\nHaplotype: {hap_name}",
        fontsize=15, fontweight="bold", pad=30,
    )
    plt.tight_layout()

    out_dir = os.path.dirname(os.path.abspath(output_path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f">>> Done. Plot saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Plot AsmSV COMPLEX SV from VCF *_SUMMARY / SUB_SVS records.",
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
