#!/usr/bin/env python3
"""
过滤 AsmSV VCF，用于 benchmark 前的变异类型筛选。

保留 PASS 的 BND/DEL/INS/DUP/INV 等单类型记录。
默认丢弃 CSV；加 --keep-csv 时额外保留 SUBTYPE 为 INS_DEL / INV-INS / INV-DEL / INV-INS-DEL 的 CSV。

Usage:
  python filter_sv.py -i results.vcf
  python filter_sv.py -i results.vcf -o clean.vcf --keep-csv
"""
import argparse
import os
import sys


def default_output(input_path, tag="filtered"):
    dirname = os.path.dirname(os.path.abspath(input_path))
    name = os.path.basename(input_path)
    if name.endswith(".vcf.gz"):
        stem = name[:-7]
    elif name.endswith(".vcf"):
        stem = name[:-4]
    else:
        stem = os.path.splitext(name)[0]
    return os.path.join(dirname, f"{stem}.{tag}.vcf")


def parse_info(info_str):
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, val = item.split('=', 1)
            info_dict[key] = val
        else:
            info_dict[item] = True
    return info_dict


def flatten_for_benchmark(input_vcf, output_vcf, retain=False):
    with open(input_vcf, 'r') as fin, open(output_vcf, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue

            fields = line.strip('\n').split('\t')
            if len(fields) < 8:
                continue

            if fields[6] != 'PASS':
                continue

            info = parse_info(fields[7])
            svtype = info.get('SVTYPE', '')

            if svtype != 'CSV':
                fout.write(line)
                continue

            if retain:
                subtype = info.get('SUBTYPE', '')
                if subtype in ['INS_DEL', 'INV-INS', 'INV-DEL', 'INV-INS-DEL']:
                    fout.write(line)


def main():
    parser = argparse.ArgumentParser(description="Filter AsmSV VCF for benchmarking.")
    parser.add_argument("-i", "--input", required=True, help="Input VCF")
    parser.add_argument("-o", "--output", help="Output VCF (default: <input>.filtered.vcf)")
    parser.add_argument("--keep-csv", action="store_true",
                        help="Also keep selected CSV subtypes (INS_DEL, INV-*)")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"Error: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    output = args.output or default_output(args.input)
    out_dir = os.path.dirname(os.path.abspath(output))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    flatten_for_benchmark(args.input, output, retain=args.keep_csv)
    print(f">>> Done. Filtered VCF saved: {output}")


if __name__ == "__main__":
    main()
