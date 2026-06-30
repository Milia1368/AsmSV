#!/usr/bin/env python3
"""
将 INV 相关 CSV 复合变异拆分为独立 INV / DEL / INS 记录，用于 benchmark。

- 原生 INV：原样保留
- CSV + SUBTYPE 为 INV-INS / INV-DEL / INV-INS-DEL：拆分为 CORE_INV 及两侧 DEL/INS
- 其他 CSV：原样保留

Usage:
  python trans_inv.py -i results.vcf
  python trans_inv.py -i results.vcf -o out.vcf --min-sv 50
"""
import argparse
import os
import sys


def default_output(input_path, tag="inv_split"):
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


def dict_to_info(info_dict):
    items = []
    for k, v in info_dict.items():
        if v is True:
            items.append(k)
        else:
            items.append(f"{k}={v}")
    return ';'.join(items)


def parse_bp_details(bp_str):
    details = {}
    for item in bp_str.split(','):
        if ':' in item:
            k, v = item.split(':')
            try:
                details[k] = int(v)
            except ValueError:
                pass
    return details


def convert_vcf(input_vcf, output_vcf, minsvlen=30):
    with open(input_vcf, 'r') as fin, open(output_vcf, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue

            fields = line.strip('\n').split('\t')
            if len(fields) < 8:
                continue

            pos = int(fields[1])
            vid = fields[2]
            info = parse_info(fields[7])
            svtype = info.get('SVTYPE', '')

            if svtype == 'INV':
                fout.write(line)
                continue

            if svtype == 'CSV':
                subtype = info.get('SUBTYPE', '')
                if subtype not in ['INV-INS', 'INV-DEL', 'INV-INS-DEL']:
                    fout.write(line)
                    continue

                try:
                    inv_span = int(info.get('INV_SPAN', 0))
                except ValueError:
                    fout.write(line)
                    continue

                bp_details = parse_bp_details(info.get('BP_DETAILS', 'LC:0,RC:0,LF:0,RF:0'))
                q_raw_left = bp_details.get('LC', 0)
                q_raw_right = bp_details.get('RC', 0)
                raw_gap_left = bp_details.get('LF', 0)
                raw_gap_right = bp_details.get('RF', 0)

                delta_left = q_raw_left - raw_gap_left
                delta_right = q_raw_right - raw_gap_right

                true_ins_left = max(0, delta_left)
                true_del_left = max(0, -delta_left)
                true_ins_right = max(0, delta_right)
                true_del_right = max(0, -delta_right)

                inv_end = pos + inv_span

                inv_info = info.copy()
                inv_info['SVTYPE'] = 'INV'
                inv_info['SVLEN'] = inv_span
                inv_info['END'] = inv_end
                inv_info.pop('SUBTYPE', None)

                inv_fields = fields.copy()
                inv_fields[2] = f"{vid}_CORE_INV"
                inv_fields[4] = "<INV>"
                inv_fields[7] = dict_to_info(inv_info)
                fout.write('\t'.join(inv_fields) + '\n')

                if true_del_left > minsvlen:
                    del_info = info.copy()
                    del_info.update({'SVTYPE': 'DEL', 'END': pos, 'SVLEN': -true_del_left})
                    del_fields = fields.copy()
                    del_fields[1] = str(pos - true_del_left)
                    del_fields[2] = f"{vid}_LEFT_DEL"
                    del_fields[4] = "<DEL>"
                    del_fields[7] = dict_to_info(del_info)
                    fout.write('\t'.join(del_fields) + '\n')

                if true_del_right > minsvlen:
                    del_info = info.copy()
                    del_info.update({'SVTYPE': 'DEL', 'END': inv_end + true_del_right, 'SVLEN': -true_del_right})
                    del_fields = fields.copy()
                    del_fields[1] = str(inv_end)
                    del_fields[2] = f"{vid}_RIGHT_DEL"
                    del_fields[4] = "<DEL>"
                    del_fields[7] = dict_to_info(del_info)
                    fout.write('\t'.join(del_fields) + '\n')

                if true_ins_left > minsvlen:
                    ins_info = info.copy()
                    ins_info.update({'SVTYPE': 'INS', 'END': pos, 'SVLEN': true_ins_left})
                    ins_fields = fields.copy()
                    ins_fields[1] = str(pos)
                    ins_fields[2] = f"{vid}_LEFT_INS"
                    ins_fields[4] = "<INS>"
                    ins_fields[7] = dict_to_info(ins_info)
                    fout.write('\t'.join(ins_fields) + '\n')

                if true_ins_right > minsvlen:
                    ins_info = info.copy()
                    ins_info.update({'SVTYPE': 'INS', 'END': inv_end, 'SVLEN': true_ins_right})
                    ins_fields = fields.copy()
                    ins_fields[1] = str(inv_end)
                    ins_fields[2] = f"{vid}_RIGHT_INS"
                    ins_fields[4] = "<INS>"
                    ins_fields[7] = dict_to_info(ins_info)
                    fout.write('\t'.join(ins_fields) + '\n')
                continue

            fout.write(line)


def main():
    parser = argparse.ArgumentParser(description="Split INV-related CSV records into simple SV types.")
    parser.add_argument("-i", "--input", required=True, help="Input VCF")
    parser.add_argument("-o", "--output", help="Output VCF (default: <input>.inv_split.vcf)")
    parser.add_argument("--min-sv", type=int, default=30, metavar="BP",
                        help="Minimum length for split DEL/INS (default: 30)")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"Error: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    if args.min_sv < 0:
        print("Error: --min-sv must be >= 0", file=sys.stderr)
        sys.exit(1)

    output = args.output or default_output(args.input)
    out_dir = os.path.dirname(os.path.abspath(output))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    convert_vcf(args.input, output, minsvlen=args.min_sv)
    print(f">>> Done. Converted VCF saved: {output}")


if __name__ == "__main__":
    main()
