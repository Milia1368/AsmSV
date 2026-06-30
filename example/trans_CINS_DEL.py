#!/usr/bin/env python3
"""
将 CSV;SUBTYPE=INS_DEL 记录拆分为独立的 DEL 和 INS，用于 benchmark。

其他记录原样保留。

Usage:
  python trans_CINS_DEL.py -i results.vcf
  python trans_CINS_DEL.py -i results.vcf -o split.vcf
"""
import argparse
import os
import sys


def default_output(input_path, tag="ins_del_split"):
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
    for entry in info_str.split(';'):
        if '=' in entry:
            key, value = entry.split('=', 1)
            info_dict[key] = value
    return info_dict


def build_info(info_dict):
    info_parts = []
    for key, value in info_dict.items():
        if value is None or value == '':
            info_parts.append(key)
        else:
            info_parts.append(f"{key}={value}")
    return ';'.join(info_parts)


def process_vcf_line(line, id_counter):
    if not line.strip():
        return [line.strip()]
    if line.startswith('#'):
        return [line.strip()]

    fields = line.strip().split('\t')
    if len(fields) < 8:
        return [line.strip()]

    pos, info = fields[1], fields[7]
    info_dict = parse_info(info)

    if info_dict.get('SVTYPE') != 'CSV' or info_dict.get('SUBTYPE') != 'INS_DEL':
        return [line.strip()]

    try:
        ins_len = int(info_dict['INS_LEN'])
        del_len = int(info_dict['DEL_LEN'])
        original_pos = int(pos)
    except (KeyError, ValueError):
        return [line.strip()]

    del_fields = fields.copy()
    del_info = info_dict.copy()
    del_info['SVTYPE'] = 'DEL'
    del_info['SVLEN'] = -del_len
    del_info['END'] = original_pos + del_len
    del_info.pop('SUBTYPE', None)
    del_info.pop('INS_LEN', None)
    del_info.pop('CSV_LEN', None)
    del_fields[2] = f"SV.DEL.{id_counter:06d}"
    del_fields[4] = "<DEL>"
    del_fields[7] = build_info(del_info)
    del_line = '\t'.join(del_fields)

    ins_fields = fields.copy()
    ins_info = info_dict.copy()
    ins_info['SVTYPE'] = 'INS'
    ins_info['SVLEN'] = ins_len
    ins_info['END'] = original_pos + 1
    ins_info.pop('SUBTYPE', None)
    ins_info.pop('DEL_LEN', None)
    ins_info.pop('CSV_LEN', None)
    ins_fields[2] = f"SV.INS.{id_counter:06d}"
    ins_fields[1] = str(original_pos + 1)
    ins_fields[4] = "<INS>"
    ins_fields[7] = build_info(ins_info)
    ins_line = '\t'.join(ins_fields)

    return [del_line, ins_line]


def main():
    parser = argparse.ArgumentParser(description="Split INS_DEL CSV records into DEL + INS.")
    parser.add_argument("-i", "--input", required=True, help="Input VCF")
    parser.add_argument("-o", "--output", help="Output VCF (default: <input>.ins_del_split.vcf)")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"Error: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    output = args.output or default_output(args.input)
    out_dir = os.path.dirname(os.path.abspath(output))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    variant_id_counter = 1
    with open(args.input, 'r', encoding='utf-8') as in_f, \
         open(output, 'w', encoding='utf-8') as out_f:
        for line in in_f:
            processed_lines = process_vcf_line(line, variant_id_counter)
            for pline in processed_lines:
                out_f.write(pline + '\n')
            if len(processed_lines) == 2:
                variant_id_counter += 1

    print(f">>> Done. Split VCF saved: {output}")


if __name__ == "__main__":
    main()
