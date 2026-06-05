#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
将 VCF 中 SVTYPE=CSV; SUBTYPE=INS_DEL 的记录拆分为独立的 INS 和 DEL 记录
非 CSV 类型记录原样保留
"""

import sys
import re

def parse_info(info_str):
    """解析VCF INFO字段为字典"""
    info_dict = {}
    for entry in info_str.split(';'):
        if '=' in entry:
            key, value = entry.split('=', 1)
            info_dict[key] = value
    return info_dict

def build_info(info_dict):
    """将字典重新构建为VCF INFO字段"""
    info_parts = []
    for key, value in info_dict.items():
        if value is None or value == '':
            info_parts.append(key)
        else:
            info_parts.append(f"{key}={value}")
    return ';'.join(info_parts)

def process_vcf_line(line, id_counter):
    """处理单行VCF记录，返回处理后的行列表"""
    # 跳过空行
    if not line.strip():
        return [line.strip()]
    
    # 注释行/标题行直接返回
    if line.startswith('#'):
        return [line.strip()]
    
    fields = line.strip().split('\t')
    if len(fields) < 8:
        return [line.strip()]  # 格式错误直接保留
    
    chrom, pos, vid, ref, alt, qual, filter_, info = fields[:8]
    info_dict = parse_info(info)
    
    # 仅处理 SVTYPE=CSV 且 SUBTYPE=INS_DEL 的记录
    if (info_dict.get('SVTYPE') == 'CSV' and 
        info_dict.get('SUBTYPE') == 'INS_DEL'):
        
        # 提取必需字段
        try:
            original_end = int(info_dict['END'])
            ins_len = int(info_dict['INS_LEN'])
            del_len = int(info_dict['DEL_LEN'])
            original_pos = int(pos)
        except (KeyError, ValueError):
            return [line.strip()]  # 字段缺失，不处理
        
        # ====================== 生成 DEL 记录 ======================
        del_fields = fields.copy()
        del_info = info_dict.copy()
        
        # 更新DEL核心信息
        del_info['SVTYPE'] = 'DEL'
        del_info['SVLEN'] = -del_len  # 缺失为负数
        del_info['END'] = original_pos + del_len
        del_info.pop('SUBTYPE', None)  # 移除子类型
        del_info.pop('INS_LEN', None)  # 移除插入长度
        del_info.pop('CSV_LEN', None)  # 清理CSV相关字段
        
        del_fields[2] = f"SV.DEL.{id_counter:06d}"  # 新ID
        del_fields[4] = "<DEL>"
        del_fields[7] = build_info(del_info)
        del_line = '\t'.join(del_fields)
        
        # ====================== 生成 INS 记录 ======================
        ins_fields = fields.copy()
        ins_info = info_dict.copy()
        
        # 更新INS核心信息
        ins_info['SVTYPE'] = 'INS'
        ins_info['SVLEN'] = ins_len
        ins_info['END'] = original_pos+1  # 插入END=POS
        ins_info.pop('SUBTYPE', None)
        ins_info.pop('DEL_LEN', None)
        ins_info.pop('CSV_LEN', None)
        
        ins_fields[2] = f"SV.INS.{id_counter:06d}"
        ins_fields[1] = str(original_pos+1)  # 插入位点=原始CSV的END
        ins_fields[4] = "<INS>"
        ins_fields[7] = build_info(ins_info)
        ins_line = '\t'.join(ins_fields)
        
        return [del_line, ins_line]
    
    # 非目标记录，原样返回
    else:
        return [line.strip()]

def main():
    if len(sys.argv) != 3:
        print(f"用法: {sys.argv[0]} <输入VCF> <输出VCF>")
        sys.exit(1)
    
    input_vcf = sys.argv[1]
    output_vcf = sys.argv[2]
    variant_id_counter = 1
    
    with open(input_vcf, 'r', encoding='utf-8') as in_f, \
         open(output_vcf, 'w', encoding='utf-8') as out_f:
        
        for line in in_f:
            processed_lines = process_vcf_line(line, variant_id_counter)
            for pline in processed_lines:
                out_f.write(pline + '\n')
            
            # 每拆分一个CSV变异，ID计数器+1
            if len(processed_lines) == 2:
                variant_id_counter += 1

if __name__ == "__main__":
    main()

