#!/usr/bin/env python3
import sys
"""
保留所有BND\DEL\INS\DUP\INV等单类型
对于CSV，只保留INS_DEL及INV-INS/DEL复合变异
TRA通过BND可以指示
retain_choice可选是否保留CSV;默认为否，即等价于bcftools view -i 'SVTYPE!="CSV"' input.vcf -Oz -o out.vcf.gz
"""
def parse_info(info_str):
    """解析 VCF INFO 字段为字典"""
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, val = item.split('=', 1)
            info_dict[key] = val
        else:
            info_dict[item] = True
    return info_dict
def flatten_for_benchmark(input_vcf, output_vcf,retain=False):
    with open(input_vcf, 'r') as fin, open(output_vcf, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue
            
            fields = line.strip('\n').split('\t')
            if len(fields) < 8:
                continue
                
            # 1. 严格过滤：非 PASS 记录直接丢弃（高置信度筛选）
            filter_status = fields[6]
            if filter_status != 'PASS':
                continue
                
            orig_pos = int(fields[1])
            vid = fields[2]
            info_str = fields[7]
            info = parse_info(info_str)
            svtype = info.get('SVTYPE', '')
            
            # 2. 如果是标准变异（非 CSV），直接保留输出
            if svtype != 'CSV':
                fout.write(line)
                continue
            if retain:  
                # 3. 处理 CSV 复杂变异
                subtype = info.get('SUBTYPE', '')
                
                # 严格筛选：只保留 INS_DEL 和 INV-相关子类型，其余直接过滤
                
                if subtype in ['INS_DEL', 'INV-INS', 'INV-DEL', 'INV-INS-DEL']:
                    fout.write(line)
                    continue

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python clean_simpleSV.py <input.vcf> [y] ")
        sys.exit(1)
    retain_choice=False
    if len(sys.argv) ==3:
        retain_choice=True
    flatten_for_benchmark(sys.argv[1], "clean_results.vcf",retain_choice)
