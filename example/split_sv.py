#!/usr/bin/env python3
import sys

def flatten_for_benchmark(input_vcf, output_vcf):
    with open(input_vcf, 'r') as fin, open(output_vcf, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue
                
            parts = line.strip().split('\t')
            info_col = parts[7]
            
            # 简单解析 INFO 字典提取关键类型
            info_dict = {kv.split('=')[0]: kv.split('=')[1] if '=' in kv else True 
                         for kv in info_col.split(';')}
            
            svtype = info_dict.get("SVTYPE", "")
            subtype = info_dict.get("SUBTYPE", "")
            
            # =======================================================
            # 逻辑 1: 保留所有断点连线 (BND)
            # 它们天然负责向评测软件解释所有的 TRA (易位) 和复杂跳跃
            # =======================================================
            if svtype == "BND":
                fout.write(line)
                continue
                
            # =======================================================
            # 逻辑 2: 过滤掉宏观包裹事件 (COMPLEX)
            # 评测软件无法理解宏观标签，直接删掉，防止引发假阳性
            # =======================================================
            if svtype == "COMPLEX":
                continue
                
            # =======================================================
            # 逻辑 3: 处理复合倒位实体 (如 INV_TRA_DUP, INV_DUP:ISP)
            # 既然 BND 解释了跳跃，实体只需转换为纯粹的 <INV>，并删掉原 DUP 标签
            # =======================================================
            if "Inv" in subtype and parts[4] != "<INV>":
                new_parts = parts.copy()
                new_parts[4] = "<INV>"  # 核心替换：直接变身为 INV
                
                # 在 INFO 里用 SVTYPE=INV 替换掉原有的 SVTYPE=DUP
                new_info = info_col.replace(f"SVTYPE={svtype}", "SVTYPE=INV")
                new_parts[7] = new_info
                
                # 输出改头换面的 INV 记录，【原有的 DUP 记录被彻底抛弃】
                fout.write('\t'.join(new_parts) + '\n')
                continue
                
            # =======================================================
            # 逻辑 4: 处理纯粹的跨染色体转移实体 (如 TRA_DUP, 且无 INV)
            # 既然是纯易位，且无倒位，BND 已经完美交代了整个事件
            # 输出这个实体可能会被评测软件误判为原位置发生了散布重复 (假阳性)
            # 所以我们果断删掉/跳过这类纯易位实体！
            # =======================================================
            if "TRA" in subtype and "Inv" not in subtype:
                continue 
                
            # =======================================================
            # 逻辑 5: 保留其他原汁原味的经典变异
            # 比如单纯的 DEL, INS, 以及本地的 DUP:TANDEM
            # =======================================================
            fout.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python benchmark_flattener.py <input.vcf> ")
        sys.exit(1)
    flatten_for_benchmark(sys.argv[1], "results_split.vcf")
