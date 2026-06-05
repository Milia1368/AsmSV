import sys
import re
"""只保留INV\INV-DEL\INV-INS\INV-INS-DEL并将复合变异拆分为单独记录"""
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

def dict_to_info(info_dict):
    """将字典转换回 VCF INFO 字符串"""
    items = []
    for k, v in info_dict.items():
        if v is True:
            items.append(k)
        else:
            items.append(f"{k}={v}")
    return ';'.join(items)

def parse_bp_details(bp_str):
    """从 BP_DETAILS=LC:0,RC:19,LF:0,RF:4220 中提取数值"""
    details = {}
    for item in bp_str.split(','):
        if ':' in item:
            k, v = item.split(':')
            details[k] = int(v)
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
                
            chrom = fields[0]
            pos = int(fields[1])
            vid = fields[2]
            info_str = fields[7]
            
            info = parse_info(info_str)
            svtype = info.get('SVTYPE', '')
            
            # 1. 直接保留原生 INV
            if svtype == 'INV':
                fout.write(line)
                continue
                
            # 2. 处理 CSV 类型中的倒位复合体，其余类型过滤
            if svtype == 'CSV':
                subtype = info.get('SUBTYPE', '')
                if subtype not in ['INV-INS', 'INV-DEL', 'INV-INS-DEL']:
                    #可选其他输出否
                    fout.write(line)
                    continue
                
                # 提取 INV_SPAN 和 BP_DETAILS
                inv_span = int(info.get('INV_SPAN', 0))
                bp_details = parse_bp_details(info.get('BP_DETAILS', 'LC:0,RC:0,LF:0,RF:0'))
                
                # LC/RC 为 Query gap (插入)，LF/RF 为 Ref gap (缺失)
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
                
                # ==========================================
                # 输出核心 INV 记录
                # ==========================================
                inv_info = info.copy()
                inv_info['SVTYPE'] = 'INV'
                inv_info['SVLEN'] = inv_span
                inv_info['END'] = inv_end
                if 'SUBTYPE' in inv_info:
                    del inv_info['SUBTYPE']
                
                inv_fields = fields.copy()
                inv_fields[2] = f"{vid}_CORE_INV"
                inv_fields[4] = "<INV>"
                inv_fields[7] = dict_to_info(inv_info)
                fout.write('\t'.join(inv_fields) + '\n')
                
                # ==========================================
                # 输出缺失 (DEL)
                # ==========================================
                if true_del_left > minsvlen:
                    del_pos = pos - true_del_left  # 修正：左侧缺失应当减去 true_del_left
                    del_end = pos
                    del_info = info.copy()
                    del_info.update({'SVTYPE': 'DEL', 'END': del_end, 'SVLEN': -true_del_left})
                    
                    del_fields = fields.copy()
                    del_fields[1] = str(del_pos)
                    del_fields[2] = f"{vid}_LEFT_DEL"
                    del_fields[4] = "<DEL>"
                    del_fields[7] = dict_to_info(del_info)
                    fout.write('\t'.join(del_fields) + '\n')

                if true_del_right > minsvlen:
                    del_pos = inv_end          # 右侧缺失紧接在 INV 之后
                    del_end = inv_end + true_del_right
                    del_info = info.copy()
                    del_info.update({'SVTYPE': 'DEL', 'END': del_end, 'SVLEN': -true_del_right})
                    
                    del_fields = fields.copy()
                    del_fields[1] = str(del_pos)
                    del_fields[2] = f"{vid}_RIGHT_DEL"
                    del_fields[4] = "<DEL>"
                    del_fields[7] = dict_to_info(del_info)
                    fout.write('\t'.join(del_fields) + '\n')

                # ==========================================
                # 输出插入 (INS)
                # ==========================================
                if true_ins_left > minsvlen:
                    ins_pos = pos
                    ins_info = info.copy()
                    ins_info.update({'SVTYPE': 'INS', 'END': ins_pos, 'SVLEN': true_ins_left})
                    
                    ins_fields = fields.copy()
                    ins_fields[1] = str(ins_pos)
                    ins_fields[2] = f"{vid}_LEFT_INS"
                    ins_fields[4] = "<INS>"
                    ins_fields[7] = dict_to_info(ins_info)
                    fout.write('\t'.join(ins_fields) + '\n')

                if true_ins_right > minsvlen:
                    ins_pos = inv_end
                    ins_info = info.copy()
                    ins_info.update({'SVTYPE': 'INS', 'END': ins_pos, 'SVLEN': true_ins_right})
                    
                    ins_fields = fields.copy()
                    ins_fields[1] = str(ins_pos)
                    ins_fields[2] = f"{vid}_RIGHT_INS"
                    ins_fields[4] = "<INS>"
                    ins_fields[7] = dict_to_info(ins_info)
                    fout.write('\t'.join(ins_fields) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("用法: python convert_asmsv_csv.py <input.vcf> <output.vcf>")
        sys.exit(1)
        
    in_vcf = sys.argv[1]
    out_vcf = sys.argv[2]
    convert_vcf(in_vcf, out_vcf, minsvlen=30)
