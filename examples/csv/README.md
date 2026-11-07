# CSV example data / CSV 示例数据

[English walkthrough](../../docs/CSV_EXAMPLES_EN.md) · [中文逐例说明](../../docs/CSV_EXAMPLES_ZH.md)

This package contains five selected real AsmSV 1.3.0 report units and one separate incomplete region from frozen HG002-to-CHM13v2.0 output. It is an output-reading example, not a truth set or a detection benchmark. The caller was not rerun to create the subset.

本包保留五个真实报告单元及一个独立未解析区域，用于结果阅读。它不是生物学真值或检测性能基准；提取时没有重新运行检测器。

| File | Contents / 内容 |
| --- | --- |
| `data/calls.vcf` | 5 main report rows / 五条主报告 |
| `data/calls.vcf.components.vcf` | 6 linked components / 六条关联组件 |
| `data/calls.vcf.indels.vcf` | 5 reference INS/DEL projections / 五条参考投影 |
| `data/calls.vcf.junctions.vcf` | 20 BND rows, 10 reciprocal mate pairs / 二十条BND、十对mate |
| `data/calls.vcf.events.jsonl` | Selected record indices, events, contexts, atoms and referenced blocks / 选择后的结构证据 |
| `data/calls.vcf.regions.jsonl` | Metadata and one PARTIAL region / 元数据及一个PARTIAL区域 |
| `data/calls.vcf.regions.bed` | Reference projections of that region / 区域的参考投影 |
| `data/provenance.json` | Original source/row hashes, exact IDs, selection and subset hashes / 来源与原始行哈希、ID、摘录规则 |
| `inspect_csv.py` | Read-only event/region inspection / 只读查询脚本 |
| `figures/` | REF/ALT cartoons for E1–E3, not to scale / 非按比例结构示意 |

VCF data lines remain unchanged; the headers add a subset note. JSONL is an excerpt: original block indices are retained, unrelated objects and the full-contig genotyping index are omitted. It cannot be used to rerun detection or infer assembly sequence from missing objects. The frozen source VCF header date is 20261102; the source bundle was retained in the 20261103 region-validation collection. The current maintenance version remains 1.3.0; provenance retains the source scoring preset and implementation hash.

VCF 数据行未改写，文件头增加摘录说明。JSONL只保留关联对象，原始block索引保留，未附完整contig分型遍历索引；不能将未收录对象视为不存在或用于重跑检测。原始VCF头日期为20261102，来源为20261103区域验证保存包。软件仍为1.3.0，元数据保留来源配置及实现哈希。

Run from the repository root / 从仓库根目录运行：

~~~bash
asmSV-cpx --vcf examples/csv/data/calls.vcf \
  --scope csv --out-prefix /tmp/asmsv-csv-demo/calls.csv

python examples/csv/inspect_csv.py --vcf examples/csv/data/calls.vcf \
  --event EVENT.47c90f8bd1c8324ff795918f

python examples/csv/inspect_csv.py --vcf examples/csv/data/calls.vcf \
  --region REGION.440fcc21c1282d2a732619c3
~~~

Expected export: five broad CSV reports, versus three with `--scope cpx`. The region is not a CSV report. Generated HTML/JSONL/VERIX files are written outside the source fixtures; use a new prefix or `--force` to repeat export.

预期广义CSV导出5个，CPX-only为3个；区域不计作CSV。派生文件写到源摘录之外，重复导出时请换前缀或使用`--force`。
