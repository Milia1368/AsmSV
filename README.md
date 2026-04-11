
````markdown
# AsmSV

**AsmSV** is an advanced Structural Variant (SV) detection tool designed specifically for genome assembly contigs. By analyzing long-read alignments (e.g., from PacBio HiFi or Oxford Nanopore), AsmSV accurately identifies both simple and complex structural variants and genotypes them across haplotypes.

## ✨ Key Features

* **Comprehensive SV Detection**: Supports standard SV types including Deletions (`DEL`), Insertions (`INS`), Duplications (`DUP`), Inversions (`INV`), Translocations (`TRA`), Breakends (`BND`), and `COMPLEX` structural variants.
* **Haplotype-Aware Genotyping**: Natively accepts phased assembly BAMs (Hap1 & Hap2) to output phased genotypes (e.g., `1|0`, `0|1`, `1|1`).
* **Advanced Refinement**: Utilizes dynamic programming chaining and segment splitting to resolve complex assembly breakpoints.
* **Standardized Output**: Generates strictly formatted VCFv4.2 files ready for downstream analysis.
* **High Performance**: Optimized with memory buffers and fast edit-distance sequence similarity calculations using `edlib`.

## 🛠️ Installation

### Prerequisites
* **Python**: >= 3.8
* **Dependencies**: `pysam`, `edlib`

### Install from Source

Clone the repository and install it using `pip`. This will automatically install the required dependencies and set up the `asmSV` command-line tool.

```bash
git clone [https://github.com/Milia1368/AsmSV.git](https://github.com/Milia1368/AsmSV.git)
cd AsmSV
pip install -e .
````

*Note: The `-e` flag installs the package in editable mode, which is highly recommended during active development.*

## 🚀 Quick Start

AsmSV requires two sorted and indexed BAM files (representing the two assembled haplotypes aligned to a reference genome) and the reference genome FASTA file.

```bash
asmSV \
  --bam1 hap1_aligned.bam \
  --bam2 hap2_aligned.bam \
  --ref reference.fasta \
  --output ./results/sv_calls.vcf \
  --threads 8
```

## ⚙️ Command-Line Arguments

### Input / Output

| Argument | Description | Required |
| :--- | :--- | :---: |
| `--bam1` | Input BAM for Hap1 (must be sorted & indexed). | **Yes** |
| `--bam2` | Input BAM for Hap2 (must be sorted & indexed). | **Yes** |
| `--ref` | Reference Genome FASTA file (indexed). | **Yes** |
| `-o`, `--output` | Path to the output VCF file. Default: `./results.vcf`. | No |

### Filtering Options

| Argument | Description | Default |
| :--- | :--- | :--- |
| `--min-sv` | Minimum SV size to report (bp). | `30` |
| `--max-sv` | Maximum SV size. Set to `-1` to disable upper bound filtering. | `1500000` |
| `--min-mapq` | Minimum mapping quality (MAPQ) score for alignment blocks. | `10` |
| `--min-tra-size`| Threshold size to classify a large gap as a Translocation. | `150000` |

### Clustering & Merging Options

| Argument | Description | Default |
| :--- | :--- | :--- |
| `--cluster-pos-tol` | Breakpoint position distance tolerance for clustering SVs. | `3` |
| `--cluster-len-tol` | Length tolerance adjustment for clustering DEL/INS. | `0` |
| `--min-seq-sim` | Minimum sequence similarity required for merging INS events. | `1.0` |

### VCF Formatting

  * `--max-alleles-len <INT>`: Maximum sequence length to write in the ALT column. Longer SVs use symbolic alleles (e.g., `<DEL>`). (Default: `5000`)
  * `--no-symbolic`: Disables symbolic alleles entirely and forces full sequence output. *(Warning: May create massive VCF files).*

### Advanced Options

  * `-t`, `--threads`: Number of threads for parallel processing. (Default: `8`)
  * `--enable-density-filter`: Enable high-density noise filtering to remove false positives in dense alignment regions.
  * `--split-segments`: Enable aggressive segment splitting and refinement before DP chaining.

## 📄 License

This project is licensed under the [MIT License](https://www.google.com/search?q=LICENSE).

## ✉️ Contact

**Author**: WuChuanmin  
**Email**: 25S103213@stu.hit.edu.cn  
**Issues**: Please report bugs or request features on the [GitHub Issues](https://www.google.com/search?q=https://github.com/Milia1368/AsmSV/issues) page.

```
```