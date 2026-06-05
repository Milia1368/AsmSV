# coding=utf-8
import os
from setuptools import setup, find_packages

this_directory = os.path.abspath(os.path.dirname(__file__))
readme_path = os.path.join(this_directory, 'README.md')
if os.path.exists(readme_path):
    with open(readme_path, 'r', encoding='utf-8') as f:
        long_description = f.read()
else:
    long_description = "A structural variant (SV) detection tool for genome assembly contigs."

setup(
    name="AsmSV",
    version="1.2.1",
    description="A structural variant (SV) detection tool for genome assembly contigs, supporting analysis of complex structural variants.",
    long_description=long_description,
    long_description_content_type='text/markdown',
    
    author="WuChuanmin",
    author_email="25S103213@stu.hit.edu.cn",
    url="https://github.com/Milia1368/AsmSV",
    license="MIT",
    
    # 指定源码目录
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    
    # 包含的额外文件
    data_files=[("", ["LICENSE"])],

    entry_points={
        'console_scripts': [
            'asmSV=AsmSV.asmSV:main',  # 
        ],
    },
    
    zip_safe=False,
    
    # =========================================================================
    # 核心依赖 (Dependencies)
    # =========================================================================
    install_requires=[
        "pysam>=0.19.0",   # 用于 BAM/FASTA 文件的解析
        "edlib>=1.3.9"     # 用于核心的序列相似度和编辑距离计算
    ],
    
    # 限制 Python 最小版本 (推荐 Python 3.8+)
    python_requires=">=3.8",
    
    # PyPI / GitHub 标识符，提升项目的规范性和可搜索性
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Operating System :: POSIX :: Linux",
    ],
    
    # 关键词
    keywords="bioinformatics structural-variant genome-assembly long-reads",
)
