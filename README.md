# 磷循环功能基因测序数据分析

分析安徽（AH）土壤样本中三个磷循环功能基因（phoC、phoD、pqqC）的测序数据。

## 样本信息

- **地点**：安徽
- **处理类型**（5种，每处理3重复）：
  - CK：对照组
  - F：化肥
  - FM：化肥+有机肥
  - FS：化肥+秸秆
  - FMS：化肥+秸秆+有机肥

## 功能基因

| 基因 | 功能 |
|------|------|
| phoC | 碱性磷酸酶C基因，有机磷矿化 |
| phoD | 酸性/碱性磷酸酶D基因，有机磷矿化 |
| pqqC | 吡咯喹啉醌合成酶C基因，磷酸溶解 |

## 脚本说明

| 脚本 | 功能 |
|-------|------|
| `01_data_preprocessing.py` | 数据预处理和质控评估 |
| `02_gene_abundance_analysis.py` | 基因丰度比较分析 |
| `03_alpha_diversity_analysis.py` | α多样性分析 |
| `04_beta_diversity_analysis.py` | β多样性分析 |
| `08_summary_analysis.py` | 综合总结分析（简化版） |

## 安装依赖

```bash
pip install -r requirements.txt
```

## 使用方法

按顺序运行脚本：

```bash
python 01_data_preprocessing.py    # 数据预处理（首次运行）
python 02_gene_abundance_analysis.py
python 03_alpha_diversity_analysis.py
python 04_beta_diversity_analysis.py
python 08_summary_analysis.py
```

## 输出结果

分析结果保存在 `output/` 目录：
- PNG格式图表
- Excel格式统计表
- TXT格式总结报告

## 目录结构

```
analysis-MicrobialGene/
├── *.py                   # 分析脚本
├── requirements.txt         # Python依赖
├── README.md              # 本文件
├── output/               # 分析结果（.gitignore）
├── processed_data/        # 处理后数据（.gitignore）
└── 微生物功能基因测序/    # 原始数据（.gitignore）
```
