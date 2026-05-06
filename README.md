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
| `01_data_preprocessing.py` | 数据预处理 |
| `02_gene_abundance_analysis.py` | 基因丰度比较分析 |
| `03_alpha_diversity_analysis.py` | α多样性分析 |
| `04_beta_diversity_analysis.py` | β多样性分析 |
| `08_summary_analysis.py` | 综合分析 |

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

### 基因丰度分析 (02)
- `gene_abundance_comparison.png` - 丰度对比图
- `gene_abundance_by_treatment.png` - 各处理丰度分布
- `gene_abundance_rank.png` - 丰度排名
- `gene_abundance_statistical_tests.xlsx` - 统计检验结果
- `gene_abundance_summary.xlsx` - 汇总统计表

### α多样性分析 (03)
- `alpha_diversity_shannon.png` - Shannon多样性箱线图
- `alpha_diversity_observed.png` - 物种数箱线图
- `alpha_diversity_chao1.png` - Chao1箱线图
- `alpha_diversity_ace.png` - ACE箱线图
- `alpha_diversity_simpson.png` - Simpson箱线图
- `alpha_diversity_coverage.png` - 覆盖度箱线图
- `alpha_diversity_phoC_all_indices.png` - phoC多样性综合图
- `alpha_diversity_phoD_all_indices.png` - phoD多样性综合图
- `alpha_diversity_pqqC_all_indices.png` - pqqC多样性综合图
- `alpha_diversity_heatmap.png` - 多样性热图
- `alpha_diversity_correlation.png` - 指标相关性图
- `alpha_diversity_statistical_tests.xlsx` - 统计检验结果
- `alpha_diversity_summary.xlsx` - 汇总统计表

### β多样性分析 (04)
- `beta_diversity_pca_phoC.png` - phoC的PCA图
- `beta_diversity_pca_phoD.png` - phoD的PCA图
- `beta_diversity_pca_pqqC.png` - pqqC的PCA图
- `beta_diversity_nmds_phoC.png` - phoC的NMDS图
- `beta_diversity_nmds_phoD.png` - phoD的NMDS图
- `beta_diversity_nmds_pqqC.png` - pqqC的NMDS图
- `beta_diversity_heatmap_phoC.png` - phoC距离热图
- `beta_diversity_heatmap_phoD.png` - phoD距离热图
- `beta_diversity_heatmap_pqqC.png` - pqqC距离热图
- `beta_diversity_comparison.png` - 三基因β多样性对比
- `beta_diversity_anosim_results.xlsx` - ANOSIM检验结果

### 综合总结分析 (08)
- `summary_community_heatmap.png` - 群落组成热图
- `summary_gene_comparison.png` - 三基因对比
- `summary_report.txt` - 总结报告

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
