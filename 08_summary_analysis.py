#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
磷循环功能基因综合总结分析
功能：合并群落组成、差异分析和三基因对比的简化版本
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 设置路径
BASE_DIR = Path(r'C:\github\analysis-MicrobialGene')
PROCESSED_DIR = BASE_DIR / 'processed_data'
OUTPUT_DIR = BASE_DIR / 'output'
OUTPUT_DIR.mkdir(exist_ok=True)

# 基因列表
GENES = ['phoC', 'phoD', 'pqqC']

# 安徽样本信息
AH_SAMPLES = {
    'CK': ['AH_CK_1', 'AH_CK_2', 'AH_CK_3'],
    'F': ['AH_F_1', 'AH_F_2', 'AH_F_3'],
    'FM': ['AH_FM_1', 'AH_FM_2', 'AH_FM_3'],
    'FS': ['AH_FS_1', 'AH_FS_2', 'AH_FS_3'],
    'FMS': ['AH_FMS_1', 'AH_FMS_2', 'AH_FMS_3']
}

TREATMENT_NAMES = {
    'CK': 'CK (对照)',
    'F': 'F (化肥)',
    'FM': 'FM (化肥+有机肥)',
    'FS': 'FS (化肥+秸秆)',
    'FMS': 'FMS (化肥+秸秆+有机肥)'
}

AH_TREATMENTS = ['AH_CK', 'AH_F', 'AH_FM', 'AH_FS', 'AH_FMS']


def load_taxon_abundance(gene_name, taxon_level):
    """加载分类水平丰度数据"""
    file_path = PROCESSED_DIR / f'{gene_name}_{taxon_level}_relative_abundance.xlsx'

    if not file_path.exists():
        return None

    df = pd.read_excel(file_path, engine='openpyxl')
    return df


def plot_community_heatmap():
    """绘制群落组成热图（简化版）"""
    print("  群落组成热图...")

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    for idx, gene in enumerate(GENES):
        ax = axes[idx]
        taxon_df = load_taxon_abundance(gene, 'genus')

        if taxon_df is None:
            continue

        treatment_cols = [col for col in taxon_df.columns if col in AH_TREATMENTS]
        mean_abundance = taxon_df[treatment_cols].mean(axis=1).sort_values(ascending=False)

        # 取Top 10
        top_taxa = mean_abundance.head(10)
        plot_data = taxon_df.loc[top_taxa.index, treatment_cols]

        # 转置
        plot_data = plot_data.T

        im = ax.imshow(plot_data.values, cmap='YlOrRd', aspect='auto')
        ax.set_xticks(range(len(plot_data.columns)))
        ax.set_yticks(range(len(plot_data.index)))
        ax.set_xticklabels(plot_data.columns, rotation=90, fontsize=7)
        ax.set_yticklabels([TREATMENT_NAMES.get(t, t) for t in plot_data.index], fontsize=8)
        ax.set_title(f'{gene} - 属水平', fontsize=11, fontweight='bold')

    plt.suptitle('群落组成热图', fontsize=14, fontweight='bold')
    plt.tight_layout()

    output_path = OUTPUT_DIR / 'summary_community_heatmap.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_gene_comparison():
    """绘制三基因对比（简化版）"""
    print("  三基因对比...")

    # 加载α多样性数据
    alpha_data = {}
    for gene in GENES:
        file_path = PROCESSED_DIR / f'{gene}_alpha_diversity.xlsx'
        if file_path.exists():
            df = pd.read_excel(file_path, engine='openpyxl')
            alpha_data[gene] = df

    if not alpha_data:
        return

    # 提取Shannon多样性
    comparison_data = []
    for gene, df in alpha_data.items():
        sample_cols = [col for col in df.columns if col in sum(AH_SAMPLES.values(), [])]
        if not sample_cols:
            continue

        shannon_row = df[df['Samples'] == 'Shannon']
        if shannon_row.empty:
            continue

        for treatment, samples in AH_SAMPLES.items():
            treatment_samples = [s for s in samples if s in df.columns]
            if treatment_samples:
                mean_shannon = shannon_row.iloc[0][treatment_samples].mean()
                comparison_data.append({
                    'Gene': gene,
                    'Treatment': treatment,
                    'Shannon': mean_shannon
                })

    comparison_df = pd.DataFrame(comparison_data)

    # 绘制
    fig, ax = plt.subplots(figsize=(12, 6))

    for gene in GENES:
        gene_data = comparison_df[comparison_df['Gene'] == gene]
        ax.plot(gene_data['Treatment'], gene_data['Shannon'],
                marker='o', label=gene, linewidth=2, markersize=8)

    ax.set_xlabel('处理', fontsize=12)
    ax.set_ylabel('Shannon多样性', fontsize=12)
    ax.set_title('三基因Shannon多样性对比', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3)

    plt.tight_layout()

    output_path = OUTPUT_DIR / 'summary_gene_comparison.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def generate_summary_report():
    """生成简化报告"""
    print("  生成报告...")

    report = [
        "=" * 50,
        "磷循环功能基因综合分析报告",
        "=" * 50,
        "",
        "分析内容：",
        "  1. 群落组成热图（属水平）",
        "  2. 三基因多样性对比",
        "  3. 差异分析（基础）",
        "",
        "结论：",
        "  详见生成的图表文件",
        "=" * 50
    ]

    report_path = OUTPUT_DIR / 'summary_report.txt'
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(report))

    # 简单的差异分析（只做top类群）
    for gene in GENES:
        taxon_df = load_taxon_abundance(gene, 'genus')
        if taxon_df is None:
            continue

        treatment_cols = [col for col in taxon_df.columns if col in AH_TREATMENTS]
        mean_abundance = taxon_df[treatment_cols].mean(axis=1).sort_values(ascending=False)

        report.append(f"\n{gene} 优势菌属 (Top 5):")
        for i, (taxon, abundance) in enumerate(mean_abundance.head(5).items(), 1):
            report.append(f"  {i}. {taxon}: {abundance:.2%}")

    with open(report_path, 'a', encoding='utf-8') as f:
        f.write('\n'.join(report))


def main():
    """主函数"""
    print("=" * 50)
    print("综合总结分析")
    print("=" * 50)

    plot_community_heatmap()
    plot_gene_comparison()
    generate_summary_report()

    print("\n完成 →", OUTPUT_DIR)


if __name__ == '__main__':
    main()
