#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
磷循环功能基因丰度比较分析脚本
功能：比较三个功能基因在不同处理下的绝对丰度差异
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
from scipy import stats

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

# 处理显示名称
TREATMENT_NAMES = {
    'CK': 'CK (对照)',
    'F': 'F (化肥)',
    'FM': 'FM (化肥+有机肥)',
    'FS': 'FS (化肥+秸秆)',
    'FMS': 'FMS (化肥+秸秆+有机肥)'
}


def load_gene_abundance(gene_name):
    """加载指定基因的OTU丰度数据"""
    file_path = PROCESSED_DIR / f'{gene_name}_otu_processed.xlsx'

    if not file_path.exists():
        print(f"  [错误] {gene_name}: 文件不存在")
        return None

    try:
        df = pd.read_excel(file_path, engine='openpyxl')
        return df
    except Exception as e:
        print(f"  [错误] {gene_name}: {e}")
        return None


def calculate_total_abundance(df, gene_name):
    """
    计算每个样本的总丰度

    参数:
        df: OTU数据框
        gene_name: 基因名称

    返回:
        DataFrame: 每个样本的总丰度
    """
    # 提取样本列（绝对丰度）
    sample_columns = [col for col in df.columns if col in sum(AH_SAMPLES.values(), [])]

    # 计算每列的总和
    abundance_data = []
    for treatment, samples in AH_SAMPLES.items():
        for sample in samples:
            if sample in df.columns:
                total_abundance = df[sample].sum()
                abundance_data.append({
                    'Gene': gene_name,
                    'Treatment': treatment,
                    'Treatment_Name': TREATMENT_NAMES[treatment],
                    'Sample': sample,
                    'Total_Abundance': total_abundance
                })

    return pd.DataFrame(abundance_data)


def plot_gene_abundance_comparison(abundance_df):
    """
    绘制基因丰度比较图

    参数:
        abundance_df: 丰度数据框
    """
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # 图1: 箱线图 - 每个基因在不同处理下的丰度
    ax1 = axes[0]
    sns.boxplot(data=abundance_df, x='Treatment', y='Total_Abundance',
                hue='Gene', ax=ax1, palette='Set2')
    ax1.set_xlabel('处理方式', fontsize=12)
    ax1.set_ylabel('基因总丰度 (reads数)', fontsize=12)
    ax1.set_title('不同处理下磷循环功能基因丰度比较', fontsize=14, fontweight='bold')
    # 修复x轴标签
    treatment_order = sorted(abundance_df['Treatment'].unique())
    ax1.set_xticklabels([TREATMENT_NAMES.get(t, t) for t in treatment_order], rotation=45, ha='right')
    ax1.legend(title='基因', fontsize=10)
    ax1.grid(axis='y', alpha=0.3)

    # 图2: 柱状图 - 每个基因在不同处理下的平均丰度
    ax2 = axes[1]
    mean_abundance = abundance_df.groupby(['Gene', 'Treatment'])['Total_Abundance'].mean().reset_index()
    mean_abundance_pivot = mean_abundance.pivot(index='Treatment', columns='Gene', values='Total_Abundance')

    mean_abundance_pivot.plot(kind='bar', ax=ax2, width=0.8, colormap='Set2')
    ax2.set_xlabel('处理方式', fontsize=12)
    ax2.set_ylabel('平均基因丰度 (reads数)', fontsize=12)
    ax2.set_title('不同处理下磷循环功能基因平均丰度', fontsize=14, fontweight='bold')
    ax2.set_xticklabels([TREATMENT_NAMES[t] for t in mean_abundance_pivot.index], rotation=45, ha='right')
    ax2.legend(title='基因', fontsize=10)
    ax2.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    output_path = OUTPUT_DIR / 'gene_abundance_comparison.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def perform_statistical_tests(abundance_df):
    """进行统计检验"""
    results = []

    for gene in GENES:
        gene_data = abundance_df[abundance_df['Gene'] == gene]
        groups = [group['Total_Abundance'].values for name, group in gene_data.groupby('Treatment')]

        stat, p_value = stats.kruskal(*groups)

        results.append({
            'Gene': gene,
            'Test': 'Kruskal-Wallis',
            'Statistic': stat,
            'P_Value': p_value,
            'Significant': p_value < 0.05
        })

        # 如果有显著差异，进行两两比较（Wilcoxon检验）
        if p_value < 0.05:
            treatments = list(AH_SAMPLES.keys())
            for i in range(len(treatments)):
                for j in range(i + 1, len(treatments)):
                    t1, t2 = treatments[i], treatments[j]
                    data1 = gene_data[gene_data['Treatment'] == t1]['Total_Abundance'].values
                    data2 = gene_data[gene_data['Treatment'] == t2]['Total_Abundance'].values

                    stat, p_val = stats.mannwhitneyu(data1, data2, alternative='two-sided')

                    results.append({
                        'Gene': gene,
                        'Test': f'Wilcoxon: {t1} vs {t2}',
                        'Statistic': stat,
                        'P_Value': p_val,
                        'Significant': p_val < 0.05
                    })

    return pd.DataFrame(results)


def plot_gene_abundance_by_treatment(abundance_df):
    """
    分别绘制每个基因在不同处理下的丰度

    参数:
        abundance_df: 丰度数据框
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    for idx, gene in enumerate(GENES):
        ax = axes[idx]
        gene_data = abundance_df[abundance_df['Gene'] == gene]

        # 箱线图 + 散点图
        sns.boxplot(data=gene_data, x='Treatment', y='Total_Abundance',
                   ax=ax, color='lightgray')
        sns.stripplot(data=gene_data, x='Treatment', y='Total_Abundance',
                     ax=ax, color='darkblue', size=6, alpha=0.6)

        ax.set_xlabel('处理方式', fontsize=11)
        ax.set_ylabel('基因丰度 (reads数)', fontsize=11)
        ax.set_title(f'{gene} 基因丰度', fontsize=12, fontweight='bold')
        treatment_order = sorted(gene_data['Treatment'].unique())
        ax.set_xticklabels([TREATMENT_NAMES.get(t, t) for t in treatment_order], rotation=45, ha='right')
        ax.grid(axis='y', alpha=0.3)

    plt.suptitle('磷循环功能基因在不同处理下的丰度分布', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    output_path = OUTPUT_DIR / 'gene_abundance_by_treatment.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_gene_rank(abundance_df):
    """
    绘制基因丰度排名图

    参数:
        abundance_df: 丰度数据框
    """
    # 计算每个基因在每个处理下的平均丰度
    mean_abundance = abundance_df.groupby(['Gene', 'Treatment'])['Total_Abundance'].mean().reset_index()

    fig, axes = plt.subplots(1, 5, figsize=(20, 5))

    for idx, treatment in enumerate(AH_SAMPLES.keys()):
        ax = axes[idx]
        treatment_data = mean_abundance[mean_abundance['Treatment'] == treatment]

        # 按丰度排序
        treatment_data = treatment_data.sort_values('Total_Abundance', ascending=True)

        # 柱状图
        colors = ['lightcoral', 'lightblue', 'lightgreen']
        bars = ax.barh(treatment_data['Gene'], treatment_data['Total_Abundance'],
                      color=colors)

        # 添加数值标签
        for i, (bar, value) in enumerate(zip(bars, treatment_data['Total_Abundance'])):
            ax.text(value, bar.get_y() + bar.get_height()/2,
                   f'{int(value):,}',
                   va='center', ha='left', fontsize=10)

        ax.set_xlabel('平均基因丰度 (reads数)', fontsize=11)
        ax.set_title(TREATMENT_NAMES[treatment], fontsize=12, fontweight='bold')
        ax.grid(axis='x', alpha=0.3)

    plt.suptitle('不同处理下磷循环功能基因平均丰度排名', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    output_path = OUTPUT_DIR / 'gene_abundance_rank.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def generate_summary_table(abundance_df, stats_df):
    """
    生成汇总表格

    参数:
        abundance_df: 丰度数据框
        stats_df: 统计检验结果

    返回:
        DataFrame: 汇总表格
    """
    # 计算每个基因在每个处理下的平均丰度和标准差
    summary_stats = abundance_df.groupby(['Gene', 'Treatment'])['Total_Abundance'].agg(
        ['mean', 'std', 'min', 'max']
    ).reset_index()

    # 重命名列
    summary_stats.columns = ['Gene', 'Treatment', 'Mean', 'Std', 'Min', 'Max']

    # 添加处理显示名称
    summary_stats['Treatment_Name'] = summary_stats['Treatment'].map(TREATMENT_NAMES)

    # 计算CV（变异系数）
    summary_stats['CV'] = (summary_stats['Std'] / summary_stats['Mean'] * 100).round(2)

    # 重新排列列
    summary_stats = summary_stats[['Gene', 'Treatment', 'Treatment_Name', 'Mean',
                                   'Std', 'CV', 'Min', 'Max']]

    return summary_stats


def main():
    """主函数"""
    print("=" * 50)
    print("基因丰度分析")
    print("=" * 50)

    all_abundance_data = []
    for gene in GENES:
        df = load_gene_abundance(gene)
        if df is not None:
            gene_abundance = calculate_total_abundance(df, gene)
            all_abundance_data.append(gene_abundance)

    if not all_abundance_data:
        print("  [错误] 没有数据")
        return

    abundance_df = pd.concat(all_abundance_data, ignore_index=True)

    plot_gene_abundance_comparison(abundance_df)
    plot_gene_abundance_by_treatment(abundance_df)
    plot_gene_rank(abundance_df)

    stats_df = perform_statistical_tests(abundance_df)
    stats_output = OUTPUT_DIR / 'gene_abundance_statistical_tests.xlsx'
    stats_df.to_excel(stats_output, index=False, engine='openpyxl')

    summary_table = generate_summary_table(abundance_df, stats_df)
    summary_output = OUTPUT_DIR / 'gene_abundance_summary.xlsx'
    summary_table.to_excel(summary_output, index=False, engine='openpyxl')

    print("\n完成 →", OUTPUT_DIR)


if __name__ == '__main__':
    main()
