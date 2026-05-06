#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
磷循环功能基因α多样性分析脚本
功能：评估不同处理对磷循环微生物群落多样性的影响
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
from scipy import stats
from scipy.stats import f_oneway, kruskal

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

# 多样性指标
DIVERSITY_INDICES = ['Observed', 'Chao1', 'ACE', 'Shannon', 'Simpson', 'Coverage']


def load_alpha_diversity(gene_name):
    """加载指定基因的α多样性数据"""
    file_path = PROCESSED_DIR / f'{gene_name}_alpha_diversity.xlsx'

    if not file_path.exists():
        print(f"  [错误] {gene_name}: 文件不存在")
        return None

    df = pd.read_excel(file_path, engine='openpyxl')
    return df


def restructure_alpha_diversity(alpha_df, gene_name):
    """
    重构α多样性数据为长格式

    参数:
        alpha_df: α多样性数据框
        gene_name: 基因名称

    返回:
        DataFrame: 重构后的长格式数据框
    """
    if alpha_df is None:
        return None

    data = []

    for treatment, samples in AH_SAMPLES.items():
        for sample in samples:
            if sample in alpha_df.columns:
                for index in alpha_df.index:
                    diversity_index = alpha_df.loc[index, 'Samples']
                    value = alpha_df.loc[index, sample]

                    data.append({
                        'Gene': gene_name,
                        'Treatment': treatment,
                        'Treatment_Name': TREATMENT_NAMES[treatment],
                        'Sample': sample,
                        'Diversity_Index': diversity_index,
                        'Value': value
                    })

    return pd.DataFrame(data)


def plot_alpha_diversity_boxplot(alpha_long_df, diversity_index):
    """
    绘制α多样性箱线图

    参数:
        alpha_long_df: 长格式α多样性数据
        diversity_index: 多样性指标名称
    """
    # 筛选指定多样性指标的数据
    data = alpha_long_df[alpha_long_df['Diversity_Index'] == diversity_index].copy()

    # 检查是否有数据
    if data.empty:
        print(f"  警告: 没有找到 {diversity_index} 的数据")
        return None

    # 创建图形
    fig, ax = plt.subplots(figsize=(10, 6))

    # 绘制箱线图
    sns.boxplot(data=data, x='Treatment', y='Value', hue='Gene',
                ax=ax, palette='Set2', width=0.7)

    # 添加散点
    sns.stripplot(data=data, x='Treatment', y='Value', hue='Gene',
                 ax=ax, color='black', size=4, alpha=0.5, dodge=True,
                 legend=False)

    ax.set_xlabel('处理方式', fontsize=12)
    ax.set_ylabel(f'{diversity_index} 指数', fontsize=12)
    ax.set_title(f'{diversity_index} 多样性指数比较', fontsize=14, fontweight='bold')
    treatment_order = sorted(data['Treatment'].unique())
    ax.set_xticklabels([TREATMENT_NAMES.get(t, t) for t in treatment_order], rotation=45, ha='right')
    ax.legend(title='基因', fontsize=10)
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()

    output_path = OUTPUT_DIR / f'alpha_diversity_{diversity_index.lower()}.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    return output_path


def plot_alpha_diversity_by_gene(alpha_long_df):
    """
    分别绘制每个基因的α多样性

    参数:
        alpha_long_df: 长格式α多样性数据
    """
    for gene in GENES:
        gene_data = alpha_long_df[alpha_long_df['Gene'] == gene]

        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.flatten()

        for idx, diversity_index in enumerate(DIVERSITY_INDICES):
            ax = axes[idx]
            idx_data = gene_data[gene_data['Diversity_Index'] == diversity_index]

            if not idx_data.empty:
                sns.boxplot(data=idx_data, x='Treatment', y='Value',
                           ax=ax, color='lightgray')
                sns.stripplot(data=idx_data, x='Treatment', y='Value',
                             ax=ax, color='darkblue', size=6, alpha=0.6)

                ax.set_xlabel('处理方式', fontsize=10)
                ax.set_ylabel(f'{diversity_index}', fontsize=10)
                ax.set_title(f'{gene} - {diversity_index}', fontsize=11, fontweight='bold')
                treatment_order = sorted(idx_data['Treatment'].unique())
                ax.set_xticklabels([TREATMENT_NAMES.get(t, t) for t in treatment_order],
                                  rotation=45, ha='right', fontsize=8)
                ax.grid(axis='y', alpha=0.3)

        plt.suptitle(f'{gene} 基因α多样性分析', fontsize=14, fontweight='bold')
        plt.tight_layout()

        output_path = OUTPUT_DIR / f'alpha_diversity_{gene}_all_indices.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()


def perform_diversity_statistical_tests(alpha_long_df):
    """进行α多样性统计检验"""
    results = []

    # 对每个基因的每个多样性指标进行检验
    for gene in GENES:
        for diversity_index in DIVERSITY_INDICES:
            # 筛选数据
            data = alpha_long_df[(alpha_long_df['Gene'] == gene) &
                                (alpha_long_df['Diversity_Index'] == diversity_index)]

            if data.empty:
                continue

            # 准备各组数据
            groups = [group['Value'].values for name, group in data.groupby('Treatment')]

            # 进行Kruskal-Wallis检验（非参数检验）
            try:
                stat, p_value = kruskal(*groups)

                results.append({
                    'Gene': gene,
                    'Diversity_Index': diversity_index,
                    'Test': 'Kruskal-Wallis',
                    'Statistic': stat,
                    'P_Value': p_value,
                    'Significant': p_value < 0.05
                })

                # 如果有显著差异，进行两两比较
                if p_value < 0.05:
                    treatments = list(AH_SAMPLES.keys())
                    for i in range(len(treatments)):
                        for j in range(i + 1, len(treatments)):
                            t1, t2 = treatments[i], treatments[j]
                            data1 = data[data['Treatment'] == t1]['Value'].values
                            data2 = data[data['Treatment'] == t2]['Value'].values

                            stat, p_val = stats.mannwhitneyu(data1, data2, alternative='two-sided')

                            results.append({
                                'Gene': gene,
                                'Diversity_Index': diversity_index,
                                'Test': f'Wilcoxon: {t1} vs {t2}',
                                'Statistic': stat,
                                'P_Value': p_val,
                                'Significant': p_val < 0.05
                            })
            except Exception as e:
                pass

    return pd.DataFrame(results)


def calculate_diversity_summary(alpha_long_df):
    """
    计算α多样性汇总统计

    参数:
        alpha_long_df: 长格式α多样性数据

    返回:
        DataFrame: 汇总统计数据
    """
    # 按基因、处理、多样性指标分组计算统计量
    summary = alpha_long_df.groupby(['Gene', 'Treatment', 'Diversity_Index'])['Value'].agg(
        ['mean', 'std', 'min', 'max', 'count']
    ).reset_index()

    # 添加处理显示名称
    summary['Treatment_Name'] = summary['Treatment'].map(TREATMENT_NAMES)

    # 计算CV
    summary['CV'] = (summary['std'] / summary['mean'] * 100).round(2)

    # 重新排列列
    summary = summary[['Gene', 'Treatment', 'Treatment_Name', 'Diversity_Index',
                      'mean', 'std', 'CV', 'min', 'max', 'count']]

    return summary


def plot_diversity_heatmap(alpha_long_df):
    """
    绘制α多样性热图

    参数:
        alpha_long_df: 长格式α多样性数据
    """
    # 计算每个基因在每个处理下的平均多样性
    diversity_pivot = alpha_long_df.groupby(['Gene', 'Treatment', 'Diversity_Index'])['Value'].mean().reset_index()
    diversity_pivot = diversity_pivot.pivot(index=['Gene', 'Diversity_Index'],
                                           columns='Treatment', values='Value')

    # 标准化（按行）
    diversity_normalized = diversity_pivot.div(diversity_pivot.mean(axis=1), axis=0)

    # 创建图形
    fig, ax = plt.subplots(figsize=(12, 8))

    # 绘制热图
    sns.heatmap(diversity_normalized, annot=True, fmt='.2f', cmap='RdYlGn',
                cbar_kws={'label': '相对值（行标准化）'}, ax=ax, linewidths=0.5)

    ax.set_title('α多样性指数热图（行标准化）', fontsize=14, fontweight='bold')
    ax.set_xlabel('处理方式', fontsize=12)
    ax.set_ylabel('基因 - 多样性指标', fontsize=12)

    plt.tight_layout()

    output_path = OUTPUT_DIR / 'alpha_diversity_heatmap.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_diversity_correlation(alpha_long_df):
    """
    绘制多样性指标相关性分析图

    参数:
        alpha_long_df: 长格式α多样性数据
    """
    # 准备数据：每个样本的所有多样性指标
    diversity_matrix = []

    for gene in GENES:
        gene_data = alpha_long_df[alpha_long_df['Gene'] == gene]

        for sample in sum(AH_SAMPLES.values(), []):
            sample_data = gene_data[gene_data['Sample'] == sample]
            if not sample_data.empty:
                row = {'Gene': gene, 'Sample': sample}
                for diversity_index in DIVERSITY_INDICES:
                    idx_data = sample_data[sample_data['Diversity_Index'] == diversity_index]
                    if not idx_data.empty:
                        row[diversity_index] = idx_data['Value'].values[0]
                diversity_matrix.append(row)

    diversity_df = pd.DataFrame(diversity_matrix)

    # 为每个基因计算相关性
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    for idx, gene in enumerate(GENES):
        ax = axes[idx]
        gene_diversity = diversity_df[diversity_df['Gene'] == gene][DIVERSITY_INDICES]

        # 计算相关性矩阵
        corr_matrix = gene_diversity.corr()

        # 绘制热图
        sns.heatmap(corr_matrix, annot=True, fmt='.2f', cmap='coolwarm',
                   center=0, vmin=-1, vmax=1, ax=ax, square=True,
                   cbar_kws={'label': '相关系数'}, linewidths=0.5)

        ax.set_title(f'{gene} - 多样性指标相关性', fontsize=12, fontweight='bold')

    plt.tight_layout()

    output_path = OUTPUT_DIR / 'alpha_diversity_correlation.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def main():
    """主函数"""
    print("=" * 50)
    print("α多样性分析")
    print("=" * 50)

    all_alpha_data = []
    for gene in GENES:
        alpha_df = load_alpha_diversity(gene)
        if alpha_df is not None:
            alpha_long = restructure_alpha_diversity(alpha_df, gene)
            if alpha_long is not None:
                all_alpha_data.append(alpha_long)

    if not all_alpha_data:
        print("  [错误] 没有数据")
        return

    alpha_long_df = pd.concat(all_alpha_data, ignore_index=True)

    for diversity_index in DIVERSITY_INDICES:
        plot_alpha_diversity_boxplot(alpha_long_df, diversity_index)

    plot_alpha_diversity_by_gene(alpha_long_df)
    plot_diversity_heatmap(alpha_long_df)
    plot_diversity_correlation(alpha_long_df)

    stats_df = perform_diversity_statistical_tests(alpha_long_df)
    stats_output = OUTPUT_DIR / 'alpha_diversity_statistical_tests.xlsx'
    stats_df.to_excel(stats_output, index=False, engine='openpyxl')

    summary_df = calculate_diversity_summary(alpha_long_df)
    summary_output = OUTPUT_DIR / 'alpha_diversity_summary.xlsx'
    summary_df.to_excel(summary_output, index=False, engine='openpyxl')

    print("\n完成 →", OUTPUT_DIR)


if __name__ == '__main__':
    main()
