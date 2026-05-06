#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
磷循环功能基因β多样性分析脚本
功能：评估不同处理间群落组成的相似性/差异性
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import MDS
from scipy.spatial.distance import braycurtis, euclidean
from scipy.stats import permutation_test

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

# 处理颜色
TREATMENT_COLORS = {
    'CK': '#1f77b4',
    'F': '#ff7f0e',
    'FM': '#2ca02c',
    'FS': '#d62728',
    'FMS': '#9467bd'
}


def load_otu_data(gene_name):
    """加载指定基因的OTU数据"""
    file_path = PROCESSED_DIR / f'{gene_name}_otu_processed.xlsx'

    if not file_path.exists():
        print(f"  [错误] {gene_name}: 文件不存在")
        return None

    df = pd.read_excel(file_path, engine='openpyxl')

    # 提取样本列（绝对丰度）
    sample_columns = [col for col in df.columns if col in sum(AH_SAMPLES.values(), [])]
    otu_data = df[sample_columns].copy()

    return otu_data


def calculate_relative_abundance(otu_df):
    """
    计算相对丰度

    参数:
        otu_df: OTU数据框（绝对丰度）

    返回:
        DataFrame: 相对丰度数据框
    """
    # 计算每列的总和
    col_sums = otu_df.sum(axis=0)

    # 计算相对丰度
    rel_abundance = otu_df.div(col_sums, axis=1)

    return rel_abundance


def calculate_bray_curtis_distance(rel_abundance_df):
    """
    计算Bray-Curtis距离矩阵

    参数:
        rel_abundance_df: 相对丰度数据框

    返回:
        numpy.ndarray: 距离矩阵
    """
    samples = rel_abundance_df.columns
    n_samples = len(samples)
    distance_matrix = np.zeros((n_samples, n_samples))

    for i in range(n_samples):
        for j in range(i, n_samples):
            dist = braycurtis(rel_abundance_df.iloc[:, i], rel_abundance_df.iloc[:, j])
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist

    return distance_matrix, samples


def perform_pca(rel_abundance_df, gene_name):
    """
    执行PCA分析

    参数:
        rel_abundance_df: 相对丰度数据框
        gene_name: 基因名称

    返回:
        tuple: (PCA结果, 样本信息)
    """
    # 转置：行=样本，列=OTU
    data = rel_abundance_df.T

    # 标准化
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data)

    # 执行PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(data_scaled)

    # 创建结果数据框
    pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'])
    pca_df['Sample'] = data.index

    # 添加处理信息
    sample_info = []
    for sample in pca_df['Sample']:
        for treatment, samples in AH_SAMPLES.items():
            if sample in samples:
                sample_info.append(treatment)
                break

    pca_df['Treatment'] = sample_info
    pca_df['Treatment_Name'] = pca_df['Treatment'].map(TREATMENT_NAMES)

    # 计算解释方差
    explained_variance = pca.explained_variance_ratio_ * 100

    return pca_df, explained_variance


def perform_mds(rel_abundance_df, gene_name):
    """
    执行MDS分析

    参数:
        rel_abundance_df: 相对丰度数据框
        gene_name: 基因名称

    返回:
        tuple: (MDS结果, 样本信息, 应力值)
    """
    # 计算距离矩阵
    distance_matrix, samples = calculate_bray_curtis_distance(rel_abundance_df)

    # 执行MDS
    mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
    mds_result = mds.fit_transform(distance_matrix)

    # 创建结果数据框
    mds_df = pd.DataFrame(mds_result, columns=['MDS1', 'MDS2'])
    mds_df['Sample'] = samples

    # 添加处理信息
    sample_info = []
    for sample in mds_df['Sample']:
        for treatment, t_samples in AH_SAMPLES.items():
            if sample in t_samples:
                sample_info.append(treatment)
                break

    mds_df['Treatment'] = sample_info
    mds_df['Treatment_Name'] = mds_df['Treatment'].map(TREATMENT_NAMES)

    stress = mds.stress_

    return mds_df, stress


def plot_pca(pca_df, explained_variance, gene_name):
    """
    绘制PCA图

    参数:
        pca_df: PCA结果数据框
        explained_variance: 解释方差
        gene_name: 基因名称
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    # 按处理分组绘制
    for treatment in AH_SAMPLES.keys():
        treatment_data = pca_df[pca_df['Treatment'] == treatment]
        ax.scatter(treatment_data['PC1'], treatment_data['PC2'],
                  label=TREATMENT_NAMES[treatment],
                  color=TREATMENT_COLORS[treatment],
                  s=100, alpha=0.7, edgecolors='black', linewidth=1.5)

        # 添加样本标签
        for _, row in treatment_data.iterrows():
            ax.annotate(row['Sample'], (row['PC1'], row['PC2']),
                       fontsize=8, alpha=0.6)

    ax.set_xlabel(f'PC1 ({explained_variance[0]:.2f}%)', fontsize=12)
    ax.set_ylabel(f'PC2 ({explained_variance[1]:.2f}%)', fontsize=12)
    ax.set_title(f'{gene_name} 基因 - PCA分析', fontsize=14, fontweight='bold')
    ax.legend(title='处理', fontsize=10, loc='best')
    ax.grid(alpha=0.3)

    plt.tight_layout()

    output_path = OUTPUT_DIR / f'beta_diversity_pca_{gene_name}.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_mds(mds_df, stress, gene_name):
    """
    绘制MDS图

    参数:
        mds_df: MDS结果数据框
        stress: 应力值
        gene_name: 基因名称
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    # 按处理分组绘制
    for treatment in AH_SAMPLES.keys():
        treatment_data = mds_df[mds_df['Treatment'] == treatment]
        ax.scatter(treatment_data['MDS1'], treatment_data['MDS2'],
                  label=TREATMENT_NAMES[treatment],
                  color=TREATMENT_COLORS[treatment],
                  s=100, alpha=0.7, edgecolors='black', linewidth=1.5)

        # 添加样本标签
        for _, row in treatment_data.iterrows():
            ax.annotate(row['Sample'], (row['MDS1'], row['MDS2']),
                       fontsize=8, alpha=0.6)

    ax.set_xlabel('MDS1', fontsize=12)
    ax.set_ylabel('MDS2', fontsize=12)
    ax.set_title(f'{gene_name} 基因 - NMDS分析 (Stress = {stress:.4f})', fontsize=14, fontweight='bold')
    ax.legend(title='处理', fontsize=10, loc='best')
    ax.grid(alpha=0.3)

    plt.tight_layout()

    output_path = OUTPUT_DIR / f'beta_diversity_nmds_{gene_name}.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_distance_heatmap(distance_matrix, samples, gene_name):
    """
    绘制Bray-Curtis距离热图

    参数:
        distance_matrix: 距离矩阵
        samples: 样本名称列表
        gene_name: 基因名称
    """
    # 提取处理信息用于分组
    sample_treatments = []
    for sample in samples:
        for treatment, t_samples in AH_SAMPLES.items():
            if sample in t_samples:
                sample_treatments.append(treatment)
                break

    # 创建分组颜色
    treatment_to_num = {t: i for i, t in enumerate(AH_SAMPLES.keys())}
    treatment_nums = [treatment_to_num[t] for t in sample_treatments]

    # 绘制热图
    fig, ax = plt.subplots(figsize=(12, 10))

    im = ax.imshow(distance_matrix, cmap='YlOrRd', aspect='auto')

    # 设置刻度
    ax.set_xticks(range(len(samples)))
    ax.set_yticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=90, fontsize=8)
    ax.set_yticklabels(samples, fontsize=8)

    # 添加颜色条
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Bray-Curtis 距离', fontsize=12)

    # 添加分组线
    current_treatment = None
    for i, treatment in enumerate(sample_treatments):
        if i > 0 and treatment != current_treatment:
            ax.axhline(i, color='white', linewidth=2)
            ax.axvline(i, color='white', linewidth=2)
        current_treatment = treatment

    ax.set_title(f'{gene_name} 基因 - Bray-Curtis 距离热图', fontsize=14, fontweight='bold')

    plt.tight_layout()

    output_path = OUTPUT_DIR / f'beta_diversity_heatmap_{gene_name}.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def perform_anosim(distance_matrix, treatments):
    """
    执行ANOSIM分析（Analysis of Similarities）

    参数:
        distance_matrix: 距离矩阵
        treatments: 处理分组列表

    返回:
        dict: ANOSIM结果
    """
    # 计算组内和组间距离
    n = len(treatments)
    within_distances = []
    between_distances = []

    for i in range(n):
        for j in range(i+1, n):
            if treatments[i] == treatments[j]:
                within_distances.append(distance_matrix[i, j])
            else:
                between_distances.append(distance_matrix[i, j])

    # 计算R值
    r_between = np.mean(between_distances)
    r_within = np.mean(within_distances)
    r_value = (r_between - r_within) / (r_between + r_within) if (r_between + r_within) > 0 else 0

    # 置换检验
    n_permutations = 999
    observed_r = r_value
    greater_or_equal = 0

    for _ in range(n_permutations):
        permuted_treatments = np.random.permutation(treatments)

        perm_within = []
        perm_between = []

        for i in range(n):
            for j in range(i+1, n):
                if permuted_treatments[i] == permuted_treatments[j]:
                    perm_within.append(distance_matrix[i, j])
                else:
                    perm_between.append(distance_matrix[i, j])

        perm_r_between = np.mean(perm_between)
        perm_r_within = np.mean(perm_within)
        perm_r = (perm_r_between - perm_r_within) / (perm_r_between + perm_r_within) if (perm_r_between + perm_r_within) > 0 else 0

        if perm_r >= observed_r:
            greater_or_equal += 1

    p_value = (greater_or_equal + 1) / (n_permutations + 1)

    return {
        'R_value': r_value,
        'P_value': p_value,
        'Within_mean': r_within,
        'Between_mean': r_between
    }


def plot_beta_diversity_comparison(pca_results, mds_results):
    """
    绘制三个基因的β多样性对比图

    参数:
        pca_results: PCA结果字典
        mds_results: MDS结果字典
    """
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    # 第一行：PCA
    for idx, gene in enumerate(GENES):
        ax = axes[0, idx]
        pca_df = pca_results[gene]['data']

        for treatment in AH_SAMPLES.keys():
            treatment_data = pca_df[pca_df['Treatment'] == treatment]
            ax.scatter(treatment_data['PC1'], treatment_data['PC2'],
                      label=TREATMENT_NAMES[treatment],
                      color=TREATMENT_COLORS[treatment],
                      s=80, alpha=0.7, edgecolors='black', linewidth=1)

        ev = pca_results[gene]['explained_variance']
        ax.set_xlabel(f'PC1 ({ev[0]:.1f}%)', fontsize=10)
        ax.set_ylabel(f'PC2 ({ev[1]:.1f}%)', fontsize=10)
        ax.set_title(f'{gene} - PCA', fontsize=11, fontweight='bold')
        ax.grid(alpha=0.3)

    # 第二行：NMDS
    for idx, gene in enumerate(GENES):
        ax = axes[1, idx]
        mds_df = mds_results[gene]['data']

        for treatment in AH_SAMPLES.keys():
            treatment_data = mds_df[mds_df['Treatment'] == treatment]
            ax.scatter(treatment_data['MDS1'], treatment_data['MDS2'],
                      label=TREATMENT_NAMES[treatment],
                      color=TREATMENT_COLORS[treatment],
                      s=80, alpha=0.7, edgecolors='black', linewidth=1)

        stress = mds_results[gene]['stress']
        ax.set_xlabel('MDS1', fontsize=10)
        ax.set_ylabel('MDS2', fontsize=10)
        ax.set_title(f'{gene} - NMDS (Stress={stress:.3f})', fontsize=11, fontweight='bold')
        ax.grid(alpha=0.3)

    # 添加图例
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, title='处理', loc='center right',
              bbox_to_anchor=(0.98, 0.5), fontsize=10)

    plt.suptitle('磷循环功能基因β多样性分析对比', fontsize=14, fontweight='bold')
    plt.tight_layout()

    output_path = OUTPUT_DIR / 'beta_diversity_comparison.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def main():
    """主函数"""
    print("=" * 50)
    print("β多样性分析")
    print("=" * 50)

    pca_results = {}
    mds_results = {}
    anosim_results = []

    for gene in GENES:
        otu_data = load_otu_data(gene)
        if otu_data is None:
            continue

        rel_abundance = calculate_relative_abundance(otu_data)

        pca_df, explained_variance = perform_pca(rel_abundance, gene)
        pca_results[gene] = {'data': pca_df, 'explained_variance': explained_variance}

        plot_pca(pca_df, explained_variance, gene)

        mds_df, stress = perform_mds(rel_abundance, gene)
        mds_results[gene] = {'data': mds_df, 'stress': stress}

        plot_mds(mds_df, stress, gene)

        distance_matrix, samples = calculate_bray_curtis_distance(rel_abundance)
        plot_distance_heatmap(distance_matrix, samples, gene)

        treatments = [t for s in samples for t, t_s in AH_SAMPLES.items() if s in t_s]
        anosim_result = perform_anosim(distance_matrix, treatments)
        anosim_result['Gene'] = gene
        anosim_results.append(anosim_result)

    plot_beta_diversity_comparison(pca_results, mds_results)

    anosim_df = pd.DataFrame(anosim_results)
    anosim_output = OUTPUT_DIR / 'beta_diversity_anosim_results.xlsx'
    anosim_df.to_excel(anosim_output, index=False, engine='openpyxl')

    print("\n完成 →", OUTPUT_DIR)


if __name__ == '__main__':
    main()
