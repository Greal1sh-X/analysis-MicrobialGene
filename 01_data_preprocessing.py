#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
磷循环功能基因测序数据预处理脚本
功能：读取原始数据，进行数据清洗和质控评估
"""

import pandas as pd
import numpy as np
import os
import warnings
from pathlib import Path

warnings.filterwarnings('ignore')

# 设置路径
BASE_DIR = Path(r'C:\github\analysis-MicrobialGene')
DATA_DIR = BASE_DIR / '微生物功能基因测序'

# 三个基因的文件夹名称
GENE_FOLDERS = {
    'phoC': '23B1201B_PHOC_1_完整分析_20240114',
    'phoD': '23B1201B_PHOD_1_完整分析_20240115',
    'pqqC': '23B1201B_PQQC_1_完整分析_20240111'
}

# 安徽样本名称（只分析AH样本）
AH_SAMPLES = [
    'AH_CK_1', 'AH_CK_2', 'AH_CK_3',
    'AH_F_1', 'AH_F_2', 'AH_F_3',
    'AH_FM_1', 'AH_FM_2', 'AH_FM_3',
    'AH_FS_1', 'AH_FS_2', 'AH_FS_3',
    'AH_FMS_1', 'AH_FMS_2', 'AH_FMS_3'
]

# 处理映射
TREATMENT_MAP = {
    'CK': 'CK (对照)',
    'F': 'F (化肥)',
    'FM': 'FM (化肥+有机肥)',
    'FS': 'FS (化肥+秸秆)',
    'FMS': 'FMS (化肥+秸秆+有机肥)'
}

# 安徽样本对应的数字列名（从Sample_Group文件中获取）
AH_SAMPLE_IDS = {
    'AH_CK_1': '28',
    'AH_CK_2': '29',
    'AH_CK_3': '30',
    'AH_F_1': '25',
    'AH_F_2': '26',
    'AH_F_3': '27',
    'AH_FM_1': '19',
    'AH_FM_2': '20',
    'AH_FM_3': '21',
    'AH_FS_1': '22',
    'AH_FS_2': '23',
    'AH_FS_3': '24',
    'AH_FMS_1': '16',
    'AH_FMS_2': '17',
    'AH_FMS_3': '18'
}

# 反向映射：数字列名 -> 样本名称
ID_TO_SAMPLE = {v: k for k, v in AH_SAMPLE_IDS.items()}

# 安徽处理名称（用于分类水平数据）
AH_TREATMENTS = ['AH_CK', 'AH_F', 'AH_FM', 'AH_FS', 'AH_FMS']


def read_otu_table(otu_file_path):
    """读取OTU丰度表"""
    try:
        df = pd.read_excel(otu_file_path, engine='openpyxl')
        return df
    except:
        try:
            df = pd.read_csv(otu_file_path, sep='\t', encoding='utf-8')
            return df
        except UnicodeDecodeError:
            df = pd.read_csv(otu_file_path, sep='\t', encoding='gbk')
            return df
        except Exception as e:
            print(f"  [错误] 读取OTU表失败: {e}")
            return None


def extract_sample_columns(otu_df, sample_list):
    """
    从OTU表中提取指定样本的列

    参数:
        otu_df: OTU数据框
        sample_list: 样本名称列表

    返回:
        DataFrame: 只包含指定样本的OTU表
    """
    # 将样本名称转换为对应的数字列名
    sample_ids = [AH_SAMPLE_IDS.get(sample) for sample in sample_list]
    sample_ids = [sid for sid in sample_ids if sid is not None]

    # 查找包含样本ID的列
    sample_columns = [col for col in otu_df.columns if str(col) in sample_ids]

    if not sample_columns:
        return None

    # 重命名列：将数字ID替换为样本名称
    rename_dict = {col: ID_TO_SAMPLE.get(str(col), col) for col in sample_columns}

    # 提取样本列和分类信息
    tax_cols = ['#OTU ID', 'Taxonomy', 'superkingdom', 'phylum', 'class',
                'order', 'family', 'genus', 'species', 'OTUsize']
    existing_tax_cols = [col for col in tax_cols if col in otu_df.columns]

    result_df = otu_df[sample_columns + existing_tax_cols].copy()
    result_df = result_df.rename(columns=rename_dict)

    return result_df


def clean_taxonomy(taxonomy_str):
    """
    清理分类学信息字符串

    参数:
        taxonomy_str: 原始分类学字符串

    返回:
        dict: 清理后的分类学信息字典
    """
    if pd.isna(taxonomy_str) or taxonomy_str == '':
        return {
            'superkingdom': 'Unclassified',
            'phylum': 'Unclassified',
            'class': 'Unclassified',
            'order': 'Unclassified',
            'family': 'Unclassified',
            'genus': 'Unclassified',
            'species': 'Unclassified'
        }

    result = {
        'superkingdom': 'Unclassified',
        'phylum': 'Unclassified',
        'class': 'Unclassified',
        'order': 'Unclassified',
        'family': 'Unclassified',
        'genus': 'Unclassified',
        'species': 'Unclassified'
    }

    # 解析分类学字符串
    parts = taxonomy_str.split(';')
    for part in parts:
        part = part.strip()
        if '{superkingdom}' in part:
            result['superkingdom'] = part.split('{superkingdom}')[0].strip()
        elif '{phylum}' in part:
            result['phylum'] = part.split('{phylum}')[0].strip()
        elif '{class}' in part:
            result['class'] = part.split('{class}')[0].strip()
        elif '{order}' in part:
            result['order'] = part.split('{order}')[0].strip()
        elif '{family}' in part:
            result['family'] = part.split('{family}')[0].strip()
        elif '{genus}' in part:
            result['genus'] = part.split('{genus}')[0].strip()
        elif '{species}' in part:
            result['species'] = part.split('{species}')[0].strip()

    return result


def parse_taxonomy_column(df):
    """
    解析Taxonomy列，提取各分类水平信息

    参数:
        df: 包含Taxonomy列的数据框

    返回:
        DataFrame: 添加了各分类水平列的数据框
    """
    if 'Taxonomy' not in df.columns:
        print("  警告: 未找到Taxonomy列")
        return df

    # 解析每个OTU的分类学信息
    taxonomy_data = df['Taxonomy'].apply(clean_taxonomy)
    taxonomy_df = pd.DataFrame(taxonomy_data.tolist())

    # 将解析结果添加到原数据框
    result_df = pd.concat([df, taxonomy_df], axis=1)

    return result_df


def calculate_relative_abundance(otu_df, sample_columns):
    """
    计算相对丰度

    参数:
        otu_df: OTU数据框
        sample_columns: 样本列名列表

    返回:
        DataFrame: 添加了相对丰度列的数据框
    """
    result_df = otu_df.copy()

    # 计算每个样本的总丰度
    for sample in sample_columns:
        total = otu_df[sample].sum()
        if total > 0:
            result_df[f'{sample}_rel'] = otu_df[sample] / total
        else:
            result_df[f'{sample}_rel'] = 0

    return result_df


def read_alpha_diversity(alpha_div_file_path):
    """读取α多样性指数数据"""
    try:
        df = pd.read_excel(alpha_div_file_path, engine='openpyxl')
        return df
    except:
        try:
            df = pd.read_csv(alpha_div_file_path, sep='\t', encoding='utf-8')
            return df
        except UnicodeDecodeError:
            df = pd.read_csv(alpha_div_file_path, sep='\t', encoding='gbk')
            return df
        except Exception as e:
            print(f"  [错误] 读取α多样性数据失败: {e}")
            return None


def extract_ah_alpha_diversity(alpha_df, sample_list):
    """
    从α多样性数据中提取安徽样本

    参数:
        alpha_df: α多样性数据框
        sample_list: 样本名称列表

    返回:
        DataFrame: 只包含AH样本的α多样性数据
    """
    if alpha_df is None:
        return None

    # 查找包含样本名称的列
    sample_columns = [col for col in alpha_df.columns if col in sample_list]

    if not sample_columns:
        return None

    # 提取样本列和Samples列
    cols_to_extract = ['Samples'] + sample_columns if 'Samples' in alpha_df.columns else sample_columns
    result_df = alpha_df[cols_to_extract].copy()

    return result_df


def read_taxon_abundance(taxon_file_path):
    """读取分类水平丰度数据"""
    try:
        df = pd.read_excel(taxon_file_path, engine='openpyxl')
        return df
    except:
        try:
            df = pd.read_csv(taxon_file_path, sep='\t', encoding='utf-8')
            return df
        except UnicodeDecodeError:
            df = pd.read_csv(taxon_file_path, sep='\t', encoding='gbk')
            return df
        except Exception as e:
            print(f"  [错误] 读取分类水平丰度数据失败: {e}")
            return None


def process_gene_data(gene_name, folder_name, output_dir):
    """处理单个基因的数据"""
    gene_dir = DATA_DIR / folder_name / 'result'

    # 1. 读取OTU表
    otu_file = gene_dir / 'Total_OTU' / 'otu.xls'
    if not otu_file.exists():
        print(f"  [错误] {gene_name}: OTU表不存在")
        return

    otu_df = read_otu_table(otu_file)
    if otu_df is None:
        return

    # 2. 提取安徽样本
    otu_ah = extract_sample_columns(otu_df, AH_SAMPLES)
    if otu_ah is None:
        return

    # 3. 解析分类学信息
    otu_ah_parsed = parse_taxonomy_column(otu_ah)

    # 4. 计算相对丰度
    sample_columns = [col for col in otu_ah_parsed.columns if col in AH_SAMPLES]
    otu_with_rel = calculate_relative_abundance(otu_ah_parsed, sample_columns)

    # 5. 保存处理后的OTU数据
    output_file = output_dir / f'{gene_name}_otu_processed.xlsx'
    otu_with_rel.to_excel(output_file, index=False, engine='openpyxl')

    # 6. 读取α多样性数据
    alpha_div_file = (gene_dir / 'Analysis' / 'G1-Compare_1' /
                     'AlphaDiversity' / '1_DiversityIndex' / 'alpha.diversity.index.xls')

    if alpha_div_file.exists():
        alpha_df = read_alpha_diversity(alpha_div_file)
        alpha_ah = extract_ah_alpha_diversity(alpha_df, AH_SAMPLES)
        if alpha_ah is not None:
            alpha_output = output_dir / f'{gene_name}_alpha_diversity.xlsx'
            alpha_ah.to_excel(alpha_output, index=False, engine='openpyxl')

    # 7. 读取分类水平丰度数据
    taxon_dir = (gene_dir / 'Analysis' / 'G1-Compare_1' /
                 'Community' / '1_TaxonLevel' / 'group')

    taxon_levels = ['phylum', 'class', 'order', 'family', 'genus', 'species']

    for level in taxon_levels:
        rel_file = taxon_dir / f'{level}.taxon.RelativeAbundance.xls'
        if rel_file.exists():
            taxon_df = read_taxon_abundance(rel_file)
            if taxon_df is not None:
                ah_cols = [col for col in taxon_df.columns if col in AH_TREATMENTS]
                taxon_cols = [col for col in taxon_df.columns if col not in AH_TREATMENTS and col != 'superkingdom']

                if ah_cols:
                    taxon_ah = taxon_df[ah_cols + taxon_cols].copy()
                    taxon_output = output_dir / f'{gene_name}_{level}_relative_abundance.xlsx'
                    taxon_ah.to_excel(taxon_output, index=False, engine='openpyxl')


def main():
    """主函数"""
    print("=" * 50)
    print("数据预处理")
    print("=" * 50)

    output_dir = BASE_DIR / 'processed_data'
    output_dir.mkdir(exist_ok=True)

    for gene_name, folder_name in GENE_FOLDERS.items():
        process_gene_data(gene_name, folder_name, output_dir)
        print(f"[{gene_name}] OK")

    print(f"\n完成 → {output_dir}")


if __name__ == '__main__':
    main()
