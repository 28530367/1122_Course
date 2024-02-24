import os
import pandas as pd
import numpy as np

# 取得當前工作目錄
current_directory = os.getcwd()

def bestHit(id, left, right):
    
    # 構建相對路徑
    file_name = f"blast{id}_{left}_inDB_{right}.csv"
    relative_path = f"input/{file_name}"
    file_path = os.path.join(current_directory, relative_path)

    column_names = [left, 'E_value', 'Score', 'query cover', 'identity', right, 'unknown1', 'unknown2']
    input_df = pd.read_csv(file_path, sep='\t', names=column_names)

    input_df.drop_duplicates(keep="first")

    # 移除 'unknown1' 和 'unknown2' 列
    input_df = input_df.drop(['unknown1', 'unknown2'], axis=1)

    # 按照 'E_value', 'Score', 'query cover', 'identity' 做篩選
    # 篩出每個 left 分組下 'E_value' 最小的值
    # 使用 transform 找出每個 left 分組下 'E_value' 最小的值
    min_evalue = input_df.groupby(left)['E_value'].transform('min')
    filter1_df = input_df[input_df['E_value'] == min_evalue]

    # 從上一個結果篩出每個 left 分組下 'Score' 最大的值
    max_score = filter1_df.groupby(left)['Score'].transform('max')
    filter2_df = filter1_df[filter1_df['Score'] == max_score]

    # 從上一個結果篩出每個 left 分組下 'query cover' 最大的值
    max_qcover = filter2_df.groupby(left)['query cover'].transform('max')
    result3_df = filter2_df[filter2_df['query cover'] == max_qcover]

    # 從上一個結果篩出每個 left 分組下 'identity' 最大的值
    max_identity = result3_df.groupby(left)['identity'].transform('max')
    result_df = result3_df[result3_df['identity'] == max_identity]

    # 最後篩選結果將 right 用 ', ' 合併 
    relay_df = result_df.groupby(left).agg({'E_value': 'first', 'Score': 'first', 'query cover': 'first', 'identity': 'first', right: ', '.join}).reset_index()

    # 構建相對路徑
    output_name = f"blast{id}_{left}to{right}_besthit.csv"
    output_relative_path = f"bestHit/{output_name}"
    output_path = os.path.join(current_directory, output_relative_path)

    # 輸出成csv
    relay_df.to_csv(output_path, index=False)

    return result_df

def readTable(id):
    file_name = f"Dsim_table.csv"
    dsim_file_path = os.path.join(current_directory, f"input/{file_name}")
    dsim_input_df = pd.read_csv(dsim_file_path, sep=',')

    file_name = f"Dmel_table.csv"
    dmel_file_path = os.path.join(current_directory, f"input/{file_name}")
    dmel_input_df = pd.read_csv(dmel_file_path, sep=',')

    if id == 'p':
        geneType_df = dsim_input_df.drop(['Symbol', 'Description', 'Transcripts', 'Transcript_Accession', 'Nomenclature ID', 'Nomenclature ID2', 'Chromosomes', 'Transcript_Genomic_Accession', 'Transcript_Genomic_Start', 'Transcript_Genomic_Stop', 'Orientation', 'Proteins', 'Synonyms'], axis=1)
        dsim_df = dsim_input_df.drop(['Symbol', 'Description', 'Gene_Type', 'Transcripts', 'Transcript_Accession', 'Nomenclature ID', 'Nomenclature ID2', 'Chromosomes', 'Transcript_Genomic_Accession', 'Transcript_Genomic_Start', 'Transcript_Genomic_Stop', 'Orientation', 'Proteins', 'Synonyms'], axis=1)
        dmel_df = dmel_input_df.drop(['Symbol', 'Description', 'Gene_Type', 'Transcripts', 'Transcript_Accession', 'Nomenclature ID', 'Nomenclature ID2', 'Chromosomes', 'Transcript_Genomic_Accession', 'Transcript_Genomic_Start', 'Transcript_Genomic_Stop', 'Orientation', 'Proteins', 'Synonyms'], axis=1)

        dsim_df = dsim_df[dsim_df['Protein_Accession'] != '-']
        dmel_df = dmel_df[dmel_df['Protein_Accession'] != '-']

        geneType_df.rename(columns={'NCBI_GeneID': 'Gene name', 'Gene_Type': 'Gene Type', 'Protein_Accession': 'Dsim_x'}, inplace=True)
        dsim_df.rename(columns={'NCBI_GeneID': 'Hit back name', 'Protein_Accession': 'Dsim_y'}, inplace=True)
        dmel_df.rename(columns={'NCBI_GeneID': 'Dmel name', 'Protein_Accession': 'Dmel'}, inplace=True)

    elif id == 'n':
        geneType_df = dsim_input_df.drop(['Symbol', 'Description', 'Transcripts', 'Nomenclature ID', 'Nomenclature ID2', 'Chromosomes', 'Transcript_Genomic_Accession', 'Transcript_Genomic_Start', 'Transcript_Genomic_Stop', 'Orientation', 'Proteins', 'Protein_Accession', 'Synonyms'], axis=1)
        dsim_df = dsim_input_df.drop(['Symbol', 'Description', 'Gene_Type', 'Transcripts', 'Nomenclature ID', 'Nomenclature ID2', 'Chromosomes', 'Transcript_Genomic_Accession', 'Transcript_Genomic_Start', 'Transcript_Genomic_Stop', 'Orientation', 'Proteins', 'Protein_Accession', 'Synonyms'], axis=1)
        dmel_df = dmel_input_df.drop(['Symbol', 'Description', 'Gene_Type', 'Transcripts', 'Nomenclature ID', 'Nomenclature ID2', 'Chromosomes', 'Transcript_Genomic_Accession', 'Transcript_Genomic_Start', 'Transcript_Genomic_Stop', 'Orientation', 'Proteins', 'Protein_Accession', 'Synonyms'], axis=1)

        dsim_df = dsim_df[dsim_df['Transcript_Accession'] != '-']
        dmel_df = dmel_df[dmel_df['Transcript_Accession'] != '-']
    
        geneType_df.rename(columns={'NCBI_GeneID': 'Gene name', 'Gene_Type': 'Gene Type', 'Transcript_Accession': 'Dsim_x'}, inplace=True)
        dsim_df.rename(columns={'NCBI_GeneID': 'Hit back name', 'Transcript_Accession': 'Dsim_y'}, inplace=True)
        dmel_df.rename(columns={'NCBI_GeneID': 'Dmel name', 'Transcript_Accession': 'Dmel'}, inplace=True)
    
    return geneType_df, dsim_df, dmel_df

def dsim_check_csv(id, merge_df, dsimGeneType_df, dmelTable_df, dsimTable_df):
    dismCheck_df = pd.merge(merge_df, dsimGeneType_df, on='Dsim_x', how='left')
    dismCheck_df = pd.merge(dismCheck_df, dmelTable_df, on='Dmel', how='left')
    dismCheck_df = pd.merge(dismCheck_df, dsimTable_df, on='Dsim_y', how='left')

    dismCheck_df.rename(columns={'Dsim_x': 'Dsim', 'E_value_x': 'evalue1', 'Score_x': 'score1', 'query cover_x': 'qcover1', 'identity_x': 'identity1', 'E_value_y': 'evalue2', 'Score_y': 'score2', 'query cover_y': 'qcover2', 'identity_y': 'identity2', 'Dsim_y': 'hit_back', 'type': 'Type'}, inplace=True)

    # 輸出csv
    if id == 'p':
        output_file_name = 'Dsim_check_protein.csv'
    elif id == 'n':
        output_file_name = 'Dsim_check_transcript.csv'
    output_path = os.path.join(current_directory, f'output/{output_file_name}')
    dismCheck_df.to_csv(output_path, columns=['Gene name', 'Dsim', 'Gene Type', 'evalue1', 'score1', 'qcover1', 'identity1', 'Dmel', 'Dmel name', 'evalue2', 'score2', 'qcover2', 'identity2', 'hit_back', 'Hit back name', 'Type'] ,index=False)


if __name__ == "__main__":
    # 輸入
    id = 'n'

    Dsim_to_Dmel = bestHit(id, 'Dsim', 'Dmel')
    Dmel_to_Dsim = bestHit(id, 'Dmel', 'Dsim')

    merge_result_df = pd.merge(Dsim_to_Dmel, Dmel_to_Dsim, on='Dmel', how='left')

    # type
    # 沒有type4啦幹
    # 沒有type4啦幹
    # 沒有type4啦幹
    # 沒有type4啦幹
    # 沒有type4啦幹
    # 沒有type4啦幹
    # 沒有type4啦幹
    # 沒有type4啦幹
    # 沒有type4啦幹
    # 沒有type4啦幹
    # 沒有type4啦幹
    # 沒有type4啦幹
    merge_result_df['type'] = np.where((merge_result_df['Dsim_x'] == merge_result_df['Dsim_y']) & (~merge_result_df['Dsim_y'].isna()), 1, 2)
    merge_result_df.loc[merge_result_df['Dsim_y'].isna(), 'type'] = 3
    


    # # 輸出csv
    # output_path = os.path.join(current_directory, 'merge_type_result.csv')
    # merge_result_df.to_csv(output_path, index=False)

    dsimGeneType_df, dsimTable_df, dmelTable_df = readTable(id)
    # print(dsimGeneType_df)
    # print(dsimTable_df)
    # print(dmelTable_df)

    dsim_check_csv(id, merge_result_df, dsimGeneType_df, dmelTable_df, dsimTable_df)