import os
import pandas as pd
import numpy as np
import time

# 取得當前工作目錄
current_directory = os.path.dirname(os.path.abspath(__file__))

def readDsimTable():
    file_name = f"Dsim_table.csv"
    dsim_file_path = os.path.join(current_directory, f"input/{file_name}")
    dsim_input_df = pd.read_csv(dsim_file_path, sep=',')

    dsim_df = dsim_input_df.drop(['Symbol', 'Description', 'Transcripts', 'Transcript_Accession', 'Nomenclature ID', 'Nomenclature ID2', 'Chromosomes', 'Transcript_Genomic_Accession', 'Transcript_Genomic_Start', 'Transcript_Genomic_Stop', 'Orientation', 'Proteins', 'Protein_Accession', 'Synonyms'], axis=1)

    return dsim_df

def readDsimCheck_csv(type):
    file_name = f"Dsim_check_{type}.csv"
    file_path = os.path.join(current_directory, f"output/{file_name}")
    check_df = pd.read_csv(file_path, sep=',')
    check_df.drop(['Dsim', 'Gene Type', 'evalue1', 'score1', 'qcover1', 'identity1', 'Dmel', 'evalue2', 'score2', 'qcover2', 'identity2', 'hit_back', 'Hit back name'], axis=1, inplace=True)  
    check_df.rename(columns={'Gene name': 'NCBI_GeneID', 'Dmel name': 'Dmel ID'}, inplace=True)
    check_df.replace('-', None, inplace=True)
    check_type1_df = check_df[check_df['Type'] == 1]
    check_type1_df = check_type1_df.groupby('NCBI_GeneID').agg({
        'Dmel ID': lambda x: [value for value in x if value is not None] if any(x) else None
    }).reset_index()
    check_df = check_df.groupby('NCBI_GeneID').agg({
        'Type': lambda x: [value for value in x if value is not None] if any(x) else None
    }).reset_index()
    check_df = pd.merge(check_df, check_type1_df, on='NCBI_GeneID', how='outer')

    return check_df

if __name__ == "__main__":
    # 記錄開始時間
    start_time = time.time()

    check_protein_df = readDsimCheck_csv('protein')
    check_transcript_df = readDsimCheck_csv('transcript')

    dsimTable_df = readDsimTable()

    p_dsimTable_df = dsimTable_df[dsimTable_df['Gene_Type'] == 'PROTEIN_CODING']
    t_dsimTable_df = dsimTable_df[dsimTable_df['Gene_Type'] != 'PROTEIN_CODING']

    p_dsimTable_df = pd.merge(p_dsimTable_df, check_protein_df, on='NCBI_GeneID', how='left').drop_duplicates(subset=['NCBI_GeneID'])
    t_dsimTable_df = pd.merge(t_dsimTable_df, check_transcript_df, on='NCBI_GeneID', how='left').drop_duplicates(subset=['NCBI_GeneID'])
    
    dsimClassification_df = pd.concat([p_dsimTable_df, t_dsimTable_df], ignore_index=True, sort=False)
    dsimClassification_df = dsimClassification_df.sort_values(by='NCBI_GeneID', ascending=True)
    dsimClassification_df['Ortholog_Num'] = dsimClassification_df['Type'].apply(lambda x: x.count(1) if isinstance(x, list) else 0)
    
    # classification.csv
    file_name = f"Dsim_classification.csv"
    dsim_file_path = os.path.join(current_directory, f"output/{file_name}")
    dsimClassification_df.to_csv(dsim_file_path, columns=['NCBI_GeneID', 'Gene_Type', 'Ortholog_Num', 'Dmel ID', 'Type'] , index=False)

    # Protein-coding and Non-coding
    p_dsimClassification_df = dsimClassification_df[dsimClassification_df['Gene_Type'] == 'PROTEIN_CODING']
    t_dsimClassification_df = dsimClassification_df[dsimClassification_df['Gene_Type'] != 'PROTEIN_CODING']

    # unique
    dsimOrtholog_unique_df = dsimClassification_df[dsimClassification_df['Ortholog_Num'] == 1]
    dsimOrtholog_unique_df['Result Type'] = 1
    # print(dsimOrtholog_unique_df)
    dsimOrtholog_unique_df.to_csv(f"{current_directory}/output/Dsim_ortholog_unique.csv", columns=['NCBI_GeneID', 'Gene_Type', 'Ortholog_Num', 'Dmel ID', 'Type', 'Result Type'] , index=False)
    
    # unique Protein-coding and Non-coding
    p_dsimOrtholog_unique_df = dsimOrtholog_unique_df[dsimOrtholog_unique_df['Gene_Type'] == 'PROTEIN_CODING']
    t_dsimOrtholog_unique_df = dsimOrtholog_unique_df[dsimOrtholog_unique_df['Gene_Type'] != 'PROTEIN_CODING']

    # none
    dsimOrtholog_none_df = dsimClassification_df[dsimClassification_df['Ortholog_Num'] == 0]
    dsimOrtholog_none_df['Result Type'] = 3
    dsimOrtholog_none_df = dsimOrtholog_none_df.fillna('-')
    # print(dsimOrtholog_none_df)
    dsimOrtholog_none_df.to_csv(f"{current_directory}/output/Dsim_ortholog_none.csv", columns=['NCBI_GeneID', 'Gene_Type', 'Ortholog_Num', 'Dmel ID', 'Type', 'Result Type'] , index=False)

    # none Protein-coding and Non-coding
    p_dsimOrtholog_none_df = dsimOrtholog_none_df[dsimOrtholog_none_df['Gene_Type'] == 'PROTEIN_CODING']
    t_dsimOrtholog_none_df = dsimOrtholog_none_df[dsimOrtholog_none_df['Gene_Type'] != 'PROTEIN_CODING']

    # multiple
    dsimOrtholog_mulitple_df = dsimClassification_df[dsimClassification_df['Ortholog_Num'] > 1]
    dsimOrtholog_mulitple_df['Result Type'] = 2
    # print(dsimOrtholog_mulitple_df)
    dsimOrtholog_mulitple_df.to_csv(f"{current_directory}/output/Dsim_ortholog_multiple.csv", columns=['NCBI_GeneID', 'Gene_Type', 'Ortholog_Num', 'Dmel ID', 'Type', 'Result Type'] , index=False)

    # multiple Protein-coding and Non-coding
    p_dsimOrtholog_mulitple_df = dsimOrtholog_mulitple_df[dsimOrtholog_mulitple_df['Gene_Type'] == 'PROTEIN_CODING']
    t_dsimOrtholog_mulitple_df = dsimOrtholog_mulitple_df[dsimOrtholog_mulitple_df['Gene_Type'] != 'PROTEIN_CODING']

    # Dsim.Readme
    with open(f'{current_directory}/output/Dsim.Readme', 'w') as file:
        file.write(f'在 Dsim 中共有 {len(dsimClassification_df)} 個基因\n')
        file.write(f'class 1.具有唯一答案(unique)\n')
        file.write(f'    在這個分類共有 {len(dsimOrtholog_unique_df)} 個基因，佔全部的 {len(dsimOrtholog_unique_df)/len(dsimClassification_df)*100:.2f}%\n')
        file.write(f'class 2.有多個答案(multiple)\n')
        file.write(f'    在這個分類共有 {len(dsimOrtholog_mulitple_df)} 個基因，佔全部的 {len(dsimOrtholog_mulitple_df)/len(dsimClassification_df)*100:.2f}%\n')
        file.write(f'class 3.沒答案\n')
        file.write(f'    在這個分類共有 {len(dsimOrtholog_none_df)} 個基因，佔全部的 {len(dsimOrtholog_none_df)/len(dsimClassification_df)*100:.2f}%\n\n')

        file.write(f'在 Dsim 中共有 {len(dsimClassification_df)} 個基因，其中 Protein-coding 有 {len(p_dsimClassification_df)} 個，佔全部的 {len(p_dsimClassification_df)/len(dsimClassification_df)*100:.2f}%\n')
        file.write(f'class 1.具有唯一答案(unique)\n')
        file.write(f'    在這個分類共有 {len(p_dsimOrtholog_unique_df)} 個基因，佔全部的 {len(p_dsimOrtholog_unique_df)/len(dsimClassification_df)*100:.2f}%\n')
        file.write(f'class 2.有多個答案(multiple)\n')
        file.write(f'    在這個分類共有 {len(p_dsimOrtholog_mulitple_df)} 個基因，佔全部的 {len(p_dsimOrtholog_mulitple_df)/len(dsimClassification_df)*100:.2f}%\n')
        file.write(f'class 3.沒答案\n')
        file.write(f'    在這個分類共有 {len(p_dsimOrtholog_none_df)} 個基因，佔全部的 {len(p_dsimOrtholog_none_df)/len(dsimClassification_df)*100:.2f}%\n\n')

        file.write(f'在 Dsim 中共有 {len(dsimClassification_df)} 個基因，其中 Non-coding 有 {len(t_dsimClassification_df)} 個，佔全部的 {len(t_dsimClassification_df)/len(dsimClassification_df)*100:.2f}%\n')
        file.write(f'class 1.具有唯一答案(unique)\n')
        file.write(f'    在這個分類共有 {len(t_dsimOrtholog_unique_df)} 個基因，佔全部的 {len(t_dsimOrtholog_unique_df)/len(dsimClassification_df)*100:.2f}%\n')
        file.write(f'class 2.有多個答案(multiple)\n')
        file.write(f'    在這個分類共有 {len(t_dsimOrtholog_mulitple_df)} 個基因，佔全部的 {len(t_dsimOrtholog_mulitple_df)/len(dsimClassification_df)*100:.2f}%\n')
        file.write(f'class 3.沒答案\n')
        file.write(f'    在這個分類共有 {len(t_dsimOrtholog_none_df)} 個基因，佔全部的 {len(t_dsimOrtholog_none_df)/len(dsimClassification_df)*100:.2f}%\n\n')

    # 記錄結束時間
    end_time = time.time()
    # 計算執行時間
    execution_time = end_time - start_time
    print(f"執行時間：{execution_time} 秒")