import os
import pandas as pd
import numpy as np
import time

# 取得當前工作目錄
current_directory = os.path.dirname(os.path.abspath(__file__))

def readDismTable():
    file_name = f"Dsim_table.csv"
    dsim_file_path = os.path.join(current_directory, f"input/{file_name}")
    dsim_input_df = pd.read_csv(dsim_file_path, sep=',')

    dsim_df = dsim_input_df.drop(['Symbol', 'Description', 'Transcripts', 'Transcript_Accession', 'Nomenclature ID', 'Nomenclature ID2', 'Chromosomes', 'Transcript_Genomic_Accession', 'Transcript_Genomic_Start', 'Transcript_Genomic_Stop', 'Orientation', 'Proteins', 'Protein_Accession', 'Synonyms'], axis=1)

    return dsim_df

if __name__ == "__main__":
    # 記錄開始時間
    start_time = time.time()

    file_name = f"Dsim_check_protein.csv"
    protein_file_path = os.path.join(current_directory, f"output/{file_name}")
    check_protein_df = pd.read_csv(protein_file_path, sep=',')
    check_protein_df.drop(['Dsim', 'Gene Type', 'evalue1', 'score1', 'qcover1', 'identity1', 'Dmel', 'evalue2', 'score2', 'qcover2', 'identity2', 'hit_back', 'Hit back name'], axis=1, inplace=True)  
    check_protein_df.rename(columns={'Gene name': 'NCBI_GeneID', 'Dmel name': 'Dmel ID'}, inplace=True)
    check_protein_df.replace('-', None, inplace=True)
    check_protein_type1_df = check_protein_df[check_protein_df['Type'] == 1]
    check_protein_type1_df = check_protein_type1_df.groupby('NCBI_GeneID').agg({
        'Dmel ID': lambda x: [value for value in x if value is not None] if any(x) else None
    }).reset_index()
    check_protein_df = check_protein_df.groupby('NCBI_GeneID').agg({
        'Type': lambda x: [value for value in x if value is not None] if any(x) else None
    }).reset_index()
    check_protein_df = pd.merge(check_protein_df, check_protein_type1_df, on='NCBI_GeneID', how='outer')

    file_name = f"Dsim_check_transcript.csv"
    transcript_file_path = os.path.join(current_directory, f"output/{file_name}")
    check_transcript_df = pd.read_csv(transcript_file_path, sep=',')
    check_transcript_df.drop(['Dsim', 'Gene Type', 'evalue1', 'score1', 'qcover1', 'identity1', 'Dmel', 'evalue2', 'score2', 'qcover2', 'identity2', 'hit_back', 'Hit back name'], axis=1, inplace=True)
    check_transcript_df.rename(columns={'Gene name': 'NCBI_GeneID', 'Dmel name': 'Dmel ID'}, inplace=True)
    check_transcript_df.replace('-', None, inplace=True)
    check_transcript_type1_df = check_transcript_df[check_transcript_df['Type'] == 1]
    check_transcript_type1_df = check_transcript_type1_df.groupby('NCBI_GeneID').agg({
        'Dmel ID': lambda x: [value for value in x if value is not None] if any(x) else None
    }).reset_index()
    check_transcript_df = check_transcript_df.groupby('NCBI_GeneID').agg({
        'Type': lambda x: [value for value in x if value is not None] if any(x) else None
    }).reset_index()
    check_transcript_df = pd.merge(check_transcript_df, check_transcript_type1_df, on='NCBI_GeneID', how='outer')

    dsimTable_df = readDismTable()

    p_dsimTable_df = dsimTable_df[dsimTable_df['Gene_Type'] == 'PROTEIN_CODING']
    t_dsimTable_df = dsimTable_df[dsimTable_df['Gene_Type'] != 'PROTEIN_CODING']

    p_dsimTable_df = pd.merge(p_dsimTable_df, check_protein_df, on='NCBI_GeneID', how='left').drop_duplicates(subset=['NCBI_GeneID'])
    t_dsimTable_df = pd.merge(t_dsimTable_df, check_transcript_df, on='NCBI_GeneID', how='left').drop_duplicates(subset=['NCBI_GeneID'])
    
    dsimClassification_df = pd.concat([p_dsimTable_df, t_dsimTable_df], ignore_index=True, sort=False)
    dsimClassification_df = dsimClassification_df.sort_values(by='NCBI_GeneID', ascending=True)

    dsimClassification_df['Ortholog_Num'] = dsimClassification_df['Type'].apply(lambda x: x.count(1) if isinstance(x, list) else 0)
    
    file_name = f"Dsim_classification.csv"
    dsim_file_path = os.path.join(current_directory, f"output/{file_name}")
    dsimClassification_df.to_csv(dsim_file_path, columns=['NCBI_GeneID', 'Gene_Type', 'Ortholog_Num', 'Dmel ID', 'Type'] , index=False)

    
    # filtered_dismCheck_df = drop_dismCheck_df[drop_dismCheck_df['Type'] == 1]

    # grouped_dismCheck_df = filtered_dismCheck_df.groupby(['NCBI_GeneID']).agg({'Dmel ID': list_of_dmels}).reset_index()
    # grouped_dismCheck_df['Ortholog_Num'] = grouped_dismCheck_df['Dmel ID'].apply(len)
    # # print(grouped_dismCheck_df)

    # dismTable = readDismTable()
    # # print(dismTable)

    # dsim_classification_df = pd.merge(dismTable, grouped_dismCheck_df, on='NCBI_GeneID', how='left')
    # dsim_classification_df['Ortholog_Num'].fillna(0, inplace=True)
    # dsim_classification_df = dsim_classification_df.sort_values(by='NCBI_GeneID')
    # print(dsim_classification_df)





    # 記錄結束時間
    end_time = time.time()
    # 計算執行時間
    execution_time = end_time - start_time
    print(f"執行時間：{execution_time} 秒")