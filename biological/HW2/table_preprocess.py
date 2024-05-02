from Bio import SeqIO
import os
import re
import pandas as pd

current_directory = os.path.dirname(os.path.abspath(__file__))

def protein_length_table(target, csv=False):
    # protein_length_table.csv
    fasta_file = "protein.faa"
    fasta_file_path = f"{current_directory}/data/{target}/{fasta_file}"

    protein_length_dict = {}
    protein_length_dict['NCBI GeneID'] = []
    protein_length_dict['Protein Accession'] = []
    protein_length_dict['Protein_length'] = []

    for record in SeqIO.parse(fasta_file_path, "fasta"):
        name = record.name
        sequence = record.seq
        description = record.description
        gene_id_match = re.search(r'GeneID=(\d+)', description)
        gene_id = gene_id_match.group(1)

        protein_length_dict['NCBI GeneID'].append(gene_id)
        protein_length_dict['Protein Accession'].append(name)
        protein_length_dict['Protein_length'].append(len(sequence))

    protein_length_df = pd.DataFrame(protein_length_dict)
    protein_length_df.drop_duplicates(keep='first', inplace=True)
    # print(protein_length_df)

    if csv == True:
        # 輸出中繼檔 protein_length_table.csv
        protein_length_table_path = f"{current_directory}/tmp/{target}/protein_length_table.csv"
        protein_length_df.to_csv(protein_length_table_path, index=False)

    return protein_length_df

def transcript_length_table(target, csv=False):
    # transcript_length_table.csv
    fasta_file = "rna.fna"
    fasta_file_path = f"{current_directory}/data/{target}/{fasta_file}"

    transcript_length_dict = {}
    transcript_length_dict['NCBI GeneID'] = []
    transcript_length_dict['Transcript Accession'] = []
    transcript_length_dict['Transcript_length'] = []

    for record in SeqIO.parse(fasta_file_path, "fasta"):
        name = record.name
        sequence = record.seq
        description = record.description
        gene_id_match = re.search(r'GeneID=(\d+)', description)
        gene_id = gene_id_match.group(1)

        transcript_length_dict['NCBI GeneID'].append(gene_id)
        transcript_length_dict['Transcript Accession'].append(name)
        transcript_length_dict['Transcript_length'].append(len(sequence))

    transcript_length_df = pd.DataFrame(transcript_length_dict)
    transcript_length_df.drop_duplicates(keep='first', inplace=True)
    # print(transcript_length_df)

    if csv == True:
        # 輸出中繼檔 transcript_length_table.csv
        transcript_length_table_path = f"{current_directory}/tmp/{target}/transcript_length_table.csv"
        transcript_length_df.to_csv(transcript_length_table_path, index=False)

    return transcript_length_df

def protein_transcript_pair(target, csv=False):

    def extract_id(text):
        start_index = text.find("ID=cds-") + len("ID=cds-")
        end_index = text.find(";", start_index)
        return text[start_index:end_index]
    
    def extract_parent(text):
        start_index = text.find("Parent=rna-") + len("Parent=rna-")
        end_index = text.find(";", start_index)
        return text[start_index:end_index]

    file_name = "genomic.gff"
    file_path = f"{current_directory}/data/{target}/{file_name}"

    genomic_df = pd.read_csv(file_path, sep='\t', comment='#', header=None)
    # print(genomic_df)

    cds_df = genomic_df[genomic_df[2] == 'CDS']
    # print(cds_df)

    protein_transcript_pair_df = pd.DataFrame()
    protein_transcript_pair_df['Protein Accession'] = cds_df[8].apply(extract_id)
    protein_transcript_pair_df['Transcript Accession'] = cds_df[8].apply(extract_parent)
    # print(protein_transcript_pair_df)
    protein_transcript_pair_df.drop_duplicates(keep='first', inplace=True)

    if csv == True:
        # 輸出中繼檔 protein_transcript_pair.csv
        protein_transcript_pair_path = f"{current_directory}/tmp/{target}/protein_transcript_pair.csv"
        protein_transcript_pair_df.to_csv(protein_transcript_pair_path, index=False)

    return protein_transcript_pair_df

def pr_tr_pair_all(target, protein_length_df, transcript_length_df, protein_transcript_pair_df, csv = False):
    pr_tr_pair_all_df = pd.merge(protein_transcript_pair_df, protein_length_df, on='Protein Accession')
    pr_tr_pair_all_df = pd.merge(pr_tr_pair_all_df, transcript_length_df[['Transcript Accession', 'Transcript_length']], on='Transcript Accession')
    pr_tr_pair_all_df = pr_tr_pair_all_df[['NCBI GeneID', 'Protein Accession', 'Protein_length', 'Transcript Accession', 'Transcript_length']]
    pr_tr_pair_all_df.drop_duplicates(keep='first', inplace=True)
    # print(pr_tr_pair_all_df)

    if csv == True:
        # 輸出中繼檔 pr_tr_pair_all.csv
        pr_tr_pair_all_df_path = f"{current_directory}/tmp/{target}/pr_tr_pair_all.csv"
        pr_tr_pair_all_df.to_csv(pr_tr_pair_all_df_path, index=False)
    return pr_tr_pair_all_df

def pr_tr_pair_represent(target, pr_tr_pair_all_df, csv = False):
    def extract_float_from_accession(accession):
        match = re.search(r'\d+\.\d+', accession)
        return float(match.group())
    
    pr_tr_pair_all_df['Protein Accession Float'] = pr_tr_pair_all_df['Protein Accession'].apply(extract_float_from_accession)
    pr_tr_pair_all_df['Transcript Accession Float'] = pr_tr_pair_all_df['Transcript Accession'].apply(extract_float_from_accession)
    # print(pr_tr_pair_all_df)

    # 按照 'Protein_length', 'Transcript_length', 'Protein Accession Float' 做篩選
    max_Protein_length = pr_tr_pair_all_df.groupby('NCBI GeneID')['Protein_length'].transform('max')
    filter1_df = pr_tr_pair_all_df[pr_tr_pair_all_df['Protein_length'] == max_Protein_length]

    max_Transcript_length = filter1_df.groupby('NCBI GeneID')['Transcript_length'].transform('max')
    filter2_df = filter1_df[filter1_df['Transcript_length'] == max_Transcript_length]

    min_Protein_Accession_Float = filter2_df.groupby('NCBI GeneID')['Protein Accession Float'].transform('min')
    pr_tr_pair_represent_df = filter2_df[filter2_df['Protein Accession Float'] == min_Protein_Accession_Float]
    pr_tr_pair_represent_df = pr_tr_pair_represent_df[['NCBI GeneID', 'Protein Accession', 'Protein_length', 'Transcript Accession', 'Transcript_length']]
    # print(pr_tr_pair_represent_df)

    if csv == True:
        # 輸出中繼檔 pr_tr_pair_represent.csv
        pr_tr_pair_represent_df_path = f"{current_directory}/tmp/{target}/pr_tr_pair_represent.csv"
        pr_tr_pair_represent_df.to_csv(pr_tr_pair_represent_df_path, index=False)
    return pr_tr_pair_represent_df

def transcript_represent(target, transcript_length_df, csv = False):
    def extract_float_from_accession(accession):
        match = re.search(r'\d+\.\d+', accession)
        return float(match.group())
    
    transcript_length_df['Transcript Accession Float'] = transcript_length_df['Transcript Accession'].apply(extract_float_from_accession)

    # 按照 'Transcript_length', 'Transcript_Accession_Float' 做篩選
    max_Transcript_length = transcript_length_df.groupby('NCBI GeneID')['Transcript_length'].transform('max')
    filter1_df = transcript_length_df[transcript_length_df['Transcript_length'] == max_Transcript_length]

    min_Transcript_Accession_Float = filter1_df.groupby('NCBI GeneID')['Transcript Accession Float'].transform('min')
    transcript_represent_df = filter1_df[filter1_df['Transcript Accession Float'] == min_Transcript_Accession_Float]
    transcript_represent_df = transcript_represent_df[['NCBI GeneID', 'Transcript Accession', 'Transcript_length']]

    if csv == True:
        # 輸出中繼檔 transcript_represent.csv
        transcript_represent_df_path = f"{current_directory}/tmp/{target}/transcript_represent.csv"
        transcript_represent_df.to_csv(transcript_represent_df_path, index=False)
    return transcript_represent_df

def protein_length_duplicate(target, pr_tr_pair_all_df, csv = False):
    # 按照 'Protein_length', 'Transcript_length', 'Protein Accession Float' 做篩選
    max_Protein_length = pr_tr_pair_all_df.groupby('NCBI GeneID')['Protein_length'].transform('max')
    filter1_df = pr_tr_pair_all_df[pr_tr_pair_all_df['Protein_length'] == max_Protein_length]

    max_Transcript_length = filter1_df.groupby('NCBI GeneID')['Transcript_length'].transform('max')
    filter2_df = filter1_df[filter1_df['Transcript_length'] == max_Transcript_length]

    duplicate_df = filter2_df.groupby('NCBI GeneID').filter(lambda x: len(x) > 1)
    protein_length_duplicate_df = duplicate_df.groupby('NCBI GeneID').agg({'Protein Accession': lambda x: ', '.join(x), 'Protein_length': 'first', 'Transcript_length': 'first'}).reset_index()
    
    if csv == True:
        # 輸出中繼檔 protein_length_duplicate.csv
        protein_length_duplicate_df_path = f"{current_directory}/tmp/{target}/protein_length_duplicate.csv"
        protein_length_duplicate_df.to_csv(protein_length_duplicate_df_path, index=False)

    return protein_length_duplicate_df

def protein_no_Accession(target, pr_tr_pair_all_df, csv = False):
    file_name = "ncbi_dataset.tsv"
    file_path = f"{current_directory}/data/{target}/{file_name}"

    ncbi_dataset_df = pd.read_csv(file_path, sep='\t')
    filter_ncbi_dataset_df = ncbi_dataset_df[['NCBI GeneID', 'Symbol', 'Gene Type']]
    filter_ncbi_dataset_df = filter_ncbi_dataset_df[filter_ncbi_dataset_df['Gene Type'] == 'PROTEIN_CODING']

    pr_tr_pair_all_df['NCBI GeneID'] = pr_tr_pair_all_df['NCBI GeneID'].astype(int)

    merge_df = pd.merge(filter_ncbi_dataset_df, pr_tr_pair_all_df, on='NCBI GeneID', how='left')
    
    nan_rows = merge_df[merge_df['Protein Accession'].isna()]
    print(nan_rows)

def result_table(target, pr_tr_pair_represent_df, transcript_represent_df, csv = False):
    file_name = "ncbi_dataset.tsv"
    file_path = f"{current_directory}/data/{target}/{file_name}"

    ncbi_dataset_df = pd.read_csv(file_path, sep='\t')
    coding_ncbi_dataset_df = ncbi_dataset_df[ncbi_dataset_df['Gene Type'] == 'PROTEIN_CODING']
    noncoding_dataset_df = ncbi_dataset_df[ncbi_dataset_df['Gene Type'] != 'PROTEIN_CODING']

    pr_tr_pair_represent_df = pr_tr_pair_represent_df.drop(columns=['Protein_length', 'Transcript_length'])
    pr_tr_pair_represent_df['NCBI GeneID'] = pr_tr_pair_represent_df['NCBI GeneID'].astype(int)

    transcript_represent_df = transcript_represent_df.drop(columns=['Transcript_length'])
    transcript_represent_df['NCBI GeneID'] = transcript_represent_df['NCBI GeneID'].astype(int)
   
    coding_merge_df = pd.merge(coding_ncbi_dataset_df, pr_tr_pair_represent_df, on='NCBI GeneID', how='left')
    noncoding_merge_df = pd.merge(noncoding_dataset_df, transcript_represent_df, on='NCBI GeneID', how='left')

    result_table_df = pd.concat([coding_merge_df, noncoding_merge_df])
    result_table_df = result_table_df.fillna('-')
    
    result_table_df = result_table_df[['NCBI GeneID', 'Symbol', 'Description', 'Gene Type', 'Transcripts', 'Nomenclature ID', 'Chromosomes', 'Annotation Genomic Range Accession ', 'Annotation Genomic Range Start', 'Annotation Genomic Range Stop', 'Orientation', 'Proteins', 'Protein Accession', 'Transcript Accession', 'Synonyms']]
    
    result_table_df = result_table_df.rename(columns={
        'NCBI GeneID': 'NCBI_GeneID',
        'Gene Type': 'Gene_Type',
        'Annotation Genomic Range Accession ': 'Transcript_Genomic_Accession',
        'Annotation Genomic Range Start': 'Transcript_Genomic_Start',
        'Annotation Genomic Range Stop': 'Transcript_Genomic_Stop',
        'Protein Accession': 'Protein_Accession',
        'Transcript Accession': 'Transcript_Accession'
        })

    if csv == True:
        # 輸出output {target}_table.csv
        result_table_df_path = f"{current_directory}/output/{target}/{target}_table.csv"
        result_table_df.to_csv(result_table_df_path, index=False)

    return coding_merge_df, noncoding_merge_df

def fa_file_output(target, coding_df, noncoding_df):
    noncoding_symbol_list = noncoding_df['Transcript Accession'].tolist()
    coding_symbol_list = coding_df['Protein Accession'].tolist()

    rna_file = "rna.fna"
    rna_file_path = f"{current_directory}/data/{target}/{rna_file}"
    filtered_records = []
    with open(rna_file_path, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            # print(record.id)
            if record.id in noncoding_symbol_list:
                filtered_records.append(record)
                

    blastn_file_path = f"{current_directory}/output/{target}/{target}_for_blastn"
    with open(blastn_file_path, 'w') as f:
        SeqIO.write(filtered_records, f, 'fasta')
        


if __name__ == "__main__":
    target = "Dsim"
    protein_length_df = protein_length_table(target)
    transcript_length_df = transcript_length_table(target)
    protein_transcript_pair_df = protein_transcript_pair(target) 
    pr_tr_pair_all_df = pr_tr_pair_all(target, protein_length_df, transcript_length_df, protein_transcript_pair_df)

    pr_tr_pair_represent_df = pr_tr_pair_represent(target, pr_tr_pair_all_df)

    transcript_represent_df = transcript_represent(target, transcript_length_df)

    protein_length_duplicate_df = protein_length_duplicate(target, pr_tr_pair_all_df)

    coding_df, noncoding_df = result_table(target, pr_tr_pair_represent_df, transcript_represent_df, True)

    fa_file_output(target, coding_df, noncoding_df)
    print('done')
    