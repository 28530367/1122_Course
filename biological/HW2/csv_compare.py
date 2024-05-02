import pandas as pd
import os

current_directory = os.path.dirname(os.path.abspath(__file__))

# 讀取csv
df1 = pd.read_csv(f'{current_directory}/Downloads/output/Dsim/Dsim_table.csv')
df2 = pd.read_csv(f'{current_directory}/output/Dsim/Dsim_table.csv')

# print(df1[['NCBI_GeneID', 'Symbol', 'Description', 'Gene_Type', 'Transcripts',
#        'Nomenclature ID', 'Chromosomes', 'Transcript_Genomic_Accession',
#        'Transcript_Genomic_Start', 'Transcript_Genomic_Stop', 'Orientation',
#        'Proteins', 'Protein_Accession', 'Transcript_Accession', 'Synonyms']].dtypes)
# print(df2[['NCBI_GeneID', 'Symbol', 'Description', 'Gene_Type', 'Transcripts',
#        'Nomenclature ID', 'Chromosomes', 'Transcript_Genomic_Accession',
#        'Transcript_Genomic_Start', 'Transcript_Genomic_Stop', 'Orientation',
#        'Proteins', 'Protein_Accession', 'Transcript_Accession', 'Synonyms']].dtypes)

# 找出不同行
different_rows = pd.concat([df1, df2]).drop_duplicates(subset=['NCBI_GeneID', 'Protein_Accession', 'Transcript_Accession'], keep=False)

print("\n不同行:")
print(different_rows)

