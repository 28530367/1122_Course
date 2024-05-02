import pandas as pd
import os

current_directory = os.path.dirname(os.path.abspath(__file__))

dsim_table_df = pd.read_csv(f'{current_directory}/Dsim_table.csv')
dsim_table_df = dsim_table_df[['NCBI_GeneID', 'Symbol', 'Transcript_Genomic_Accession', 'Transcript_Genomic_Start', 'Transcript_Genomic_Stop', 'Orientation']]

dsim_table_df['first_gene_location'] = dsim_table_df.apply(lambda row: row['Transcript_Genomic_Start'] if row['Orientation'] == 'plus' else row['Transcript_Genomic_Stop'], axis=1)
print(dsim_table_df)