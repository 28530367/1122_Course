import pandas as pd
import os

current_directory = os.path.dirname(os.path.abspath(__file__))

def level(group):
    num_rows = group.shape[0]
    row_value = [-1] * num_rows
    for index, row in group.iterrows():
        for y_index, y in enumerate(row_value):
            if row['Transcript_Genomic_Start'] > y:
                group.loc[index, 'level'] = y_index + 1
                break
            else:
                pass
        row_value[y_index] = row['Transcript_Genomic_Stop']

    return group

dsim_table_df = pd.read_csv(f'{current_directory}/Dmel_table.csv')
# print(dsim_table_df)
dsim_table_df = dsim_table_df[['NCBI_GeneID', 'Chromosomes', 'Transcript_Genomic_Accession', 'Transcript_Genomic_Start', 'Transcript_Genomic_Stop', 'Orientation']]

dsim_sort_df = dsim_table_df.sort_values(by=['Transcript_Genomic_Accession', 'Transcript_Genomic_Start', 'Transcript_Genomic_Stop'], ascending=True)
# print(dsim_sort_df)
dsim_sort_df.to_csv(f'{current_directory}/Dmel_sort.csv', index=False)

level_df = dsim_sort_df.groupby('Transcript_Genomic_Accession').apply(level)

level_df = level_df.rename(columns={
    'NCBI_GeneID': 'Dmel',
    'Chromosomes': 'Chromosomes',
    'Transcript_Genomic_Accession': 'scaffold',
    'Transcript_Genomic_Start': 'start',
    'Transcript_Genomic_Stop': 'end',
    'Orientation': 'frame',
})
print(level_df)

level_df.to_csv(f'{current_directory}/Dmel_GeneOrder_level.csv', index=False)