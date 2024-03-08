import pandas as pd
import os

# 取得當前工作目錄
current_directory = os.path.dirname(os.path.abspath(__file__))

def readCSV(file_name):
    filePath = f"{current_directory}/{file_name}"
    fileData_df = pd.read_csv(filePath)

    return fileData_df

def outputCSV(output_df, file_name):
    filePath = f"{current_directory}/{file_name}"
    output_df.to_csv(filePath, index=False)
    print('CSV SUCCESS!!')

if __name__ == "__main__":
    pathwayData_df = readCSV('humanmine_pathway_AllData.csv')
    
    column_rename_dict = {'Gene.primaryIdentifier': 'Gene Primary Identifier',
                          'Gene.symbol': 'Gene Symbol',
                          'Gene.proteins.primaryAccession': 'Proteins Primary Accession',
                          'Gene.proteins.primaryIdentifier': 'Proteins Primary Identifier',
                          'Gene.proteins.pathways.identifier': 'Pathways Identifier',
                          'Gene.proteins.pathways.name': 'Pathways Name',
                          'Gene.proteins.pathways.dataSets.name': 'Data Sets Name'
                          }
    
    pathwayData_df.rename(columns=column_rename_dict, inplace=True)

    pathwayData_df['Species'] = pathwayData_df['Proteins Primary Identifier'].str.split('_', expand=True)[1]

    ## 分離 HUMAN 與其他物種
    human_df = pathwayData_df[pathwayData_df['Species'].str.contains('HUMAN')]
    # others_df = pathwayData_df[~pathwayData_df['Species'].str.contains('HUMAN')]

    drop_column_list = ['Gene Primary Identifier', 'Proteins Primary Accession', 'Proteins Primary Identifier', 'Pathways Identifier', 'Data Sets Name', 'Species']
    human_df = human_df.drop(columns=drop_column_list)

    print('before drop_duplicates:', len(human_df))
    human_df = human_df.drop_duplicates()
    print('after drop_duplicates:', len(human_df))

    # groupby
    quantity_df = human_df.groupby('Pathways Name').size().reset_index(name='Quantity')
    humanArrangement_df = human_df.groupby('Pathways Name').agg({'Gene Symbol': ', '.join}).reset_index()
    
    # merge
    humanArrangement_df = humanArrangement_df.merge(quantity_df, on='Pathways Name')

    target_column_list = ['Pathways Name', 'Quantity', 'Gene Symbol']
    humanArrangement_df = humanArrangement_df[target_column_list]

    # 輸出csv
    outputCSV(humanArrangement_df, 'humanmine_pathway_arrangement.csv')
