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
    readData_df = readCSV('humanmine_Human-Disease_AllData.csv')
    
    column_rename_dict = {'Disease.primaryIdentifier': 'Disease Primary Identifier',
                          'Disease.name': 'Disease Name',
                          'Disease.genes.primaryIdentifier': 'Genes Primary Identifier',
                          'Disease.genes.symbol': 'Gene Symbol',
                          'Disease.genes.homologues.homologue.primaryIdentifier': 'Homologue Primary Identifier',
                          'Disease.genes.homologues.homologue.symbol': 'Homologue Symbol'
                          }
    
    data_df = readData_df.rename(columns=column_rename_dict)

    drop_column_list = ['Disease Primary Identifier', 'Genes Primary Identifier', 'Homologue Primary Identifier', 'Homologue Symbol']
    data_df = data_df.drop(columns=drop_column_list)
    
    print('before drop_duplicates:', len(data_df))
    data_df = data_df.drop_duplicates()
    print('after drop_duplicates:', len(data_df))

    # groupby
    quantity_df = data_df.groupby('Disease Name').size().reset_index(name='Quantity')
    arrangement_df = data_df.groupby('Disease Name').agg({'Gene Symbol': ', '.join}).reset_index()
    
    # merge
    arrangement_df = arrangement_df.merge(quantity_df, on='Disease Name')
 
    target_column_list = ['Disease Name', 'Quantity', 'Gene Symbol']
    arrangement_df = arrangement_df[target_column_list]

    # 輸出csv
    outputCSV(arrangement_df, 'humanmine_Human-Disease_arrangement.csv')