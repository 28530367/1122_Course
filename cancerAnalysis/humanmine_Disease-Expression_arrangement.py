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
    readData_df = readCSV('humanmine_Disease-Expression_AllData.csv')
    
    column_rename_dict = {'AtlasExpression.condition': 'Atlas Expression Condition',
                          'AtlasExpression.gene.primaryIdentifier': 'Gene Primary Identifier',
                          'AtlasExpression.gene.symbol': 'Gene Symbol',
                          'AtlasExpression.gene.name': 'Gene Name',
                          'AtlasExpression.expression': 'Atlas Expression Expression',
                          'AtlasExpression.pValue': 'Atlas Expression P Value',
                          'AtlasExpression.tStatistic': 'Atlas Expression T Statistic',
                          'AtlasExpression.dataSets.name': 'Data Sets Name'
                          }
    
    data_df = readData_df.rename(columns=column_rename_dict)

    drop_column_list = ['Gene Primary Identifier', 'Gene Name', 'Atlas Expression Expression', 'Atlas Expression P Value', 'Atlas Expression T Statistic', 'Data Sets Name']
    data_df = data_df.drop(columns=drop_column_list)
    
    print('before drop_duplicates:', len(data_df))
    data_df = data_df.drop_duplicates()
    print('after drop_duplicates:', len(data_df))

    # groupby
    quantity_df = data_df.groupby('Atlas Expression Condition').size().reset_index(name='Quantity')
    arrangement_df = data_df.groupby('Atlas Expression Condition').agg({'Gene Symbol': ', '.join}).reset_index()
    
    # merge
    arrangement_df = arrangement_df.merge(quantity_df, on='Atlas Expression Condition')

    target_column_list = ['Atlas Expression Condition', 'Quantity', 'Gene Symbol']
    arrangement_df = arrangement_df[target_column_list]

    # 輸出csv
    outputCSV(arrangement_df, 'humanmine_Disease-Expression_arrangement.csv')