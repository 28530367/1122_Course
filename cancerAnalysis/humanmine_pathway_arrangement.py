import pandas as pd
import os

# 取得當前工作目錄
current_directory = os.path.dirname(os.path.abspath(__file__))

def readCSV(file_type):
    filePath = f"{current_directory}/humanmine_{file_type}_AllData.csv"
    fileData_df = pd.read_csv(filePath)

    return fileData_df

def outputCSV(df, file_name):
    filePath = f"{current_directory}/{file_name}"
    df.to_csv(filePath, index=False)
    print('CSV SUCCESS!!')

if __name__ == "__main__":
    pathwayData_df = readCSV('pathway')

    pathwayData_df[['Gene.proteins.primaryIdentifier', 'Species']] = pathwayData_df['Gene.proteins.primaryIdentifier'].str.split('_', expand=True)

    # 分離 HUMAN 與其他物種
    human_df = pathwayData_df[pathwayData_df['Species'].str.contains('HUMAN')]
    others_df = pathwayData_df[~pathwayData_df['Species'].str.contains('HUMAN')]

    # groupby
    quantity_df = human_df.groupby('Gene.proteins.pathways.name').size().reset_index(name='quantity')
    humanArrangement_df = human_df.groupby('Gene.proteins.pathways.name').agg({'Gene.primaryIdentifier': ', '.join, 'Gene.symbol': ', '.join, 'Gene.proteins.primaryAccession': ', '.join, 'Gene.proteins.primaryIdentifier': ', '.join, 'Gene.proteins.pathways.identifier': ', '.join}).reset_index()
    
    # merge
    humanArrangement_df = humanArrangement_df.merge(quantity_df, on='Gene.proteins.pathways.name')

    # 輸出csv
    outputCSV(humanArrangement_df, 'human_pathwayArrangement.csv')
