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
    readData_df = readCSV('humanmine_GO-terms_AllData.csv')
    
    column_rename_dict = {'Gene.primaryIdentifier': 'Gene Primary Identifier',
                          'Gene.symbol': 'Gene Symbol',
                          'Gene.name': 'Gene Name',
                          'Gene.goAnnotation.ontologyTerm.identifier': 'Ontology Term Identifier',
                          'Gene.goAnnotation.ontologyTerm.name': 'Ontology Term Name',
                          'Gene.organism.shortName': 'Organism Short Name',
                          'Gene.goAnnotation.ontologyTerm.parents.name': 'Parents Name'
                          }
    
    data_df = readData_df.rename(columns=column_rename_dict)

    drop_column_list = ['Gene Primary Identifier', 'Gene Name', 'Ontology Term Identifier', 'Ontology Term Name', 'Organism Short Name',]
    data_df = data_df.drop(columns=drop_column_list)
    data_df = data_df.drop_duplicates()

    # groupby
    quantity_df = data_df.groupby('Parents Name').size().reset_index(name='Quantity')
    arrangement_df = data_df.groupby('Parents Name').agg({'Gene Symbol': ', '.join}).reset_index()
    
    # merge
    arrangement_df = arrangement_df.merge(quantity_df, on='Parents Name')

    target_column_list = ['Parents Name', 'Quantity', 'Gene Symbol']
    arrangement_df = arrangement_df[target_column_list]

    # 輸出csv
    outputCSV(arrangement_df, 'humanmine_GO-terms_arrangement.csv')
