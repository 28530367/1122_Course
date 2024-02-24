import pandas as pd
import os

# 取得當前.py路徑
script_path = os.path.dirname(os.path.realpath(__file__))

# 讀檔
protein_links_file_path = os.path.join(script_path, f'../Download/9606.protein.links.v12.0.txt')

df_protein_links = pd.read_csv(protein_links_file_path, sep=' ')

protein_info_file_path = os.path.join(script_path, f'../Download/9606.protein.info.v12.0.txt')

df_protein_info = pd.read_csv(protein_info_file_path, sep='\t')

# 輸入
min_score = 700

# 篩選
df_pass_protein_links_id = df_protein_links[df_protein_links['combined_score'] >= min_score]

# 相同protein合併
df_combined_pass_protein_links_id = df_pass_protein_links_id.groupby('protein2', as_index=False).agg({'protein1': lambda x: ','.join(x)})

# 計算數量
df_combined_pass_protein_links_id['count'] = df_pass_protein_links_id.groupby('protein2', as_index=False).agg({'protein1': 'count'})['protein1']

# 合併兩datafeame
df_combined_pass_protein_links_id = pd.merge(df_combined_pass_protein_links_id, df_protein_info, left_on='protein2', right_on='#string_protein_id', how='outer')

# 將count為NaN的改為0
df_combined_pass_protein_links_id['count'].fillna(0, inplace=True)

# 移除protein2
df_combined_pass_protein_links_id = df_combined_pass_protein_links_id.drop('protein2', axis=1)

# 將#string_protein_id改名為protein2
df_combined_pass_protein_links_id = df_combined_pass_protein_links_id.rename(columns={'#string_protein_id': 'protein2'})

# 輸出csv
desired_columns_order = ['protein2', 'count', 'protein1']

output_csv_path = os.path.join(script_path, f'output_csv/protein_map_{min_score}.csv')

df_combined_pass_protein_links_id.to_csv(output_csv_path, index=False, columns=desired_columns_order)

