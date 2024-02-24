import pandas as pd
import os

# 取得當前.py路徑
script_path = os.path.dirname(os.path.realpath(__file__))

# 輸入
min_score = 700

# 讀檔
protein_map_file_path = os.path.join(script_path, f'output_csv/protein_map_{min_score}.csv')

df_protein_map = pd.read_csv(protein_map_file_path)

protein_info_file_path = os.path.join(script_path, f'../Download/9606.protein.info.v12.0.txt')

df_protein_info = pd.read_csv(protein_info_file_path, sep='\t')

# 將含有NaN的row刪除
df_protein_map = df_protein_map.dropna()

# 將資料轉成原本分開的樣子
df_split_protein_map = df_protein_map.assign(protein1=df_protein_map['protein1'].str.split(',')).explode('protein1')
df_split_protein_map = df_split_protein_map.reset_index(drop=True)

# protein1 id to name
df_split_protein_map['protein1'] = df_split_protein_map['protein1'].map(df_protein_info.set_index('#string_protein_id')['preferred_name'])

# 相同protein合併
df_combined_protein_map_name = df_split_protein_map.groupby('protein2', as_index=False).agg({'protein1': lambda x: ','.join(x)})

# protein2 id to name
df_combined_protein_map_name['protein2'] = df_combined_protein_map_name['protein2'].map(df_protein_info.set_index('#string_protein_id')['preferred_name'])

# 計算數量
df_combined_protein_map_name['count'] = df_split_protein_map.groupby('protein2', as_index=False).agg({'protein1': 'count'})['protein1']

# 合併兩datafeame
df_combined_protein_map_name = pd.merge(df_combined_protein_map_name, df_protein_info, left_on='protein2', right_on='preferred_name', how='outer')

# 將count為NaN的改為0
df_combined_protein_map_name['count'].fillna(0, inplace=True)

# 移除protein2
df_combined_protein_map_name = df_combined_protein_map_name.drop('protein2', axis=1)

# 將preferred_name改名為protein2
df_combined_protein_map_name = df_combined_protein_map_name.rename(columns={'preferred_name': 'protein2'})

# 輸出csv
desired_columns_order = ['protein2', 'count', 'protein1']

output_csv_path = os.path.join(script_path, f'output_csv/protein_map_name_{min_score}.csv')

df_combined_protein_map_name.to_csv(output_csv_path, index=False, columns=desired_columns_order)