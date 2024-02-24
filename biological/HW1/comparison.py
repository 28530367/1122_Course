import pandas as pd
import os

# 取得當前工作目錄
current_directory = os.getcwd()

# 讀取csv
df1 = pd.read_csv(f'{current_directory}/output/Dsim_check_transcript.csv')
df2 = pd.read_csv(f'{current_directory}/outputAns/Dsim_check_transcript.csv')

df2 = df2[df2['score2'] != '-']
df1 = df1[df1['score2'].notna()]

df1.drop(columns=['evalue1', 'evalue2'], inplace=True)
df2.drop(columns=['evalue1', 'evalue2'], inplace=True)

convert_dict = {'score1': float, 'qcover1': float, 'identity1': float, 'score2': float, 'qcover2': float, 'identity2': float}
df1 = df1.astype(convert_dict)
df2 = df2.astype(convert_dict)

# 找出不同行
different_rows = pd.concat([df1, df2]).drop_duplicates(keep=False)

print("\n不同行:")
print(different_rows)

