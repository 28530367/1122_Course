import pandas as pd

# 輸入
min_score = 700

# 讀取csv
df1 = pd.read_csv(f'/home/huangshouwei/1122_Course/HW3/ans_csv/ans_protein_map_{min_score}.csv')
df2 = pd.read_csv(f'/home/huangshouwei/1122_Course/HW3/output_csv/protein_map_{min_score}.csv')

# 找出共同行
common_rows = pd.merge(df1, df2, how='inner')

# 找出不同行
different_rows = pd.concat([df1, df2]).drop_duplicates(keep=False)

# print结果
print("共同行:")
print(common_rows)

print("\n不同行:")
print(different_rows)
