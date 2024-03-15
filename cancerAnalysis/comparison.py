import pandas as pd
import os

# 取得當前工作目錄
current_directory = os.path.dirname(os.path.abspath(__file__))
print(current_directory)


# 讀取csv
df1 = pd.read_csv(f'{current_directory}/humanmine_Human-Disease_arrangement.csv')
df2 = pd.read_csv(f'{current_directory}/humanmine_Disease_arrangement.csv')

# 找出不同行
different_rows = pd.concat([df1, df2]).drop_duplicates(keep=False)

print("\n不同行:")
print(different_rows)