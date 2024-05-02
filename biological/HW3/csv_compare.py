import pandas as pd
import os

current_directory = os.path.dirname(os.path.abspath(__file__))

# 讀取csv
df1 = pd.read_csv(f'{current_directory}/Dsim_GeneOrder_level.csv')
df2 = pd.read_csv(f'{current_directory}/Ans/Dsim_GeneOrder_level.csv')

# 找出不同行
different_rows = pd.concat([df1, df2]).drop_duplicates(keep=False)

print("\n不同行:")
print(different_rows)

