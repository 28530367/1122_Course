import pandas as pd

# 创建一个示例DataFrame
data = {'column1': ['A', 'A', 'B', 'B', 'C', 'C', 'D', 'F', 'F'],
        'column2': [1, 2, 1, 2, 1, 2, 5, 2 ,3]}
df = pd.DataFrame(data)

# 根据 'column1' 进行分组，并筛选出每个分组中大小大于1的行
result = df.groupby('column1').filter(lambda x: len(x) > 1)

print(result)