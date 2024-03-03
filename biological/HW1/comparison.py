import pandas as pd
import os

# 取得當前工作目錄
current_directory = os.path.dirname(os.path.abspath(__file__))
print(current_directory)

def bestHit_comparsion():
    # 'n' or 'p'
    target = 'p'

    left = 'Dmel'
    right = 'Dsim'
    # left = 'Dsim'
    # right = 'Dmel'

    # 讀取csv
    df1 = pd.read_csv(f'{current_directory}/bestHit/blast{target}_{left}to{right}_besthit.csv')
    df2 = pd.read_csv(f'{current_directory}/bestHitAns/blast{target}_{left}to{right}_besthit.csv')

    # 找出不同行
    different_rows = pd.concat([df1, df2]).drop_duplicates(keep=False)

    print("\n不同行:")
    print(different_rows)


    # # 使用 merge 方法，指定 indicator=True
    # merged_df = pd.merge(df1, df2, how='outer', indicator=True)

    # # 找出不同的行
    # different_rows_df1 = merged_df[merged_df['_merge'] == 'left_only'].drop('_merge', axis=1)
    # different_rows_df2 = merged_df[merged_df['_merge'] == 'right_only'].drop('_merge', axis=1)

    # print("\n不同行 Mine:")
    # print(different_rows_df1)

    # print("\n不同行 Ans:")
    # print(different_rows_df2)

def dsimCheck_comarsion():
    # 'protein' 'transcript'
    target = 'transcript'

    if target == 'transcript':
        # 讀取csv for transcript
        df1 = pd.read_csv(f'{current_directory}/output/Dsim_check_{target}.csv', dtype=str)
        df2 = pd.read_csv(f'{current_directory}/outputAns/Dsim_check_{target}.csv', dtype=str)
    elif target == 'protein':
        # 讀取csv for protein
        df1 = pd.read_csv(f'{current_directory}/output/Dsim_check_{target}.csv')
        df2 = pd.read_csv(f'{current_directory}/outputAns/Dsim_check_{target}.csv')

        df1['qcover2'] = list(map(lambda x: float(x) if x != '-' else x, df1['qcover2']))
        df2['qcover2'] = list(map(lambda x: float(x) if x != '-' else x, df2['qcover2']))

    # 找出不同行
    different_rows = pd.concat([df1, df2]).drop_duplicates(keep=False)

    print("\n不同行:")
    print(different_rows)

    # # 使用 merge 方法，指定 indicator=True
    # merged_df = pd.merge(df1, df2, how='outer', indicator=True)

    # # 找出不同的行
    # different_rows_df1 = merged_df[merged_df['_merge'] == 'left_only'].drop('_merge', axis=1)
    # different_rows_df2 = merged_df[merged_df['_merge'] == 'right_only'].drop('_merge', axis=1)

    # print("\n不同行 Mine:")
    # print(different_rows_df1)

    # print("\n不同行 Ans:")
    # print(different_rows_df2)

def dsimClassification():
    # 讀取csv
    df1 = pd.read_csv(f'{current_directory}/output/Dsim_classification.csv')
    df2 = pd.read_csv(f'{current_directory}/outputAns/Dsim_classification.csv')

    # # 找出不同行
    # different_rows = pd.concat([df1, df2]).drop_duplicates(keep=False)

    # print("\n不同行:")
    # print(different_rows)

    # 使用 merge 方法，指定 indicator=True
    merged_df = pd.merge(df1, df2, how='outer', indicator=True)

    # 找出不同的行
    different_rows_df1 = merged_df[merged_df['_merge'] == 'left_only'].drop('_merge', axis=1)
    different_rows_df2 = merged_df[merged_df['_merge'] == 'right_only'].drop('_merge', axis=1)

    print("\n不同行 Mine:")
    print(different_rows_df1)

    print("\n不同行 Ans:")
    print(different_rows_df2)

def dsimOrtholog():
    # 'none' 'unique' 'multiple'
    type = 'multiple'

    # 讀取csv
    df1 = pd.read_csv(f'{current_directory}/output/Dsim_ortholog_{type}.csv')
    df2 = pd.read_csv(f'{current_directory}/outputAns/Dsim_ortholog_{type}.csv')

    # # 找出不同行
    # different_rows = pd.concat([df1, df2]).drop_duplicates(keep=False)

    # print("\n不同行:")
    # print(different_rows)

    # 使用 merge 方法，指定 indicator=True
    merged_df = pd.merge(df1, df2, how='outer', indicator=True)

    # 找出不同的行
    different_rows_df1 = merged_df[merged_df['_merge'] == 'left_only'].drop('_merge', axis=1)
    different_rows_df2 = merged_df[merged_df['_merge'] == 'right_only'].drop('_merge', axis=1)

    print("\n不同行 Mine:")
    print(different_rows_df1)

    print("\n不同行 Ans:")
    print(different_rows_df2)

if __name__ == "__main__":
    dsimOrtholog()