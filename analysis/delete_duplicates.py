import pandas as pd

df = pd.read_csv('../results/result.csv')

df = df.drop_duplicates()

df.to_csv('../results/new_result.csv', index=False)