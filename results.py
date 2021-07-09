import pandas as pd
import numpy as np

df = pd.read_csv("Fourier.out", delimiter=' ')
table=df.groupby('#file').median()
df2=table[['vsini','vsini_err']]
df2.to_csv('results.csv')
print(df2)
