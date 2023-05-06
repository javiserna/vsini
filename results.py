import pandas as pd
import numpy as np

df = pd.read_csv("Fourier.out", delimiter=' ', comment='!')
table=df.groupby('#file')[['vsini', 'vsini_err']].median()
df2=table
df2.to_csv('results.csv')
print(table)
