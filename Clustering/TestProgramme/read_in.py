import pandas as pd

df = pd.read_csv("Clustering/cluster-data/bananas-1-4d.csv")
data_list = df.to_numpy().tolist()

print(data_list)