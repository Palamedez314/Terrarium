import pandas as pd

# import csv

# with open('Clustering/Daten/bananas-1-4d.csv', newline='') as csvfile:
#     reader = csv.reader(csvfile, delimiter=' ', quotechar='|', #quoting=csv.QUOTE_NONNUMERIC
#                         )

#     datapoints = list(reader)

# print(datapoints)

df = pd.read_csv('Clustering/Daten/bananas-1-4d.csv')
data_list = df.to_numpy().tolist()

print(data_list)

