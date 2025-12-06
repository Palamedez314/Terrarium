from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import pandas as pd

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("datasetname", help="tbd")

# funktioniert iwi nur mit "--" statt "<" ">" :/
parser.add_argument("-d", "--delta", type=float, default=0.05, help="Parameter, that determines the fineness of the underlying lattice")
parser.add_argument("-e", "--eps-factor", type=float, default=3.0, help="tbd")
parser.add_argument("-t", "--tau-factor", type=float, default=2.000001, help="tbd")

args = vars(parser.parse_args())

datasetname = args["datasetname"]
delta = args["delta"]
eps = args["eps_factor"]
tau = args["tau_factor"]

dataPath = "cluster-data/" + datasetname + ".csv"
resultPath = "cluster-results/team-12-" + datasetname + ".result.csv"
logPath = "cluster-results/team-12-" + datasetname + ".log"
resultlogPath = "cluster-results/team-12-" + datasetname + ".result.log"

df = pd.read_csv(dataPath)
data_list = df.to_numpy().tolist()

# print(data_list)

def cluster_analysis(data_list : list[list[float]]) -> list[list[float]]:
    return [[]]
    # return [[1,3,4],[1,3,7],[2,8,9,],[3,2,4],[0,2,1]]

clustered_data = pd.DataFrame(cluster_analysis(data_list))

clustered_data.to_csv(resultPath, index=False)