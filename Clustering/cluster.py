from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import pandas as pd
import matplotlib.pyplot as plt

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("datasetname", help="tbd")

# funktioniert iwi nur mit "--" statt "<" ">" :/
parser.add_argument("-d", "--delta", type=float, default=0.05, help="Parameter, that determines the fineness of the underlying lattice")
parser.add_argument("-e", "--eps-factor", type=float, default=3.0, help="tbd")
parser.add_argument("-t", "--tau-factor", type=float, default=2.000001, help="tbd")

args = vars(parser.parse_args())

datasetname = args["datasetname"]
delta = args["delta"]
eps_factor = args["eps_factor"]
tau_factor = args["tau_factor"]

data_path = "cluster-data/" + datasetname + ".csv"
result_path = "cluster-results/team-12-" + datasetname + ".result.csv"
log_path = "cluster-results/team-12-" + datasetname + ".log"
result_log_path = "cluster-results/team-12-" + datasetname + ".result.log"

df = pd.read_csv(data_path)
data_list : list[list[float]] = df.to_numpy().tolist()

####################################################################################################
# Verarbeitung der Daten (ohne Verwendung der Module!)
####################################################################################################

# def cluster_analysis(data_list : list[list[float]]) -> list[list[float]]:
#     return [[]]
# clustered_data = pd.DataFrame(cluster_analysis(data_list))
# clustered_data.to_csv(result_path, index=False)

data_lattice_list : list[tuple] = [tuple([int((component + 1) // delta)
                                          for component in point])
                                    for point in data_list]

density_dict = {}

# def dictMaybeInside(dic,key):
#     if key in dic:
#         return dic[key]
#     else:
#         return 0

def getDensity(lattice_point: tuple[int]) -> int:
    if lattice_point in density_dict:
        return density_dict[lattice_point]
    else:
        return 0

for lattice_point in data_lattice_list:
    density_dict[lattice_point] = getDensity(lattice_point) + 1

#rho_bar ist einfach nur step_count + 1, rho ist nur rho_bar mit Vorfaktor
def survivingLatticePoints(step_count: int) -> set[tuple]:
    return {lattice_point for lattice_point in data_lattice_list 
                if density_dict[lattice_point] > step_count}

# Am besten oben definieren:
dim = len(data_list[0])
# n_data = len(data_list)
# density_factor = n_data * (2 ** dim) * (delta ** dim)
# tau = tau_factor * delta

def cartesian_product(X : list[list]):
    dim = len(X)
    prod = [()]
    for l in range(dim):
        prod = [ tup + (item,) for tup in prod for item in X[l] ]
    return prod

def cartesian_potentiation(lst : list, dim : int):
    return cartesian_product([lst for _ in range(dim)])

# Intervall [-(int(tau_factor)+1), ... , int(tau_factor)+1] Reicht das ??
eff_tau_interval_list = list(range(-int(tau_factor)-1,int(tau_factor)+2))

# dim-dimensionale Box, die alle Gitterpunkte enthält, für die eff_origin_distance<=tau_factor gelten könnte
tau_box = cartesian_potentiation(eff_tau_interval_list, dim)

def eff_origin_distance(lattice_point: tuple):
    eff_lattice_point_list = [coord-1 if coord > 0 else coord+1 if coord < 0 else coord for coord in lattice_point]
    eff_coord_squares = [coord**2 for coord in eff_lattice_point_list]
    return sum(eff_coord_squares)**.5

# dist(box1, box2) <= tau <=> eff_origin_distance() der Differenz der assoziierten lattice points <= tau_factor (= tau/delta) 
tau_distance_set = {lattice_point for lattice_point in tau_box if eff_origin_distance(lattice_point) <= tau_factor}



####################################################################################################
# Visualisierung
####################################################################################################

# print(tau_factor)
# lst = list(tau_distance_set)
# x = [point[0] for point in lst]
# y = [point[1] for point in lst]
# plt.scatter(x, y)
# plt.show()

# if dim == 2:
#     for i in [5]:
#         tuple_set = survivingLatticePoints(i)
#         survived = [data_list[i] for i in range(len(data_list)) if data_lattice_list[i] in tuple_set]
#         xdata = [item[0] for item in survived]
#         ydata = [item[1] for item in survived]
#         # print(xdata)
#         plt.scatter(xdata, ydata, s=10, alpha=0.5)
#         x = [(item[0]+.5)*delta - 1 for item in tuple_set]
#         y = [(item[1]+.5)*delta - 1 for item in tuple_set]
#         plt.scatter(x, y, s=50, alpha=0.3)
#         plt.show()
