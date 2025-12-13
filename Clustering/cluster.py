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

data_lattice_list : list[tuple] = [tuple([int((component + 1) // delta)
                                          for component in point])
                                    for point in data_list]

density_dict = {}

def getDensity(lattice_point: tuple[int]) -> int:
    return density_dict.get(lattice_point, 0)

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

def epsDensityTest(component: set[tuple], step_count: int) -> bool:
    for lattice_point in component:
        if density_dict[lattice_point] > step_count: #+ eps_bar: #ToDo : richtige Bedingung einfügen!!!
            return True
    return False

def connectedComponents(step_count: int) -> list[set[tuple]]:
    connected_component_list : list[set[tuple]] = []

    for lattice_point in survivingLatticePoints(step_count):

        # if step_count in [0]:
        #     print(f"\nmerge {lattice_point} into\n{connected_component_list}")

        neighbors = {tuple(map(sum, zip(lattice_point, diff))) for diff in tau_box}
        # ist das schneller:? = [tuple([sum(components) for components in zip(latticePoint, diff)])
        #                         for diff in tau_box]
        neighboring_component_indices = [i for i in range(len(connected_component_list)) if bool(neighbors & connected_component_list[i])]

        # if step_count in [0]:
        #     print(bool(neighboring_component_indices))

        if neighboring_component_indices: # aka wenn ... nichtleer

            # if step_count in [0]:
            #     print(neighboring_component_indices)

            min_index = neighboring_component_indices.pop(0)
            final_component = connected_component_list[min_index]
            # # richtig rum funktioniert nicht weil Rest der Liste sich verschiebt:
            # for i in neighboring_component_indices:
            #     final_component.update(connected_component_list.pop(i))
            for i in reversed(neighboring_component_indices):
                final_component.update(connected_component_list.pop(i))
            final_component.add(lattice_point)
            
            # # Alternative: (geht das deleten noch besser?)
            # merge_components = [connected_component_list[i] for i in neighboring_component_indices]
            # final_component.update(*merge_components)
            # for i in reversed(neighboring_component_indices):
            #     del connected_component_list[i]
        else:
            connected_component_list.append({lattice_point})

        
    # Wäre es besser ein dictionary statt der Liste "connected_component_list" zu benutzen?
    
    surviving_components = [component for component in connected_component_list if epsDensityTest(component,step_count)]
    return surviving_components

# woanders hin?
step_count_limit = 10000

step_count = 0
connected_component_count = 1
connected_component_list = []
while connected_component_count == 1:
    if step_count > step_count_limit:
        raise RuntimeError(f"Clustering nach {step_count_limit} Versuchen abgebrochen")
    
    connected_component_list = connectedComponents(step_count)
    connected_component_count = len(connected_component_list)
    step_count += 1

if connected_component_count == 0:
    clustered_data = [[1] + point for point in data_list]

print(step_count)
print(connected_component_list)

# else:
    # ToDo:
    # dictionary, dass Gitterpunkten ihre Zusammenhangskomponente zuweist
    # Datenpunkte über zugeordnete Gitterpunkte mit Zusammenhangskomponente markieren
    


####################################################################################################
# Visualisierung
####################################################################################################

# print(tau_factor)
# lst = list(tau_distance_set)
# x = [point[0] for point in lst]
# y = [point[1] for point in lst]
# plt.scatter(x, y)
# plt.show()

if dim == 2:
    for i in [0]:
        tuple_set = survivingLatticePoints(i)
        survived = [data_list[i] for i in range(len(data_list)) if data_lattice_list[i] in tuple_set]
        xdata = [item[0] for item in survived]
        ydata = [item[1] for item in survived]
        # print(xdata)
        plt.scatter(xdata, ydata, s=10, alpha=0.5)
        x = [(item[0]+.5)*delta - 1 for item in tuple_set]
        y = [(item[1]+.5)*delta - 1 for item in tuple_set]
        plt.scatter(x, y, s=50, alpha=0.3)
        # first_cluster = connectedComponents(0)
        # xfc = [(item[0]+.5)*delta - 1 for item in first_cluster[0]]
        # yfc = [(item[1]+.5)*delta for item in first_cluster[0]]
        # plt.scatter(xfc, yfc, s=20, alpha=0.7)
        a_cluster = list(connected_component_list[0])
        b_cluster = list(connected_component_list[1])
        xa = [(item[0]+.5)*delta - 1 for item in a_cluster]
        ya = [(item[1]+.5)*delta - 1 for item in a_cluster]
        xb = [(item[0]+.5)*delta - 1 for item in b_cluster]
        yb = [(item[1]+.5)*delta - 1 for item in b_cluster]
        plt.scatter(xa, ya, s=40, alpha=0.7)
        plt.scatter(xb, yb, s=40, alpha=0.7)

        plt.show()


####################################################################################################
# Schreiben der Daten in .csv/.png-Dateien
####################################################################################################

# result_data_frame = pd.DataFrame(clustered_data(data_list))
# result_data_frame.to_csv(result_path, index=False)



# ToDo (global):
# markieren der Daten mit Cluster-Indices durch else-Fall vervollständigen
# epsDensityTest richtig definieren
# log-Datei Daten implementieren
# richtige Sachen in die richtigen Dateien schreiben
# finale Datenvisualisierung für 2D
# ...
