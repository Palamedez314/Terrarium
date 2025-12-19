from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import pandas as pd
import matplotlib.pyplot as plt

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("datasetname", help="tbd")

# gleiches tau_distance_set aber verschiedenes clustering???????

# funktioniert iwi nur mit "--" statt "<" ">" :/
parser.add_argument("-d", "--delta", type=float, default=0.05, help="Parameter, that determines the fineness of the underlying lattice")
parser.add_argument("-e", "--eps-factor", type=float, default=3.0, help="tbd")
parser.add_argument("-t", "--tau-factor", type=float, default=2.0, help="tbd")
parser.add_argument("-n", "--norm", type=float, default="inf", help="tbd")

args = vars(parser.parse_args())

datasetname = args["datasetname"]
delta = args["delta"]
eps_factor = args["eps_factor"]
tau_factor = args["tau_factor"]
norm = args["norm"]

if norm < 1:
    raise ValueError("norm must be a positive float >= 1")    

data_path = "cluster-data/" + datasetname + ".csv"
result_path = "cluster-results/team-12-" + datasetname + ".result.csv"
log_path = "cluster-results/team-12-" + datasetname + ".log"
result_log_path = "cluster-results/team-12-" + datasetname + ".result.log"

df = pd.read_csv(data_path)
data_list : list[list[float]] = df.to_numpy().tolist()

####################################################################################################
# Verarbeitung der Daten (ohne Verwendung der Module!)
####################################################################################################

n_data = len(data_list)
dim = len(data_list[0])

#Liste an Kugeln in Reihenfolge der Datenpunkte
data_lattice_list : list[tuple[int,...]] = [tuple([int((component + 1) // delta)
                                                   for component in point])
                                            for point in data_list]

density_dict = {}

for lattice_point in data_lattice_list:
    density_dict[lattice_point] = density_dict.get(lattice_point, 0) + 1

h_max_bar = max(density_dict.values())
eps_bar = eps_factor * ((h_max_bar) ** .5)


#rho_bar ist einfach nur step_count + 1, rho ist nur rho_bar mit Vorfaktor
def survivingLatticePoints(step_count: int) -> set[tuple[int,...]]:
    return {lattice_point for lattice_point in data_lattice_list 
                if density_dict[lattice_point] > step_count}

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
tau_box : list[tuple[int,...]] = cartesian_potentiation(eff_tau_interval_list, dim)

def eff_origin_distance(lattice_point: tuple[int,...]):
    # TODO: wie schlimm ist die Code-Dopplung?
    if norm == float("inf"):
        eff_abs_lattice_point_list = [coord-1 if coord > 0 else -coord-1 if coord < 0 else 0 for coord in lattice_point]
        return max(eff_abs_lattice_point_list)
    elif norm == 1:
        eff_abs_lattice_point_list = [-coord+1 if coord > 0 else -coord-1 if coord < 0 else 0 for coord in lattice_point]
        return sum(eff_abs_lattice_point_list)
    else:
        eff_lattice_point_list = [coord-1 if coord > 0 else coord+1 if coord < 0 else 0 for coord in lattice_point]
        eff_coord_squares = [coord**norm for coord in eff_lattice_point_list]
        return sum(eff_coord_squares)**(1/norm)

# dist(box1, box2) < tau <=> eff_origin_distance() der Differenz der assoziierten lattice points < tau_factor (= tau/delta) 
tau_distance_set = {lattice_point for lattice_point in tau_box if eff_origin_distance(lattice_point) < tau_factor}

def epsDensityTest(component: set[tuple[int,...]], step_count : int) -> bool:

    #TODO Was ist schneller?
    return True in [density_dict[lattice_point] > step_count + eps_bar for lattice_point in component]

    # for lattice_point in component:
    #     if density_dict[lattice_point] > step_count: + eps_bar:
    #         return True
    # return False

def connectedComponents(step_count: int) -> list[set[tuple[int,...]]]:
    connected_component_list : list[set[tuple[int,...]]] = []

    for lattice_point in survivingLatticePoints(step_count):
        
        # if step_count in [0]:
        #     print(f"\nmerge {lattice_point} into\n{connected_component_list}")

        neighbors = {tuple(map(sum, zip(lattice_point, diff))) for diff in tau_distance_set}
        # ist das schneller:? = [tuple([sum(components) for components in zip(latticePoint, diff)])
        #                         for diff in tau_box]
        # TODO Schaut sich der Schnitt alle Elemente an -> ja! geht das schneller?
        # TODO Schnitt besser, nicht die ganze Menge berechnen -> educated guess welche Methode besser ist
        neighboring_component_indices = [i for i in range(len(connected_component_list)) if 
                                         bool(neighbors & connected_component_list[i])]

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
    surviving_list : list[bool] = [epsDensityTest(component,step_count) 
                                   for component in connected_component_list]
    surviving_components = [component for component, survived in zip(connected_component_list,surviving_list) if survived]
    dead_components : set[tuple[int,...]] = set.union(*[component for component, survived 
                                                        in zip(connected_component_list,surviving_list) 
                                                        if not survived]) if False in surviving_list else set()
    
    surviving_components.insert(0,dead_components)

    return surviving_components

# woanders hin?
step_count_limit = 10000

step_count = 0
connected_component_count = 1
connected_component_list : list[set[tuple[int,...]]] = []
while connected_component_count == 1:
    if step_count > step_count_limit:
        raise RuntimeError(f"Clustering nach {step_count_limit} Versuchen abgebrochen")
    
    connected_component_list = connectedComponents(step_count)
    connected_component_count = len(connected_component_list) - 1
    step_count += 1

if connected_component_count == 0:
    clustered_data = [[1] + point for point in data_list]

else:
    clustered_data = []
    cluster_lattice_points = set.union(*connected_component_list)
    clustered_data = [[connected_component_list.index(component)] + data_point 
                        for component,data_point in zip(connected_component_list,data_list) 
                        if tuple(data_point) in component]
    
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
    for i in [0,step_count-1]:
        tuple_set = survivingLatticePoints(i)
        # print(len(tuple_set))
        survived = [data_list[i] for i in range(n_data) if data_lattice_list[i] in tuple_set]
        xdata = [item[0] for item in survived]
        ydata = [item[1] for item in survived]
        # print(xdata)
        plt.scatter(xdata, ydata, s=10, alpha=0.5)
        x = [(item[0]+.5)*delta - 1 for item in tuple_set]
        y = [(item[1]+.5)*delta - 1 for item in tuple_set]
        plt.scatter(x, y, s=50, alpha=0.3)
        # first_cluster = connectedComponents(i)
        # xfc = [(item[0]+.5)*delta - 1 for item in first_cluster[1]]
        # yfc = [(item[1]+.5)*delta -1  for item in first_cluster[1]]
        # plt.scatter(xfc, yfc, s=20, alpha=0.7)

        if i == step_count-1:
            for cluster_set in connected_component_list[1:]:
                cluster_list = list(cluster_set)
                cluster_x = [(item[0]+.5)*delta - 1 for item in cluster_list]
                cluster_y = [(item[1]+.5)*delta - 1 for item in cluster_list]
                plt.scatter(cluster_x,cluster_y, s=40, alpha=0.7)

        plt.show()


####################################################################################################
# Schreiben der Daten in .csv/.png-Dateien
####################################################################################################

# result_data_frame = pd.DataFrame(clustered_data(data_list))
# result_data_frame.to_csv(result_path, index=False)



# TODO (global):
# log-Datei Daten implementieren
# richtige Sachen in die richtigen Dateien schreiben
# finale Datenvisualisierung für 2D
# ...
