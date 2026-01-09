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

############################################################################################
# Provisorische Timer-Klasse (am Ende entfernen!)
############################################################################################

import time
class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""
class Timer:
    def __init__(self):
        self._start_time = None
        self._timing_name : str

    def start(self, name="Elapsed time"):
        """Start a new timer"""
        if self._start_time is not None:
            raise TimerError(f"Timer is running. Use .stop() to stop it")
        self._start_time = time.perf_counter()
        assert type(name) == str
        self._timing_name = name

    def stop(self):
        """Stop the timer, and report the elapsed time"""
        if self._start_time is None:
            raise TimerError(f"Timer is not running. Use .start() to start it")
        elapsed_time = time.perf_counter() - self._start_time
        self._start_time = None
        print(f"{self._timing_name}: {elapsed_time:0.4f} seconds")

timer = Timer()

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

timer.start("tau_distance_set")
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
timer.stop()

def epsDensityTest(component: set[tuple[int,...]], step_count : int) -> bool:

    #TODO Was ist schneller?
    return True in [density_dict[lattice_point] > step_count + eps_bar for lattice_point in component]

    # for lattice_point in component:
    #     if density_dict[lattice_point] > step_count: + eps_bar:
    #         return True
    # return False


# Bestimmen der tau-Zusammenhangskomponenten von surviving_lattice_points durch Aufstellen eines Nähe-Graphen und Tiefensuche
def connectedComponents(step_count : int) -> list[set[tuple[int,...]]]:
    
    print("")
    surviving_lattice_points = survivingLatticePoints(step_count)

    print(surviving_lattice_points)

    # lattice_large_box_dict = {point : tuple(map(lambda x: x // int(tau_factor), point)) for point in surviving_lattice_points}
    large_boxes = {tuple(map(lambda x: x // int(tau_factor), point)) for point in surviving_lattice_points}
    
    # langsam aber allgemeineres Konzept:
    # def singletonPreimageDict(d : dict) -> dict:
    #   return {target_val : {key for key, val in d.items() if val == target_val} for target_val in d.values()}
    # large_box_lattice_dict = singletonPreimageDict(lattice_large_box_dict)

    # schneller:
    timer.start("inverse_dict")
    large_box_lattice_dict = {box : set(cartesian_product([list(range(int(tau_factor)*comp,int(tau_factor)*(comp+1))) for comp in box])) & surviving_lattice_points 
                              for box in large_boxes}
    timer.stop()

    tau_connection_graph = {}
    coord_shift = lambda coord : coord-1 if coord > 0 else -coord-1 if coord < 0 else 0
    sub = lambda x, y : x - y
    # large_boxes = large_box_lattice_dict.keys()
    relative_neighbor_boxes = cartesian_potentiation(list(range(-1,1+1)), dim)
    timer.start("box-loop")
    # t1 = 0
    # t2 = 0
    # print(len(large_boxes))
    for box in large_boxes:
        # t2_start = time.perf_counter()
        neighboring_boxes = {tuple(map(sum, zip(box, diff))) for diff in relative_neighbor_boxes} & large_boxes
        # t2_stop = time.perf_counter()
        # t2 += t2_stop- t2_start
        # t1_start = time.perf_counter()
        neighboring_lattice_points = set.union(*[large_box_lattice_dict[neighboring_box] for neighboring_box in neighboring_boxes])
        # t1_stop = time.perf_counter()
        # t1 += t1_stop - t1_start
        # print(len(neighboring_lattice_points))
        for new_point in large_box_lattice_dict[box]:
            connected_points = {point for point in neighboring_lattice_points - {new_point}
                                if max([coord_shift(coord) for coord in map(sub, point, new_point)]) < tau_factor}
            tau_connection_graph[new_point] = connected_points
    timer.stop()
    # print(f"neigboring boxes: {t2:0.4f}")
    # print(f"neighboring_lattice_points: {t1:0.4f}")


    # tau_connection_graph2 : dict[tuple[int,...],set[tuple[int,...]]] = {}

    # distance_calculations = 0
    # # distance_calc_time = 0
    # # set_creation_time = 0
    # # graph_add_time = 0

    # coord_shift = lambda coord : coord-1 if coord > 0 else -coord-1 if coord < 0 else 0
    # sub = lambda x, y : x - y

    # timer.start("alt")
    # for new_point in surviving_lattice_points:

    #     # set_creation_time_start = time.perf_counter()

    #     old_points = list(tau_connection_graph2.keys())
    #     # rel_points = [tuple(map(lambda x, y : x - y, point, new_point)) for point in old_points]

    #     # distance_calc_time_start = time.perf_counter()
    #     connected_points = [point for point in old_points if max([coord_shift(coord) for coord in map(sub, point, new_point)]) < tau_factor]
    #     # distance_calc_time_stop = time.perf_counter()
    #     # distance_calc_time += distance_calc_time_stop - distance_calc_time_start


    #     tau_connection_graph2[new_point] = set(connected_points)

    #     # set_creation_time_stop = time.perf_counter()
    #     # set_creation_time += set_creation_time_stop  - set_creation_time_start

    #     distance_calculations += len(old_points)

    #     # graph_add_time_start = time.perf_counter()

    #     for point in connected_points:
    #         tau_connection_graph2[point].add(new_point)

    #     # graph_add_time_stop = time.perf_counter()
    #     # graph_add_time += graph_add_time_stop - graph_add_time_start

    # # timer.stop()
    # # print(f"with {distance_calculations} distance calculations from {len(surviving_lattice_points)} lattice points")
    # # print(f"set creation time: {set_creation_time:0.4f}")
    # # print(f"Of That: distance calculation time:{distance_calc_time:0.4f}\n")
    # # print(f"time to add edges to the graph: {graph_add_time:0.4f}\n")
    # timer.stop()

    timer.start("finding components with dfs-alghorithm")

    # Bestimmen der Komponenten des Graphen
    component_list : list[set[tuple[int,...]]] = [] 
    visited_vertices = set()
    for vertex in surviving_lattice_points:
        # Depth-first-search im Graphen tau_connection_graph ausgehend vom Knoten vertex (siehe Wikipedia-Pseudocode)
        if vertex not in visited_vertices:
            component = set()
            stack = []
            stack.append(vertex)
            while bool(stack):
                queued_vertex = stack.pop()
                if queued_vertex not in component:
                    component.add(queued_vertex)
                    stack.extend(tau_connection_graph[queued_vertex])
            visited_vertices.update(component)
            component_list.append(component)
    timer.stop()

    timer.start("packing together the data")
    surviving_list : list[bool] = [epsDensityTest(component,step_count) 
                                   for component in component_list]
    surviving_components = [component for component, survived in zip(component_list,surviving_list) if survived]
    dead_components : set[tuple[int,...]] = set.union(*[component for component, survived 
                                                        in zip(component_list,surviving_list) 
                                                        if not survived]) if False in surviving_list else set()
    # 0-te Komponente beinhaltet Kästchen, die zu keiner überlebenden Komponente gehören
    surviving_components.insert(0,dead_components)
    timer.stop()

    return surviving_components




# woanders hin?
step_count_limit = 10000

step_count = 0
connected_component_count = 1
connected_component_list : list[set[tuple[int,...]]] = []
while connected_component_count == 1:
    if step_count > step_count_limit:
        raise RuntimeError(f"Clustering nach {step_count_limit} Versuchen abgebrochen")
    
    # timer.start("connectedComponents("+ str(step_count) + ")")
    connected_component_list = connectedComponents(step_count)
    connected_component_count = len(connected_component_list) - 1
    step_count += 1
    # timer.stop()

if connected_component_count == 0:
    clustered_data = [[1] + point for point in data_list]

else:
    clustered_data = []
    cluster_lattice_points = set.union(*connected_component_list)
    clustered_data = [[connected_component_list.index(component)] + data_point 
                        for component,data_point in zip(connected_component_list,data_list) 
                        if tuple(data_point) in component]
    

############################################################################################ 
# Visualisierung
############################################################################################ 

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
