#!/usr/bin/env python

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("datasetname", help="Name of the dataset out of the folder 'cluster-data' (without .csv extension)")

# gleiches tau_distance_set aber verschiedenes clustering???????

# funktioniert iwi nur mit "--" statt "<" ">" :/ Ja Windows erlaubt keine <> in Dateinamen

def prase_args():
    parser = ArgumentParser(prog="cluster", description= "Team 12, Aufgabe 2, Cluster Algorithmus", formatter_class=ArgumentDefaultsHelpFormatter) 

    parser.add_argument("datasetname",help="Name of the dataset in 'cluster-data' folder (without .csv extension)")
    parser.add_argument("-d", "--delta", type=float, default=0.05, help="Parameter, that determines the fineness of the underlying lattice")
    parser.add_argument("-e", "--eps-factor", type=float, default=3.0, help="tbd")
    parser.add_argument("-t", "--tau-factor", type=float, default=2.0, help="tbd")
    parser.add_argument("-n", "--norm", type=float, default=float("inf"), help="tbd")
    return parser.parse_args()

############################################################################################
# Provisorische Timer-Klasse (am Ende entfernen!)
############################################################################################

import time
class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""
class Timer:
    def __init__(self):
        self._start_times : dict[str,float] = {}
        self._timing_descriptions : dict[str,str] = {}
        self._cumul_times : dict[str,float] = {}

    def start_variable(self, varname:str, description:str=""):
        """Start timer with name varname"""
        assert type(varname) == str
        if varname in self._start_times.keys():
            raise TimerError(f"Timer is running. Use .stop()/.stop_variable to stop it")
        self._start_times[varname] = time.perf_counter()
        assert type(description) == str
        if description == "":
            description = varname
        self._timing_descriptions[varname] = description

    def print_variable(self, varname, val, cumul:bool=False):
        assert type(varname) == str
        if cumul:
            print(f"{self._timing_descriptions[varname]}: {val:0.4f} seconds (cumul)")
        else:
            print(f"{self._timing_descriptions[varname]}: {val:0.4f} seconds")

    def stop_variable(self, varname:str, print_single:bool=True, print_cumul:bool=False):
        """Stop timer with name varname"""
        assert type(varname) == str
        if varname not in self._start_times.keys():
            raise TimerError(f"Timer is not running. Use .start()/.start_variable() to start it")
        elapsed_time = time.perf_counter() - self._start_times.pop(varname)
        if print_single:
            self.print_variable(varname, elapsed_time, cumul=False)
        if print_cumul:
            cumul_time = self._cumul_times.pop(varname, 0) + elapsed_time
            self.print_variable(varname, cumul_time, cumul=False)
        else:
            self._cumul_times.pop(varname, 0)

    def print_cumul_variable(self, varname:str):
        if varname not in self._cumul_times.keys():
            raise TimerError(f"Timer has no recorded time yet. Use .start()/.start_variable() and .pause()/.pause_variable() to record time values")
        self.print_variable(varname, self._cumul_times[varname], cumul=True)

    def pause_variable(self, varname:str, print_single:bool=False, print_cumul:bool=False) :
        assert type(varname) == str
        if varname not in self._start_times.keys():
            raise TimerError(f"Timer is not running. Use .start()/.start_variable() to start it")
        elapsed_time = time.perf_counter() - self._start_times.pop(varname)
        if varname not in self._cumul_times.keys():
            self._cumul_times[varname] = elapsed_time
        else:
            self._cumul_times[varname] += elapsed_time
        if print_single:
            self.print_variable(varname, elapsed_time, cumul=False)
        if print_cumul:
            self.print_cumul_variable(varname)

    def start(self, description="Elapsed time"):
        """Start a new timer"""
        self.start_variable("standard", description)

    def stop(self, print_single:bool=True, print_cumul:bool=False):
        """Stop the timer, and report the elapsed time"""
        self.stop_variable("standard", print_single=print_single, print_cumul=print_cumul)

    def pause(self, print_single:bool=False, print_cumul:bool=False):
        self.pause_variable("standard", print_single=print_single, print_cumul=print_cumul)
    
    def print_cumul(self):
        self.print_cumul_variable("standard")
    
timer = Timer()

############################################################################################
# Haupt-Methoden
############################################################################################


#rho_bar ist einfach nur step_count + 1, rho ist nur rho_bar mit Vorfaktor
def survivingLatticePoints(step_count, data_lattice_list, density_dict):
    print(len(data_lattice_list))
    l = {lp for lp in data_lattice_list
        if density_dict[lp] > step_count}
    print(len(list(l)))
    return l

def cartesian_product(X : list[list]):
    dim = len(X)
    prod = [()]
    for l in range(dim):
        prod = [ tup + (item,) for tup in prod for item in X[l] ]
    return prod

def cartesian_potentiation(lst : list, dim : int):
    return cartesian_product([lst for _ in range(dim)])

def eff_origin_distance(lattice_point: tuple[int,...], norm):
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

def epsDensityTest(component: set[tuple[int,...]], step_count : int, density_dict, eps_bar) -> bool:

    #TODO Was ist schneller?
    return True in [density_dict[lattice_point] > step_count + eps_bar for lattice_point in component]

    # for lattice_point in component:
    #     if density_dict[lattice_point] > step_count: + eps_bar:
    #         return True
    # return False


# Bestimmen der tau-Zusammenhangskomponenten von surviving_lattice_points durch Aufstellen eines Nähe-Graphen und Tiefensuche
def connectedComponents(data_lattice_list, density_dict, tau_factor, step_count, norm, eps_bar) -> list[set[tuple[int,...]]]:
    
    timer.start_variable("surviving_lattice_points")
    surviving_lattice_points = survivingLatticePoints(step_count, data_lattice_list, density_dict)
    timer.pause_variable("surviving_lattice_points")
    dim = len(next(iter(surviving_lattice_points)))
    timer.start_variable("tau_distance_sets")

    eff_tau_interval_list = list(range(-int(tau_factor)-1, int(tau_factor)+2))
    tau_box = cartesian_potentiation(eff_tau_interval_list, dim)

    tau_distance_set = { lattice_point for lattice_point in tau_box if eff_origin_distance(lattice_point, norm) < tau_factor}

    timer.pause_variable("tau_distance_sets")
    # lattice_large_box_dict = {point : tuple(map(lambda x: x // int(tau_factor), point)) for point in surviving_lattice_points}
    large_boxes = {tuple(map(lambda x: x // int(tau_factor), point)) for point in surviving_lattice_points}
    
    # langsam aber allgemeineres Konzept:
    # def singletonPreimageDict(d : dict) -> dict:
    #   return {target_val : {key for key, val in d.items() if val == target_val} for target_val in d.values()}
    # large_box_lattice_dict = singletonPreimageDict(lattice_large_box_dict)

    # schneller:
    timer.start_variable("inverse_dicts")
    large_box_lattice_dict = {box : set(cartesian_product([list(range(int(tau_factor)*comp,int(tau_factor)*(comp+1))) for comp in box])) & surviving_lattice_points 
                              for box in large_boxes}
    timer.pause_variable("inverse_dicts")

    tau_connection_graph = {}
    coord_shift = lambda coord : coord-1 if coord > 0 else -coord-1 if coord < 0 else 0
    sub = lambda x, y : x - y
    # large_boxes = large_box_lattice_dict.keys()
    relative_neighbor_boxes = cartesian_potentiation(list(range(-1,1+1)), dim)
    timer.start_variable("box_loops")
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
            if tuple(map(sub, point, new_point)) in tau_distance_set}
            tau_connection_graph[new_point] = connected_points
    timer.pause_variable("box_loops")
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

    timer.start_variable("dfs_algorithm")

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
    timer.pause_variable("dfs_algorithm")

    surviving_list : list[bool] = [epsDensityTest(component, step_count, density_dict, eps_bar) 
                                   for component in component_list]
    surviving_components = [component for component, survived in zip(component_list,surviving_list) if survived]
    dead_components : set[tuple[int,...]] = set.union(*[component for component, survived 
                                                        in zip(component_list,surviving_list) 
                                                        if not survived]) if False in surviving_list else set()
    # 0-te Komponente beinhaltet Kästchen, die zu keiner überlebenden Komponente gehören
    surviving_components.insert(0,dead_components)

    return surviving_components




# woanders hin?
# step_count_limit = 10000

# step_count = 0
# connected_component_count = 1
# connected_component_list : list[set[tuple[int,...]]] = []
# while connected_component_count == 1:
#     if step_count > step_count_limit:
#         raise RuntimeError(f"Clustering nach {step_count_limit} Versuchen abgebrochen")
    
#     # timer.start("connectedComponents("+ str(step_count) + ")")
#     connected_component_list = connectedComponents(step_count)
#     connected_component_count = len(connected_component_list) - 1
#     step_count += 1
#     # timer.stop()

# if connected_component_count == 0:
#     clustered_data = [[1] + point for point in data_list]

# else:
#     clustered_data = []
#     cluster_lattice_points = set.union(*connected_component_list)
#     clustered_data = [[connected_component_list.index(component)] + data_point 
#                         for component,data_point in zip(connected_component_list,data_list) 
#                         if tuple(data_point) in component]
    

############################################################################################ 
# Visualisierung
############################################################################################ 


def plot_dataset(data: list[tuple[float,...]], path):
    xs = [p[0] for p in data]
    ys = [p[1] for p in data]

    plt.figure(figsize=(6, 6))
    plt.scatter(xs, ys, s=40, color="blue")
    plt.gca().set_aspect("equal", adjustable="box")
    
    plt.xlim(min(xs)-0.05, max(xs)+0.05)
    plt.ylim(min(ys)-0.05, max(ys)+0.05)

    plt.savefig(path)
    plt.close()

def plot_clusters(clustered_data: list[list], path):
    if not clustered_data:
        print("Keine Daten zum Plotten!")
        return

    xs = [p[1] for p in clustered_data]
    ys = [p[2] for p in clustered_data]
    labels = [p[0] for p in clustered_data]

    plt.figure(figsize=(6, 6))
    plt.gca().set_aspect("equal", adjustable="box")

    plt.xlim(min(xs)-0.05, max(xs)+0.05)
    plt.ylim(min(ys)-0.05, max(ys)+0.05)

    # Punkte nach Cluster: 0 grau, andere farbig
    for cluster_id in set(labels):
        xs_c = [x for x, l in zip(xs, labels) if l == cluster_id]
        ys_c = [y for y, l in zip(ys, labels) if l == cluster_id]
        if cluster_id == 0:
            color = "lightgray"
            plt.scatter(xs_c, ys_c, s=5, c=color)
        else:
            plt.scatter(xs_c, ys_c, s=5)

    plt.savefig(path, dpi=800)
    plt.close()



####################################################################################################
# Schreiben der Daten in .csv/.png-Dateien
####################################################################################################

####################################################################################################
# Main-Funktion damit Benutzer Code über Konsole nutzen kann
####################################################################################################

def main():
    args = prase_args()

    datasetname = args.datasetname
    delta = args.delta
    eps_factor = args.eps_factor
    tau_factor = args.tau_factor
    norm = args.norm

    if norm < 1:
        raise ValueError("norm must be >= 1")

    # Pfade
    data_folder = Path("cluster-data")
    result_folder = Path("cluster-results")
    result_folder.mkdir(exist_ok=True)

    data_path = data_folder / f"{datasetname}.csv"
    result_picture_data = result_folder / f"team-12-{datasetname}.train.png"
    result_picture_clusters = result_folder / f"team-12-{datasetname}.result.png"

    if not data_path.exists():
        raise FileNotFoundError(f"Dataset '{datasetname}' nicht gefunden in {data_folder}")

    # Daten einlesen
    df = pd.read_csv(data_path)
    data_list = df.to_numpy().tolist()
    n_data = len(data_list)
    dim = len(data_list[0])

    # Lattice vorbereiten
    data_lattice_list = [tuple(int((c + 1) // (2 * delta)) for c in pt) for pt in data_list]
    density_dict = {}
    for pt in data_lattice_list:
        density_dict[pt] = density_dict.get(pt, 0) + 1

    h_max_bar = max(density_dict.values())
    eps_bar = eps_factor * (h_max_bar ** 0.5)

    # Clustering-Schleife
    step_count = 0
    connected_component_count = 1
    connected_component_list = []
    step_count_limit = 10000

    while connected_component_count == 1:
        if step_count > step_count_limit:
            raise RuntimeError(f"Clustering nach {step_count_limit} Versuchen abgebrochen")
        timer.start_variable("connectedComponents")
        connected_component_list = connectedComponents(
            data_lattice_list, density_dict, tau_factor, step_count, norm, eps_bar
        )
        timer.pause_variable("connectedComponents")
        connected_component_count = len(connected_component_list) - 1
        step_count += 1

    if connected_component_count == 0:
        clustered_data = [[1] + pt for pt in data_list]
    else:
        clustered_data = []
        for pt in data_list:
            lattice_pt = tuple(int((c + 1) // (2 * delta)) for c in pt)
            cluster_id = 0  # Standard Cluster 0
            for i, component in enumerate(connected_component_list):
                if lattice_pt in component:
                    cluster_id = i
                    break
            clustered_data.append([cluster_id] + pt)

    timer.start_variable("plotting")

    # 2D Plots
    if dim == 2:
        plot_dataset([tuple(item) for item in data_list] , result_picture_data)
        plot_clusters(clustered_data, result_picture_clusters)
    
    timer.pause_variable("plotting")


# result_data_frame = pd.DataFrame(clustered_data(data_list))
# result_data_frame.to_csv(result_path, index=False)


# TODO (global):
# log-Datei Daten implementieren
# richtige Sachen in die richtigen Dateien schreiben
# finale Datenvisualisierung für 2D
# ...


if __name__ == "__main__":

    timer.start()
    main()
    timer.stop()
    timer.print_cumul_variable("surviving_lattice_points")
    timer.print_cumul_variable("tau_distance_sets")
    timer.print_cumul_variable("inverse_dicts")
    timer.print_cumul_variable("box_loops")
    timer.print_cumul_variable("dfs_algorithm")
    timer.print_cumul_variable("connectedComponents")
    timer.print_cumul_variable("plotting")

# Berechtigung zum ausführen (u oder a)
#!/usr/bin/python3???
# aufpassen mit Windows/Linus (ein /r was Probleme macht?)