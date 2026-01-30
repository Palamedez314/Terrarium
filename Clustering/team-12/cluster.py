#!/usr/bin/env python3

import pandas as pd
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from pathlib import Path
from customTimer import CustomTimer
from plotting import plot_dataset, plot_clusters

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("datasetname", help="Name of the dataset out of the folder 'cluster-data' (without .csv extension)")

def parse_args():
    parser = ArgumentParser(prog="cluster", description= "Team 12, Aufgabe 2, Cluster Algorithmus", formatter_class=ArgumentDefaultsHelpFormatter) 

    parser.add_argument("datasetname",help="Name of the dataset in 'cluster-data' folder (without .csv extension)")
    parser.add_argument("-d", "--delta", type=float, default=0.05, help="Parameter, that determines the fineness of the underlying lattice")
    parser.add_argument("-e", "--eps-factor", type=float, default=3.0, help="tbd")
    parser.add_argument("-t", "--tau-factor", type=float, default=2.0, help="tbd")
    parser.add_argument("-n", "--norm", type=float, default=float("inf"), help="tbd")
    return parser.parse_args()

############################################################################################
# Methoden
############################################################################################

def round_up(x : float) -> int:
    return int(x) if float(int(x)) == x else int(x)+1

def cartesian_product(X : list[list]):
    dim = len(X)
    prod = [()]
    for l in range(dim):
        prod = [ tup + (item,) for tup in prod for item in X[l] ]
    return prod

def cartesian_potentiation(lst : list, dim : int):
    return cartesian_product([lst for _ in range(dim)])

def neighborTest(b1, b2):
    for coord1, coord2 in zip(b1, b2):
        if abs(coord1 - coord2) > 1:
            return False
    return True

def isBoundaryPoint(pt, mod, tau_eff):
    for coord in pt:
        if coord % mod < tau_eff or coord % mod >= mod - tau_eff:
            return True
    return False

####################################################################################################
# Ausgeführter Code
####################################################################################################

def main():

    ############### Übergeben der Argumente aus Interface ###############
    args = parse_args()

    datasetname = args.datasetname
    delta = args.delta
    eps_factor = args.eps_factor
    tau_factor = args.tau_factor
    norm = args.norm

    if norm < 1:
        raise ValueError("norm must be >= 1")

    # Pfade
    data_folder = Path("../cluster-data")
    result_folder = Path("../cluster-results")
    result_folder.mkdir(exist_ok=True)

    data_path = data_folder / f"{datasetname}.csv"
    result_dataset_path = result_folder / f"team-12-{datasetname}.result.csv"
    result_picture_data_path = result_folder / f"team-12-{datasetname}.train.png"
    result_picture_clusters_path = result_folder / f"team-12-{datasetname}.result.png"
    log_path = result_folder / f"team-12-{datasetname}.log"
    result_log_path = result_folder / f"team-12-{datasetname}.result.log"

    if not data_path.exists():
        raise FileNotFoundError(f"Dataset '{datasetname}' nicht gefunden in {data_folder}")


    ############### Einlesen der Daten ###############

    df = pd.read_csv(data_path, header=None)
    data_list = df.to_numpy().tolist()

    ############### Datenverarbeitung (ohne Verwendung von externen Modulen) ###############
    timer = CustomTimer()
    timer.start()

    n_data = len(data_list)
    dim = len(data_list[0])
    tau_eff = round_up(tau_factor/2)
    rho_step = 1 / (n_data * (2 * delta) ** dim)

    timer.start_variable("density-dict")

    # Lattice vorbereiten
    timer.start_variable("data_lattice_list")
    data_lattice_list = [tuple(int((c + 1) // (2 * delta)) for c in pt) for pt in data_list]
    timer.pause_variable("data_lattice_list")
    density_dict = {}
    inverse_density_dict = {}
    for pt in data_lattice_list:
        old_density = density_dict.get(pt, 0)
        density_dict[pt] = old_density + 1
        if old_density > 0:
            inverse_density_dict[old_density].remove(pt)
        if old_density + 1 not in inverse_density_dict.keys():
            inverse_density_dict[old_density + 1] = set()
        inverse_density_dict[old_density + 1].add(pt)
    timer.pause_variable("density-dict")

    density_values = inverse_density_dict.keys()
    
    h_max_bar = max(density_dict.values())
    eps_bar = eps_factor * (h_max_bar ** 0.5)

    # Clustering-Schleife
    step_count = 0
    connected_component_count = 1
    connected_component_list = []
    surviving_lattice_points = set(data_lattice_list)

    rho_list = [0.0]
    cluster_sizes_1 = [n_data]
    cluster_sizes_2 = [0]

    while connected_component_count == 1:
        if step_count > h_max_bar:
            raise RuntimeError(f"Clustering nach {h_max_bar} Versuchen abgebrochen weil finden von Clustern ab jetzt unmöglich ist")
        
        # timer.start_variable("connected components")

        if step_count in density_values:
            surviving_lattice_points -= inverse_density_dict[step_count]

        # if step_count == 0:
        #     print(surviving_lattice_points)

        def tauDistanceTest(pt1, pt2):
            for coord1, coord2 in zip(pt1, pt2):
                if abs(coord1 - coord2) > tau_eff:
                    return False
            return True

        timer.start_variable("boxes_same_method")

        box_size = 2
        mod = box_size * tau_eff

        lattice_large_boxe_dict = {point : tuple(map(lambda x: x // mod, point)) for point in surviving_lattice_points}
        large_boxes = set(lattice_large_boxe_dict.values())
        large_box_lattice_dict = {box : set(cartesian_product([list(range(mod*comp,mod*(comp+1))) for comp in box])) & surviving_lattice_points 
                                        for box in large_boxes}

        large_box_neighbor_graph = {}
        old_boxes = set()
        for new_box in large_boxes:
            neighboring_boxes = {box for box in old_boxes if neighborTest(box, new_box)}
            large_box_neighbor_graph[new_box] = neighboring_boxes
            for box in neighboring_boxes:
                large_box_neighbor_graph[box].add(new_box)
            old_boxes.add(new_box)

        tau_connection_graph = {}
        old_points = set()
        for box in large_boxes:
            
            box_lattice_points = large_box_lattice_dict[box]
            neighboring_boxes_lattice_points = set.union(*[large_box_lattice_dict[neighboring_box]
                                                           for neighboring_box 
                                                           in large_box_neighbor_graph[box]]) if bool(large_box_neighbor_graph[box]) else set()
            neighboring_boxes_lattice_points.update(box_lattice_points)
            for new_point in large_box_lattice_dict[box]:
                # lieber absolut (ohne sowas wie old_points) laufen und dann ohne den Schnitt?
                if isBoundaryPoint(new_point, mod, tau_eff):
                    relevant_points = neighboring_boxes_lattice_points & old_points
                else:
                    relevant_points = box_lattice_points & old_points
                connected_points = {point for point in relevant_points if tauDistanceTest(point, new_point)}
                tau_connection_graph[new_point] = connected_points
                for point in connected_points:
                    tau_connection_graph[point].add(new_point)
                old_points.add(new_point)

        timer.pause_variable("boxes_same_method")


        timer.start_variable("dfs algorithm")
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
        timer.pause_variable("dfs algorithm")

        def epsDensityTest(component: set[tuple[int,...]]) -> bool:
            for lattice_point in component:
                if density_dict[lattice_point] > step_count + eps_bar:
                    return True
            return False

        surviving_list : list[bool] = [epsDensityTest(component) for component in component_list]
        
        connected_component_list = [component for component, survived in zip(component_list,surviving_list) if survived]
        dead_components : set[tuple[int,...]] = set.union(*[component for component, survived 
                                                            in zip(component_list,surviving_list) 
                                                            if not survived]) if False in surviving_list else set()
        # 0-te Komponente beinhaltet Kästchen, die zu keiner überlebenden Komponente gehören
        connected_component_list.insert(0,dead_components)
        # timer.pause_variable("connected components")

        connected_component_count = len(connected_component_list) - 1
        step_count += 1
        
        rho_list.append(step_count * rho_step)
        cluster_sizes_1.append(len(connected_component_list[1]) if connected_component_count > 0 else 0)
        cluster_sizes_2.append(len(connected_component_list[2]) if connected_component_count > 1 else 0)

    timer.start_variable("packing data")
    # geht das effizienter?
    if connected_component_count == 0:
        clustered_data = [[1] + pt for pt in data_list]
    else:
        clustered_data = []
        for i, pt in enumerate(data_list):
            lattice_pt = data_lattice_list[i]
            id = 0
            for j, component in enumerate(connected_component_list):
                if lattice_pt in component:
                    id = j
                    break
            clustered_data.append([id] + pt)
    timer.pause_variable("packing data")

    log_data = list(map(list, zip(rho_list, cluster_sizes_1, cluster_sizes_2)))

    timer.stop(print_single=False)

    runtime = timer.get_cumul()
    result_log_data = [runtime, cluster_sizes_1[-1], cluster_sizes_2[-1], rho_list[-1]]

    ############### Schreiben der Daten in .csv/.log-Dateien ###############

    result_data_frame = pd.DataFrame(clustered_data)
    result_data_frame.to_csv(result_dataset_path, index=False, header=False)
    log_data_frame = pd.DataFrame(log_data)
    log_data_frame.to_csv(log_path, index=False, header=False)
    result_log_data_frame = pd.DataFrame(result_log_data)
    result_log_data_frame.to_csv(result_log_path, index=False, header=False)

    ############### Plotten der Daten ###############

    timer.start_variable("plotting")
    # # 2D Plots
    if dim == 2:
        plot_dataset([tuple(item) for item in data_list] , result_picture_data_path)
        plot_clusters(clustered_data, result_picture_clusters_path)
    timer.pause_variable("plotting")

    
    # Zum Optimieren, manche Zeiten doppeln sich (z.B. "box-loops" und "def algorithm" in connected components enthalten)
    timer.print_cumul()
    timer.print_cumul_variable("data_lattice_list")
    timer.print_cumul_variable("density-dict")
    # timer.print_cumul_variable("inverse-dict")
    # timer.print_cumul_variable("box-loops")
    timer.print_cumul_variable("boxes_same_method")
    timer.print_cumul_variable("dfs algorithm")
    # timer.print_cumul_variable("connected components")
    timer.print_cumul_variable("packing data")
    timer.print_cumul_variable("plotting")


if __name__ == "__main__":
    main()

# TODO:

# Performance verbessern (besserer Alghorithmus für "box-loop"?)

# alt:
# sed 's/\r$//' cluster.py > cluster1.py
# chmod +rwx cluster1.py
# <Umbenennen in cluster.py>
# tar -czf team-12.tar.gz cluster.py customTimer.py plotting.py

# neu:
# Einstellung LF
# chmod +rwx cluster.py
# tar -czf team-12.tar.gz cluster.py customTimer.py plotting.py