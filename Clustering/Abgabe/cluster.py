#!/usr/bin/env python


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

def epsDensityTest(component: set[tuple[int,...]], step_count : int, density_dict, eps_bar) -> bool:
    for lattice_point in component:
        if density_dict[lattice_point] > step_count + eps_bar:
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
    data_folder = Path("cluster-data")
    result_folder = Path("cluster-results")
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

    timer.start_variable("read in")

    df = pd.read_csv(data_path)
    data_list = df.to_numpy().tolist()
    
    timer.pause_variable("read in")

    ############### Datenverarbeitung (ohne Verwendung von externen Modulen) ###############

    dim = len(data_list[0])
    tau_eff = round_up(tau_factor/2)

    timer.start_variable("density-dict")
    # Lattice vorbereiten
    data_lattice_list = [tuple(int((c + 1) // (2 * delta)) for c in pt) for pt in data_list]
    density_dict = {}
    for pt in data_lattice_list:
        density_dict[pt] = density_dict.get(pt, 0) + 1
    timer.pause_variable("density-dict")

    density_values = density_dict.values()
    timer.start_variable("inverse-dict")
    inverse_density_dict = {k : {point for point, density in density_dict.items() if density == k} for k in density_values}
    timer.pause_variable("inverse-dict")
    
    h_max_bar = max(density_dict.values())
    eps_bar = eps_factor * (h_max_bar ** 0.5)

    # Clustering-Schleife
    step_count = 0
    connected_component_count = 1
    connected_component_list = []
    surviving_lattice_points = set(data_lattice_list)

    while connected_component_count == 1:
        if step_count > h_max_bar:
            raise RuntimeError(f"Clustering nach {h_max_bar} Versuchen abgebrochen weil finden von Clustern ab jetzt unmöglich ist")
        
        timer.start_variable("connected components")

        if step_count in density_values:
            surviving_lattice_points -= inverse_density_dict[step_count]

        timer.start_variable("box-loops")
        tau_connection_graph = {}
        coord_shift = lambda coord : coord-1 if coord > 0 else -coord-1 if coord < 0 else 0
        sub = lambda x, y : x - y
        old_points = set()
        for new_point in surviving_lattice_points:
            # rel_points = [tuple(map(lambda x, y : x - y, point, new_point)) for point in old_points]
            connected_points = {point for point in old_points if max(map(coord_shift, map(sub, point, new_point))) < tau_eff}
            tau_connection_graph[new_point] = connected_points
            for point in connected_points:
                tau_connection_graph[point].add(new_point)
            old_points.add(new_point)
        timer.pause_variable("box-loops")

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

        surviving_list : list[bool] = [epsDensityTest(component, step_count, density_dict, eps_bar) 
                                       for component in component_list]
        
        connected_component_list = [component for component, survived in zip(component_list,surviving_list) if survived]
        dead_components : set[tuple[int,...]] = set.union(*[component for component, survived 
                                                            in zip(component_list,surviving_list) 
                                                            if not survived]) if False in surviving_list else set()
        # 0-te Komponente beinhaltet Kästchen, die zu keiner überlebenden Komponente gehören
        connected_component_list.insert(0,dead_components)
        timer.pause_variable("connected components")

        connected_component_count = len(connected_component_list) - 1
        step_count += 1

    timer.start_variable("packing data")

    if connected_component_count == 0:
        clustered_data = [[1] + pt for pt in data_list]
    else:
        clustered_data = []
        for i, pt in enumerate(data_list):
            lattice_pt = data_lattice_list[i]
            id = 0  # Standard Cluster 0
            for j, component in enumerate(connected_component_list):
                if lattice_pt in component:
                    id = j
                    break
            clustered_data.append([id] + pt)

    timer.pause_variable("packing data")

    ############### Schreiben der Daten in .csv/.log-Dateien ###############

    result_data_frame = pd.DataFrame(clustered_data)
    result_data_frame.to_csv(result_dataset_path, index=False)

    ############### Plotten der Daten ###############

    timer.start_variable("plotting")
    # 2D Plots
    if dim == 2:
        plot_dataset([tuple(item) for item in data_list] , result_picture_data_path)
        plot_clusters(clustered_data, result_picture_clusters_path)
    timer.pause_variable("plotting")

if __name__ == "__main__":
    timer = CustomTimer() 
    timer.start()
    main()
    timer.stop()

    # Zum Optimieren, manche Zeiten doppeln sich (z.B. "box-loops" und "def algorithm" in connected components enthalten)
    timer.print_cumul_variable("read in")
    timer.print_cumul_variable("density-dict")
    timer.print_cumul_variable("inverse-dict")
    timer.print_cumul_variable("box-loops")
    timer.print_cumul_variable("dfs algorithm")
    timer.print_cumul_variable("connected components")
    timer.print_cumul_variable("packing data")
    timer.print_cumul_variable("plotting")


# TODO:

# log-Datei Dinge implementieren

# Performance verbessern (besserer Alghorithmus für "box-loop"?)

# "optionale weitere Aspekte"

# Genaue Anforderungen auf seiner Seite erfüllen:
# Berechtigung zum ausführen (u oder a)?
# aufpassen mit Windows/Linux (ein /r was Probleme macht?)