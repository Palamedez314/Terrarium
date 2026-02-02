#!/usr/bin/env python3

import pandas as pd
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from pathlib import Path
from customTimer import CustomTimer
from plotting import plot_dataset, plot_clusters
from classter import Cluster, point

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
    max_density = max(density_values)
    eps_bar = eps_factor * (max_density ** 0.5)

    def tauDistanceTest(pt1, pt2) -> bool:
        for coord1, coord2 in zip(pt1, pt2):
            if abs(coord1 - coord2) > tau_eff:
                return False
        return True

    def epsDensityTest(ct: Cluster) -> bool:
        return ct.max_density > rho_bar + eps_bar
    
    timer.start_variable("cluster-loop")

    # Länge einer Box in jede Dimension
    small_box_size = 2 * tau_eff
    large_box_size = 4 * tau_eff

    # Weist Punkt die große Box (bei tau_eff=1: 4 lang) zu, in der alle tau entfernten Punkte auch liegen
    lattice_large_box_dict : dict[point,point] = {}

    # Weist kleiner Box (bei tau_eff=1: 2 lang) Punkte in ihr zu
    small_box_lattice_dict : dict[point,set[point]] = {}

    # Weist großer Box (bei tau_eff=1: 4 lang) kleine Boxen in ihr zu
    large_box_small_boxes_dict : dict[point,set[point]] = {}

    #Gleichzeitiges Befüllen der dictionaries
    for pt in data_lattice_list:

        corr_small_box = tuple(map(lambda x: x // small_box_size, pt))
        # Punkt zu corr_small_box ins dictionary
        if corr_small_box not in small_box_lattice_dict.keys():
            small_box_lattice_dict[corr_small_box] = set()
        small_box_lattice_dict[corr_small_box].add(pt)

        # Relative Position innerhalb der kleinen Box
        relative_large_box_offset = tuple(map(lambda coord: -1 if coord % small_box_size < tau_eff else 0, pt))
        # Große Box, in der alle tau entfernten Punkte auch liegen
        corr_large_box = tuple(map(int.__add__, corr_small_box, relative_large_box_offset))
        lattice_large_box_dict[pt] = corr_large_box

    small_boxes = set(small_box_lattice_dict.keys())
    large_boxes = set(lattice_large_box_dict.values())

    def smallBoxinLargeBoxTest (sb : point, lb : point) -> bool:
        for coord_sb, coord_lb in zip(sb, lb):
            diff = coord_sb - coord_lb
            if diff != 0 and diff != 1:
                return False
        return True
    
    timer.start_variable("large_box_shit")

    timer.start_variable("cart_pot")
    neg_one_zero_list = cartesian_potentiation([-1,0],dim)
    timer.pause_variable("cart_pot")

    large_box_count = len(large_boxes)

    if large_box_count > 2**dim:
        for small_box in small_boxes:
            potential_large_boxes : list[point] = [tuple(map(int.__add__, small_box, diff)) for diff in neg_one_zero_list]
            corr_large_boxes : set[point] = set(potential_large_boxes) & large_boxes
            for large_box in corr_large_boxes:
                if large_box not in large_box_small_boxes_dict.keys():
                    large_box_small_boxes_dict[large_box] = set()
                large_box_small_boxes_dict[large_box].add(small_box)
    else:
        for small_box in small_boxes:

            # small_box zu jeder Großen Box ins dictionary, in der sie liegt
            for large_box in large_boxes:

                if smallBoxinLargeBoxTest(small_box, large_box):
                    
                    # small_box zu large_box ins dictionary
                    if large_box not in large_box_small_boxes_dict.keys():
                        large_box_small_boxes_dict[large_box] = set()
                    large_box_small_boxes_dict[large_box].add(small_box)

    # print("used cart_pot method" if large_box_count > 2**dim else "used large_box iteration method")

    timer.pause_variable("large_box_shit")

    clusters : set[Cluster] = set()
    iterated_points = set()
    point_cluster_dict : dict[point, Cluster] = {}

    result_cluster_sets : list[set[point]] = []
    result_cluster_count : int = 1
    result_rho_bar = 0

    first_cluster_sizes_rev = []
    second_cluster_sizes_rev = []

    for rho_bar in range(max_density, 0, -1):
        
        new_layer = inverse_density_dict[rho_bar]

        if not bool(new_layer):
            continue
        
        # TODO: vllt noch tau_connection_graph Ansatz mit Pointern

        # tau_connection_graph = {}
        # old_points = set()
        # for box in large_boxes:
            
            # box_lattice_points = large_box_lattice_dict[box]
            # neighboring_boxes_lattice_points = set.union(*[large_box_lattice_dict[neighboring_box]
            #                                                for neighboring_box 
            #                                                in large_box_neighbor_graph[box]]) if bool(large_box_neighbor_graph[box]) else set()
            # neighboring_boxes_lattice_points.update(box_lattice_points)
            # for new_point in large_box_lattice_dict[box]:
            #     # lieber absolut (ohne sowas wie old_points) laufen und dann ohne den Schnitt?
            #     if isBoundaryPoint(new_point, box_size, tau_eff):
            #         relevant_points = neighboring_boxes_lattice_points & old_points
            #     else:
            #         relevant_points = box_lattice_points & old_points
            #     connected_points = {pt for pt in relevant_points if tauDistanceTest(pt, new_point)}
            #     tau_connection_graph[new_point] = connected_points
            #     for pt in connected_points:
            #         tau_connection_graph[pt].add(new_point)
            #     old_points.add(new_point)

        # for new_point in new_layer:

        #     associated_box = lattice_large_box_dict[new_point]

        #     neighboring_boxes = large_box_neighbor_graph[associated_box]

        #     neighboring_boxes_lattice_points

        for new_point in new_layer:

            corr_large_box = lattice_large_box_dict[new_point]

            corr_small_boxes = large_box_small_boxes_dict[corr_large_box]

            large_box_points = set.union(*[small_box_lattice_dict[sb] for sb in corr_small_boxes])

            connected_points = {pt for pt in iterated_points & large_box_points if tauDistanceTest(pt, new_point)}
            connected_clusters = {point_cluster_dict[pt] for pt in connected_points}

            if len(connected_clusters) == 0:
                new_cluster = Cluster({new_point}, {density_dict[new_point]})
                point_cluster_dict[new_point] = new_cluster
                clusters.add(new_cluster)

            else:
                first_cluster = connected_clusters.pop()

                if len(connected_clusters) > 0:
                    first_cluster.merge(*connected_clusters)
                    clusters -= connected_clusters
                    connected_cluster_points = set.union(*[ct.points for ct in connected_clusters])
                    for pt in connected_cluster_points:
                        point_cluster_dict[pt] = first_cluster

                first_cluster.add_point(new_point, density_dict[new_point])
                point_cluster_dict[new_point] = first_cluster

            iterated_points.add(new_point)


        for cluster in clusters:
            cluster.update_visibility(epsDensityTest(cluster))

        visible_clusters = {ct for ct in clusters if ct.visible}

        visible_cluster_count = len(visible_clusters)

        if visible_cluster_count != 1:
            result_cluster_sets = [ct.points.copy() for ct in visible_clusters]
            hidden_clusters = clusters - visible_clusters
            hidden_cluster_points = set.union(*[ct.points for ct in hidden_clusters]) if bool(hidden_clusters) else set()
            result_cluster_sets.insert(0, hidden_cluster_points)
            result_cluster_count = visible_cluster_count
            result_rho_bar = rho_bar

        visible_clusters_list = list(visible_clusters)
        first_cluster_sizes_rev.append(len(visible_clusters_list[0].points) if visible_cluster_count > 0 else 0)
        second_cluster_sizes_rev.append(len(visible_clusters_list[1].points) if visible_cluster_count > 1 else 0)

    assert(result_cluster_count != 1)

    timer.pause_variable("cluster-loop")

    first_cluster_sizes_rev.append(n_data)
    second_cluster_sizes_rev.append(0)

    first_cluster_sizes = list(reversed(first_cluster_sizes_rev))[:result_rho_bar + 1]
    second_cluster_sizes = list(reversed(second_cluster_sizes_rev))[:result_rho_bar + 1]
    rho_list = [i * rho_step for i in range(result_rho_bar + 1)]

    timer.start_variable("packing data")
    # geht das effizienter?
    if result_cluster_count == 0:
        clustered_data = [[1] + pt for pt in data_list]
    else:
        clustered_data = []
        for i, pt in enumerate(data_list):
            lattice_pt = data_lattice_list[i]
            id = 0
            for j, cluster_set in enumerate(result_cluster_sets):
                if lattice_pt in cluster_set:
                    id = j
                    break
            clustered_data.append([id] + pt)
    timer.pause_variable("packing data")

    log_data = list(map(list, zip(rho_list, first_cluster_sizes, second_cluster_sizes)))

    timer.stop(print_single=False)

    runtime = timer.get_cumul()
    result_log_data = [runtime, first_cluster_sizes[-1], second_cluster_sizes[-1], rho_list[-1]]

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
    # timer.print_cumul()
    # timer.print_cumul_variable("data_lattice_list")
    # timer.print_cumul_variable("density-dict")
    # timer.print_cumul_variable("large_box_shit")
    # timer.print_cumul_variable("cart_pot")
    # timer.print_cumul_variable("cluster-loop")
    # timer.print_cumul_variable("packing data")
    # timer.print_cumul_variable("plotting")


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