#!/usr/bin/env python3

############################################################################################
# Importe
############################################################################################

############### Externe Module ###############

import pandas as pd
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from pathlib import Path

############### Eigene Klassen und Methoden (und Typabkürzung) ###############

from customTimer import CustomTimer
from plotting import plot_dataset, plot_clusters
from classter import Cluster, point

############################################################################################
# Interface zum Ausführen im Terminal
############################################################################################

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("datasetname", help="Name of the dataset out of the folder 'cluster-data' (without .csv extension)")

def parse_args():
    parser = ArgumentParser(prog="cluster", description= "Team 12, Aufgabe 2, Cluster Algorithmus", formatter_class=ArgumentDefaultsHelpFormatter) 

    parser.add_argument("datasetname",help="Name of the dataset in 'cluster-data' folder (without .csv extension)")
    parser.add_argument("-d", "--delta", type=float, default=0.05, help="Parameter, that determines the fineness of the underlying lattice")
    parser.add_argument("-e", "--eps-factor", type=float, default=3.0, help="Additional density threshold for a connected component (in respect to tau) to be considered a cluster.")
    parser.add_argument("-t", "--tau-factor", type=float, default=2.0, help="Parameter, that determines the minimum distance of different clusters")
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

def smallBoxinLargeBoxTest (sb : point, lb : point) -> bool:
    for coord_sb, coord_lb in zip(sb, lb):
        diff = coord_sb - coord_lb
        if diff != 0 and diff != 1:
            return False
    return True

def epsDensityTest(ct: Cluster, rho_bar : int, eps_bar: float) -> bool:
    return ct.max_density > rho_bar + eps_bar

def tauDistanceTest(pt1 : point, pt2 : point, tau_eff : int) -> bool:
    for coord1, coord2 in zip(pt1, pt2):
        if abs(coord1 - coord2) > tau_eff:
            return False
    return True

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


    ############### Bau der Ergebnis Pfade ###############
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

    ############### Starten des Timers ###############
    timer = CustomTimer()
    timer.start()

    ############### Datenverarbeitung (ohne Verwendung von externen Modulen) ###############

    n_data : int = len(data_list)
    dim : int = len(data_list[0])
    tau_eff : int = round_up(tau_factor/2)
    rho_step : float = 1 / (n_data * (2 * delta) ** dim)
    
    # Zuweisung der Datenpunkte auf nächstgelegene Eckpunkte im Gitter 
    data_lattice_list : list[point] = [tuple(int((c + 1) // (2 * delta)) for c in pt) for pt in data_list]

    # timer.start_variable("density-dict")

    # Weist Gittenpunkten ihre Dichte zu
    density_dict : dict[point, int] = {}

    # und umgekehrt Dichten die Gitterpunkte mit entsprechender Dichte
    inverse_density_dict : dict[int, set[point]] = {}

    # Befüllen der dictionaries, indem die einzelnen Datenpunkte durchgezählt werden
    for pt in data_lattice_list:
        old_density = density_dict.get(pt, 0)
        density_dict[pt] = old_density + 1
        if old_density > 0:
            inverse_density_dict[old_density].remove(pt)
        if old_density + 1 not in inverse_density_dict.keys():
            inverse_density_dict[old_density + 1] = set()
        inverse_density_dict[old_density + 1].add(pt)

    density_values = inverse_density_dict.keys()
    max_density = max(density_values)

    # timer.pause_variable("density-dict")
    
    eps_bar = eps_factor * (max_density ** 0.5)

    # timer.start_variable("preparing boxes")

    # Länge einer kleinen Box in jede Dimension
    small_box_size = 2 * tau_eff

    # Weist Punkt die große Box (bei tau_eff=1: 4 lang) zu, in der alle tau entfernten Punkte auch liegen
    lattice_large_box_dict : dict[point,point] = {}

    # Weist große Box die Punkte zu, die in lattice_large_box_dict auf sie gemapped werden ("inverse Multifunktion")
    large_box_lattice_dict : dict[point,set[point]] = {}

    # Weist kleiner Box (bei tau_eff=1: 2 lang) Punkte in ihr zu
    small_box_lattice_dict : dict[point,set[point]] = {}

    # Weist großer Box (bei tau_eff=1: 4 lang) kleine Boxen in ihr zu
    large_box_small_boxes_dict : dict[point,set[point]] = {}

    # Gleichzeitiges Befüllen der dictionaries
    for pt in data_lattice_list:

        # Kleine Box, in der pt liegt
        corr_small_box = tuple(map(lambda x: x // small_box_size, pt))
        
        if corr_small_box not in small_box_lattice_dict.keys():
            small_box_lattice_dict[corr_small_box] = set()
        small_box_lattice_dict[corr_small_box].add(pt)

        # Relative Position innerhalb der kleinen Box
        relative_large_box_offset = tuple(map(lambda coord: -1 if coord % small_box_size < tau_eff else 0, pt))

        # Große Box, in der alle tau entfernten Punkte auch liegen
        corr_large_box = tuple(map(int.__add__, corr_small_box, relative_large_box_offset))

        lattice_large_box_dict[pt] = corr_large_box
        if corr_large_box not in large_box_lattice_dict.keys():
            large_box_lattice_dict[corr_large_box] = set()
        large_box_lattice_dict[corr_large_box].add(pt)

    # Bewohnte kleine Boxen
    small_boxes = set(small_box_lattice_dict.keys())

    # Große Boxen, die benutzt werden
    large_boxes = set(lattice_large_box_dict.values())

    large_box_count = len(large_boxes)

    # timer.initialize_variable("cart_pot")

    # Abschätzung welche Methode schneller ist
    if large_box_count > 2**dim:

        # timer.start_variable("cart_pot")

        # Alle möglichen dim-Tupel bestehend aus -1 und 0
        neg_one_zero_list = cartesian_potentiation([-1,0],dim)

        # timer.pause_variable("cart_pot")

        for small_box in small_boxes:

            # timer.start_variable("cart_pot")

            potential_large_boxes : set[point] = {tuple(map(int.__add__, small_box, diff)) for diff in neg_one_zero_list}
            
            # timer.pause_variable("cart_pot")

            # Große Boxen, in denen small_box liegt
            corr_large_boxes : set[point] = potential_large_boxes & large_boxes

            # small_box zu allen corr_large_boxes ins dictionary
            for large_box in corr_large_boxes:
                if large_box not in large_box_small_boxes_dict.keys():
                    large_box_small_boxes_dict[large_box] = set()
                large_box_small_boxes_dict[large_box].add(small_box)
    else:
        for small_box in small_boxes:

            # small_box zu jeder Großen Box ins dictionary, in der sie liegt
            for large_box in large_boxes:
                
                # Abfrage, ob small_box in large_box liegt
                if smallBoxinLargeBoxTest(small_box, large_box):
                    
                    # small_box zu large_box ins dictionary
                    if large_box not in large_box_small_boxes_dict.keys():
                        large_box_small_boxes_dict[large_box] = set()
                    large_box_small_boxes_dict[large_box].add(small_box)

    # print("used cart_pot method" if large_box_count > 2**dim else "used large_box iteration method")

    # timer.pause_variable("preparing boxes")
    
    # timer.start_variable("cluster-loop")

    clusters : set[Cluster] = set()
    visible_clusters = set()
    iterated_points = set()
    point_cluster_dict : dict[point, Cluster] = {}

    result_cluster_sets : list[set[point]] = []
    result_cluster_count : int = 1
    result_rho_bar = 0

    first_cluster_sizes_rev = []
    second_cluster_sizes_rev = []

    # Iteration von rho von maximaler Dichte bis 0
    for rho_bar in range(max_density, 0, -1):

        old_visible_clusters = visible_clusters.copy()
        
        # Punkte mit Dichte rho_bar
        new_layer = inverse_density_dict[rho_bar]
        
        # TODO: vllt noch tau_connection_graph Ansatz mit Pointern, sonst auskommentierten Code entfernen

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

        iterating_set = new_layer.copy()

        while bool(iterating_set):

            layer_point = iterating_set.pop()

            corr_large_box = lattice_large_box_dict[layer_point]
            
            corr_small_boxes = large_box_small_boxes_dict[corr_large_box]

            large_box_points = set.union(*[small_box_lattice_dict[sb] for sb in corr_small_boxes])

            # Punkte, die auch auf corr_large_box mappen aus new_layer
            corr_large_box_new_layer_points = large_box_lattice_dict[corr_large_box] & new_layer

            for new_point in corr_large_box_new_layer_points:
            
                # Punkte, die weniger als tau Abstand von new_point haben
                connected_points = {pt for pt in iterated_points & large_box_points if tauDistanceTest(pt, new_point, tau_eff)}

                connected_clusters = {point_cluster_dict[pt] for pt in connected_points}

                if bool(connected_clusters):
                    first_cluster = connected_clusters.pop()

                    if bool(connected_clusters):
                        # Mergen aller anderen Cluster in first_cluster
                        first_cluster.merge(*connected_clusters)
                        clusters -= connected_clusters
                        connected_cluster_points = set.union(*[ct.points for ct in connected_clusters])
                        for pt in connected_cluster_points:
                            point_cluster_dict[pt] = first_cluster

                    # Hinzufügen von new_point in first_cluster
                    first_cluster.add_point(new_point, density_dict[new_point])
                    point_cluster_dict[new_point] = first_cluster

                else:
                    # Erstellen eines neuen Clusters, das nur aus dem Punkt new_point besteht
                    new_cluster = Cluster({new_point}, {density_dict[new_point]})
                    point_cluster_dict[new_point] = new_cluster
                    clusters.add(new_cluster)

                iterated_points.add(new_point)

            iterating_set -= corr_large_box_new_layer_points

        for cluster in clusters:
            cluster.update_visibility(epsDensityTest(cluster, rho_bar, eps_bar))

        visible_clusters = {ct for ct in clusters if ct.visible}

        visible_cluster_count = len(visible_clusters)
        visible_clusters_list = list(visible_clusters)

        if visible_cluster_count != 1:

            if visible_clusters != old_visible_clusters:

                # Sammeln der möglichen finalen Cluster
                result_cluster_sets = [ct.points.copy() for ct in visible_clusters_list]
                hidden_clusters = clusters - visible_clusters
                hidden_cluster_points = set.union(*[ct.points for ct in hidden_clusters]) if bool(hidden_clusters) else set()
                result_cluster_sets.insert(0, hidden_cluster_points)
                result_cluster_count = visible_cluster_count

            # Mögliches finales rho_bar
            result_rho_bar = rho_bar

        # Cluster-Längen für log-Dateien
        first_cluster_sizes_rev.append(len(visible_clusters_list[0].points) if visible_cluster_count > 0 else 0)
        second_cluster_sizes_rev.append(len(visible_clusters_list[1].points) if visible_cluster_count > 1 else 0)

    # TODO: wollen wir so nen assert hier drin haben oder ist das zu gefährlich :)
    # Überprüfen, ob das Höchzählende Verfahren terminiert wäre (also ob das Ergebnis korrekt ist)
    assert(result_cluster_count != 1)

    # timer.pause_variable("cluster-loop")

    # Cluster-Längen für rho_bar = 0
    first_cluster_sizes_rev.append(n_data)
    second_cluster_sizes_rev.append(0)

    # timer.start_variable("packing data")

    # TODO: geht das effizienter?
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

    # timer.pause_variable("packing data")

    # Daten für .log Datei
    first_cluster_sizes = list(reversed(first_cluster_sizes_rev))[:result_rho_bar + 1]
    second_cluster_sizes = list(reversed(second_cluster_sizes_rev))[:result_rho_bar + 1]
    rho_list = [i * rho_step for i in range(result_rho_bar + 1)]
    log_data = list(map(list, zip(rho_list, first_cluster_sizes, second_cluster_sizes)))

    # Daten für .result.log Datei
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

    ############### Plotten der Daten für dim == 2 ###############

    # timer.start_variable("plotting")

    if dim == 2:
        # Plotten des Datensatzes
        plot_dataset([tuple(item) for item in data_list] , result_picture_data_path)

        # Plotten des geclusterten Datensatzes
        plot_clusters(clustered_data, result_picture_clusters_path)

    # timer.pause_variable("plotting")

    ############### Laufzeiten einzelner Teile des Programms zum testen ###############

    # timer.print_cumul()
    # timer.print_cumul_variable("density-dict")
    # timer.print_cumul_variable("cart_pot")
    # timer.print_cumul_variable("preparing boxes")
    # timer.print_cumul_variable("cluster-loop")
    # timer.print_cumul_variable("packing data")
    # timer.print_cumul_variable("plotting")

if __name__ == "__main__":
    main()