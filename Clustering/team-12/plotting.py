############################################################################################ 
# Visualisierung
############################################################################################ 

import matplotlib.pyplot as plt

def plot_dataset(data: list[tuple[float,...]], path):
    xs = [p[0] for p in data]
    ys = [p[1] for p in data]

    plt.figure(figsize=(6, 6))
    plt.scatter(xs, ys, s=5)
    plt.gca().set_aspect("equal", adjustable="box")
    
    plt.xlim(min(xs)-0.05, max(xs)+0.05)
    plt.ylim(min(ys)-0.05, max(ys)+0.05)

    plt.savefig(path, dpi=400)
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

    plt.savefig(path, dpi=400)
    plt.close()