import matplotlib.pyplot as plt

def plot_dataset(data, datasetname: str):
    points = [point[0] for point in data]

    plt.figure(figsize=(6, 6))
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    plt.gca().set_aspect("equal", adjustable="box")

    xs = [p[0] for p in points]
    ys = [p[1] for p in points]

    plt.scatter(xs, ys, c="violet", s=40)
    plt.savefig(f"team-12-{datasetname}.train.png")
    plt.show()

def plot_clusters(data, datasetname: str):
    points = [point[0] for point in data]
    labels = [point[1] for point in data]

    plt.figure(figsize=(6, 6))
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    plt.gca().set_aspect("equal", adjustable="box")

    xs = [p[0] for p in points]
    ys = [p[1] for p in points]

    plt.scatter(xs, ys, c=labels, s=40, cmap="tab10")
    plt.savefig(f"team-12-{datasetname}.result.png")
    plt.show()


example = [
    [(0.8, 0.5), 1],
    [(0.5, 1), 1],
    [(-0.7, 0), 2],
    [(-0.6, -0.4), 2],
    [(0, 0), 0],       
    [(1, -1), 0]        ]

plot_clusters(example, "example_dataset")