############################################################################################
# Cluster-Klasse 
############################################################################################

type point = tuple[int,...]

class Cluster:
    def __init__(self, points : set[point], densities : set[int], visible : bool = False) -> None:
        self.points : set[point] = points
        self.max_density : int = max(*list(densities)) if len(densities) > 1 else next(iter(densities))
        self.visible : bool = visible

    def __str__(self) -> str:
        return f"Cluster-Object with the following fields\npoints: {self.points},\nmax_density: {self.max_density},\nvisible: {self.visible}\n"

    def add_point(self, pt : point, density : int) -> None:
        """Adding a single Point to the cluster"""
        self.points.add(pt)
        self.max_density = max(self.max_density, density)

    def merge(self, *merge_clusters : 'Cluster') -> None:
        """Merging an iterable of Clustern into this Cluster"""
        self.points.update(*[cluster.points for cluster in merge_clusters])
        self.max_density = max(*[cluster.max_density for cluster in merge_clusters], self.max_density)

    def update_visibility(self, visible) -> None:
        """changing the visibility flag of the Cluster"""
        self.visible = visible