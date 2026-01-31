point_component_dict = {"a" : 1, "b" : 1, "c" : 2, "d" : 1}
component_point_dict = {1 : {"a", "b", "d"}, 2 : {"c"}}

e = "e"
assert type(e) == str
comp_a = point_component_dict["a"]
# print(comp_a)
point_component_dict[e] = comp_a
component_point_dict[comp_a].add(e)

# print(point_component_dict)
# print(component_point_dict)


# new_point:
# compare with each old point
# save all clusters that are near
# merge clusters and 

type point = str

density_dict = {"a" : 2, "b" : 1, "c" : 1}
inverse_density_dict = {1 : {"b" , "c"}, 2 : {"a"}}
densities = set(inverse_density_dict.keys())
# triviale Fälle ausschließen?
max_density = max(densities)

class Cluster:
    def __init__(self, points : set[point], hidden : bool = True) -> None:
        self._points : set[point] = points
        self._max_density : int = max(*[density_dict[pt] for pt in points]) if len(points) > 1 else density_dict[next(iter(points))]
        self._hidden : bool = hidden

    def __str__(self) -> str:
        return f"Cluster with \npoints: {self._points},\nmax_density: {self._max_density},\nhidden: {self._hidden}"

    def get_points(self) -> set[point]:
        return self._points
    
    def get_max_density(self) -> int:
        return self._max_density

    def is_hidden(self) -> bool:
        return self._hidden

    def add_point(self, pt : point) -> None:
        self._points.add(pt)

    def merge(self, *merge_clusters : 'Cluster', single_point : point | None = None) -> None:
        self._points.update(*[cluster.get_points() for cluster in merge_clusters])
        max_density = max(*[cluster.get_max_density() for cluster in merge_clusters], self._max_density)
        if single_point != None:
            self._points.add(single_point)
            max_density = max(max_density, density_dict[single_point])
        self._max_density = max_density
    
    # unnnötig?
    def set_hidden(self) -> None:
        self._hidden = True
    
    def set_visible(self) -> None:
        self._hidden = False

tau_eff = 1

def tauDistanceTest(pt1, pt2) -> bool:
    for coord1, coord2 in zip(pt1, pt2):
        if abs(coord1 - coord2) > tau_eff:
            return False
    return True

new_point = (1,2)

# lohnt das sich?:
# point_component_dict = {}

clusters = set()
used_points = set()

for rho_bar in range(max_density, 0, -1):
    
    new_layer = inverse_density_dict[rho_bar]

    modified_clusters = []

    for new_point in new_layer:
        
        # near_points =
        # near_clusters = 
        # first_cluster = 

        # alter, zu modifizierender Code:
        tau_connection_graph = {}
        old_points = set()
        for new_point in surviving_lattice_points:
            connected_points = {point for point in old_points if tauDistanceTest(point, new_point)}
            tau_connection_graph[new_point] = connected_points
            for point in connected_points:
                tau_connection_graph[point].add(new_point)
            old_points.add(new_point)


        # Falls mehr als eins: first_cluster.merge(near_clusters, single_point=new_point)
        # Falls genau eins: first_cluster.add_point(new_point)
        # Sonst: new_cluster = Cluster({n})
        # 

        used_points.add(new_point)





class Pointer:
    def __init__(self, target):
        self.__target = target

    @property
    def target(self):
        return self.__target

    def change_target(self, target):
        self.__target = target