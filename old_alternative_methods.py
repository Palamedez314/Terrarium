
# alter Ansatz von connectedComponents ohne Graph direkt durch 

def connectedComponents(step_count: int) -> list[set[tuple[int,...]]]:
    connected_component_list : list[set[tuple[int,...]]] = []

    for lattice_point in survivingLatticePoints(step_count):

        neighbors = {tuple(map(sum, zip(lattice_point, diff))) for diff in tau_distance_set}
        # TODO Schaut sich der Schnitt alle Elemente an -> ja! geht das schneller?
        # TODO Schnitt besser, nicht die ganze Menge berechnen -> educated guess welche Methode besser ist
        neighboring_component_indices = [i for i in range(len(connected_component_list)) if 
                                         bool(neighbors & connected_component_list[i])]

        if neighboring_component_indices:
            min_index = neighboring_component_indices.pop(0)
            final_component = connected_component_list[min_index]

            for i in reversed(neighboring_component_indices):
                final_component.update(connected_component_list.pop(i))

            final_component.add(lattice_point)
            
        else:
            connected_component_list.append({lattice_point})

    surviving_list : list[bool] = [epsDensityTest(component,step_count) 
                                   for component in connected_component_list]
    surviving_components = [component for component, survived in zip(connected_component_list,surviving_list) if survived]
    dead_components : set[tuple[int,...]] = set.union(*[component for component, survived 
                                                        in zip(connected_component_list,surviving_list) 
                                                        if not survived]) if False in surviving_list else set()
    
    surviving_components.insert(0,dead_components)

    return surviving_components