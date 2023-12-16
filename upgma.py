import numpy as np
import string
from typing import List
from pylocluster import *
from newick import loads

pure_distance_matrix = [[0.0000000, 2.2266620, 3.3353841, 3.0926941, 1.8896296, 4.9589915, 4.3943885, 5.4122139],
                   [2.2266620, 0.0000000, 1.1567076, 0.9140176, 1.3807486, 4.4501105, 3.8855074, 4.9033329],
                   [3.3353841, 1.1567076, 0.0000000, 1.8198296, 2.4894707, 5.5588327, 4.9942296, 6.0120551],
                   [3.0926941, 0.9140176, 1.8198296, 0.0000000, 2.2467807, 5.3161427, 4.7515396, 5.7693651],
                   [1.8896296, 1.3807486, 2.4894707, 2.2467807, 0.0000000, 3.3471233, 2.7825202, 3.8003457],
                   [4.9589915, 4.4501105, 5.5588327, 5.3161427, 3.3471233, 0.0000000, 1.1966232, 2.2144487],
                   [4.3943885, 3.8855074, 4.9942296, 4.7515396, 2.7825202, 1.1966232, 0.0000000, 1.2176200],
                   [5.4122139, 4.9033329, 6.0120551, 5.7693651, 3.8003457, 2.2144487, 1.2176200, 0.0000000]]

class Coordinate:
    def __init__(self, label: str, value: int) -> None:
        self.value = value
        self.label = label

    def __eq__(self, other):
        if isinstance(other, Coordinate):
            return self.value == other.get_value()
        return False

    def get_value(self) -> int:
        return self.value
    
    def get_label(self) -> str:
        return self.label
    
    def set_value(self, value: int) -> None:
        self.value = value

    def from_new_label(self, new_label) -> 'Coordinate':
        self.label = new_label
        return self

class Distance:
    def __init__(self, value: float, deep: int, coordinate_x: Coordinate, coordinate_y: Coordinate) -> None:
        self.value = value
        self.deep = deep
        self.coordinate_x = coordinate_x
        self.coordinate_y = coordinate_y

    def __eq__(self, other):
        if isinstance(other, Distance):
            return (self.coordinate_x == other.get_coordinate_x() and self.coordinate_y == other.get_coordinate_y()) or (
                self.coordinate_y == other.get_coordinate_x() and self.coordinate_x == other.get_coordinate_y()
            )
        return False

    def get_value(self) -> float:
        return self.value
    
    def get_deep(self) -> int:
        return self.deep
    
    def get_coordinate_x(self) -> Coordinate:
        return self.coordinate_x
    
    def get_coordinate_y(self) -> Coordinate:
        return self.coordinate_y

class DistanceRow:
    def __init__(self, coordinate_x: Coordinate) -> None:
        self.row: List[Distance] = []
        self.size = 0
        self.coordinate_x = coordinate_x
            
    def append(self, item: Distance) -> None:
        self.row.append(item)
        self.size = self.size + 1

    def get_as_list(self) -> List[Distance]:
        return self.row

    def get_item(self, y: int) -> Distance:
        return self.row[y]
    
    def get_size(self) -> int:
        return self.size
    
    def get_coordinate_x(self) -> Coordinate:
        return self.coordinate_x
    
    def is_empty(self) -> bool:
        return self.size == 0


class DistanceMatrix:
    def __init__(self) -> None:
        self.matrix: List[DistanceRow] = []
        self.size = 0
    
    def append(self, row: DistanceRow) -> None:
        self.matrix.append(row)
        self.size = self.size + 1

    def get_as_list(self) -> List[DistanceRow]:
        return self.matrix
    
    def get_row(self, x: int) -> DistanceRow:
        return self.matrix[x]
    
    def get_item(self, x: int, y: int) -> Distance:
        return self.matrix[x].get_item(y)
    
    def get_distance_by_coordenate_value(self, x: int, y: int) -> Distance:
        if x < y:
            return self.get_row(x).get_item(y - (x + 1))
        else:
            return self.get_row(y).get_item(x - (y + 1))
    
    def get_distance_by_coordenate(self, coordinate_x: Coordinate, coordinate_y: Coordinate) -> Distance:
        return self.get_distance_by_coordenate_value(coordinate_x.get_value(), coordinate_y.get_value())
    
    def get_size(self) -> int:
        return self.size
    
    def get_min_distance(self) -> Distance:
        invalid_coordenate = Coordinate('invalid', 0)
        min_distance = Distance(float('inf'), 1, invalid_coordenate, invalid_coordenate)

        for distance_row in self.get_as_list():
            for distance in distance_row.get_as_list():
                if distance.get_value() < min_distance.get_value():
                    min_distance = distance

        return min_distance
    
    def print(self) -> None:
        matrix_size = self.get_size() + 1
        readable_distance_matrix = np.zeros((matrix_size, matrix_size))
        
        for distance_row in self.get_as_list():
            for distance in distance_row.get_as_list():
                readable_distance_matrix[distance.get_coordinate_x().get_value()][distance.get_coordinate_y().get_value()] = distance.get_value()
                readable_distance_matrix[distance.get_coordinate_y().get_value()][distance.get_coordinate_x().get_value()] = distance.get_value()

        print(readable_distance_matrix) 
    
    def is_final(self):
        return self.size == 1
    
def get_label(index: int):
    alphabetical_order = list(string.ascii_lowercase)
    index += 1
    label = ''
    while index > 0:
        index -= 1
        current_index = index % len(alphabetical_order)
        label = alphabetical_order[current_index] + label
        index //= len(alphabetical_order)

    return label.capitalize()

def get_label_from_combination(distance: Distance) -> str:
    labels = [distance.get_coordinate_x().get_label(), distance.get_coordinate_y().get_label()]
    labels.sort()
    return '(' + labels[0] + ',' + labels[1] + ')'

def validate_pure_distance_matrix(pure_distance_matrix) -> None:
    size = len(pure_distance_matrix)
    for x in range(size):
        for y in range(size):
            if pure_distance_matrix[x][y] != pure_distance_matrix[y][x]:
                raise Exception("Matrix is not symmetric")

def initilize_distance_matrix(pure_distance_matrix) -> DistanceMatrix:
    validate_pure_distance_matrix(pure_distance_matrix)

    matrix_size = len(pure_distance_matrix)
    distance_matrix = DistanceMatrix()

    coordinates = []

    for x in range(0, matrix_size):
        coordinates.append(Coordinate(get_label(x), x))

    first_column = 1
    for x in range(0, matrix_size - 1):
        distance_row = DistanceRow(coordinates[x])
        for y in range(first_column, matrix_size):
            if x != y:
                distance_row.append(Distance(pure_distance_matrix[x][y], 1, coordinates[x], coordinates[y]))
        
        distance_matrix.append(distance_row)
        first_column = first_column + 1

    return distance_matrix

def are_distances_neighbor(distanceOne: Distance, distanceTwo: Distance) -> bool:
    return (
        distanceOne.get_coordinate_x().get_value() == distanceTwo.get_coordinate_x().get_value() or distanceOne.get_coordinate_x().get_value() == distanceTwo.get_coordinate_y().get_value()
    ) or ( 
        distanceOne.get_coordinate_y().get_value() == distanceTwo.get_coordinate_x().get_value() or distanceOne.get_coordinate_y().get_value() == distanceTwo.get_coordinate_y().get_value()
    )

def get_not_shared_coordinate(distanceOne: Distance, distanceTwo: Distance) -> Coordinate:
    is_first_x_not_equal_second_x = distanceOne.get_coordinate_x().get_value() != distanceTwo.get_coordinate_x().get_value()
    is_first_x_not_equal_second_y = distanceOne.get_coordinate_x().get_value() != distanceTwo.get_coordinate_y().get_value()
    is_coordinate_x_not_shared = is_first_x_not_equal_second_x and is_first_x_not_equal_second_y
    if is_coordinate_x_not_shared:
        return distanceOne.get_coordinate_x()
    else:
        return distanceOne.get_coordinate_y()

def update_distance_matrix(distance_matrix: DistanceMatrix) -> DistanceMatrix:
    min_distance = distance_matrix.get_min_distance()

    min_distance_label = get_label_from_combination(min_distance)

    new_distance_matrix = DistanceMatrix()

    row_and_column_to_ignore = min_distance.get_coordinate_y().get_value()
    
    new_x = 0

    for distance_row in distance_matrix.get_as_list():
        new_distance_row = DistanceRow(Coordinate(distance_row.get_coordinate_x().get_label(), new_x))
        new_y = new_x + 1

        if distance_row.get_coordinate_x().get_value() == row_and_column_to_ignore:
            continue

        for distance in distance_row.get_as_list():
            if distance.get_coordinate_y().get_value() == row_and_column_to_ignore:
                continue

            if are_distances_neighbor(distance, min_distance):
                neighbor_coordinate_from_min_distance = get_not_shared_coordinate(distance, min_distance)

                distance_from_x = distance_matrix.get_distance_by_coordenate(neighbor_coordinate_from_min_distance, min_distance.get_coordinate_x())
                distance_from_y = distance_matrix.get_distance_by_coordenate(neighbor_coordinate_from_min_distance, min_distance.get_coordinate_y())

                distance_with_weight_from_x = distance_from_x.get_value() * distance_from_x.get_deep()
                distance_with_weight_from_y = distance_from_y.get_value() * distance_from_y.get_deep()

                distance_divider = distance_from_x.get_deep() + distance_from_y.get_deep();

                new_deep = max(distance_from_x.get_deep(), distance_from_y.get_deep()) + 1
                new_distance_value = (distance_with_weight_from_x + distance_with_weight_from_y) / distance_divider

                if distance.get_coordinate_x().get_value() < min_distance.get_coordinate_x().get_value():
                    new_coordinate_y = Coordinate(min_distance_label, new_y)
                    new_coordinate_x = Coordinate(distance.get_coordinate_x().get_label(), new_x)
                else:
                    new_coordinate_y = Coordinate(distance.get_coordinate_y().get_label(), new_y)
                    new_coordinate_x = Coordinate(min_distance_label, new_x)
                
                new_distance_row.append(Distance(new_distance_value, new_deep, new_coordinate_x, new_coordinate_y))

                new_y = new_y + 1
            else:
                new_coordinate_x = Coordinate(distance.get_coordinate_x().get_label(), new_x)
                new_coordinate_y = Coordinate(distance.get_coordinate_y().get_label(), new_y)
                new_distance_row.append(Distance(distance.get_value(), distance.get_deep(), new_coordinate_x, new_coordinate_y))
                new_y = new_y + 1
                
        if not new_distance_row.is_empty():
            new_x = new_x + 1
            new_distance_matrix.append(new_distance_row)
    
    return new_distance_matrix

def get_result(distance_matrix: DistanceMatrix) -> str:
    if not distance_matrix.is_final():
        raise Exception('Extracting result from non final matrix')
    
    final_distance = distance_matrix.get_distance_by_coordenate_value(0, 1)

    return get_label_from_combination(final_distance)

def upgma(pure_distance_matrix: list[list[float]]):
    distance_matrix = initilize_distance_matrix(pure_distance_matrix)
    distance_matrix.print()
    print('\n')

    while not distance_matrix.is_final():
        distance_matrix = update_distance_matrix(distance_matrix)
        distance_matrix.print()
        print('\n')

    print(get_result(distance_matrix))

upgma(pure_distance_matrix)

nwk = linkage(pure_distance_matrix, taxa=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'], method='upgma')
'((H:0.86,(F:0.60,G:0.60):0.26):1.49,((C:0.74,(B:0.46,D:0.46):0.29):0.49,(A:0.94,E:0.94):0.29):1.11);'
print(loads(nwk).ascii_art())