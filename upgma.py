import numpy as np
from typing import List

pure_distance_matrix = [[0, 17, 21, 31, 23],
                        [17, 0, 30, 34, 21],
                        [21, 30, 0, 28, 39],
                        [31, 34, 28, 0, 43],
                        [23, 21, 39, 43, 0]]

class Distance:
    def __init__(self, value: float, deep: int, x: int, y: int) -> None:
        self.value = value
        self.deep = deep
        self.pure_x = x
        self.pure_y = y

    def get_value(self) -> float:
        return self.value
    
    def get_deep(self) -> int:
        return self.deep
    
    def get_pure_x(self) -> int:
        return self.pure_x
    
    def get_pure_y(self) -> int:
        return self.pure_y
    
    def get_min_pure_coordenate(self) -> int:
        if self.pure_x < self.pure_y:
            return self.pure_x
        else:
            return self.pure_y

class DistanceRow:
    def __init__(self, x: int) -> None:
        self.row: List[Distance] = []
        self.size = 0
        self.pure_x = x
            
    def append(self, item: Distance) -> None:
        self.row.append(item)
        self.size = self.size + 1

    def get_as_list(self) -> List[Distance]:
        return self.row

    def get_item(self, y: int) -> Distance:
        return self.row[y]
    
    def get_size(self) -> int:
        return self.size
    
    def get_pure_x(self) -> int:
        return self.pure_x

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
    
    def get_item_by_pure_coordenate(self, x: int, y: int) -> Distance:
        if x < y:
            return self.get_row(x).get_item(y - (x + 1))
        else:
            return self.get_row(y).get_item(x - (y + 1))
    
    def get_size(self) -> int:
        return self.size
    
def validate_pure_distance_matrix(pure_distance_matrix) -> None:
    size = len(pure_distance_matrix)
    for x in range(size):
        for y in range(size):
            if pure_distance_matrix[x][y] != pure_distance_matrix[y][x]:
                raise Exception("Matrix is not symmetric")

def initilize_distance_matrix(pure_distance_matrix) -> DistanceMatrix:
    matrix_size = len(pure_distance_matrix)
    distance_matrix = DistanceMatrix()

    first_column = 0
    for x in range(0, matrix_size):
        distance_row = DistanceRow(x)
        for y in range(first_column, matrix_size):
            if x != y:
                distance_row.append(Distance(pure_distance_matrix[x][y], 0, x, y))
        
        distance_matrix.append(distance_row)
        first_column = first_column + 1

    return distance_matrix

def print_distance_matrix(distance_matrix: DistanceMatrix) -> None:
    matrix_size = distance_matrix.get_size()
    readable_distance_matrix = np.zeros((matrix_size, matrix_size))
    
    for distance_row in distance_matrix.get_as_list():
        for distance in distance_row.get_as_list():
            readable_distance_matrix[distance.get_pure_x()][distance.get_pure_y()] = distance.get_value()
            readable_distance_matrix[distance.get_pure_y()][distance.get_pure_x()] = distance.get_value()

    print(readable_distance_matrix)

def find_min_distance(distance_matrix: DistanceMatrix) -> Distance:
    min_distance = Distance(float('inf'), 0, 0, 0)

    for distance_row in distance_matrix.get_as_list():
        for distance in distance_row.get_as_list():
            if distance.get_value() < min_distance.get_value():
                min_distance = distance

    return min_distance

def update_distance_matrix(distance_matrix: DistanceMatrix, min_distance: Distance) -> DistanceMatrix:
    new_distance_matrix = DistanceMatrix()
    new_x = 0

    for distance_row in distance_matrix.get_as_list():
        new_distance_row = DistanceRow(new_x)
        new_y = 0
        first_distance_from_min_distance: Distance = None
        for distance in distance_row.get_as_list():
            if distance.get_pure_x() != min_distance.get_pure_x() and distance.get_pure_y() != min_distance.get_pure_y():
                new_distance_row.append(Distance(distance.get_value(), distance.get_deep(), new_x, new_y))
                new_y = new_y + 1
            elif first_distance_from_min_distance == None:
                first_distance_from_min_distance = distance
            else:
                new_distance_value = 0 # TODO Calculate new distance based on first_distance_from_min_distance and distance / 2
                new_distance_row.append(Distance(new_distance_value, distance.get_deep(), new_x, new_y))
                new_y = new_y + 1
        new_x = new_x + 1
        new_distance_matrix.append(new_distance_row)
    
    return new_distance_matrix

validate_pure_distance_matrix(pure_distance_matrix)

distance_matrix = initilize_distance_matrix(pure_distance_matrix)

print_distance_matrix(distance_matrix)

min_distance = find_min_distance(distance_matrix)

distance_matrix = update_distance_matrix(distance_matrix, min_distance)

# print(min_distance_pair)

# distance_matrix = calculate_new_distance_matrix(min_distance_pair)

# print(distance_matrix)

# min_distance_pair = find_min_distance_pair()

# print(min_distance_pair)

# distance_matrix = calculate_new_distance_matrix(min_distance_pair)

# print(distance_matrix)