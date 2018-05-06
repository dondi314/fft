def transpose(matrix):
    new_matrix = []
    for i in range(len(matrix[0])):
        row = [j[i] for j in matrix]
        new_matrix.append(row)
    return new_matrix
