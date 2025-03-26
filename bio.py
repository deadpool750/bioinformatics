import numpy as np
from numpy.matrixlib.defmatrix import matrix


def maxvalue_needleman(matrix,x_value,y_value):

    #matrix x value = x_value
    #matrix y value = y value

    match_score = 1
    mismatch_penalty = -1
    gap_penalty = -1

    if matrix[x_value][0] == matrix[0][y_value]:
        diagonal = matrix[x_value - 1][y_value - 1] + match_score
    else: diagonal = matrix[x_value - 1][y_value - 1] + mismatch_penalty

    left = matrix[x_value - 1][y_value] + gap_penalty
    top = matrix[x_value][y_value + 1] + gap_penalty

    matrix[x_value][y_value] = max(left, top, diagonal)

strand1 = input("Enter your strand 1: ").upper()
strand2 = input("Enter your strand 2: ").upper()

row_1_length = len(strand1) + 2
row_2_length = len(strand2) + 2

bases = {'A', 'C', 'G', 'T'}
if set(strand1).issubset(bases) and set(strand2).issubset(bases):
    print("Your strands are valid.")
else:
    print("Invalid bases detected!")

matrix = np.zeros((row_1_length, row_2_length), dtype='object')

matrix[0, 2:] = list(strand2)
matrix[2:, 0] = list(strand1)

matrix[0, 0] = '-'
matrix[0, 1] = '-'
matrix[1, 0] = '-'
matrix[1, 1] = 0

n = -1
for i in range(2, row_1_length):
    matrix[i, 1] = n
    n -= 1

n = -1
for i in range(2, row_2_length):
    matrix[1, i] = n
    n -= 1

#maxvalue_needleman(matrix, 2,2)

for i in range(2, row_1_length - 1):
    for j in range(2, row_2_length - 1):
        maxvalue_needleman(matrix,i,j)

print(matrix)
