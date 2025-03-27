from array import array

import numpy as np

def maxvalue_needleman(matrix, x_pos, y_pos, match_score=1, mismatch_penalty=-1, gap_penalty=-1):

    if matrix[x_pos][0] == matrix[0][y_pos]:
        diagonal = matrix[x_pos - 1][y_pos - 1] + match_score
    else:
        diagonal = matrix[x_pos - 1][y_pos - 1] + mismatch_penalty

    left = matrix[x_pos - 1][y_pos] + gap_penalty
    top = matrix[x_pos][y_pos - 1] + gap_penalty

    matrix[x_pos][y_pos] = max(left, top, diagonal)

strand1 = input("Enter your strand 1: ").upper()
strand2 = input("Enter your strand 2: ").upper()

row_length = len(strand1) + 2
column_length = len(strand2) + 2

bases = {'A', 'C', 'G', 'T'}
if set(strand1).issubset(bases) and set(strand2).issubset(bases):
    print("Your strands are valid.")
else:
    print("Invalid bases detected!")

matrix = np.zeros((row_length, column_length), dtype='object')

matrix[1, 0] = ''
matrix[0, 0] = ''
matrix[0, 1] = ''
matrix[1, 0] = ''
matrix[1, 1] = 0


matrix[0, 2:2+len(strand2)] = list(strand2)
matrix[2:2+len(strand1), 0] = list(strand1)

n = -1
for i in range(2, row_length):
    matrix[i, 1] = n
    n -= 1

n = -1
for i in range(2, column_length):
    matrix[1, i] = n
    n -= 1

for i in range(2, row_length):
    for j in range(2, column_length):
        maxvalue_needleman(matrix, i, j)

print(matrix)
