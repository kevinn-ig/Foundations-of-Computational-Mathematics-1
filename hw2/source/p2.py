import random as rnd
import numpy as np


import random
import math

def one_norm(matrix):
    norm = 0
    for col in range(len(matrix[0])):
        column_sum = sum(abs(row[col]) for row in matrix)
        norm = max(norm, column_sum)
    return norm

def frobenius_norm(matrix):
    return math.sqrt(sum(sum(x**2 for x in row) for row in matrix))

def vector_one_norm(vector):
    return sum(abs(x) for x in vector)

def two_norm(vector):
    return math.sqrt(sum(x**2 for x in vector))

def lu_decomposition_no_pivot(A):
    n = len(A)
    P = list(range(n))
    Q = list(range(n))

    for k in range(n):
        for i in range(k + 1, n):
            A[i][k] /= A[k][k]
            for j in range(k + 1, n):
                A[i][j] -= A[i][k] * A[k][j]

    return P, Q, A

def lu_decomposition_partial_pivot(A):
    n = len(A)
    P = list(range(n))
    Q = list(range(n))

    for k in range(n):
        pivot_row = max(range(k, n), key=lambda i: abs(A[i][k]))
        if pivot_row != k:
            A[k], A[pivot_row] = A[pivot_row], A[k]
            P[k], P[pivot_row] = P[pivot_row], P[k]

        for i in range(k + 1, n):
            A[i][k] /= A[k][k]
            for j in range(k +1 , n):
                A[i][j] -= A[i][k] * A[k][j]

    return P, Q, A

def lu_decomposition_complete_pivot(A):
    n = len(A)
    P = list(range(n))
    Q = list(range(n))

    for k in range(n):
        pivot_element = max((abs(A[i][j]), i, j) for i in range(k, n) for j in range(k, n))
        i, j = pivot_element[1], pivot_element[2]

        if i != k:
            A[k], A[i] = A[i], A[k]
            P[k], P[i] = P[i], P[k]

        if j != k:
            Q[k], Q[j] = Q[j], Q[k]

        for i in range(k + 1, n):
            A[i][k] /= A[k][k]
            for j in range(k + 1, n):
                A[i][j] -= A[i][k] * A[k][j]

    return P, Q, A

def forward_substitution(LU, b, P):
    n = len(LU)
    y = [0.0] * n

    b = [b[P[i]] for i in range(n)]

    for i in range(n):
        y[i] = b[i]
        for j in range(i):
            y[i] -= LU[i][j] * y[j]

    return y

def backward_substitution(LU, y):
    n = len(LU)
    x = [0.0] * n

    for i in range(n - 1, -1, -1):
        x[i] = y[i]
        for j in range(i + 1, n):
            x[i] -= LU[i][j] * x[j]
        x[i] /= LU[i][i]

    return x

def perform_lu_decomposition(A, pivot_type):
    if pivot_type == "partial":
        P, Q, LU = lu_decomposition_partial_pivot(A)
    elif pivot_type == "complete":
        P, Q, LU = lu_decomposition_complete_pivot(A)
    elif pivot_type == "none":
        P, Q, LU = lu_decomposition_no_pivot(A)
    else:
        raise ValueError("Invalid pivot_type. Use 'partial', 'complete', or 'none'.")

    return P, Q, LU

# User input for pivot type
pivot_type = "partial"


# Define the size of the matrix
n = 3

A = [[2,1,0],[-4,0,4],[2,5,-10]]
b = [3,0,17]

P, Q, LU = perform_lu_decomposition(A, pivot_type)

PAQ_minus_LU = [[0.0] * n for _ in range(n)]
for i in range(n):
    for j in range(n):
        PAQ_minus_LU[i][j] = A[P[i]][Q[j]] - LU[i][j]

y = forward_substitution(LU, b, P)
x_hat = backward_substitution(LU, y)

Ax_hat = [sum(A[i][j] * x_hat[j] for j in range(n)) for i in range(n)]

print("Partial pivoting\n")
print(A,P,x_hat)

pivot_type = "complete"

A = [[2,1,0],[-4,0,4],[2,5,-10]]
b = [3,0,17]

P, Q, LU = perform_lu_decomposition(A, pivot_type)

PAQ_minus_LU = [[0.0] * n for _ in range(n)]
for i in range(n):
    for j in range(n):
        PAQ_minus_LU[i][j] = A[P[i]][Q[j]] - LU[i][j]

y = forward_substitution(LU, b, P)
x_hat = backward_substitution(LU, y)

Ax_hat = [sum(A[i][j] * x_hat[j] for j in range(n)) for i in range(n)]
print("\nComplete pivoting\n")
print(A,P,Q,x_hat)

