import random
import math

def generate_spd_matrix(n):
    A = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            A[i][j] = random.uniform(-100, 100)
            A[j][i] = A[i][j]  # Make the matrix symmetric
    for i in range(n):
        A[i][i] += n  # Ensure diagonal dominance for positive definiteness
    return A

def generate_diagonally_dominant_matrix(n):
    A = [[random.uniform(-100, 100) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        diagonal_sum = sum(abs(A[i][j]) for j in range(n) if j != i)
        A[i][i] = diagonal_sum + 1  # Adjust diagonal elements for dominance
    return A

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
pivot_type = "complete"

# User input for matrix type
matrix_type = "spd"

# Define the size of the matrix
n = 20

# Initialize lists to store errors
factorization_error_1norm_list = []
factorization_error_frobenius_list = []
relative_error_one_norm_list = []
relative_error_two_norm_list = []
relative_error_one_norm_residual_list = []
relative_error_two_norm_residual_list = []

# Repeat the process 10 times
for _ in range(10):
    if matrix_type == 'spd':
        A = generate_spd_matrix(n)
    elif matrix_type == 'diagonally_dominant':
        A = generate_diagonally_dominant_matrix(n)
    else:
        raise ValueError("Invalid matrix_type. Use 'spd' or 'diagonally_dominant'.")

    x = [random.uniform(-100, 100) for _ in range(n)]

    P, Q, LU = perform_lu_decomposition(A, pivot_type)

    PAQ_minus_LU = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            PAQ_minus_LU[i][j] = A[P[i]][Q[j]] - LU[i][j]

    factorization_error_1norm = one_norm(PAQ_minus_LU) / one_norm(A)
    factorization_error_frobenius = frobenius_norm(PAQ_minus_LU) / frobenius_norm(A)

    b = [sum(A[i][j] * x[j] for j in range(n)) for i in range(n)]
    y = forward_substitution(LU, b, P)
    x_hat = backward_substitution(LU, y)

    one_norm_x = vector_one_norm(x)
    two_norm_x = two_norm(x)

    one_norm_x_hat = vector_one_norm(x_hat)
    two_norm_x_hat = two_norm(x_hat)

    relative_error_one_norm = one_norm_x_hat / one_norm_x
    relative_error_two_norm = two_norm_x_hat / two_norm_x

    Ax_hat = [sum(A[i][j] * x_hat[j] for j in range(n)) for i in range(n)]
    r = [b[i] - Ax_hat[i] for i in range(n)]

    one_norm_r = vector_one_norm(r)
    two_norm_r = two_norm(r)

    one_norm_b = vector_one_norm(b)
    two_norm_b = two_norm(b)

    relative_error_one_norm_residual = one_norm_r / one_norm_b
    relative_error_two_norm_residual = two_norm_r / two_norm_b

    factorization_error_1norm_list.append(factorization_error_1norm)
    factorization_error_frobenius_list.append(factorization_error_frobenius)
    relative_error_one_norm_list.append(relative_error_one_norm)
    relative_error_two_norm_list.append(relative_error_two_norm)
    relative_error_one_norm_residual_list.append(relative_error_one_norm_residual)
    relative_error_two_norm_residual_list.append(relative_error_two_norm_residual)

# Calculate the averages
average_factorization_error_1norm = sum(factorization_error_1norm_list) / len(factorization_error_1norm_list)
average_factorization_error_frobenius = sum(factorization_error_frobenius_list) / len(factorization_error_frobenius_list)
average_relative_error_one_norm = sum(relative_error_one_norm_list) / len(relative_error_one_norm_list)
average_relative_error_two_norm = sum(relative_error_two_norm_list) / len(relative_error_two_norm_list)
average_relative_error_one_norm_residual = sum(relative_error_one_norm_residual_list) / len(relative_error_one_norm_residual_list)
average_relative_error_two_norm_residual = sum(relative_error_two_norm_residual_list) / len(relative_error_two_norm_residual_list)

print("\nAverage Factorization Error (1-norm):")
print(average_factorization_error_1norm)

print("\nAverage Factorization Error (Frobenius norm):")
print(average_factorization_error_frobenius)

print("\nAverage Relative Error using One-Norm:")
print(average_relative_error_one_norm)

print("\nAverage Relative Error using Two-Norm:")
print(average_relative_error_two_norm)

print("\nAverage Relative Error via residual using One-Norm:")
print(average_relative_error_one_norm_residual)

print("\nAverage Relative Error via residual using Two-Norm:")
print(average_relative_error_two_norm_residual)
