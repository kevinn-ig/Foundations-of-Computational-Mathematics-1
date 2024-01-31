import random
import math

def calculate_l2_norm(vector):
    max_value = max(abs(x) for x in vector)
    normalized_vector = [x / max_value for x in vector]
    return max_value * math.sqrt(sum(x**2 for x in normalized_vector))

def norm_ratio(r, b):
    assert len(r) == len(b)
    l2_norm_r = calculate_l2_norm(r)
    l2_norm_b = calculate_l2_norm(b)
    ratio = l2_norm_r / l2_norm_b
    return ratio

def identity_matrix(n):
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]

def extract_diagonal(matrix):
    n = len(matrix)
    return [[matrix[i][i] if i == j else 0 for j in range(n)] for i in range(n)]


def solve_system(P, r):
    P_inverse = invert_matrix(P)
    z = matrix_vector_multiply(P_inverse, r)
    return z

def matrix_vector_multiply(matrix, vector):
    return [sum(a * b for a, b in zip(row, vector)) for row in matrix]

def invert_matrix(matrix):
    n = len(matrix)
    identity = identity_matrix(n)

    for i in range(n):
        pivot = matrix[i][i]
        if pivot == 0:
            raise ValueError("Matrix is singular and cannot be inverted")

        for j in range(n):
            matrix[i][j] /= pivot
            identity[i][j] /= pivot

        for k in range(n):
            if k != i:
                factor = matrix[k][i]
                for j in range(n):
                    matrix[k][j] -= factor * matrix[i][j]
                    identity[k][j] -= factor * identity[i][j]

    return identity

def lower_triangular_matrix(matrix):
    n = len(matrix)
    return [[matrix[i][j] if i >= j else 0 for j in range(n)] for i in range(n)]

def multiply_lower_triangular_transpose(lower_triangular_matrix):
    n = len(lower_triangular_matrix)
    result = [[0] * n for _ in range(n)]

    for i in range(n):
        for j in range(n):
            for k in range(j + 1):
                result[i][j] += lower_triangular_matrix[i][k] * lower_triangular_matrix[j][k]

    return result

def steepest_descent(P, A, b, x):
    Ax = matrix_vector_multiply(A, x)
    r = [bi - ai for ai, bi in zip(Ax, b)]
    z = solve_system(P, r)

    w = []
    k = -1
    Ax = matrix_vector_multiply(A, x)
    r_norm = [bi - ai for ai, bi in zip(Ax, b)]


    while calculate_l2_norm(r_norm) > 1e-6:
       
        w = matrix_vector_multiply(A, z)

        z_transpose = [z]
        

        # Handle the case where the denominator is close to zero
        alpha = sum(a * b for a, b in zip(z_transpose[0], r)) / sum(a * b for a, b in zip(z_transpose[0], w))
        
        x = [a + alpha * b for a, b in zip(x, z)]
        r = [a - alpha * b for a, b in zip(r, w)]
        z = solve_system(P, r)
        k += 1
        Ax = matrix_vector_multiply(A, x)
        r_norm = [bi - ai for ai, bi in zip(Ax, b)]

    return x, k


def cg(P, A, b, x):
    Ax = matrix_vector_multiply(A, x)
    r = [bi - ai for ai, bi in zip(Ax, b)]
    z = solve_system(P, r)
    p = z
    v = []
    k = -1
    Ax = matrix_vector_multiply(A, x)
    r_norm = [bi - ai for ai, bi in zip(Ax, b)]

    while calculate_l2_norm(r_norm) > 1e-6:
 
        v = matrix_vector_multiply(A, p)
        r_transpose = [r]
        p_transpose = [p]
        alpha = sum(a * b for a, b in zip(r_transpose[0], z)) / sum(a * b for a, b in zip(p_transpose[0], v))

        x_new = [a + alpha * b for a, b in zip(x, p)]
        r_new = [a - alpha * b for a, b in zip(r, v)]
        z_new = solve_system(P, r_new)
        r_new_transpose = [r_new]
        
        beta = sum(a * b for a, b in zip(r_new_transpose[0], z_new)) / sum(a * b for a, b in zip(r_transpose[0], z))
        p = [a + beta * b for a, b in zip(z_new, p)]
        r = r_new
        z = z_new
        x = x_new
        k += 1
        Ax = matrix_vector_multiply(A, x_new)
        r_norm = [bi - ai for ai, bi in zip(Ax, b)]

    return x, k

# Define the size of the matrix
n = 4
A = [[5, 7, 6, 5], [7, 10, 8, 7], [6, 8, 10, 9], [5, 7, 9, 10]]
x = [0 for _ in range(n)]
b = [57, 79, 88, 86]

precond_types = ["identity", "jacobi", "sgs"]

for precond_type in precond_types:
    if precond_type == "identity":
        P = identity_matrix(n)
    elif precond_type == "jacobi":
        P = extract_diagonal(A)
    elif precond_type == "sgs":
        C = lower_triangular_matrix(A)
        P = multiply_lower_triangular_transpose(C)

    x_sd, k_sd = steepest_descent(P, A, b, x.copy())
    x_cg, k_cg = cg(P, A, b, x.copy())

    print(f"{precond_type.capitalize()} - Steepest Descent: x =", x_sd, "Iterations =", k_sd)
    print(f"{precond_type.capitalize()} - Conjugate Gradient: x =", x_cg, "Iterations =", k_cg)
    print("\n")
