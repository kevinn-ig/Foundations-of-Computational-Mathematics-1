import random
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def scale_matrix(matrix, factor):
    return [[element * factor for element in row] for row in matrix]

def scale_vector(vector, factor):
    return [element * factor for element in vector]

def generate_spd_matrix(n):
    A = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            A[i][j] = random.uniform(-10, 10)
            A[j][i] = A[i][j]  # Make the matrix symmetric
    for i in range(n):
        A[i][i] += 100  # Ensure diagonal dominance for positive definiteness
    return A

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

def relative_error(x, x_true):
    # Calculate the relative error
    x_norm = calculate_l2_norm(x)
    x_true_norm = calculate_l2_norm(x_true)
    return x_norm / x_true_norm

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
    r_norm_list = []

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
        r_norm_list.append(r_norm)
        r_norm = [bi - ai for ai, bi in zip(Ax, b)]
        


    return x, k, r_norm_list


def cg(P, A, b, x):
    Ax = matrix_vector_multiply(A, x)
    r = [bi - ai for ai, bi in zip(Ax, b)]
    z = solve_system(P, r)
    p = z
    v = []
    k = -1
    Ax = matrix_vector_multiply(A, x)
    r_norm = [bi - ai for ai, bi in zip(Ax, b)]
    r_norm_list = []

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
        r_norm_list.append(r_norm)
        r_norm = [bi - ai for ai, bi in zip(Ax, b)]
        

    return x, k, r_norm_list
    


# Define the size of the matrix
n = 20
num_trials = 100
precond_types = ["identity", "jacobi", "sgs"]

total_iterations_steepest = {precond_type: 0 for precond_type in precond_types}
total_iterations_cg = {precond_type: 0 for precond_type in precond_types}

fig, axes = plt.subplots(2, 1, figsize=(10, 8))

errors_steepest = []
errors_cg = []
condition_numbers = []
iterations_steepest = []
iterations_cg = []

for _ in range(num_trials):
    A = generate_spd_matrix(n)
    b = scale_vector([random.uniform(-100, 100) for _ in range(n)], 1e-6)
    A = scale_matrix(A, 1e-2)

    condition_number = np.linalg.cond(A)
    print("Condition Number:", condition_number)

    for precond_type in precond_types:
        x = [0 for _ in range(n)]  # Reset x to zero vector for each trial
        if precond_type == "identity":
            P = identity_matrix(n)
        elif precond_type == "jacobi":
            P = extract_diagonal(A)
        elif precond_type == "sgs":
            C = lower_triangular_matrix(A)
            P = multiply_lower_triangular_transpose(C)

        x1, k_steepest, sd_norms = steepest_descent(P, A, b, x)
        x2, k_cg, cg_norms = cg(P, A, b, x)

        total_iterations_steepest[precond_type] += k_steepest
        total_iterations_cg[precond_type] += k_cg

        # Calculate errors
        error_steepest = norm_ratio(sd_norms[-1], b)
        error_cg = norm_ratio(cg_norms[-1], b)

        # Calculate relative error of x
        relative_error_steepest = relative_error(x1, np.linalg.solve(A, b))
        relative_error_cg = relative_error(x2, np.linalg.solve(A, b))

        errors_steepest.append(error_steepest)
        errors_cg.append(error_cg)

        # Store condition numbers, iteration counts, and relative errors
        condition_numbers.append(condition_number)
        iterations_steepest.append(k_steepest)
        iterations_cg.append(k_cg)

        # Plot residual norms
        plt.figure()
        plt.plot(sd_norms, label="Steepest Descent")
        plt.title(f"Residual Norms SD - {precond_type}")
        plt.xlabel("Iteration")
        plt.ylabel("Residual Norm")
        plt.savefig(f"/Users/kevin_smith/Desktop/FSU_Relevant_Stuff/fall_2023/FCM1/homework/programming/hw3/figures/n20/residual_norms_sd_{precond_type}.png")
        plt.close()
        plt.figure()
        plt.plot(cg_norms, label="Conjugate Gradient")
        plt.title(f"Residual Norms CG - {precond_type}")
        plt.xlabel("Iteration")
        plt.ylabel("Residual Norm")
        plt.savefig(f"/Users/kevin_smith/Desktop/FSU_Relevant_Stuff/fall_2023/FCM1/homework/programming/hw3/figures/n20/residual_norms_cg_{precond_type}.png")
        plt.close()
    print(f"iteration {_ +1 } done")

# Calculate average iterations
average_iterations_steepest = {precond_type: total_iterations_steepest[precond_type] / num_trials for precond_type in precond_types}
average_iterations_cg = {precond_type: total_iterations_cg[precond_type] / num_trials for precond_type in precond_types}

# Print average iterations
print("Average Iterations - Steepest Descent:", average_iterations_steepest)
print("Average Iterations - Conjugate Gradient:", average_iterations_cg)

# Create a DataFrame
data = {
    "Error Steepest": errors_steepest,
    "Error CG": errors_cg,
    "Relative Error Steepest": relative_error_steepest,
    "Relative Error CG": relative_error_cg,
    "Condition Number": condition_numbers,
    "Iterations Steepest": iterations_steepest,
    "Iterations CG": iterations_cg,
    "Preconditioner Type": [precond_type for _ in range(len(errors_steepest))],
}

df = pd.DataFrame(data)

# Display the table
print(df)

# Save the DataFrame to a CSV file
df.to_csv("/Users/kevin_smith/Desktop/FSU_Relevant_Stuff/fall_2023/FCM1/homework/programming/hw3/data/n20.csv", index=False)
