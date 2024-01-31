import numpy as np
import random as rd
import matplotlib.pyplot as plt
from sympy import symbols, lambdify, exp
import sympy as sp
import os

x = symbols('x')

def bisection(f, x, a_range=(-100, 100)):
    fx_list = []
    x_list = []
    dist_list = []
    a, b = a_range

    while abs(a - b) > 1e-6:
        x_value = (a + b) / 2
        x_list.append(x_value)
        fx_list.append(f.subs(x, x_value))

        if (f.subs(x, a) * f.subs(x, b)) < 0:
            b = x_value
        else:
            a = x_value

        error = abs(a - b)
        dist_list.append(error)

    # Plot Bisection
    plt.figure(figsize=(8, 6))
    plt.plot(x_list, fx_list, marker='o', linestyle='-', color='b', label='Bisection')
    plt.xlabel('x_k')
    plt.ylabel('f(x_k)')
    plt.title("Bisection Method - x_k vs f(x_k)")
    plt.legend()
    plt.grid(True)

    # Save the plot
    plt.savefig('/Users/kevin_smith/Desktop/FSU_Relevant_Stuff/fall_2023/FCM1/homework/programming/hw4/figures/prob3/bisection_plot.png')
    

def newton(f, x, y):
    # Convert the symbolic function and its derivative to numerical functions
    f_numeric = lambdify(x, f, 'numpy')
    f_prime_numeric = lambdify(x, f.diff(x), 'numpy')

    fx_list = []
    x_list = [y]
    fx = f_numeric(y)
    fx_list.append(fx)
    fpx = f_prime_numeric(y)

    while abs(fx) > 1e-4:
        if fpx == 0:
            print("Derivative is 0, exiting")
            y = "error"
            break
        else:
            y = y - fx / fpx
            fx = f_numeric(y)
            fpx = f_prime_numeric(y)
            fx_list.append(fx)
            x_list.append(y)

    # Plot Newton's Method
    plt.figure(figsize=(8, 6))
    plt.plot(x_list, fx_list, marker='o', linestyle='-', label="Newton's Method")
    plt.xlabel('x_k')
    plt.ylabel('f(x_k)')
    plt.title("Newton's Method - x_k vs f(x_k)")
    plt.legend()
    plt.grid(True)

    # Save the plot
    plt.savefig('/Users/kevin_smith/Desktop/FSU_Relevant_Stuff/fall_2023/FCM1/homework/programming/hw4/figures/prob3/newton_plot.png')
    

def fixed_point_iteration(f, g, y, max_iter=100):
    fx_list = [f(y)]
    x_list = [y]

    for _ in range(max_iter):
        x_next = g(y)
        fx_next = f(x_next)
        fx_list.append(fx_next)
        x_list.append(x_next)

        if abs(fx_next) < 1e-6:
            break

        y = x_next

    return y, x_list, fx_list


output_directory = "/Users/kevin_smith/Desktop/FSU_Relevant_Stuff/fall_2023/FCM1/homework/programming/hw4/figures/prob3"

# Function Plotting
f = x * exp(-x) - 0.06064

# Bisection method
bisection(f, x, a_range=(1, 13))

# Newton's method
y_list = [2]
for i, y in enumerate(y_list):
    newton(f, x, y)

# Fixed-Point Iteration
x_init = 0
f_numeric = lambdify(x, f, 'numpy')

# Modify the iteration function to prevent overflow
g_numeric = lambda x: np.clip(x + x * np.exp(-x) - 0.06064, -1e10, 1e10)

fpx, fpx_list, fpfx_list = fixed_point_iteration(f_numeric, g_numeric, x_init)
plt.figure(figsize=(8, 6))
plt.plot(fpx_list, fpfx_list, marker='o', linestyle='-', color='g', label='Fixed-Point Iteration - Stopped')
plt.axhline(0, color='black', linestyle='--', linewidth=0.8)
plt.axvline(0, color='black', linestyle='--', linewidth=0.8)
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title("Fixed-Point Iteration with Iteration Points")
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(output_directory, 'fixed_point_iteration_plot.png'))
