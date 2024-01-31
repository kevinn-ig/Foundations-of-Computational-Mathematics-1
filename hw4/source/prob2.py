import os
import numpy as np
import random as rd
import matplotlib.pyplot as plt
from sympy import symbols, lambdify
import sympy as sp

x = symbols('x')

# Bisection
def bisection(f, x):
    fx_list = []
    x_list = []
    dist_list = []
    a = 0
    b = 0
    iteration = 0  # Track the iteration number

    while (f.subs(x, a) * f.subs(x, b)) > 0:
        a = rd.uniform(-100, 100)
        b = rd.uniform(-100, 100)

    # Create a figure outside the loop
    plt.figure(figsize=(12, 8))
    plt.subplot(2, 1, 1)  # Create subplot for bisection iterations
    plt.axhline(0, color='black', linestyle='--', linewidth=0.8)
    plt.axvline(0, color='black', linestyle='--', linewidth=0.8)

    # Define a list of colors for iterations
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

    while abs(a - b) > 1e-6:
        iteration += 1  # Increment the iteration number

        # Select color based on iteration (cycling if needed)
        color = colors[(iteration - 1) % len(colors)]

        if iteration > 10:
            break

        # Different markers for initial values and iteration points
        if len(x_list) == 0:
            plt.scatter((a + b) / 2, f.subs(x, (a + b) / 2), color=color, marker='o', label=f'Initial Value - Iteration {iteration}')
        else:
            plt.plot([(a + b) / 2], [f.subs(x, (a + b) / 2)], color=color, marker='x', linestyle='-', label=f'Iteration {iteration}')

        if (f.subs(x, a) * f.subs(x, (a + b) / 2)) < 0:
            b = (a + b) / 2
        else:
            a = (a + b) / 2

        error = abs(a - b)
        dist_list.append(error)
        x_list.append((a + b) / 2)
        fx_list.append(f.subs(x, (a + b) / 2))

    # Add legend after the loop
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.title("Bisection Method with Iteration Points")
    plt.grid(True)

    # Plot errors over iteration number
    plt.subplot(2, 1, 2)  # Create subplot for errors
    plt.plot(range(1, iteration), dist_list, marker='o', linestyle='-', color='b')
    plt.xlabel('Iteration Number')
    plt.ylabel('Error')
    plt.title('Bisection Method - Error over Iteration Number')
    plt.grid(True)

    # Save the plot
    output_directory = '/Users/kevin_smith/Desktop/FSU_Relevant_Stuff/fall_2023/FCM1/homework/programming/hw4/figures/prob2'
    plt.savefig(os.path.join(output_directory, 'bisection_plot.png'))

    return x_list

# Newton's method
# Newton's method
# Newton's method
def newton(f, x, y, max_iter=10000):
    fx_list = []
    x_list = [y]
    error_list = []  # List to store errors at each iteration
    f_prime = f.diff(x)
    fx = f.subs(x, y)
    fx_list.append(fx)
    fpx = f_prime.subs(x, y)
    
    for _ in range(max_iter):
        if fpx == 0:
            print("derivative = 0, exiting")
            y = "error"
            break
        else:
            y = y - fx / fpx
            fx = f.subs(x, y)
            fpx = f_prime.subs(x, y)
            fx_list.append(fx)
            x_list.append(y)
            error_list.append(abs(x_list[-1] - x_list[-2]))

            if abs(fx) < 1e-6:
                print("yipee")  # Exit when f(x_k) is less than 1e-6
                break

    if y != "error":
        # Plot Newton's Method
        plt.figure(figsize=(8, 6))
        plt.plot(range(1, len(error_list)+1), error_list, marker='o', linestyle='-', label="Error")
        plt.xlabel('Iteration Number')
        plt.ylabel('Error')
        plt.title("Newton's Method - Convergence Rate for x_0 = " + str(x_list[0]))
        plt.legend()
        plt.grid(True)
        output_directory = '/Users/kevin_smith/Desktop/FSU_Relevant_Stuff/fall_2023/FCM1/homework/programming/hw4/figures/prob2'
        plt.savefig(os.path.join(output_directory, 'quad_conv_' + str(x_list[0]) + ".png"))  # Save the plot after showing

    return y, x_list, fx_list

# Fixed-Point Iteration
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

# Create directory if it doesn't exist
output_directory = '/Users/kevin_smith/Desktop/FSU_Relevant_Stuff/fall_2023/FCM1/homework/programming/hw4/figures/prob2'
os.makedirs(output_directory, exist_ok=True)
f = x**3 - x - 6

# Bisection method
x_list = bisection(f, x)
n_list = []
abs_e = []
for i in range(len(x_list)):
    abs_e.append(abs(2 - x_list[i]))
    n_list.append(i)

# Create a new figure for the bisection error
plt.figure(figsize=(8, 6))
plt.plot(n_list, abs_e, marker='o', linestyle='-', color='b')
plt.xlabel('Iteration Number')
plt.ylabel('|x - xk|')
plt.title('Bisection Method - Error over Iteration Number')
plt.grid(True)
plt.savefig(os.path.join(output_directory, 'bisection_error_plot.png'))

# Newton's method
y_list = [1/(3**(1/2)), 0.57735]
for i, y in enumerate(y_list):
    nx, nx_list, nfx_list = newton(f, x, y)
    if nx != "error":
        plt.figure(figsize=(8, 6))
        plt.plot(nx_list, nfx_list, marker='o', linestyle='-', label=f"Newton's Method - Initial Value: {y} - Stopped")
        plt.axhline(0, color='black', linestyle='--', linewidth=0.8)
        plt.axvline(0, color='black', linestyle='--', linewidth=0.8)
        plt.xlabel('x')
        plt.ylabel('f(x)')
        plt.title(f"Newton's Method with Iteration Points - Initial Value: {y}")
        plt.legend()
        plt.grid(True)
        plt.savefig(os.path.join(output_directory, f'newton_plot_{i}.png'))

# Fixed-Point Iteration
x_init = 0
f = lambdify(x, f, 'numpy')
g = lambda x: (x + 6)**(1/3)  # Define the iteration function g
fpx, fpx_list, fpfx_list = fixed_point_iteration(f, g, x_init)
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