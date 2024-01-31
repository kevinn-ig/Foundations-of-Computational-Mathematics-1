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
    output_directory = '/Users/kevin_smith/Desktop/FSU_Relevant_Stuff/fall_2023/FCM1/homework/programming/hw4/figures/prob1'
    plt.savefig(os.path.join(output_directory, 'bisection_plot.png'))
    print(f.subs(x, (a + b) / 2))
    print((a + b) / 2)

# Newton's method
def newton(f, x, y):
    fx_list = []
    x_list = [y]
    f_prime = f.diff(x)
    fx = f.subs(x, y)
    fx_list.append(fx)
    fpx = f_prime.subs(x, y)
    while abs(fx) > 1e-4:
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


    return y, x_list, fx_list

# Fixed-Point Iteration
def fixed_point_iteration(f, g, y, x, max_iter=100):
    fx_list = [f.subs(x, y)]
    x_list = [y]

    for _ in range(max_iter):
        x_next = g.subs(x, y)
        fx_next = f.subs(x, x_next)
        fx_list.append(fx_next)
        x_list.append(x_next)

        if abs(fx_next) < 1e-6:
            break

        y = x_next

    return y, x_list, fx_list

# Create directory if it doesn't exist
output_directory = '/Users/kevin_smith/Desktop/FSU_Relevant_Stuff/fall_2023/FCM1/homework/programming/hw4/figures/prob1'
os.makedirs(output_directory, exist_ok=True)
f = x * sp.exp(-x) - 0.06064

# Bisection method
bisection(f, x)

# Newton's method with multiple initial guesses
initial_guesses = [1, 0.99, 0.5, 1.5]  # Add more initial guesses as needed

for i, initial_guess in enumerate(initial_guesses):
    nx, nx_list, nfx_list = newton(f, x, initial_guess)
    
    # Check for the special case when derivative is zero
    if nx == "error":
        print(f"Newton's method with initial guess {initial_guess}: Derivative is zero. Exiting.")
    else:
        plt.figure(figsize=(8, 6))
        plt.plot(nx_list, nfx_list, marker='o', linestyle='-', label=f"Newton's Method - Initial Guess: {initial_guess}")
        plt.axhline(0, color='black', linestyle='--', linewidth=0.8)
        plt.axvline(0, color='black', linestyle='--', linewidth=0.8)
        plt.xlabel('x')
        plt.ylabel('f(x)')
        plt.title(f"Newton's Method with Iteration Points - Initial Guess: {initial_guess}")
        plt.legend()
        plt.grid(True)
        plt.savefig(os.path.join(output_directory, f'newton_plot_{i}.png'))

# Fixed-Point Iteration
x_init = -0.5
fpx, fpx_list, fpfx_list = fixed_point_iteration(f, x - x * sp.exp(-x) - 0.06064, x_init, x)
plt.figure(figsize=(8, 6))
plt.plot(fpx_list, fpfx_list, marker='o', linestyle='-', color='g', label='Fixed-Point Iteration ')
plt.axhline(0, color='black', linestyle='--', linewidth=0.8)
plt.axvline(0, color='black', linestyle='--', linewidth=0.8)
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title("Fixed-Point Iteration with Iteration Points")
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(output_directory, 'fixed_point_iteration_plot.png'))

# Plot f(x)
f_numeric = lambdify(x, f, 'numpy')
x_values = np.linspace(-1, 3, 1000)
y_values = f_numeric(x_values)
plt.figure(figsize=(8, 6))
plt.plot(x_values, y_values, label='$f(x) = x \cdot e^{-x} - 0.06064$')
plt.axhline(0, color='black', linestyle='--', linewidth=0.8)
plt.axvline(0, color='black', linestyle='--', linewidth=0.8)
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('Plot of $f(x) = x \cdot e^{-x} - 0.06064$')
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(output_directory, 'fx_plot.png'))

plt.show()
