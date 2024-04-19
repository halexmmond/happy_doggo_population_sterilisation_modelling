# This file contains functions for equations describing the population growth of a stray dog population
# and the growth in the number of sterilised dogs in a population of stray dogs

import math
import numpy as np
import matplotlib.pyplot as plt

# Number of sterilised dogs
def num_sterilised_dogs(S0, s_a, m, n_s, t):
    S_t = (n_s + ((m + s_a - 1) * S0)) * t
    return S_t


# Total population of stray dogs as a differential equation
def dNdt(N, k, p_f, s_a, s_i, l, r_r, S0, n_s, t, m):
    dN_dt = (N * ((k - N) / k) * ((p_f * s_a * s_i * l * r_r * (1 - (((n_s + ((m + s_a - 1) * S0)) * t) / N))) - (1 - s_a) + m))
    return dN_dt

# Runge Kutta method of solving the differential equation
def runge_kutta(t0, N0, desired_t, h, k, p_f, s_a, s_i, l, r_r, S0, m, n_s):
    # Initialize variables
    t = t0  # initial time
    N = N0  # initial value of N
    h = h  # chosen step size
    desired_t = desired_t  # target time

    t_values = [t0]
    N_values = [N0]
    S_t_values = [S0]
    S_t_N_values = [S0/N0]

    # Iterate using fourth order Runge-Kutta
    while t < desired_t:
        # Calculate k1, k2, k3, k4
        k1 = h * dNdt(t=t, N=N, k=k, p_f=p_f, s_a=s_a, s_i=s_i, l=l, r_r=r_r, S0=S0, m=m, n_s=n_s)
        k2 = h * dNdt(t = t + 0.5 * h, N = N + 0.5 * k1, k=k, p_f=p_f, s_a=s_a, s_i=s_i, l=l, r_r=r_r, S0=S0, m=m, n_s=n_s)
        k3 = h * dNdt(t = t + 0.5 * h, N = N + 0.5 * k2, k=k, p_f=p_f, s_a=s_a, s_i=s_i, l=l, r_r=r_r, S0=S0, m=m, n_s=n_s)
        k4 = h * dNdt(t = t + h, N = N + k3, k=k, p_f=p_f, s_a=s_a, s_i=s_i, l=l, r_r=r_r, S0=S0, m=m, n_s=n_s)

        # Update N using weighted average of k's
        N = N + (k1 + 2 * k2 + 2 * k3 + k4) / 6

        # Update t
        t = t + h

        # Update S(t)
        S_t = (n_s + ((m + s_a - 1) * S0)) * t

        # Update S(t)/N(t)

        # Store values in list
        t_values.append(t)
        N_values.append(N)
        S_t_values.append(S_t)
        S_t_N_values.append(S_t/N*100)

        # Check if t is close enough to desired_t
        if abs(t - desired_t) < 0.01:
            break

        # Optionally, you can print or store the values of t and N
        #print(f"t = {t}, N = {N}")

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(t_values, N_values, label='Total Population', marker='o')
    plt.plot(t_values, S_t_values, label='Sterilised Dogs', marker='o')
    plt.xlabel('t')
    plt.ylabel('N(t)')
    plt.title('Approximation of N(t) using Fourth Order Runge-Kutta Method')
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(True)
    plt.legend()
    plt.show()

    # Show sterilisation proportion
    plt.figure(figsize=(10, 6))
    plt.plot(t_values, S_t_N_values, label='Sterilised Dogs', marker='x')
    plt.xlabel('t')
    plt.ylabel('S(t)/N(t)')
    plt.title('Approximation of S(t)/N(t) using Fourth Order Runge-Kutta Method')
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(True)
    plt.legend()
    plt.show()

    return N


