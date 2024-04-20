# This file contains functions for equations describing the population growth of a stray dog population
# and the growth in the number of sterilised dogs in a population of stray dogs

import math
import numpy as np
import matplotlib.pyplot as plt


# Number of sterilised dogs
def num_sterilised_dogs(S0, s_a, m, n_s, t):
    S_t = ((n_s + ((m + s_a - 1) * S0)) * t) + S0
    return S_t


# Total population of stray dogs as a differential equation
def dNdt(N_t, k, p_f, s_a, s_i, l, r_r, S0, n_s, t, m):
    dN_dt = (N_t * ((k - N_t) / k) * (
            (p_f * s_a * s_i * l * r_r * (1 - (((n_s + ((m + s_a - 1) * S0)) * t) / N_t))) - (1 - s_a) + m))
    return dN_dt


# Runge Kutta method of solving the differential equation
def runge_kutta(t0, N0, desired_t, h, k, p_f, s_a, s_i, l, r_r, S0, m, n_s):
    # Initialize variables
    t = t0  # initial time
    N_t = N0  # initial value of N
    h = h  # chosen step size
    desired_t = desired_t  # target time

    t_values = [t0]
    N_t_values = [N0]
    S_t_values = [S0]
    S_N_values = [(S0 / N0) * 100]

    # Iterate using fourth order Runge-Kutta
    while t < desired_t+1:
        # Calculate k1, k2, k3, k4
        k1 = h * dNdt(t=t, N_t=N_t, k=k, p_f=p_f, s_a=s_a, s_i=s_i, l=l, r_r=r_r, S0=S0, m=m, n_s=n_s)
        k2 = h * dNdt(t=t + 0.5 * h, N_t=N_t + 0.5 * k1, k=k, p_f=p_f, s_a=s_a, s_i=s_i, l=l, r_r=r_r, S0=S0, m=m,
                      n_s=n_s)
        k3 = h * dNdt(t=t + 0.5 * h, N_t=N_t + 0.5 * k2, k=k, p_f=p_f, s_a=s_a, s_i=s_i, l=l, r_r=r_r, S0=S0, m=m,
                      n_s=n_s)
        k4 = h * dNdt(t=t + h, N_t=N_t + k3, k=k, p_f=p_f, s_a=s_a, s_i=s_i, l=l, r_r=r_r, S0=S0, m=m, n_s=n_s)

        # Update N using weighted average of k's
        N_t = N_t + (k1 + 2 * k2 + 2 * k3 + k4) / 6

        #if N_t < 0:
        #    N_t = 0

        # Update t
        t = t + h

        # Update S(t)
        S_t = ((n_s + ((m + s_a - 1) * S0)) * t) + S0

        if S_t > N_t:
            S_t = N_t

        #if S_t < 0:
        #    S_t = 0

        # Update S(t)/N(t)
        #if N_t == 0:
        #    S_N_t = 0
        #else:
        #    S_N_t = (S_t / N_t) * 100

        # Store values in list
        t_values.append(t)
        N_t_values.append(N_t)
        S_t_values.append(S_t)
        S_N_values.append((S_t / N_t) * 100)

        # Check if t is close enough to desired_t
        if abs(t - desired_t) < 0.01:
            break

        # Optionally, you can print or store the values of t and N_t
        # print(f"t = {t}, N_t = {N_t}")

    return math.ceil(N_t), t_values, N_t_values, S_t_values, S_N_values


# This function uses the Runge-Kutta method to plot the population and sterilisation growth on a graph
def graph_runge_kutta(t0, N0, desired_t, h, k, p_f, s_a, s_i, l, r_r, S0, m, n_s):

    population_model = runge_kutta(t0=t0, N0=N0, desired_t=desired_t, h=h, k=k, p_f=p_f, s_a=s_a, s_i=s_i, l=l, r_r=r_r, S0=S0, m=m, n_s=n_s)
    t_values = population_model[1]
    N_t_values = population_model[2]
    S_t_values = population_model[3]
    S_N_values = population_model[4]

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(t_values, N_t_values, label='Total Population', marker='o')
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
    plt.plot(t_values, S_N_values, label='Sterilised Dogs', marker='x')
    plt.xlabel('t')
    plt.ylabel('Percentage of population sterilised')
    plt.title('Approximation of percentage of population sterilised using Fourth Order Runge-Kutta Method')
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(True)
    plt.legend()
    plt.show()


# This equation calculates the number of sterilisations required per month to reach a target sterilisation proportion
# within a specified timeframe using binary search
def growth_binary_search_n_s(t0, N0, desired_t, h, k, p_f, s_a, s_i, l, r_r, S0, m, desired_S_N):

    lower_limit = 0
    upper_limit = 2000

    while upper_limit - lower_limit > 1:
        # Try n_s at the midrange point
        n_s = math.ceil((upper_limit + lower_limit) / 2)

        # Calculate S(t) and N(t) at desired time using chosen n_s
        S_t = ((n_s + ((m + s_a - 1) * S0)) * desired_t) + S0

        N_t = runge_kutta(t0=t0, N0=N0, desired_t=desired_t, h=h, k=k, p_f=p_f, s_a=s_a, s_i=s_i, l=l, r_r=r_r, S0=S0, m=m, n_s=n_s)
        if N_t < 0:
            N_t = 1

       # Calculate sterilisation proportion
        S_N = S_t / N_t

        if S_N < desired_S_N:
            lower_limit = n_s
        else:
            upper_limit = n_s

    return n_s + 1


# This function calculates the derivative of S(t)/N(t)
def dSNdt(t0, N0, desired_t, h, k, p_f, s_a, s_i, l, r_r, S0, m, n_s):
    S_t = ((n_s + ((m + s_a - 1) * S0)) * desired_t) + S0

    dS_dt = n_s + ((m + s_a - 1) * S0)

    N_t = runge_kutta(t0=t0, N0=N0, desired_t=desired_t, h=h, k=k, p_f=p_f, s_a=s_a, s_i=s_i, l=l, r_r=r_r, S0=S0, m=m,
                      n_s=n_s)

    dN_dt = (N_t * ((k - N_t) / k) * (
                (p_f * s_a * s_i * l * r_r * (1 - (((n_s + ((m + s_a - 1) * S0)) * desired_t) / N_t))) - (1 - s_a) + m))

    dSN_dt = ((N_t * dS_dt) - (S_t * dN_dt)) / N_t ** 2

    return dSN_dt


# This function finds a value for n_s which keeps the sterilisation proportion at equilibrium
def maintenance_binary_search_n_s(t0, N0, desired_t, h, k, p_f, s_a, s_i, l, r_r, S0, m):
    lower_limit = 0
    upper_limit = 2000

    while upper_limit - lower_limit > 1:
        # Try n_s at the midrange point
        n_s = math.ceil((upper_limit + lower_limit) / 1.5)

        # Calculate the derivative of S(t)/N(t) across the desired timeframe
        dSN_dt = dSNdt(t0=t0, N0=N0, desired_t=desired_t, h=h, k=k, p_f=p_f, s_a=s_a, s_i=s_i, l=l, r_r=r_r, S0=S0, m=m,
                       n_s=n_s)
        print(dSN_dt)

        if dSN_dt < 0:
            lower_limit = n_s
        else:
            upper_limit = n_s

    return n_s + 1
