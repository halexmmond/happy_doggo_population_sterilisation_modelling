# This file is testing ground for running RK4 backwards
import math
import pandas as pd
import matplotlib.pyplot as plt

### Functions
def dN_dt_reverse_RK4(N_t, k, p_a, p_f, s_a, s_i, l, r_r, S0, m, t):
    dNdt = -(N_t * ((k - N_t) / k) * ((p_a * p_f * s_a * s_i * l * (r_r/12) * (1 - (S0/N_t))) - ((1 - s_a) * p_a) + m))
    return dNdt

def reverse_births_RK4(N_t, k, p_a, s_a, m, dNdt):
    births = (dNdt / ((k - N_t) / k)) - m + ((1 - s_a) * p_a)
    return births

# Function used to find 12 months worth of previous data (N and births) to initialise model
def runge_kutta(**kwargs):
    # Initialize variables from dictionary
    t = kwargs.get("t0")  # initial time
    N_t = kwargs.get("N0")  # initial value of N
    desired_t = kwargs.get("desired_t")  # target
    h = kwargs.get("h")  # chosen step size
    k = kwargs.get("k")
    p_a = kwargs.get("p_a")
    p_f = kwargs.get("p_f")
    s_a = kwargs.get("s_a")
    s_i = kwargs.get("s_i")
    l = kwargs.get("l")
    r_r = kwargs.get("r_r")
    m = kwargs.get("m")
    S_t = 0

    # Not adding month 0 values here because we already have them and are only interested in -1 onwards
    t_values = []
    N_t_values = []
    births_values = []

    # Iterate using fourth order Runge-Kutta
    while t < desired_t + 1:
        # Calculate k1, k2, k3, k4
        k1 = h * dN_dt_reverse_RK4(t=t, N_t=N_t, k=k, p_a=p_a, p_f=p_f, s_a=s_a, s_i=s_i, l=l, S0=S_t, r_r=r_r, m=m)
        k2 = h * dN_dt_reverse_RK4(t=t + 0.5 * h, N_t=N_t + 0.5 * k1, k=k, p_a=p_a, p_f=p_f, s_a=s_a, s_i=s_i, l=l,
                                   S0=S_t, r_r=r_r, m=m)
        k3 = h * dN_dt_reverse_RK4(t=t + 0.5 * h, N_t=N_t + 0.5 * k2, k=k, p_a=p_a, p_f=p_f, s_a=s_a, s_i=s_i, l=l,
                                   S0=S_t, r_r=r_r, m=m)
        k4 = h * dN_dt_reverse_RK4(t=t + h, N_t=N_t + k3, k=k, p_a=p_a, p_f=p_f, s_a=s_a, s_i=s_i, l=l, S0=S_t, r_r=r_r, m=m)

        # Update N using weighted average of k's
        N_t = round(N_t + (k1 + 2 * k2 + 2 * k3 + k4) / 6)

        # Update N_t in dictionary

        if N_t < 0:
            N_t = 0.1

        # Update t
        t = t + h

        # Store values in list
        t_values.append(t)
        N_t_values.append(N_t)

        # Check if t is close enough to desired_t
        if abs(t - desired_t) < 0.01:
            break

    for t in t_values[:-1]:
        pop_change = N_t_values[t-1] - N_t_values[t]

        births = round(reverse_births_RK4(N_t=N_t_values[t], k=k, p_a=p_a, s_a=s_a, m=m, dNdt=pop_change))

        births_values.append(births)

    return t_values, N_t_values, births_values

# This function uses the Runge-Kutta method to plot the population and sterilisation growth on a graph
def graph_runge_kutta(t, N):
    #months = df_sterilisation_prop["Month"]
    #sterilisation_props = df_sterilisation_prop["Sterilisation Proportion"]
    #total_population = df_population["Total Population"]
    #sterilised_population = df_population["Sterilised Population"]

    # Show sterilisation proportion
    #plt.figure(figsize=(10, 6))
    #plt.plot(months, sterilisation_props, label='', marker='x')
    #plt.xlabel('Month')
    #plt.ylabel('Sterilisation Proportion')
    #plt.title('Approximation sterilisation proportion using Fourth Order Runge-Kutta Method')
    #plt.axhline(0, color='black', linewidth=0.5)
    #plt.axvline(0, color='black', linewidth=0.5)
    #plt.axhline(lower_S_N, color='red', linestyle='dotted')
    #plt.axhline(upper_S_N, color='red', linestyle='dotted')
    #plt.grid(True)
    #plt.legend()
    #plt.show()

    # Show dog population with sterilised dogs
    plt.figure(figsize=(10, 6))
    plt.plot(t, N, label='Total Population', marker='x')
    #plt.plot(t, sterilised_population, label='Sterilised Population', marker='x')
    plt.xlabel('Month')
    plt.ylabel('Number of Dogs')
    plt.title('Approximation of Total Population and Sterilised Population using Fourth Order Runge-Kutta Method')
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(True)
    plt.legend()
    plt.show()


### Code

# Create a dictionary to keep arguments/parameters of model in
initial_parameters = {"t0": 0, "N0": 1000, "initial_S_N": 0, "desired_S_N": None, "lower_S_N": None,
                      "upper_S_N": None, "n_s": 0, "desired_t": 13, "h": 1, "k": 100000, "p_f": 0.5, "s_a": 1,
                      "s_i": 0.4, "l": 6, "r_r": 2, "m": 0, "p_a": 1}

model = runge_kutta(**initial_parameters)

print(model[0][:-1])
print(model[1][:-1])
print(model[2])

graph_runge_kutta(t=model[0], N=model[1])
