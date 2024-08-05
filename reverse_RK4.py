# This file is testing ground for running RK4 backwards
import math
import pandas as pd
import matplotlib.pyplot as plt

### Functions

# Total population of stray dogs as a differential equation
def dNdt(N_t, k, p_f, s_a, s_i, l, r_r, n_s, S0, t, m):
    dN_dt = (N_t * ((k - N_t) / k) * (
            (p_f * s_a * s_i * l * (r_r/12) * (1 - ((((n_s + ((m + s_a - 1) * S0)) * t) + S0) / N_t))) - (1 - s_a) + m))
    return dN_dt


# Runge Kutta method of solving the differential equation
def runge_kutta(**kwargs):
    # Initialize variables from dictionary
    t = kwargs.get("t0")  # initial time
    N_t = kwargs.get("N0")  # initial value of N
    initial_S_N = kwargs.get("initial_S_N")
    n_s = kwargs.get("n_s")
    desired_t = kwargs.get("desired_t")  # target
    h = kwargs.get("h")  # chosen step size
    k = kwargs.get("k")
    p_f = kwargs.get("p_f")
    s_a = kwargs.get("s_a")
    s_i = kwargs.get("s_i")
    l = kwargs.get("l")
    r_r = kwargs.get("r_r")
    m = kwargs.get("m")

    S0 = N_t * initial_S_N # N_t here will be N0

    t_values = [t]
    N_t_values = [N_t]
    S_t_values = [S0]
    S_N_values = [(S0 / N_t)]

    # Iterate using fourth order Runge-Kutta
    while t < desired_t+1:
        # Calculate k1, k2, k3, k4
        k1 = h * dNdt(t=t, N_t=N_t, k=k, p_f=p_f, s_a=s_a, s_i=s_i, l=l, S0=S0, r_r=r_r, m=m, n_s=n_s)
        k2 = h * dNdt(t=t + 0.5 * h, N_t=N_t + 0.5 * k1, k=k, p_f=p_f, s_a=s_a, s_i=s_i, l=l, S0=S0, r_r=r_r, m=m, n_s=n_s)
        k3 = h * dNdt(t=t + 0.5 * h, N_t=N_t + 0.5 * k2, k=k, p_f=p_f, s_a=s_a, s_i=s_i, l=l, S0=S0, r_r=r_r, m=m, n_s=n_s)
        k4 = h * dNdt(t=t + h, N_t=N_t + k3, k=k, p_f=p_f, s_a=s_a, s_i=s_i, l=l, S0=S0, r_r=r_r, m=m, n_s=n_s)

        # Update N using weighted average of k's
        N_t = N_t + (k1 + 2 * k2 + 2 * k3 + k4) / 6

        # Update N_t in dictionary

        if N_t < 0:
            N_t = 0.1

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
        S_N_values.append((S_t / N_t))
        sterilisation_prop_data = pd.DataFrame()
        sterilisation_prop_data["Month"] = t_values
        sterilisation_prop_data["Sterilisation Proportion"] = S_N_values
        population_data = pd.DataFrame()
        population_data["Month"] = t_values
        population_data["Total Population"] = N_t_values
        population_data["Sterilised Population"] = S_t_values

        # Check if t is close enough to desired_t
        if abs(t - desired_t) < 0.01:
            break

        # Optionally, you can print or store the values of t and N_t
        # print(f"t = {t}, N_t = {N_t}")

    return math.ceil(N_t), t_values, N_t_values, S_t_values, S_N_values


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
                      "s_i": 0.4, "l": 6, "r_r": 2, "m": 0}

model = runge_kutta(**initial_parameters)

print(model[1])
print(model[2])
print(model[3])

graph_runge_kutta(t=model[1], N=model[2])
