import matplotlib.pyplot as plt
import math
import pandas as pd

# Total population of stray dogs as a differential equation
def dNdt(N, r, t):
    dN_dt = -N * r
    return dN_dt

# Runge Kutta method of solving the differential equation
def runge_kutta(t, desired_t, h, N, r):

    t_values = [t]
    N_values = [N]

    # Iterate using fourth order Runge-Kutta
    while t < desired_t+1:
        # Calculate k1, k2, k3, k4
        k1 = h * dNdt(t=t, N=N, r=r)
        k2 = h * dNdt(t=t + 0.5 * h, N=N + 0.5 * k1, r=r)
        k3 = h * dNdt(t=t + 0.5 * h, N=N + 0.5 * k2, r=r)
        k4 = h * dNdt(t=t + h, N=N + k3, r=r)

        # Update N using weighted average of k's
        N = N + (k1 + 2 * k2 + 2 * k3 + k4) / 6

        # Update t
        t = t + h


        # Store values in list
        t_values.append(t)
        N_values.append(N)

        # Check if t is close enough to desired_t
        if abs(t - desired_t) < 0.01:
            break

        # Optionally, you can print or store the values of t and N_t
        # print(f"t = {t}, N_t = {N_t}")

    return t_values, N_values


# This function uses the Runge-Kutta method to plot the population and sterilisation growth on a graph
def graph_runge_kutta(t, N):
    # Show dog population with sterilised dogs
    plt.figure(figsize=(10, 6))
    plt.plot(t, N, label='Total Population', marker='x')
    # plt.plot(t, sterilised_population, label='Sterilised Population', marker='x')
    plt.xlabel('Month')
    plt.ylabel('Number of Dogs')
    plt.title('Approximation of Total Population and Sterilised Population using Fourth Order Runge-Kutta Method')
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(True)
    plt.legend()
    plt.show()


### Code

model = runge_kutta(t=0, desired_t=6, N=1000, r=1, h=1)

print(model[0])
print(model[1])

graph_runge_kutta(t=model[0], N=model[1])