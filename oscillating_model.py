# This script will explore creating a model where you can input initial conditions to see sterilisation
# rate growth/decline and then add "interventions" where you're sterilising a certain number of dogs each
# month, etc.
# I'm thinking there will be graph of the sterilisation rate that you can add data to e.g. this number of
# sterilisations for this many months, or work that out, and also have lines to show your upper and lower
# sterilisation rates. Let's see!

import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# This file contains functions for equations describing the population growth of a stray dog population
# and the growth in the number of sterilised dogs in a population of stray dogs

import math
import numpy as np
import matplotlib.pyplot as plt

### Functions

# Calculate number of sterilised dogs for initial calculation from initial sterilisation rate
# Number of sterilised dogs
def num_sterilised_dogs(initial_S_N, N0, s_a, m, n_s, t):
    S0 = N0 * initial_S_N
    S_t = ((n_s + ((m + s_a - 1) * S0)) * t) + S0
    return S_t


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

    return math.ceil(N_t), t_values, N_t_values, S_t_values, S_N_values, sterilisation_prop_data, population_data


# This function uses the Runge-Kutta method to plot the population and sterilisation growth on a graph
def graph_runge_kutta(df_sterilisation_prop, df_population):

    months = df_sterilisation_prop["Month"]
    sterilisation_props = df_sterilisation_prop["Sterilisation Proportion"]
    total_population = df_population["Total Population"]
    sterilised_population = df_population["Sterilised Population"]
    

    # Show sterilisation proportion
    plt.figure(figsize=(10, 6))
    plt.plot(months, sterilisation_props, label='', marker='x')
    plt.xlabel('Month')
    plt.ylabel('Sterilisation Proportion')
    plt.title('Approximation sterilisation proportion using Fourth Order Runge-Kutta Method')
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.axhline(lower_S_N, color='red', linestyle='dotted')
    plt.axhline(upper_S_N, color='red', linestyle='dotted')
    plt.grid(True)
    plt.legend()
    plt.show()

    # Show dog population with sterilised dogs
    plt.figure(figsize=(10, 6))
    plt.plot(months, total_population, label='Total Population', marker='x')
    plt.plot(months, sterilised_population, label='Sterilised Population', marker='x')
    plt.xlabel('Month')
    plt.ylabel('Number of Dogs')
    plt.title('Approximation of Total Population and Sterilised Population using Fourth Order Runge-Kutta Method')
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(True)
    plt.legend()
    plt.show()


# This equation calculates the number of sterilisations required per month to reach a target sterilisation proportion
# within a specified timeframe using binary search
def growth_binary_search_n_s(**kwargs):
    # Initialize variables from dictionary
    t0 = kwargs.get("t0")  # initial time
    N0 = kwargs.get("N0")  # initial value of N
    initial_S_N = kwargs.get("initial_S_N")
    desired_S_N = kwargs.get("desired_S_N")
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

    lower_limit = 0
    upper_limit = 2000

    S0 = N0 * initial_S_N

    while upper_limit - lower_limit > 1:
        # Try n_s at the midrange point
        n_s = math.ceil((upper_limit + lower_limit) / 2)

        # Calculate S(t) and N(t) at desired time using chosen n_s
        S_t = ((n_s + ((m + s_a - 1) * S0)) * desired_t) + S0

        N_t = runge_kutta(t0=t0, N0=N0, initial_S_N=initial_S_N, desired_t=desired_t, h=h, k=k, p_f=p_f, s_a=s_a, s_i=s_i, l=l, r_r=r_r, m=m, n_s=n_s)[0]
        if N_t < 0:
            N_t = 1

       # Calculate sterilisation proportion
        S_N = S_t / N_t

        if S_N < desired_S_N:
            lower_limit = n_s
        else:
            upper_limit = n_s

    return n_s + 1


# This equation calculates the number of months required of a provided number of sterilisation
# per month to reach a target sterilisation proportion using binary search
def growth_binary_search_t(t0, N0, initial_S_N, desired_n_s, h, k, p_f, s_a, s_i, l, r_r, m, desired_S_N):
    # Note that I haven't changed anything in here to use the dictionaries and kwargs

    lower_limit = 0
    upper_limit = 2000

    S0 = N0 * initial_S_N

    while upper_limit - lower_limit > 1:
        # Try t at the midrange point
        t = math.ceil((upper_limit + lower_limit) / 2)

        # Calculate S(t) and N(t) at desired time using chosen n_s
        S_t = ((desired_n_s + ((m + s_a - 1) * S0)) * t) + S0

        N_t = runge_kutta(t0=t0, N0=N0, initial_S_N=initial_S_N, desired_t=t, h=h, k=k, p_f=p_f, s_a=s_a, s_i=s_i, l=l, r_r=r_r, m=m, n_s=desired_n_s)[0]
        if N_t < 0:
            N_t = 1

       # Calculate sterilisation proportion
        S_N = S_t / N_t

        if S_N < desired_S_N:
            lower_limit = t
        else:
            upper_limit = t

    return t + 1


# Add intervention
# This function will allow the user to add an intervention i.e. sterilising a certain number of dogs
# per month
def add_intervention(current_df_sp, current_df_pop, **kwargs):
    # All that happens in this function is that the dataframe gets changed so that it includes
    # added interventions
    while True:
        # First we need to ask what information the user has and build the graph based on that input
        intervention_type = input("Which of the following describes your situation?\n"
                                  "1. I know the timeframe and monthly sterilisation rate and want to add this to the model.\n"
                                  "2. I know the timeframe and target sterilisation proportion, but I want to calculate "
                                  "the required monthly sterilisation rate to reach the target proportion and add this to the model.\n").strip()

        if intervention_type == "1": # Have all required data and just want to plot
            # Get the inputs and run the function for 1

            # Intervention parameter assignment
            start_month = int(input("Which month do you want to start the intervention?: ").strip())
            intervention_parameters["N0"] = current_df_pop["Total Population"].iloc[start_month]
            intervention_parameters["initial_S_N"] = current_df_sp["Sterilisation Proportion"].iloc[start_month]
            #intervention_parameters["lower_S_N"] = float(input("Please enter a value for the lower sterilisation proportion "
            #                                                 "limit as a decimal, which will be shown on the graph as a dotted line: ").strip())
            #intervention_parameters["upper_S_N"] = float(input("Please enter a value for the upper sterilisation proportion "
            #                                                 "limit as a decimal, which will be shown on the graph as a dotted line: ").strip())
            intervention_parameters["n_s"] = int(input("How many sterilisations per month will be taking place?: ").strip())
            intervention_parameters["desired_t"] = int(input("How many months would you like this intervention to run?: ").strip())
            intervention_parameters["k"] = kwargs.get("k")
            intervention_parameters["p_f"] = kwargs.get("p_f")
            intervention_parameters["s_a"] = kwargs.get("s_a")
            intervention_parameters["s_i"] = kwargs.get("s_i")
            intervention_parameters["l"] = kwargs.get("l")
            intervention_parameters["r_r"] = kwargs.get("r_r")
            intervention_parameters["m"] = kwargs.get("m")


            # Create a copy of current dataframe to modify
            intervention_df_sp = current_df_sp.copy()
            intervention_df_pop = current_df_pop.copy()

            # Calculate population for each month of intervention
            intervention_data = runge_kutta(**intervention_parameters)
            intervention_data_sp = intervention_data[5]
            intervention_data_pop = intervention_data[6]

            # Isolate sterilisation props and population data
            intervention_sterilisation_props = intervention_data_sp.iloc[1:, 1]
            intervention_population_total = intervention_data_pop.iloc[1:, 1]
            intervention_population_sterilised = intervention_data_pop.iloc[1:, 2]

            print(intervention_data_pop)
            print(intervention_population_total)
            print(intervention_population_sterilised)

            # Replace intervention timeframe with intervention data
            intervention_df_sp["Sterilisation Proportion"][start_month+1:start_month+intervention_parameters.get("desired_t")
                                                                       +1] = intervention_sterilisation_props
            intervention_df_pop["Total Population"][start_month+1:start_month+intervention_parameters.get("desired_t")
                                                                +1] = intervention_population_total
            intervention_df_pop["Sterilised Population"][start_month+1:start_month+intervention_parameters.get("desired_t")
                                                                     +1] = intervention_population_sterilised

            # Set post intervention parameters using what we've calculated
            post_intervention_parameters["N0"] = intervention_df_pop["Total Population"].iloc[start_month+intervention_parameters.get("desired_t")]
            post_intervention_parameters["initial_S_N"] = intervention_df_sp["Sterilisation Proportion"].iloc[start_month+intervention_parameters.get("desired_t")]
            post_intervention_parameters["lower_S_N"] = intervention_parameters.get("lower_S_N")
            post_intervention_parameters["upper_S_N"] = intervention_parameters.get("upper_S_N")
            post_intervention_parameters["n_s"] = kwargs.get("n_s")
            post_intervention_parameters["desired_t"] = kwargs.get("desired_t")-start_month-intervention_parameters.get("desired_t")
            post_intervention_parameters["k"] = kwargs.get("k")
            post_intervention_parameters["p_f"] = kwargs.get("p_f")
            post_intervention_parameters["s_a"] = kwargs.get("s_a")
            post_intervention_parameters["s_i"] = kwargs.get("s_i")
            post_intervention_parameters["l"] = kwargs.get("l")
            post_intervention_parameters["r_r"] = kwargs.get("r_r")
            post_intervention_parameters["m"] = kwargs.get("m")

            # Calculate population for each month post intervention
            post_intervention_data = runge_kutta(**post_intervention_parameters)
            post_intervention_data_sp = post_intervention_data[5]
            post_intervention_data_pop = post_intervention_data[6]

            post_intervention_sp = post_intervention_data_sp.iloc[1:, 1]
            post_intervention_pop_total = post_intervention_data_pop.iloc[1:, 1]
            post_intervention_pop_sterilised = post_intervention_data_pop.iloc[1:, 2]

            # Replace post-intervention timeframe with post-intervention data
            intervention_df_sp.iloc[start_month+intervention_parameters.get("desired_t")+1:kwargs.get("desired_t")+1, 1] = post_intervention_sp
            intervention_df_pop.iloc[start_month+intervention_parameters.get("desired_t")+1:kwargs.get("desired_t")+1, 1] = post_intervention_pop_total
            intervention_df_pop.iloc[start_month+intervention_parameters.get("desired_t")+1:kwargs.get("desired_t")+1, 2] = post_intervention_pop_sterilised

            return intervention_df_sp, intervention_df_pop

        elif intervention_type == "2": # Have a target sterilisation and timeframe and need monthly sterilisation rate
            # Get the inputs and run the function for 2

            # Intervention parameter assignment
            start_month = int(input("Which month do you want to start the intervention?: ").strip())
            intervention_parameters["N0"] = current_df_pop["Total Population"].iloc[start_month]
            intervention_parameters["initial_S_N"] = current_df_sp["Sterilisation Proportion"].iloc[start_month]
            intervention_parameters["desired_S_N"] = float(input("Please enter the target sterilisation proportion as a decimal: ").strip())
            intervention_parameters["lower_S_N"] = float(
                input("Please enter a value for the lower sterilisation proportion "
                      "limit as a decimal, which will be shown on the graph as a dotted line: ").strip())
            intervention_parameters["upper_S_N"] = float(
                input("Please enter a value for the upper sterilisation proportion "
                      "limit as a decimal, which will be shown on the graph as a dotted line: ").strip())
            intervention_parameters["desired_t"] = int(
                input("How many months would you like this intervention to run?: ").strip())
            intervention_parameters["k"] = kwargs.get("k")
            intervention_parameters["p_f"] = kwargs.get("p_f")
            intervention_parameters["s_a"] = kwargs.get("s_a")
            intervention_parameters["s_i"] = kwargs.get("s_i")
            intervention_parameters["l"] = kwargs.get("l")
            intervention_parameters["r_r"] = kwargs.get("r_r")
            intervention_parameters["m"] = kwargs.get("m")

            monthly_sterilisations = growth_binary_search_n_s(**intervention_parameters)
            intervention_parameters["n_s"] = monthly_sterilisations
            print("************", monthly_sterilisations, "*************")

            # Create a copy of current dataframe to modify
            intervention_df_sp = current_df_sp.copy()
            intervention_df_pop = current_df_pop.copy()

            # Calculate population for each month of intervention
            intervention_data = runge_kutta(**intervention_parameters)
            intervention_data_sp = intervention_data[5]
            intervention_data_pop = intervention_data[6]

            # Isolate sterilisation props and population data
            intervention_sterilisation_props = intervention_data_sp.iloc[1:, 1]
            intervention_population_total = intervention_data_pop.iloc[1:, 1]
            intervention_population_sterilised = intervention_data_pop.iloc[1:, 2]

            print(intervention_data_pop)
            print(intervention_population_total)
            print(intervention_population_sterilised)

            # Replace intervention timeframe with intervention data
            intervention_df_sp["Sterilisation Proportion"][
            start_month + 1:start_month + intervention_parameters.get("desired_t")
                            + 1] = intervention_sterilisation_props
            intervention_df_pop["Total Population"][
            start_month + 1:start_month + intervention_parameters.get("desired_t")
                            + 1] = intervention_population_total
            intervention_df_pop["Sterilised Population"][
            start_month + 1:start_month + intervention_parameters.get("desired_t")
                            + 1] = intervention_population_sterilised

            # Set post intervention parameters using what we've calculated
            post_intervention_parameters["N0"] = intervention_df_pop["Total Population"].iloc[
                start_month + intervention_parameters.get("desired_t")]
            post_intervention_parameters["initial_S_N"] = intervention_df_sp["Sterilisation Proportion"].iloc[
                start_month + intervention_parameters.get("desired_t")]
            post_intervention_parameters["lower_S_N"] = intervention_parameters.get("lower_S_N")
            post_intervention_parameters["upper_S_N"] = intervention_parameters.get("upper_S_N")
            post_intervention_parameters["n_s"] = kwargs.get("n_s")
            post_intervention_parameters["desired_t"] = kwargs.get(
                "desired_t") - start_month - intervention_parameters.get("desired_t")
            post_intervention_parameters["k"] = kwargs.get("k")
            post_intervention_parameters["p_f"] = kwargs.get("p_f")
            post_intervention_parameters["s_a"] = kwargs.get("s_a")
            post_intervention_parameters["s_i"] = kwargs.get("s_i")
            post_intervention_parameters["l"] = kwargs.get("l")
            post_intervention_parameters["r_r"] = kwargs.get("r_r")
            post_intervention_parameters["m"] = kwargs.get("m")

            # Calculate population for each month post intervention
            post_intervention_data = runge_kutta(**post_intervention_parameters)
            post_intervention_data_sp = post_intervention_data[5]
            post_intervention_data_pop = post_intervention_data[6]

            post_intervention_sp = post_intervention_data_sp.iloc[1:, 1]
            post_intervention_pop_total = post_intervention_data_pop.iloc[1:, 1]
            post_intervention_pop_sterilised = post_intervention_data_pop.iloc[1:, 2]

            # Replace post-intervention timeframe with post-intervention data
            intervention_df_sp.iloc[start_month+intervention_parameters.get("desired_t")+1:kwargs.get("desired_t")+1, 1] = post_intervention_sp
            intervention_df_pop.iloc[start_month+intervention_parameters.get("desired_t")+1:kwargs.get("desired_t")+1, 1] = post_intervention_pop_total
            intervention_df_pop.iloc[start_month+intervention_parameters.get("desired_t")+1:kwargs.get("desired_t")+1, 2] = post_intervention_pop_sterilised

            return intervention_df_sp, intervention_df_pop, monthly_sterilisations

        else:
            print("I didn't understand that. Could you try again?")

def undo(models):
    models = models[:-1]
    return models


### Code
# Choose values for equations to test plotting and estimation
t0 = 0                 # Initial time
N0 = 2500              # Initial total population
initial_S_N = 0.8      # Initial sterilisation rate
desired_t = 24         # Time over which to model (in months)
h = 1                  # Step size for Runge Kutta method (1 recommended)
k = 10000              # Carrying capacity of environment
p_f = 0.5              # Proportion of population that are female (0-1)
s_a = 0.99             # Adult survivability rate (0-1)
s_i = 0.4              # Infant survivability rate (0-1)
l = 6                  # Reproduction rate/average litter size
r_r = 2                # Average number of litters (yearly)
m = 0                  # Net migration rate (positive into environment, negative out) (0-1)
n_s = 0                # Number of sterilisations taking place each month
desired_S_N = 0.8      # Desired sterilisation proportion
lower_S_N = 0.5        # Lower limit for sterilisation rate at which intervention is needed
upper_S_N = 0.8        # Upper limit for sterilisation rate we want to reach with intervention


# Create a dictionary to keep arguments/parameters of model in
initial_parameters = {
    "t0": 0,
    "N0": None,
    "initial_S_N": None,
    "desired_S_N": None,
    "lower_S_N": None,
    "upper_S_N": None,
    "n_s": 0,
    "desired_t": None,
    "h": 1,
    "k": 100000,
    "p_f": None,
    "s_a": 1,
    "s_i": 0.4,
    "l": None,
    "r_r": 2,
    "m": 0
}

intervention_parameters = {
    "t0": 0,
    "N0": None,
    "initial_S_N": None,
    "desired_S_N": None,
    "lower_S_N": None,
    "upper_S_N": None,
    "n_s": 0,
    "desired_t": None,
    "h": 1,
    "k": 100000,
    "p_f": None,
    "s_a": 1,
    "s_i": 0.4,
    "l": None,
    "r_r": 2,
    "m": 0
}

post_intervention_parameters = {
    "t0": 0,
    "N0": None,
    "initial_S_N": None,
    "desired_S_N": None,
    "lower_S_N": None,
    "upper_S_N": None,
    "n_s": 0,
    "desired_t": None,
    "h": 1,
    "k": 100000,
    "p_f": None,
    "s_a": 1,
    "s_i": 0.4,
    "l": None,
    "r_r": 2,
    "m": 0
}


# NEW MAIN CODE

# Create empty lists in which we'll store dfs
model_iterations = []

# First we want to use inputs to create initial model with data stored in df
# Find initial conditions of model
#initial_parameters["N0"] = int(input("What is the initial total number of dogs in the population?: ").strip())
#initial_parameters["initial_S_N"] = float(input("What is the initial sterilisation proportion as a decimal? (e.g. 0.4 for 40% of "
#                                      "population sterilised): ").strip())
#initial_parameters["n_s"] = int(input("How many sterilisations per month are currently taking place? (If none, enter 0): ").strip())
#initial_parameters["desired_t"] = int(input("How many months would you like to model over?: ").strip())
#initial_parameters["k"] = int(input("What is the carrying capacity of the environment? (If unknown, enter 100000): ").strip())
#initial_parameters["p_f"] = float(input("What proportion of the population are female (as a decimal)?: ").strip())
#initial_parameters["s_i"] = float(input("What is the infant survivability rate (as a decimal)?: ").strip())
#initial_parameters["l"] = int(input("What is the average litter size for this area?: ").strip())
#initial_parameters["r_r"] = int(input("On average, how many litters per year will a female dog have?: ").strip())
#initial_parameters["m"] = float(input("What is the net proportion of dogs (from whole population) that migrate into the area? "
#                            "(Please note that if the net migration is out of the community, you'll need to enter a minus "
#                            "sign in front of the value): "))


# Find initial conditions of model
initial_parameters["N0"] = 1200
initial_parameters["initial_S_N"] = 0.8
initial_parameters["n_s"] = 0
initial_parameters["desired_t"] = 18
initial_parameters["k"] = 100000
initial_parameters["p_f"] = 0.5
initial_parameters["s_i"] = 0.4
initial_parameters["l"] = 6
initial_parameters["r_r"] = 2
initial_parameters["m"] = 0


# Calculate model here and add dataframe to list
initial_model = runge_kutta(**initial_parameters)
df_sp = initial_model[5]
df_pop = initial_model[6]
model_iterations.append(initial_model)

# Plot model
graph_runge_kutta(df_sp, df_pop)

# On click of add_intervention
# Run function which spits out a new dataframe and add to list
intervention_model = add_intervention(df_sp, df_pop, **initial_parameters)
df_int_sp = intervention_model[0]
print(df_int_sp)
df_int_pop = intervention_model[1]
print(df_int_pop)
model_iterations.append(intervention_model)

# Add another
#another_model = add_intervention(df_int_sp, df_int_pop, **initial_parameters)
#df_int_two_sp = another_model[0]
#df_int_two_pop = another_model[1]
#model_iterations.append(another_model)

print(model_iterations)
print(len(model_iterations))

#model_iterations = undo(model_iterations)

# Plot last model in list
graph_runge_kutta(model_iterations[-1][0], model_iterations[-1][1])

# On click of go_back
# Remove last dataframe from list and plot new last df in list