# Instead of creating a differential model and using the Runge Kutta method to solve it,
# we'll try a model that calculates the population step by step, which will make it easier to
# incorporate the two-month gestation period for dogs before puppies are born

# This file will follow the same steps as oscillating_model.py but with the step model instead
import matplotlib.pyplot as plt

### Functions

# Adult monthly survivability
def adult_monthly_survivability(lifespan):
    s_a_month = 1 - (1 / (lifespan - 12))
    return s_a_month

# Monthly births
# i.e. total pop two months ago x logistic factor x female prop x monthly adult survivability
# x infant survivability x average litter size x litters per month x adult prop 2 months ago
# x sterilisation proportion 2 months ago (assuming same prop sterilised under and over 1)
def monthly_births(k, p_f, s_a_month, s_i, l, r_r, N_total_start, prop_adult_start, S_total_start):
    births = N_total_start * ((k - N_total_start) / k) * p_f * s_a_month * s_i * l * (r_r/12) * prop_adult_start * (1 - (S_total_start/N_total_start))
    # Note that this is the formula we are calculating
    # monthly_births = N_list[t-2] * ((k - N_list[t-2]) / k) * p_f * s_a_month * s_i * l * (r_r/12) * p_adult_list[t-2] * (1 - (S_list[t-2]/N_list[t-2]))
    # We need to make sure we have the correct index for t
    return births

# Monthly deaths
# i.e. total pop one month ago x adult prop one month ago x monthly adult mortality +
# pups born lifespan-12 months ago that have survived until lifespan months, but look at deaths monthly
def monthly_deaths(N_total_start, prop_adult_start, s_a_month):
    deaths = (N_total_start * prop_adult_start * (1 - s_a_month))
    # Note that this is the formula we are calculating
    # monthly_deaths = (N_list[t-1] * p_adult_list[t-1] * (1 - s_a_month)) + (monthly_births_list[t - lifespan] * (1 - ((1 - s_a_month) * lifespan)))
    # We need to make sure we have the correct index for t
    return deaths

# Monthly net migration
def monthly_net_migration_in(N_total_start, m_net_in):
    net_migration_in = N_total_start * m_net_in
    # Note that this is the formula we are calculating
    # net_migration_in = N_list[t-1] * m_net_in # This is for m as a percentage of total population
    # We need to make sure we have the correct index for t
    return net_migration_in

# Monthly change in population
def dN_dt(births, deaths, net_migration_in):
    dNdt = births - deaths + net_migration_in
    return dNdt

# Total number of sterilised dogs in population
def num_sterilised_dogs(ster_prop_start, S_total_start, n_s, monthly_deaths):
    S_t = S_total_start - (monthly_deaths * ster_prop_start) + n_s
    return S_t # CHECK THIS FORMULA AND ADD TO NOTES

# Total puppy population
def puppy_population(N_puppy_one_month_ago, births_t, births_t_one_year_ago, net_migration_in, prop_adult_start):
    N_puppy_current = N_puppy_one_month_ago + births_t - births_t_one_year_ago + (net_migration_in * (1 - prop_adult_start))
    # Note that this is the formula we are calculating
    # N_puppy = N_puppy_list[t-1] + monthly_births[t] - monthly_births_list[t-12] + (net_migration_in[t] * (1 - p_adult_list[t-1]))
    # We need to make sure we have the correct index for t
    return N_puppy_current

# Total adult population
def adult_population(N_adult_one_month_ago, births_t_one_year_ago, deaths_t, net_migration_in, prop_adult_start):
    N_adult_current = N_adult_one_month_ago + births_t_one_year_ago - deaths_t + (net_migration_in * prop_adult_start)
    # Note that this is the formula we are calculating
    # N_adult = N_adult_list[t-1] + monthly_births_list[t-12] - monthly_deaths[t] + (net_migration_in[t] * p_adult_list[t-1])
    # We need to make sure we have the correct index for t
    return N_adult_current

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
    desired_t = 13
    h = kwargs.get("h")  # chosen step size
    k = kwargs.get("k")
    p_a = kwargs.get("p_a")
    p_f = kwargs.get("p_f")
    s_a = 1 - (1 / (kwargs.get("lifespan") - 12))
    s_i = kwargs.get("s_i")
    l = kwargs.get("l")
    r_r = kwargs.get("r_r")
    m = kwargs.get("m")
    S_t = 0

    # Results
    t_values = [t]
    N_t_values = [N_t]
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
        pop_change = N_t_values[t] - N_t_values[t+1]

        previous_births = round(reverse_births_RK4(N_t=N_t_values[t], k=k, p_a=p_a, s_a=s_a, m=m, dNdt=pop_change))

        births_values.append(previous_births)

    return t_values[:-1], N_t_values[:-1], births_values

# Plot our results
def graph_results(t_list, N_list, adult_list, puppy_list):
    months = t_list
    total_pop = N_list
    adult_pop = adult_list
    puppy_pop = puppy_list

    # Show sterilisation proportion
    plt.figure(figsize=(10, 6))
    plt.plot(months, total_pop, label='Total population', marker='*')
    plt.plot(months, adult_pop, label='Adult population', marker='x')
    plt.plot(months, puppy_pop, label='Puppy population', marker='o')
    plt.xlabel('Month')
    plt.ylabel('Population')
    plt.title('')
    plt.grid(True)
    plt.legend()
    plt.show()

### Code

# Empty list in which to store full models
model_iterations = []

# Variable dictionary to make it easier when calling a function
initial_parameters = {
    "t0": 0,
    "N0": None,
    "initial_ster_prop": None,
    "desired_ster_prop": None,
    "lower_ster_prop": None,
    "upper_ster_prop": None,
    "initial_prop_adult": None,
    "n_s": 0,
    "desired_t": None,
    "h": 1,
    "k": 100000,
    "p_f": None,
    "p_a": None,
    "s_a": 1,
    "s_i": 0.4,
    "l": None,
    "r_r": 2,
    "m": 0,
    "lifespan": None
}

# Test variables
#initial_parameters["N0"] = 1000
#initial_parameters["desired_t"] = 13
#initial_ster_prop = 0.8
#p_f = 0.5
#s_i = 0.4
#lifespan = 72 # lifespan in months
#p_lifespan = 0.24 # proportion of dogs that reach lifespan
#initial_prop_adult = 0
#l = 6
#r_r = 2
#m_net_in = 0
#k = 100000
#t0 = 0 # don't need this?
#t_duration = 5
#s_a_month = 0.99
#n_s = 0

# Calculate monthly population from t=-12 up to t=0
# Set initial variables for reverse process
initial_parameters["N0"] = 1200
initial_parameters["p_f"] = 0.5
initial_parameters["s_i"] = 0.4
initial_parameters["lifespan"] = 72
initial_parameters["p_a"] = 1
initial_parameters["l"] = 6
initial_parameters["r_r"] = 2
initial_parameters["m"] = 0
initial_parameters["k"] = 100000
initial_parameters["initial_ster_prop"] = 0
initial_parameters["initial_prop_adult"] = 1
initial_parameters["desired_t"] = 18
initial_parameters["s_a"] = 1 - (1 / (initial_parameters.get("lifespan") - 12))
initial_parameters["n_s"] = 0


t_previous_12 = [-12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0]

N_previous_12 = runge_kutta(**initial_parameters)[1]
N_previous_12 = list(reversed(N_previous_12))

births_previous_12 = runge_kutta(**initial_parameters)[2]
births_previous_12 = list(reversed(births_previous_12))

# Now want to fill these lists with our 12 months previous data
# Lists to store data we want to keep track of
N_total_list = N_previous_12
S_total_list = ([0] * 11) + [(initial_parameters.get("N0") * initial_parameters.get("initial_ster_prop")) * 2]
ster_prop_list = ([0] * 11) + [(initial_parameters.get("initial_ster_prop")) * 2] # Note we are currently saying equally likely for adult and puppy to get sterilised
N_adult_list = (["NaN"] * 11) + ([round(initial_parameters.get("N0") * initial_parameters.get("initial_prop_adult"))] * 2)
N_puppy_list = (["NaN"] * 11) + ([round(initial_parameters.get("N0") * (1 - initial_parameters.get("initial_prop_adult")))] * 2)
monthly_births_list = births_previous_12
monthly_deaths_list = [0] * 13
t_list = t_previous_12
prop_adult_list = ([1] * 12) + [initial_parameters.get("initial_prop_adult")]

#print("t_list", t_list)
#print("N", N_total_list)
#print("S", S_total_list)
#print("STER", ster_prop_list)
#print("adult", N_adult_list)
#print("puppy", N_puppy_list)
#print("births", monthly_births_list)
#print("deaths", monthly_deaths_list)
#print("propadult", prop_adult_list)


######### We have initial conditions for the model, now we want to calculate the following
# Total population next month
# Number of puppies
# Number of adults
# Adult proportion
# Sterilisation proportion

# If we want a model that calculates the population after 5 months
for i in range(12, initial_parameters.get("desired_t")+13):

    # Define initial variables
    if i == 12:
        t = 1
    n_s = initial_parameters.get("n_s")
    k = initial_parameters.get("k")
    p_f = initial_parameters.get("p_f")
    s_a = initial_parameters.get("s_a")
    s_i = initial_parameters.get("s_i")
    l = initial_parameters.get("l")
    r_r = initial_parameters.get("r_r")
    m_net_in = initial_parameters.get("m")
    lifespan = initial_parameters.get("lifespan")

    # Calculate number of births
    births = monthly_births(k=k, p_f=p_f, s_a_month=s_a, s_i=s_i, l=l, r_r=r_r, N_total_start=N_total_list[i-2],
                            prop_adult_start=prop_adult_list[i-2], S_total_start=S_total_list[i-2])

    # Calculate number of deaths
    deaths = monthly_deaths(N_total_start=N_total_list[i-1], prop_adult_start=prop_adult_list[i-1], s_a_month=s_a)

    # Calculate net migration
    net_migration_in = monthly_net_migration_in(N_total_start=N_total_list[i-1], m_net_in=m_net_in)

    # Calculate population change
    dNdt = dN_dt(births=births, deaths=deaths, net_migration_in=net_migration_in)

    # Find new total population for this month
    new_N_total = N_total_list[i-1] + dNdt

    # Calculate number of sterilised dogs
    new_S_total = num_sterilised_dogs(ster_prop_start=ster_prop_list[i-1], S_total_start=S_total_list[i-1],
                                      n_s=n_s, monthly_deaths=deaths)

    # Calculate number of adults
    num_adults = adult_population(N_adult_one_month_ago=N_adult_list[i-1], births_t_one_year_ago=monthly_births_list[i-12],
                                  deaths_t=deaths, net_migration_in=net_migration_in, prop_adult_start=prop_adult_list[i-1])

    # Calculate number of puppies
    num_puppies = puppy_population(N_puppy_one_month_ago=N_puppy_list[i-1], births_t=births, births_t_one_year_ago=monthly_births_list[i-12],
                                   net_migration_in=net_migration_in, prop_adult_start=prop_adult_list[i-1])

    if births < 0:
        births = 0

    if deaths < 0:
        deaths = 0

    if new_N_total < 0:
        new_N_total = 1

    if num_adults < 0:
        num_adults = 0

    if num_puppies < 0:
        num_puppies = 0

    # Update lists with this month's values
    N_total_list.append(round(new_N_total))
    S_total_list.append(round(new_S_total))
    ster_prop_list.append(round(new_S_total/new_N_total, 2))
    N_adult_list.append(round(num_adults))
    N_puppy_list.append(round(num_puppies))
    monthly_births_list.append(round(births))
    monthly_deaths_list.append(round(deaths))
    prop_adult_list.append(round(num_adults/new_N_total, 2))
    t_list.append(t)

    # Update month counter
    t += 1

print(t_list)
print(N_total_list)
print(S_total_list)
print(ster_prop_list)
print(N_adult_list)
print(N_puppy_list)
print(monthly_births_list)
print(monthly_deaths_list)
print(prop_adult_list)

graph_results(t_list=t_list[12:], N_list=N_total_list[12:], adult_list=N_adult_list[12:], puppy_list=N_puppy_list[12:])