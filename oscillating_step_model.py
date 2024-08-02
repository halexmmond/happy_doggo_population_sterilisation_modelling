# Instead of creating a differential model and using the Runge Kutta method to solve it,
# we'll try a model that calculates the population step by step, which will make it easier to
# incorporate the two-month gestation period for dogs before puppies are born

# This file will follow the same steps as oscillating_model.py but with the step model instead


### Functions

# Adult monthly survivability
def adult_monthly_survivability(p_lifespan, lifespan):
    s_a_month = 1 - ((1 - p_lifespan) / (lifespan - 12))
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
def monthly_deaths(N_total_start, prop_adult_start, s_a_month, monthly_births, lifespan):
    deaths = (N_total_start * prop_adult_start * (1 - s_a_month)) + (monthly_births * (1 - ((1 - s_a_month) * lifespan)))
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


### Code

# Create a dictionary to keep arguments/parameters of model in
initial_parameters = {
    "p_lifespan": None,
    "lifespan": None,
    "k": None,
    "p_f": None,
    "s_a_month": None,
    "s_i": None,
    "l": None,
    "r_r": None,
    "N_total": None, # don't think we need this in initial conditions
    "initial_prop_adult": None,
    "S_total": None, # don't need this here either
    "m_net_in": None,
    "initial_ster_prop": None,
    "N0_total": None
    "t_duration": None # don't need this here either
}


# Empty list in which to store full models
model_iterations = []

# Pick initial parameters for model
initial_parameters["p_lifespan"] = 0
initial_parameters["lifespan"] = 0
initial_parameters["k"] = 0
initial_parameters["p_f"] = 0
initial_parameters["s_a_month"] = 0
initial_parameters["s_i"] = 0
initial_parameters["l"] = 0
initial_parameters["r_r"] = 0
initial_parameters["N_total"] = 0 # do we need this in initial conditions? I don't think so
initial_parameters["initial_prop_adult"] = 0
initial_parameters["S_total"] = 0 # don't need this here either
initial_parameters["m_net_in"] = 0
initial_parameters["initial_ster_prop"] = 0
initial_parameters["N0_total"] = 0
initial_parameters["t_duration"] = 0 # don't need this here either


########

# We need to create lists I'm thinking to store the data in
# Try this way instead of dataframes and see what happens - should still be able to replace at certain index

N0_total_initial = 1000
initial_ster_prop = 0.8
p_f = 0.5
s_i = 0.4
lifespan = 72 # lifespan in months
p_lifespan = 0.24 # proportion of dogs that reach lifespan
initial_prop_adult = 0
l = 6
r_r = 2
m_net_in = 0
k = 100000
t0 = 0 # don't need this?
t_duration = 5
s_a_month = 0.99
n_s = 0

N_total_list = [N0_total_initial]
S_total_list = [N0_total_initial * initial_ster_prop]
ster_prop_list = [initial_ster_prop] # Note we are currently saying equally likely for adult and puppy to get sterilised
N_adult_list = [N0_total_initial * initial_prop_adult]
N_puppy_list = [N0_total_initial * (1 - initial_prop_adult)]
monthly_births_list = [0]
monthly_deaths_list = []
t_list = [t0]
prop_adult_list = [initial_prop_adult]


######### We have initial conditions for the model, now we want to calculate the following
# Total population next month
# Number of puppies
# Number of adults
# Adult proportion
# Sterilisation proportion

# If we want a model that calculates the population after 5 months
for t in range(1, t_duration+1):

    # Calculate number of births
    births = monthly_births(k=k, p_f=p_f, s_a_month=s_a_month, s_i=s_i, l=l, r_r=r_r, N_total_start=N_total_list[t-2],
                            prop_adult_start=prop_adult_list[t-2], S_total_start=S_total_list[t-2])

    # Calculate number of deaths
    deaths = monthly_deaths(N_total_start=N_total_list[t-1], prop_adult_start=prop_adult_list[t-1], s_a_month=s_a_month,
                            monthly_births=monthly_births_list[t-lifespan], lifespan=lifespan)

    # Calculate net migration
    net_migration_in = monthly_net_migration_in(N_total_start=N_total_list[t-1], m_net_in=m_net_in)

    # Calculate population change
    dNdt = dN_dt(births=births, deaths=deaths, net_migration_in=net_migration_in)

    # Find new total population for this month
    new_N_total = N_total_list[t-1] + dNdt

    # Calculate number of sterilised dogs
    new_S_total = num_sterilised_dogs(ster_prop_start=ster_prop_list[t-1], S_total_start=S_total_list[t-1],
                                      n_s=n_s, monthly_deaths=deaths)

    # Calculate number of adults
    num_adults = adult_population(N_adult_one_month_ago=N_adult_list[t-1], births_t_one_year_ago=monthly_births_list[t-12],
                                  deaths_t=deaths, net_migration_in=net_migration_in, prop_adult_start=prop_adult_list[t-1])

    # Calculate number of puppies
    num_puppies = puppy_population(N_puppy_one_month_ago=N_puppy_list[t-1], births_t=births, births_t_one_year_ago=monthly_births_list[t-12],
                                   net_migration_in=net_migration_in, prop_adult_start=prop_adult_list[t-1])

    # Update lists with this month's values
    N_total_list.append(new_N_total)
    S_total_list.append(new_S_total)
    ster_prop_list.append(new_S_total/new_N_total)
    N_adult_list.append(num_adults)
    N_puppy_list.append(num_puppies)
    monthly_births_list.append(births)
    monthly_deaths_list.append(deaths)
    prop_adult_list.append(num_adults/new_N_total)