# Instead of creating a differential model and using the Runge Kutta method to solve it,
# we'll try a model that calculates the population step by step, which will make it easier to
# incorporate the two-month gestation period for dogs before puppies are born

# This file will follow the same steps as oscillating_model.py but with the step model instead

### Functions

# Calculate number of sterilised dogs for initial calculation from initial sterilisation rate
# Number of sterilised dogs
def num_sterilised_dogs(initial_S_N, N0, s_a, m, n_s, t):
    S0 = N0 * initial_S_N
    S_t = ((n_s + ((m + s_a - 1) * S0)) * t) + S0
    return S_t


### Code

# We need to create lists I'm thinking to store the data in
# Try this way instead of dataframes and see what happens - should still be able to replace at certain index

N0 = 1000
initial_ster_prop = 0.8
p_f = 0.5
s_i = 0.4
lifespan = 72 # lifespan in months
p_lifespan = 0.24 # proportion of dogs that reach lifespan
initial_p_adult = 0
l = 6
r_r = 2
m_net_in = 0
k = 100000
t0 = 0

N_list = [N0]
S_list = [N0 * initial_ster_prop]
ster_prop_list = [initial_ster_prop] # Note we are currently saying equally likely for adult and puppy to get sterilised
N_adult_list = [N0 * initial_p_adult]
N_puppy_list = [N0 * (1 - initial_p_adult)]
monthly_births_list = [0]
t_list = [t0]
p_adult_list = [initial_p_adult]

### Functions

# Adult monthly survivability
def adult_monthly_survivability(p_lifespan, lifespan):
    s_a_month = 1 - ((1 - p_lifespan) / (lifespan - 12))
    return s_a_month

# Monthly births
# i.e. total pop two months ago x logistic factor x female prop x monthly adult survivability
# x infant survivability x average litter size x litters per month x adult prop 2 months ago
# x sterilisation proportion 2 months ago (assuming same prop sterilised under and over 1)
def monthly_births(k, p_f, s_a_month, s_i, l, r_r, N_total, p_adult, S_total):
    births = N_total * ((k - N_total) / k) * p_f * s_a_month * s_i * l * (r_r/12) * p_adult * (1 - (S_total/N_total))
    # Note that this is the formula we are calculating
    # monthly_births = N_list[t-2] * ((k - N_list[t-2]) / k) * p_f * s_a_month * s_i * l * (r_r/12) * p_adult_list[t-2] * (1 - (S_list[t-2]/N_list[t-2]))
    # We need to make sure we have the correct index for t
    return births

# Monthly deaths
# i.e. total pop one month ago x adult prop one month ago x monthly adult mortality +
# pups born lifespan-12 months ago that have survived until lifespan months, but look at deaths monthly
def monthly_deaths(N_total, p_adult, s_a_month, monthly_births, lifespan):
    deaths = (N_total * p_adult * (1 - s_a_month)) + (monthly_births * (1 - ((1 - s_a_month) * lifespan)))
    # Note that this is the formula we are calculating
    # monthly_deaths = (N_list[t-1] * p_adult_list[t-1] * (1 - s_a_month)) + (monthly_births_list[t - lifespan] * (1 - ((1 - s_a_month) * lifespan)))
    # We need to make sure we have the correct index for t
    return deaths

# Monthly net migration
def monthly_net_migration_in(N_total, m_net_in):
    net_migration_in = N_total * m_net_in
    # Note that this is the formula we are calculating
    # net_migration_in = N_list[t-1] * m_net_in # This is for m as a percentage of total population
    # We need to make sure we have the correct index for t
    return net_migration_in

# Monthly change in population
def dN_dt(births, deaths, net_migration_in):
    dN_dt = births - deaths + net_migration_in
    return dN_dt



## Keeping track of adult/puppy population

# Total puppy population
N_puppy = N_puppy_list[t-1] + monthly_births[t] - monthly_births_list[t-12] + (net_migration_in[t] * (1 - p_adult_list[t-1]))

# Total adult population
N_adult = N_adult_list[t-1] + monthly_births_list[t-12] - monthly_deaths[t] + (net_migration_in[t] * p_adult_list[t-1])