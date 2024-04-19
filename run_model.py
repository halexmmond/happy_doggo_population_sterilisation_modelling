from equation_functions import dNdt
from equation_functions import runge_kutta

# Choose values for equations to test plotting and estimation
t0 = 0            # Initial time
N0 = 2000         # Initial total population
desired_t = 18    # Time over which to model
h = 1             # Step size for Runge Kutta method (1 recommended)
k = 10000         # Carrying capacity of environment
p_f = 0.5         # Proportion of population that are female
s_a = 0.99        # Adult survivability rate
s_i = 0.4         # Infant survivability rate
l = 6             # Average litter size
r_r = 2/12        # Average number of litters (monthly)
S0 = 0            # Initial number of sterilised dogs
m = 0             # Net migration rate (positive into environment, negative out)
n_s = 250         # Number of sterilisations taking place each month
# N_t             # Population at time t
# S_t             # Number of sterilised dogs at time t

# add more here for sterilisation proportion when new code is added

runge_kutta(t0=t0, N0=N0, desired_t=desired_t, h=h, k=k, p_f=p_f, s_a=s_a, s_i=s_i, l=l, r_r=r_r, S0=S0, m=m, n_s=n_s)
