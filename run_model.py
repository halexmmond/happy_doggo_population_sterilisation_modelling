from equation_functions import dNdt
from equation_functions import runge_kutta
from equation_functions import growth_binary_search_n_s
from equation_functions import maintenance_binary_search_n_s
from equation_functions import graph_runge_kutta
import equation_functions

# Choose values for equations to test plotting and estimation
t0 = 0                 # Initial time
N0 = 1000              # Initial total population
desired_t = 24         # Time over which to model (in months)
h = 1                  # Step size for Runge Kutta method (1 recommended)
k = 10000              # Carrying capacity of environment
p_f = 0.5              # Proportion of population that are female (0-1)
s_a = 0.99             # Adult survivability rate (0-1)
s_i = 0.4              # Infant survivability rate (0-1)
l = 6                  # Reproduction rate/average litter size
r_r = 2/12             # Average number of litters (monthly)
S0 = 800               # Initial number of sterilised dogs
m = 0                  # Net migration rate (positive into environment, negative out) (0-1)
n_s = 100              # Number of sterilisations taking place each month
desired_S_N = 0.8      # Desired sterilisation proportion
# N_t                  # Population at time t
# S_t                  # Number of sterilised dogs at time t

# add more here for sterilisation proportion when new code is added

#print(runge_kutta(t0=t0, N0=N0, desired_t=desired_t, h=h, k=k, p_f=p_f, s_a=s_a, s_i=s_i, l=l, r_r=r_r, S0=S0, m=m, n_s=n_s))

graph_runge_kutta(t0=t0, N0=N0, desired_t=desired_t, h=h, k=k, p_f=p_f, s_a=s_a, s_i=s_i, l=l, r_r=r_r, S0=S0, m=m, n_s=n_s)

#print(growth_binary_search_n_s(t0=t0, N0=N0, desired_t=desired_t, h=h, k=k, p_f=p_f, s_a=s_a, s_i=s_i, l=l, r_r=r_r, S0=S0, m=m, desired_S_N=desired_S_N))

#print(maintenance_binary_search_n_s(t0=t0, N0=N0, desired_t=desired_t, h=h, k=k, p_f=p_f, s_a=s_a, s_i=s_i, l=l, r_r=r_r, S0=S0, m=m))
