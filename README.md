<pre>
Yo yo yo herbernd.

Here are the variables I currently have as inputs for modelling either the total population at time t, or the number of sterilised dogs at time t.
All rates are in months and here are some examples.

# Choose values for equations to test plotting and estimation
t0 = 0            # Initial time  
N0 = 2000         # Initial total population  
desired_t = 18    # Time over which to model (in months)  
h = 1             # Step size for Runge Kutta method (1 recommended)  
k = 10000         # Carrying capacity of environment  
p_f = 0.5         # Proportion of population that are female (0-1)  
s_a = 0.99        # Adult survivability rate (0-1)  
s_i = 0.4         # Infant survivability rate (0-1)  
l = 6             # Average litter size  
r_r = 2/12        # Average number of litters (monthly)  
S0 = 0            # Initial number of sterilised dogs  
m = 0             # Net migration rate (positive into environment, negative out) (0-1)  
n_s = 250         # Number of sterilisations taking place each month  
  
N_t               # Total population at time t  
S_t               # Number of sterilised dogs at time t  
</pre>
