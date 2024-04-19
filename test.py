from equation_functions import dNdt
from equation_functions import runge_kutta

# Choose values for equations to test plotting and estimation
t0 = 0
N0 = 2000
desired_t = 18
h = 1
k = 10000
p_f = 0.5
s_a = 0.99
s_i = 0.4
l = 6
r_r = 2/12
S0 = 0
m = 0
n_s = 250

runge_kutta(t0=t0, N0=N0, desired_t=desired_t, h=h, k=k, p_f=p_f, s_a=s_a, s_i=s_i, l=l, r_r=r_r, S0=S0, m=m, n_s=n_s)
