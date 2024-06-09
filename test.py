import matplotlib.pyplot as plt
# dN_dt = N * (p_f * a)

# Want to plot a against t
total_pop_list = [1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000]
adult_count = [600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600]
#adult_count = [40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40]
pup_count = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
adult_proportion_list = [0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6]
pups_born_list = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
t_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
p_f = 0.5
s_i = 0.4
s_a = 0.7777777
k = 40000

for t in range(12, 400):

    pups_born = total_pop_list[t-2] * ((k - total_pop_list[t-1])/k) * p_f * s_i * s_a * adult_proportion_list[t-2]
    pups_born_list.append(pups_born)

    total_pups = pup_count[t-1] + pups_born - pups_born_list[t-12]
    pup_count.append(total_pups)

    if t-96 < 0:
        eight_year_olds = 0
    else:
        eight_year_olds = pups_born_list[t-96]

    total_adults = adult_count[t-1] + pups_born_list[t-12] - eight_year_olds
    adult_count.append(total_adults)

    total_pop = total_pop_list[t-1] + pups_born - eight_year_olds
    total_pop_list.append(total_pop)

    adult_proportion = total_adults/total_pop
    adult_proportion_list.append(adult_proportion)

    t_list.append(t+1)

plt.figure(figsize=(10, 6))
plt.plot(t_list, adult_proportion_list, label='', marker='x')
plt.xlabel('Month')
plt.ylabel('Adult Proportion')
plt.axhline(0, color='black', linewidth=0.5)
plt.axvline(0, color='black', linewidth=0.5)
plt.grid(True)
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(t_list, total_pop_list, label='', marker='x')
plt.xlabel('Month')
plt.ylabel('Total Proportion')
plt.axhline(0, color='black', linewidth=0.5)
plt.axvline(0, color='black', linewidth=0.5)
plt.grid(True)
plt.legend()
plt.show()


print(adult_proportion_list[-1])
print(total_pop_list[-1])