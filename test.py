import pandas as pd
import matplotlib.pyplot as plt

t = [0, 1, 2, 3, 4, 5, 6, 7, 8]
S_N = [0.5, 0.6, 0.8, 0.9, 0.7, 0.4, 0.2, 0.3, 0.5]

df = pd.DataFrame()
df["Month"] = t
df["SR"] = S_N

print(df)

# Show sterilisation proportion
plt.figure(figsize=(10, 6))
plt.plot(df["Month"], df["SR"], marker='x')
plt.xlabel('t')
plt.ylabel('SR')
plt.axhline(0, color='black', linewidth=0.5)
plt.axvline(0, color='black', linewidth=0.5)
plt.grid(True)
plt.legend()
plt.show()