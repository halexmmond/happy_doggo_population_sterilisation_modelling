import pandas as pd
import matplotlib.pyplot as plt

t = [0, 1, 2, 3, 4, 5, 6, 7, 8]
S_N = [0.5, 0.6, 0.8, 0.8, 0.7, 0.4, 0.2, 0.3, 0.5]

df = pd.DataFrame()
df["Month"] = t
df["SR"] = S_N

df_sr = pd.DataFrame()
df_sr["Month"] = [2, 3, 4, 5]
df_sr["SR"] = [0.9, 0.9, 0.9, 0.9]

start_month = 2
duration = 3

print(df)
print(df_sr)

# Create a copy of current dataframe to modify
intervention_df = df.copy()

# Isolate sterilisation rates
intervention_sterilisation_rates = df_sr.iloc[:, 1]

print(intervention_sterilisation_rates)

intervention_df["SR"][start_month:start_month+duration+1] = intervention_sterilisation_rates

print(intervention_df)


params = {
    "t0": None
}

print(params)

print(params.get("t0"))

params["t0"] = 5

print(params.get("t0"))



# Find initial conditions of model
print(int(input("What is the initial total number of dogs in the population?: ").strip()))

