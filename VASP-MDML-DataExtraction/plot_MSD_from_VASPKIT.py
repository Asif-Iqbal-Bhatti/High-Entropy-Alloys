#!/usr/bin/env python3

import pandas as pd
import plotly.express as px

# Replace with your actual filename
filename = "MSD.dat"
#filename = "DIFFUSION_COEFFICIENT.dat"

# Read the file, skipping comment lines
df = pd.read_csv(filename, sep=r'\s+', header=None, skiprows = 1)
print(df)

# Optional: Rename columns if known (example for MSD data)
df.columns = ["Time_fs", "x_MSD", "y_MSD", "z_MSD", "total_MSD", "sqrt_MSD"]

# Plot MSD components
fig = px.line(df, x="Time_fs", y=["x_MSD", "y_MSD", "z_MSD", "total_MSD", "sqrt_MSD"],
              labels={"value": "MSD (Å²)", "Time_fs": "Time (fs)", "variable": "Component"},
              title="Mean Squared Displacement (MSD) vs Time")
fig.write_html("msd_plot.html")



##################

# Define the filename
filename = "DIFFUSION_COEFFICIENT.dat"

# Read the file, skipping comment lines and using whitespace as delimiter
df2 = pd.read_csv(filename, sep=r'\s+', skiprows = 2, header=None)
print(df2)
# Assign column names based on the file structure
df2.columns = ["Time(fs)", "x-D(A^2)", "y-D(cm^2/s)", "z-D(cm^2/s)", "tot-D(cm^2/s)"]


# Plot diffusion coefficients
fig = px.line(df2, x="Time(fs)", y=["x-D(A^2)", "y-D(cm^2/s)", "z-D(cm^2/s)", "tot-D(cm^2/s)"],
              labels={"value": "Diffusion Coefficient (cm²/s)", "Time_fs": "Time (fs)", "variable": "Component"},
              title="Time-dependent Diffusion Coefficients")

# Save to HTML
fig.write_html("diffusion_plot.html")

