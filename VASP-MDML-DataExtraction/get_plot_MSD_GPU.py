#!/usr/bin/env python3

# AUTHOR AsifEM2R
# get MSD from XDATCAR using NVIDIA GPU

import cupy as cp
import pandas as pd
import plotly.express as px

def read_xdatcar(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    header_lines = 7
    natoms = sum(map(int, lines[6].split()))
    positions = []

    for i, line in enumerate(lines[header_lines:], start=header_lines):
        if "Direct configuration=" in line:
            step_positions = []
            for j in range(natoms):
                coords = list(map(float, lines[i + 1 + j].split()))
                step_positions.append(coords)
            positions.append(step_positions)

    return cp.array(positions)

# Load positions to GPU
positions = read_xdatcar("XDATCAR")

# Compute MSD on GPU
initial_positions = positions[0]
displacements = positions - initial_positions
squared_displacements = cp.sum(displacements**2, axis=2)
msd = cp.mean(squared_displacements, axis=1)

# Move data back to CPU for saving/plotting
msd_cpu = cp.asnumpy(msd)
time_fs = cp.asnumpy(cp.arange(len(msd)) * 3)

# Save and plot
df = pd.DataFrame({"Time_fs": time_fs, "MSD_Angstrom2": msd_cpu})
df.to_csv("msd_data.csv", index=False)

fig = px.line(df, x="Time_fs", y="MSD_Angstrom2", title="MSD (GPU-accelerated)",
              labels={"Time_fs": "Time (fs)", "MSD_Angstrom2": "MSD (Å²)"})
fig.write_html("msd_plot.html")
