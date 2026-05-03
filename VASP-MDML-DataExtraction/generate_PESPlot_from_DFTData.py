#!/usr/bin/env python3

from ase.io import read, write
from pathlib import Path
from scipy.spatial.distance import pdist
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.interpolate import griddata
from umap import UMAP
from joblib import Parallel, delayed
from tqdm import tqdm

def pair_distance_histogram(atoms, cutoff=6.0, bins=20):
    dists = pdist(atoms.get_positions(), metric='euclidean')
    dists = dists[dists <= cutoff]
    if len(dists) == 0:
        return np.zeros(bins)
    hist, _ = np.histogram(dists, bins=bins, range=(0.0, cutoff))
    n = np.linalg.norm(hist)
    return hist / n if n > 0 else hist


base_dirs = [
    #Path("DFT_or_MLIP_data_in_EXTXYZ"),
]

index = -1  # Change to ":" or "::" as needed
images = []
files = []

for base_dir in base_dirs:
    for file in base_dir.rglob("CONTCAR_mace.extxyz"):
        files.append(file)
files = sorted(files)
for file in tqdm(files, colour='red'):
    data = read(file, index=index, format="extxyz")
    if isinstance(data, list):
        images.extend(data)
    else:
        images.append(data)
    write("train_data.extxyz", images, format="extxyz")
print(f"Found {len(files)} files and {len(images)} structures.")


# Compute descriptors in parallel
cutoff = 6.0
bins = 20
print("Computing histograms in parallel...")
descs = Parallel(n_jobs=-1, backend="loky")(
    delayed(pair_distance_histogram)(atoms, cutoff, bins) for atoms in tqdm(images)
)

# Compute energies
energies = [atoms.get_potential_energy() / len(atoms) for atoms in images]

X = np.vstack(descs)
E = np.array(energies)

## PCA reduction
#pc = PCA(n_components=2).fit_transform(X)
#pc1, pc2 = pc[:, 0], pc[:, 1]

# Apply UMAP for dimensionality reduction
print("Applying UMAP...")
umap_model = UMAP(n_components=2)
embedding = umap_model.fit_transform(X)

pc1, pc2 = embedding[:, 0], embedding[:, 1]

# Create grid for interpolation
grid_x, grid_y = np.meshgrid(
    np.linspace(pc1.min(), pc1.max(), 100),
    np.linspace(pc2.min(), pc2.max(), 100)
)

# Interpolate energies on grid
grid_E = griddata((pc1, pc2), E, (grid_x, grid_y), method='cubic')

# Plot smooth PES surface
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(grid_x, grid_y, grid_E, cmap="viridis", edgecolor='none')
ax.set_xlabel("UMAP-1 (Structural variation)")
ax.set_ylabel("UMAP-2 (Structural variation)")
ax.set_zlabel("Energy [eV/atom]")
fig.colorbar(surf, ax=ax, shrink=0.6, label="Energy")
plt.tight_layout()
plt.savefig("smooth_PES_umap.png", dpi=600, bbox_inches='tight')
plt.show()


plt.figure(figsize=(8, 6))
plt.contour(grid_x, grid_y, grid_E, levels=20, colors='black', linewidths=0.5)
plt.scatter(pc1, pc2, c=E, cmap="viridis", s=20, edgecolor='k')
plt.colorbar(label="Energy [eV/atom]")
plt.xlabel("UMAP-1")
plt.ylabel("UMAP-2")
plt.tight_layout()
plt.savefig("contour_lines_umap.png", dpi=600)


import plotly.graph_objects as go
fig = go.Figure(data=[go.Surface(x=grid_x, y=grid_y, z=grid_E, colorscale='Viridis')])
fig.update_layout(scene=dict(
    xaxis_title='UMAP-1',
    yaxis_title='UMAP-2',
    zaxis_title='Energy [eV/atom]'
))
fig.write_html("smooth_PES_interactive.html")
print("Interactive 3D plot saved as smooth_PES_interactive.html")
