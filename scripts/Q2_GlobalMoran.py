
"""
GEOG41160 Assignment 2, Question 2: Global Spatial Autocorrelation
Team: 2
ANSWERS

1. Distribution: Raw extension density was positively skewed
   (skewness = 1.39). We selected a square-root transformation,
   which reduced skewness to 0.35 and gave a more balanced
   distribution than log1p (skewness = -1.14).
   
2. Spatial weights: We used row-standardised Queen contiguity.
   EDs sharing an edge or a vertex were treated as neighbours.
   All 162 EDs had neighbours. Each ED had 2–11 neighbours,
   with a mean of 5.70.

3. Global Moran's I: Raw density I = 0.3319
   (permutation p = 0.0001); square-root-transformed density
   I = 0.3751 (permutation p = 0.0001). Both tests used
   9,999 permutations and random seed 42.

4. Interpretation: Both variables show significant positive
   spatial autocorrelation, indicating that similar extension-
   density values tend to cluster in neighbouring EDs. The
   square-root transformation slightly increased Moran's I,
   but did not change the overall conclusion.
"""

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from esda.moran import Moran
from libpysal.weights import Queen

## =========== SetupLoaded =========== ##
## ============ GOOD LUCK ============ ##

ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = ROOT / "outputs"  ## Output files path.
DENSITY_FILE = OUTPUTS / "DublinED_ExtensionDensity.gpkg"  ## Electoral division in Dublin and extension density.
FIGURE_FILE = OUTPUTS / "ExtensionDensity_Distribution.png"  ## Distribution histogram plot(subplots).

SEED = 42  ## Random seed.
PERMUTATIONS = 9999  ##Number of permutation tests.
DENSITY_COLUMN = "app_density_km2"

## Read ED zone vector dataset.
eds = gpd.read_file(DENSITY_FILE)

## Convert the density column to a numeric type. Values cannot be converted will be set as NaN.
eds[DENSITY_COLUMN] = pd.to_numeric(
    eds[DENSITY_COLUMN],
    errors="coerce"
)
valid_rows = np.isfinite(eds[DENSITY_COLUMN])
eds = eds.loc[valid_rows].copy()  ## Keep the limited and non-empty valid rows.

## Extract raw density series.
density = eds[DENSITY_COLUMN]
sqrt_density = np.sqrt(density)  ## Square root transformation.
log_density = np.log1p(density)  ## Log1p transformation.

## Print valid sample count and skewness after transformation.
print("EDs with valid density:", len(eds))
print(f"Raw density skewness: {density.skew():.2f}")
print(f"Square-root density skewness: {sqrt_density.skew():.2f}")
print(f"Log-transformed density skewness: {log_density.skew():.2f}")


## Create 1 row with 3 subplots, and side-by-side compare the three histograms.
fig, axes = plt.subplots(1, 3, figsize=(15, 4)  ## Compare the three distributions.

axes[0].hist(density, bins=20, edgecolor="black")
axes[0].set_title("Raw extension density")
axes[0].set_xlabel("Applications per km²")
axes[0].set_ylabel("Number of EDs")

axes[1].hist(sqrt_density, bins=20, edgecolor="black")
axes[1].set_title("Square-root-transformed density")
axes[1].set_xlabel("Square root of applications per km²")
axes[1].set_ylabel("Number of EDs")

axes[2].hist(log_density, bins=20, edgecolor="black")
axes[2].set_title("Log-transformed density")
axes[2].set_xlabel("log1p(applications per km²)")
axes[2].set_ylabel("Number of EDs")

## Auto-adjust subplots layout. 
fig.tight_layout()

## Export the image in following format.
fig.savefig(FIGURE_FILE, dpi=300, bbox_inches="tight")
plt.show()


## Build Qween contiguity weights.
weights = Queen.from_dataframe(eds, use_index=False)
weights.transform = "R"

## How to predicate neighbours.
neighbour_counts = [
    len(neighbours)
    for neighbours in weights.neighbors.values()
]

print("Number of EDs in weights:", weights.n)
print("Number of islands:", len(weights.islands))
print("Minimum neighbours:", min(neighbour_counts))
print("Maximum neighbours:", max(neighbour_counts))
print(f"Mean neighbours: {np.mean(neighbour_counts):.2f}")
## Summarise the statistics of spatial weight matrix.


np.random.seed(SEED)

raw_moran = Moran(
    density.to_numpy(),
    weights,
    permutations=PERMUTATIONS
)

## Compute Global Moran's I for the raw density variable.
sqrt_moran = Moran(
    sqrt_density.to_numpy(),
    weights,
    permutations=PERMUTATIONS
)

## Output Moran and permutation-based p-value.
print(f"Raw density Moran's I: {raw_moran.I:.4f}")
print(f"Raw density permutation p-value: {raw_moran.p_sim:.4f}")
print(f"Square-root density Moran's I: {sqrt_moran.I:.4f}")
print(f"Square-root density permutation p-value: {sqrt_moran.p_sim:.4f}")
print("Distribution figure saved to:", FIGURE_FILE)
