
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

#########################################
## ============== Setup ============== ##
## ============ GOOD LUCK ============ ##
#########################################

ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = ROOT / "outputs"  ## Output files path.
DENSITY_FILE = OUTPUTS / "DublinED_ExtensionDensity.gpkg"  ## Electoral division in Dublin and extension density.
FIGURE_FILE = OUTPUTS / "ExtensionDensity_Distribution.png"  ## Distribution histogram plot(subplots).

SEED = 42  ## Random seed.
PERMUTATIONS = 9999  ##Number of permutation tests.
DENSITY_COLUMN = "app_density_km2"
# Use a fixed seed and a stated number of permutations for reproducibility.


###################################################
############1. DISTRIBUTION AND TRANSFORMATION #####
####################################################

eds = gpd.read_file(DENSITY_FILE)
## Load the ED polygons and extension-density values produced by Q1.

eds[DENSITY_COLUMN] = pd.to_numeric(
    eds[DENSITY_COLUMN],
    errors="coerce"
)
# Convert the density column to numeric values.
# Values that cannot be converted become NaN (Not a Number).

valid_rows = np.isfinite(eds[DENSITY_COLUMN])
excluded_count = int((~valid_rows).sum())
eds = eds.loc[valid_rows].copy().reset_index(drop=True)
## Keep the limited and non-empty valid rows.
## Apply the filter to the whole table so geometries and values stay aligned.
## Create a copy and give the remaining rows consecutive index labels.


density = eds[DENSITY_COLUMN]
## Extract raw density series.

if len(eds) < 3:
    raise ValueError("At least three EDs with valid density are required.")

if (density < 0).any():
    raise ValueError("Extension density must be non-negative.")

if density.nunique() < 2:
    raise ValueError("Moran's I requires variation in density.")
## Check that the retained data are suitable for the calculations below.

sqrt_density = np.sqrt(density)  ## Square root transformation.
log_density = np.log1p(density)  ## Log1p transformation.

print("EDs with valid density:", len(eds))
print(f"Raw density skewness: {density.skew():.2f}")
print(f"Square-root density skewness: {sqrt_density.skew():.2f}")
print(f"Log-transformed density skewness: {log_density.skew():.2f}")
# Positive skewness indicates a longer right tail.
# Compare skewness and histograms when selecting the transformation.


fig, axes = plt.subplots(1, 3, figsize=(15, 4))
# Create one figure with three subplots arranged in one row.
# axes[0], axes[1], and axes[2] refer to the three plotting areas.

axes[0].hist(density, bins=20, edgecolor="black")
axes[0].set_title("Raw extension density")
axes[0].set_xlabel("Applications per km²")
axes[0].set_ylabel("Number of EDs")
# Group raw density values into 20 intervals.
# Each bar shows the number of EDs within that density interval.

axes[1].hist(sqrt_density, bins=20, edgecolor="black")
axes[1].set_title("Square-root-transformed density")
axes[1].set_xlabel("Square root of applications per km²")
axes[1].set_ylabel("Number of EDs")

axes[2].hist(log_density, bins=20, edgecolor="black")
axes[2].set_title("Log-transformed density")
axes[2].set_xlabel("log1p(applications per km²)")
axes[2].set_ylabel("Number of EDs")
# Compare both transformations with the raw distribution.

fig.tight_layout()
fig.savefig(FIGURE_FILE, dpi=300, bbox_inches="tight")
plt.show()
# Adjust spacing, save the figure at 300 dots per inch, and display it.


#################################################################
###################### 2. SPATIAL WEIGHTS #######################
#################################################################


weights = Queen.from_dataframe(eds, use_index=False)
# Build Queen contiguity weights from the ED polygons.
# Neighbours share an edge or a vertex.
# use_index=False assigns positional IDs in the current row order.

weights.transform = "R"
# Row-standardise so each ED's neighbour weights sum to one.
# An ED with four neighbours gives each neighbour a weight of 0.25.
# The spatial lag is then the average value of those neighbours.


neighbour_counts = [
    len(neighbours)
    for neighbours in weights.neighbors.values()
]

# weights.neighbors stores the neighbour IDs for each ED.
# Count the IDs in each list to obtain one neighbour count per ED.

print("Number of EDs in weights:", weights.n)
print("Number of islands:", len(weights.islands))
print("Minimum neighbours:", min(neighbour_counts))
print("Maximum neighbours:", max(neighbour_counts))
print(f"Mean neighbours: {np.mean(neighbour_counts):.2f}")
# An island is an ED with no neighbours under the chosen definition.
# Summarise the neighbour counts before running the statistical tests.


###################################################
############## 3. GLOBAL MORAN'S I ################
###################################################
np.random.seed(SEED)
# Initialise NumPy's random sequence once before running both tests.
# Repeating the same code with the same inputs reproduces the shuffles.

raw_moran = Moran(
    density.to_numpy(),
    weights,
    permutations=PERMUTATIONS
)
# Convert the density Series to a NumPy array.
# Calculate Global Moran's I and its permutation test using Queen weights.

sqrt_moran = Moran(
    sqrt_density.to_numpy(),
    weights,
    permutations=PERMUTATIONS
)
# Repeat the test for the square-root variable using the same ED order
# and spatial weights. The random sequence continues from the first test.

print(f"Raw density Moran's I: {raw_moran.I:.4f}")
print(f"Raw density permutation p-value: {raw_moran.p_sim:.4f}")

print(f"Square-root density Moran's I: {sqrt_moran.I:.4f}")
print(
    "Square-root density permutation p-value: "
    f"{sqrt_moran.p_sim:.4f}"
)
# I stores the observed statistic; p_sim stores the permutation p-value.
# Display four decimal places to retain the reported p-value of 0.0001.

print("Permutations per test:", PERMUTATIONS)
print("Initial random seed:", SEED)
print("Distribution figure saved to:", FIGURE_FILE)




#########################################
##### 4. INTERPRETATION #####
#########################################

# The written interpretation is in item 4 of the opening ANSWERS docstring.
# Update those answers if the input data or printed results change.
