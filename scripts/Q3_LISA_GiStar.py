"""
GEOG41160 Assignemnt 2
Questions 3: Local Spatial Autocorrelation

3.1 Local Morans I 
Local Morans I (LISA) was calculated using 
square-root-transformed extension density, 9,999 permutations 
and a fixed seed to identify significant High-High, Low-High, Low-Low, and High-Low clusters 
across Dublins EDs (p ≤ 0.05).

3.2 False Discovery Rate (FDR) 
A false discovery rate (FDR) correction was applied to the LISA p-values, 
resulting in an adjusted significance threshold of 0.007716 with 25 EDs remaining statistically significant.

3.3 LISA Cluster Map
A colour-coded map of the LISA cluster classifications was created,
including a legend, and saved as LISA_ExtensionDensity_DublinED.png.


3.4 Getis-Ord Gi*
x
3.5 Gi* Cluster Map
x
3.6 Comparison of LISA and Gi* Cluster Maps
x

"""

from pathlib import Path
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from libpysal.weights import Queen
from esda.moran import Moran_Local
from esda import fdr


ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = ROOT / "outputs"

DENSITY_FILE = OUTPUTS / "DublinED_ExtensionDensity.gpkg"

SEED = 42
PERMUTATIONS = 9999
DENSITY_COLUMN = "app_density_km2"


eds = gpd.read_file(DENSITY_FILE)
eds[DENSITY_COLUMN] = pd.to_numeric(eds[DENSITY_COLUMN], errors="coerce")
valid_rows = np.isfinite(eds[DENSITY_COLUMN])
eds = eds.loc[valid_rows].copy()

density = eds[DENSITY_COLUMN]
sqrt_density = np.sqrt(density)

weights = Queen.from_dataframe(eds, use_index=False)
weights.transform = "R"

# LISA cluster classification

lisa = Moran_Local(sqrt_density, weights, permutations=9999, seed=42)
quadrants = {
    1: "High-High",
    2: "Low-High",
    3: "Low-Low",
    4: "High-Low"
}

# False discovery rate correction
original_significant = lisa.p_sim <= 0.05
fdr_cutoff = fdr(lisa.p_sim, alpha=0.05)
fdr_significant = lisa.p_sim <= fdr_cutoff

# Classify EDs
eds["LISA_Class"] = "Not Significant"

quadrant_labels = pd.Series(
    lisa.q,
    index=eds.index
).map(quadrants)

eds.loc[original_significant, "LISA_Class"] = quadrant_labels.loc[original_significant]

print(f"FDR-adjusted p-value cutoff: {fdr_cutoff:.6f}")
print(f"Significant EDs after FDR correction: {fdr_significant.sum()}")

# Map LISA classes
colours = {
    "High-High": "red",
    "Low-High": "lightblue",
    "Low-Low": "blue",
    "High-Low": "pink",
    "Not Significant": "lightgrey"
}

fig, ax = plt.subplots(figsize=(8, 8))

for label, colour in colours.items():
    eds[eds["LISA_Class"] == label].plot(
        ax=ax,
        color=colour,
        edgecolor="black",
        linewidth=0.3,
        label=label
    )
ax.set_title("LISA Clusters of Extension Density, Dublin EDs")
ax.set_axis_off()
ax.legend(title="LISA class", loc="lower right")

plt.savefig(
    OUTPUTS / "LISA_ExtensionDensity_DublinED.png",
    dpi=300,
    bbox_inches="tight"
)
plt.show()