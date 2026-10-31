
"""
GEOG41160 Assignment 2, Question 2: Global spatial autocorrelation
Team: 2
ANSWERS

1. Distribution: Raw extension density was positively skewed
   (skewness = 1.39). We selected a square-root transformation,
   which reduced skewness to 0.35 and gave a more balanced
   distribution than log1p (skewness = -1.14).
"""

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = ROOT / "outputs"
DENSITY_FILE = OUTPUTS / "DublinED_ExtensionDensity.gpkg"

eds = gpd.read_file(DENSITY_FILE)

density = pd.to_numeric(
    eds["app_density_km2"],
    errors="coerce"
).dropna()

sqrt_density = np.sqrt(density)
log_density = np.log1p(density)

print("EDs with valid density:", len(density))
print(f"Raw density skewness: {density.skew():.2f}")
print(f"Square-root density skewness: {sqrt_density.skew():.2f}")
print(f"Log-transformed density skewness: {log_density.skew():.2f}")

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

axes[0].hist(density, bins=20, edgecolor="black")
axes[0].set_title("Raw extension density")
axes[0].set_xlabel("Extension applications per km²")
axes[0].set_ylabel("Number of EDs")

axes[1].hist(sqrt_density, bins=20, edgecolor="black")
axes[1].set_title("Square-root-transformed density")
axes[1].set_xlabel("Square root of applications per km²")
axes[1].set_ylabel("Number of EDs")

plt.tight_layout()
fig.savefig(
    OUTPUTS / "ExtensionDensity_DublinED_histogram.png",
    dpi=200,
    bbox_inches="tight"
)
plt.show()
