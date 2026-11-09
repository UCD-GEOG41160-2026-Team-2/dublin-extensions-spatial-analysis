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
The script calculates Gi* for each ED using binary Queen contiguity weights and includes the focal ED (star=1.0).
It classifies hotspots where Z > 1.96 and coldspots where Z < -1.96,
calculates two-sided p-values from the Gi* Z-scores, and applies FDR correction at alpha = 0.05.
The script prints the hotspot/coldspot counts, FDR threshold, and the EDs that remain significant after correction.

3.5 Gi* Cluster Map
Hotspot and coldspot map: The script maps the Gi* Z-scores using the diverging RdBu_r colour scheme.
EDs that remain significant after FDR are outlined in black.
The map is saved as: outputs/GiStar_ExtensionDensity_DublinED.png

3.6 Comparison of LISA and Gi* Cluster Maps
Both methods identify high extension application density clusters in south-central Dublin.
However, their results differ because LISA identifies local spatial clusters and outliers,
while Gi* identifies concentrations of high or low values across neighbouring areas.
Despite these differences, both methods indicate spatial clustering in extension application density.
Their combined results provide a better understanding of its spatial distribution across Dublin.

"""

from pathlib import Path
import esda 
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from libpysal.weights import Queen
from esda.moran import Moran_Local
from esda import fdr
from esda.getisord import G_Local
from libpysal.weights import Queen
from scipy.stats import norm 


ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = ROOT / "outputs"

DENSITY_FILE = OUTPUTS / "DublinED_ExtensionDensity.gpkg"

SEED = 42
PERMUTATIONS = 9999
DENSITY_COLUMN = "app_density_km2"
FDR_ALPHA = 0.05
Z_CUTOFF = 1.96

if not DENSITY_FILE.exists():
    raise FileNotFoundError(
        f"{DENSITY_FILE} not found. Run the earlier assignment script "
        "that creates DublinED_ExtensionDensity.gpkg first."
    )

eds = gpd.read_file(DENSITY_FILE)

if DENSITY_COLUMN not in eds.columns:
    raise KeyError(
        f"'{DENSITY_COLUMN}' not found. Available columns: {list(eds.columns)}"
    )

# Ensure density values are numeric and remove missing/infinite observations.
eds[DENSITY_COLUMN] = pd.to_numeric(eds[DENSITY_COLUMN], errors="coerce")
eds = eds.loc[np.isfinite(eds[DENSITY_COLUMN])].copy()

if eds.empty:
    raise ValueError("No valid extension-density values were found.")

if (eds[DENSITY_COLUMN] < 0).any():
    raise ValueError("Extension density contains negative values.")

print("Input layer:", DENSITY_FILE)
print("EDs included:", len(eds))
print("CRS:", eds.crs)
print("Density variable:", DENSITY_COLUMN)

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

eds.loc[fdr_significant, "LISA_Class"] = quadrant_labels.loc[original_significant]

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


# Construct Queen contiguity weights. Binary weights are used for local Gi*.
## The star=1.0 setting in Section 2 includes the focal ED itself.
w_bin = Queen.from_dataframe(eds, use_index=False)
w_bin.transform = "B"

print("\n1. Spatial weights")
print("Number of EDs in weights:", w_bin.n)
print("Number of islands:", len(w_bin.islands))

if w_bin.islands:
    raise ValueError(
        "One or more EDs have no Queen neighbours. Review the geometries "
        "before running Gi*."
    )


# Calculate local Getis-Ord Gi* Z-scores using binary Queen weights.
## The p-values below are two-sided normal-approximation p-values.
np.random.seed(SEED)

gi = G_Local(
    eds[DENSITY_COLUMN].to_numpy(),
    w_bin,
    transform="B",
    star=1.0,
    permutations=PERMUTATIONS,
    seed=SEED
)

eds["gi_z"] = gi.Zs

# Calculate two-sided p-values BEFORE the FDR correction.
eds["gi_p"] = 2 * norm.sf(np.abs(eds["gi_z"]))

eds["gi_class"] = np.select(
    [
        eds["gi_z"] > Z_CUTOFF,
        eds["gi_z"] < -Z_CUTOFF
    ],
    [
        "Hotspot",
        "Coldspot"
    ],
    default="Not significant"
)

hotspots = eds["gi_class"].eq("Hotspot")
coldspots = eds["gi_class"].eq("Coldspot")

print("\n2. Gi* results before FDR correction")
print(f"Z-score cut-off: +/-{Z_CUTOFF}")
print("Hotspots (Z > 1.96):", int(hotspots.sum()))
print("Coldspots (Z < -1.96):", int(coldspots.sum()))
print("Other EDs:", int((~(hotspots | coldspots)).sum()))


# Apply the Benjamini-Hochberg false discovery rate procedure to the two-sided p-values, controlling the FDR at 5%.

fdr_threshold = esda.fdr(eds["gi_p"].to_numpy(), FDR_ALPHA)
eds["fdr_significant"] = eds["gi_p"] <= fdr_threshold

# Classify the direction of EDs that pass FDR and the Z-score cut-off.
eds["gi_fdr_class"] = np.select(
    [
        eds["fdr_significant"] & (eds["gi_z"] > Z_CUTOFF),
        eds["fdr_significant"] & (eds["gi_z"] < -Z_CUTOFF)
    ],
    [
        "Hotspot (FDR)",
        "Coldspot (FDR)"
    ],
    default="Not significant after FDR"
)

fdr_hotspots = eds["gi_fdr_class"].eq("Hotspot (FDR)")
fdr_coldspots = eds["gi_fdr_class"].eq("Coldspot (FDR)")

print("\n3. False discovery rate correction")
print(f"FDR alpha: {FDR_ALPHA:.2f}")
print(f"FDR p-value threshold: {fdr_threshold:.8f}")
print("EDs passing FDR:", int(eds["fdr_significant"].sum()))
print("FDR-significant hotspots:", int(fdr_hotspots.sum()))
print("FDR-significant coldspots:", int(fdr_coldspots.sum()))

# Print ED names and statistics for all EDs passing FDR.
name_column = "ED_ENGLISH" if "ED_ENGLISH" in eds.columns else None
result_columns = ([name_column] if name_column else []) + [
    "gi_z", "gi_p", "gi_fdr_class"
]

print("\nEDs passing the FDR correction:")
if eds["fdr_significant"].any():
    print(
        eds.loc[eds["fdr_significant"], result_columns]
        .sort_values("gi_z")
        .to_string(index=False)
    )
else:
    print("No EDs passed the FDR correction.")


# Map the Gi* Z-scores using a diverging colour scheme.
## EDs passing FDR are outlined in black.

limit = float(np.nanmax(np.abs(eds["gi_z"])))
if not np.isfinite(limit) or limit == 0:
    limit = Z_CUTOFF

fig, ax = plt.subplots(figsize=(8, 7))

eds.plot(
    column="gi_z",
    cmap="RdBu_r",          # blue = cold; red = hot
    vmin=-limit,
    vmax=limit,
    legend=True,
    edgecolor="grey",
    linewidth=0.35,
    ax=ax,
    legend_kwds={"label": "Getis-Ord Gi* Z-score"}
)

fdr_eds = eds.loc[eds["fdr_significant"]]
if not fdr_eds.empty:
    fdr_eds.boundary.plot(
        ax=ax,
        color="black",
        linewidth=1.5
    )

ax.set_title(
    "Getis-Ord Gi* – Extension Application Density, Dublin EDs\n"
    "Red = hotspots | Blue = coldspots | "
    "Black outline = significant after FDR",
    fontsize=11
)
ax.set_axis_off()

fig.tight_layout()
MAP_FILE = OUTPUTS / "GiStar_ExtensionDensity_DublinED.png"
fig.savefig(MAP_FILE, dpi=300, bbox_inches="tight")
plt.show()

print("\nMap saved to:", MAP_FILE)


# Print a concise, assignment-ready interpretation using the results calculated above.
## The geographic pattern should also be described with reference to the saved map, rather than inferred from statistics alone.
total_eds = len(eds)
n_hotspots = int(hotspots.sum())
n_coldspots = int(coldspots.sum())
n_fdr = int(eds["fdr_significant"].sum())
n_fdr_hot = int(fdr_hotspots.sum())
n_fdr_cold = int(fdr_coldspots.sum())

print("\n" + "=" * 65)
print("QUESTION 4: GETIS-ORD GI* RESULTS")
print("=" * 65)

print(
    f"Before FDR correction, {n_hotspots} EDs were classified as "
    f"hotspots and {n_coldspots} as coldspots out of {total_eds} EDs."
)
print(
    f"After FDR correction at alpha = {FDR_ALPHA:.2f} "
    f"(p-value threshold = {fdr_threshold:.8f}), "
    f"{n_fdr} EDs remained statistically significant: "
    f"{n_fdr_hot} hotspots and {n_fdr_cold} coldspots."
)

if n_fdr_hot > 0:
    print(
        "This indicates statistically significant clustering of high "
        "extension application density in the EDs listed above."
    )
else:
    print("No statistically significant hotspots remained after FDR correction.")

if n_fdr_cold > 0:
    print(
        "The significant coldspots indicate statistically significant "
        "clustering of low extension application density."
    )
else:
    print(
        "No coldspots remained statistically significant after FDR correction. "
        "EDs with negative Gi* Z-scores should therefore not be described as "
        "FDR-significant coldspots."
    )

if n_fdr > 0:
    top_hotspots = (
        eds.loc[fdr_hotspots]
        .sort_values("gi_z", ascending=False)
        .head(5)
    )
    print("\nFive strongest FDR-significant hotspots (highest Gi* Z-scores):")
    for _, row in top_hotspots.iterrows():
        ed_name = row[name_column] if name_column else f"ED index {row.name}"
        print(f"  - {ed_name}: Gi* Z = {row['gi_z']:.3f}, p = {row['gi_p']:.6f}")

print(
    "\nMap interpretation guidance: use the saved map to identify where "
    "the black-outlined significant EDs are geographically concentrated. "
    "Red areas have positive Gi* Z-scores and blue areas have negative "
    "Gi* Z-scores; black outlines identify EDs passing the FDR correction."
)

# Export the FDR-significant EDs as a separate CSV table for the report.
SIGNIFICANT_CSV = OUTPUTS / "GiStar_FDR_Significant_EDs.csv"
if n_fdr > 0:
    eds.loc[eds["fdr_significant"], result_columns].sort_values(
        "gi_z", ascending=False
    ).to_csv(SIGNIFICANT_CSV, index=False)
    print("\nFDR-significant ED table saved to:", SIGNIFICANT_CSV)
else:
    print("\nNo significant ED table was exported because no EDs passed FDR.")
