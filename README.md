# Dublin Extensions: Spatial Analysis

## Project Overview

This project analyses planning applications for building extensions in Dublin City by Electoral Division (ED). It calculates extension-application density, investigates global spatial autocorrelation, identifies local spatial clusters using Local Moran's I (LISA), and maps hotspots and coldspots using Getis-Ord Gi*.

## Repository Structure

- `data/planning/` — National planning-application GeoPackage.
- `data/boundaries/` — National Electoral Division boundary GeoPackage.
- `scripts/` — Python scripts for Questions 1, 2 and 3.
- `outputs/` — Generated database, spatial datasets, figures and maps.

## Data Sources and Licences

1. **Planning Application Points:** Department of Housing, Local Government and Heritage. Source: https://data-housinggovie.opendata.arcgis.com/datasets/housinggovie::planning-application-points  
   Licence: CC BY 4.0. Download date: [INSERT ACTUAL DOWNLOAD DATE].

2. **Electoral Divisions:** Tailte Éireann, CSO Electoral Divisions — National Statistical Boundaries — 2022 — Ungeneralised. Source: https://data-osi.opendata.arcgis.com/datasets/osi::cso-electoral-divisions-national-statistical-boundaries-2022-ungeneralised  
   Licence: CC BY 4.0.

## Coordinate Reference System

The analysis uses Irish Transverse Mercator (EPSG:2157) for spatial calculations, including area and density calculations. Planning applications will be reprojected to the boundary layer's coordinate reference system before spatial analysis.

## Execution Order

Run the scripts from the repository root in this order:

1. `python scripts/Q1_ExtensionDensity.py`
2. `python scripts/Q2_GlobalMoran.py`
3. `python scripts/Q3_LISA_GiStar.py`

The input GeoPackages must be saved in their specified locations before running the scripts. The outputs from Question 1 are used in the subsequent spatial-statistics analyses.


## Generated Outputs
- outputs/extensions.duckdb — DuckDB database created in Question 1.
- outputs/DublinED_ExtensionDensity.gpkg — Extension density by Electoral Division.
- outputs/ExtensionDensity_DublinED.png — Extension density map.
- outputs/LISA_ExtensionDensity_DublinED.png — Local Moran’s I cluster map.
- outputs/GiStar_ExtensionDensity_DublinED.png — Getis-Ord Gi* hotspot and coldspot map.


## Sprint Log

| Sprint | Date | Work completed | Next steps |
|---|---|---|---|
| Sprint 0 | 23/10/2026 | Repository setup and initial README | Prepare and verify the input data |
| Sprint 1 | 29/10/2026 | Completed Question 1: Extension Density Calculation, DuckDB filtering, spatial join, and saved density choropleth map | Begin Question 2 with density gpkg file produced from Q1|
| Sprint 2 | 04/11/2026 | Completed Question 2: Checked data distribution and compared both transformations with the raw distribution by plotting them side by side, applied square-root transformation, built spatial weights and computed Global Moran's I | Begin Question 3 (Local spatial autocorrelation and hotspot analysis with/without FDR)
| Sprint 3 | 09/11/2026 | Completed Question 3: LISA, Gi* and maps | Final review and quality checks |
