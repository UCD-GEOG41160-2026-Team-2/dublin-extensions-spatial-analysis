"""
GEOG41160 Assignment 2, Question 1: Dublin Tree and Planning Density
Team: 2

ANSWERS
1. Dataset import and filtering
    - Dublin City EDs count: 162 features (imported from national boundaries, filtered for COUNTY_ENGLISH = 'Dublin City').
    - Filtered planning applications count: 9,914 rows (filtered for Dublin City Council, years 2019-2025, and development description containing 'extension').

2. Inspecting both layers' CRS, geometry type and feature counts and reprojecting extension apps.
    -Used GeoPandas to verify CRS projection, will explain in details why did not do it with SQL query on individual report.
    -Printing geometry types and feature count via DuckDB
    -Reprojecting extension applications from EPSG:3857 to match EDs' EPSG:2157
        # Why reproject to Irish Transverse Mercator (EPSG:2157)?
        # The raw data showed up in Web Mercator (EPSG:3857) — which is great
        # for web maps like Google Maps, but messes up local measurements.
        # Web Mercator distorts areas and distances in Ireland, so we need
        # EPSG:2157 (the official local system) to keep our distances, areas,
        # and spatial joins for Dublin totally accurate.

3. Spatial join and Boundary Check:
    -Successfully performed a spatial join between EDs and reprojected extension points, printing out the aggregated table.
    -Top EDs by application volume: Clontarf East B (304 applications), followed by Rathmines East C (199), and Pembroke West A (188) in 3rd place.
    -Applications falling outside the Dublin's ED boundaries: 14. These points fell outside because:
        a. Some fall just across the border into neighboring council areas (like Fingal or South Dublin).
        b. Some are near water or coastlines where the land shape boundary stops.
        c. Some address coordinates could be slightly off and landed just outside the official map border.
"""

from pathlib import Path

import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import duckdb
#########################################
##### SETUP #####
#########################################

ROOT = Path(__file__).resolve().parents[1]  # the repository folder
DATA = ROOT / "data"
OUTPUTS = ROOT / "outputs"
OUTPUTS.mkdir(exist_ok=True)

SEED = 42  # the same seed in every script
TREE_DENSITY = OUTPUTS / "DublinED_TreeDensity.gpkg"

# defining data file paths
EDS = DATA / "boundaries" / "eds_national.gpkg"
PLANNING = DATA / "planning" / "planning_applications_national.gpkg"

# connecting to DuckDB and installing and loading spatial extensions to read from gpkg
DATABASE = OUTPUTS / "extensions.duckdb"
con = duckdb.connect(str(DATABASE))
print("Connected to:", DATABASE)
con.install_extension("spatial")
con.load_extension("spatial")

#########################################
##### 1. Import and Filter tables on database #####
#########################################

## a. Importing Dublin City EDs and reporting row count
con.execute(f"""
    CREATE OR REPLACE TABLE eds AS
    SELECT * FROM ST_Read('{EDS}')
    WHERE COUNTY_ENGLISH = 'DUBLIN CITY';
""")
ed_count = con.execute(
    "SELECT COUNT (*) FROM eds;"
).fetchone()[
    0
]  # returns the shape as numpy array and [0] returns the row number, which is our ed count
print(f"Dublin City EDs: {ed_count}")

## b. Importing and filtering planning applications for Dublin in years 2019-2025, containing "extension"
con.execute(f"""
    CREATE OR REPLACE TABLE extensions AS
    SELECT * FROM ST_Read('{PLANNING}')
    WHERE PlanningAuthority = 'Dublin City Council' 
        AND EXTRACT(YEAR FROM ReceivedDate) BETWEEN 2019 AND 2025
        AND LOWER(DevelopmentDescription) LIKE '%extension%' -- LOWER() and LIKE ensure we catch variations like "Extension" or "extensions"
""")
ext_count = con.execute("SELECT COUNT (*) FROM extensions;").fetchone()[0]
print(f"Filtered extensions rows: {ext_count}")

#########################################
##### 2. CRS AND REPROJECTION #####
#########################################

## a. Reading a small portion to verify CRS
eds_gdf = gpd.read_file(EDS, rows=5)
print("ED Boundaries CRS:", eds_gdf.crs)

ext_gdf = gpd.read_file(PLANNING, rows=5)
print("Extension Planning applications CRS:", ext_gdf.crs)

## b. Inspecting geometry types and feature counts via DuckDB
eds_info = con.execute("""
    SELECT 
                ST_GeometryType(SHAPE) AS geom_type,
                COUNT (*) AS total_features
    FROM eds
    GROUP BY geom_type
""").fetchone()

print("EDs Info (Geometry Type, Count):", eds_info)

ext_info = con.execute("""
    SELECT 
                ST_GeometryType(SHAPE) AS geom_type,
                COUNT (*) AS total_features
    FROM extensions
    GROUP BY geom_type
""").fetchone()

print("Extension applications' info (Geometry Type, Count):", ext_info)

## c. Reprojecting extension applications
con.execute("""
    CREATE or REPLACE TABLE extensions_projected AS
    SELECT * EXCLUDE (SHAPE),
            ST_Transform(SHAPE, 'EPSG:3857', 'EPSG:2157', always_xy := true) AS SHAPE
            FROM extensions 
""")
print("Successfully reprojected extensions layer to EPSG:2157!")


#########################################
##### 3. SPATIAL JOIN AND AGGREGATION #####
#########################################

## a. Joining eds in dublin with reprojected extension applications by their spatial attributes
con.sql(""" 
    CREATE OR REPLACE TABLE eds_applications AS
    SELECT eds.ED_ENGLISH, eds.SHAPE,
        COUNT(ex_p.ApplicationNumber) AS extension_applications
    FROM eds AS eds
    LEFT JOIN extensions_projected AS ex_p ON ST_Within(ex_p.SHAPE, eds.SHAPE)
    GROUP BY eds.ED_ENGLISH, eds.SHAPE
""")
print(
    con.sql(
        "SELECT ED_ENGLISH, extension_applications FROM eds_applications ORDER BY extension_applications DESC LIMIT 5"
    )
)

## b. Finding out how many applications fall outside of dublin with WHERE NOT EXISTS query filter

outside_count = con.execute("""
    SELECT COUNT (*) AS outside_eds
    FROM extensions_projected as ex_p
    WHERE NOT EXISTS (       -- Spatial anti-join filter that only keeps the planning application if the subquery inside the parentheses returns no match (0)
            SELECT 1 FROM EDS    -- Standart SQL optimization, returns dummy value 1 when it finds a geometric match instead of column name
            WHERE ST_Within(ex_p.SHAPE, eds.SHAPE)            )    
""").fetchone()[0]

print("Applications outside Dublin City ED boundaries:", outside_count)