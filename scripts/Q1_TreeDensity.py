"""
GEOG41160 Assignment 2, Question 1: Dublin Tree and Planning Density
Team: 2

ANSWERS
1. Dataset import and filtering
    - Dublin City EDs count: 162 features (imported from national boundaries, filtered for COUNTY_ENGLISH = 'Dublin City').
    - Filtered planning applications count: 9,914 rows (filtered for Dublin City Council, years 2019-2025, and development description containing 'extension').


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


