import sys
import numpy
import pandas
import pyvo
from IPython.core.display import display

from tabulate import tabulate
import warnings
import matplotlib

import seaborn
import matplotlib.pyplot as plt
from astropy import units
from astropy import constants
from astropy.visualization import quantity_support
from matplotlib_inline.backend_inline import set_matplotlib_formats
from matplotlib import cm
import matplotlib.colors
from matplotlib.colors import ListedColormap
from collections import OrderedDict
from matplotlib.lines import Line2D
import itertools
from functools import partial
import matplotlib.ticker as mticker
from cycler import cycler
from astropy import units
from astropy import constants

if len(sys.argv) < 2:
    raise SystemExit("No star name provided")

if len(sys.argv) > 2:
    raise SystemExit("Too many arguments")

queryStar = sys.argv[1]
print(f"Processing star name: {queryStar}")
queryStarNoSpace = queryStar.strip().replace(" ", "-")

quantity_support()
set_matplotlib_formats('svg')

plt.rc('legend', frameon=False)
plt.rc('figure', figsize=(10, 10))
plt.rc('font', size=12)

mass_coeff = constants.M_jup / constants.M_earth
radius_coeff = constants.R_jup / constants.R_earth

serviceExoplanets = pyvo.dal.TAPService("http://voparis-tap-planeto.obspm.fr/tap")
serviceNASA = pyvo.dal.TAPService("https://exoplanetarchive.ipac.caltech.edu/TAP")

print("[1/3] Querying Exoplanets...")

# For individual planetary systems: data from Exoplanet.eu and cross-cheking with NASA
tableNameExoplanets = "exoplanet.epn_core"

fieldsExoplanets = [
    "star_name",
    "star_spec_type",
    "star_teff",
    "star_radius",
    "star_mass",
    "mass",
    "radius",
    "granule_uid"
]

queryExoplanets = (
    " ".join((
        "SELECT",
        ", ".join(fieldsExoplanets),
        f"FROM {tableNameExoplanets}",
        f"WHERE star_name = '{queryStar}'"
    ))
)

resultsExoplanets = serviceExoplanets.search(queryExoplanets)
planetsFound = len(resultsExoplanets)
print(f"Total planets found in {queryStar}: {planetsFound}")
if planetsFound == 0:
    print("No planets found!")
    raise SystemExit(0)

workingTableExoplanets = resultsExoplanets.to_table().to_pandas(index="granule_uid")

# NASA

tableNameNASA = "ps"

fieldsNASA = [
    "hostname",
    "pl_name",
    "pl_pubdate"
]

def lookForParameterInNASA(systemName, planetName, param):
    print(f"Got this parameter to check in NASA: {param}")
    if param in dictForCheckingstar:
        # get stellar parameter
        queryNASA = (
            " ".join((
                f"SELECT {param}",
                f"FROM {tableNameNASA}",
                f"WHERE hostname = '{systemName}' AND {param} IS NOT NULL",
                "ORDER BY pl_pubdate DESC"
            ))
        )
        resultsNASA = serviceNASA.search(queryNASA)
        if len(resultsNASA) != 0:
            return resultsNASA[0].get(param)
        else:
            return None
    else:
        # get planet parameter
        queryNASA = (
            " ".join((
                f"SELECT {param}",
                f"FROM {tableNameNASA}",
                f"WHERE pl_name = '{planetName}' AND {param} IS NOT NULL",
                "ORDER BY pl_pubdate DESC"
        ))
    )
    resultsNASA = serviceNASA.search(queryNASA)
    if len(resultsNASA) != 0:
        if planetName == "Kepler-102 b":
            print(f"All results from NASA for this parameter: {resultsNASA}")
        return resultsNASA[0].get(param)
    else:
        return None


# checking
print("[2/3] Cross-checking with NASA, adding...")

dictForCheckingplanet = {
    "mass": "pl_massj",
    "radius": "pl_radj"
}

dictForCheckingstar = {
    "star_teff": "st_teff",
    "star_radius": "st_rad",
    "star_mass": "st_mass"
}



dictForChecking = {}
dictForChecking.update(dictForCheckingplanet.items())
dictForChecking.update(dictForCheckingstar.items())

listForAdding = {
    "ra": 0,
    "sy_dist": 0
}
for valueToAdd, valueType in listForAdding.items():
    # 1 - str
    if valueType == 1:
        workingTableExoplanets[valueToAdd] = numpy.array(numpy.NaN, dtype=str)
    # 0 - float
    else:
        workingTableExoplanets[valueToAdd] = numpy.array(numpy.NaN, dtype=float)

countEnrichedPlanets = 0
for index, row in workingTableExoplanets.iterrows():
    rowAlreadyCounted = False
    systemName = row["star_name"]
    print(f"Iterating {index}, star name: {systemName}")

    for valueToCheck, valueToFind in dictForCheckingstar.items():
        if pandas.isnull(row[valueToCheck]):
            print(f"Value: {row[valueToCheck]}")
            valueFromNASA = lookForParameterInNASA(systemName, index.replace("'", ""), valueToFind)
            if valueFromNASA:
                print(f"Will save this value from NASA: {valueFromNASA}")
                workingTableExoplanets.at[index, valueToCheck] = valueFromNASA
                if not rowAlreadyCounted:
                    rowAlreadyCounted = True
                    countEnrichedPlanets += 1
    # automated adding
    for valueToAdd in listForAdding:
        valueFromNASA = lookForParameterInNASA(systemName, index.replace("'", ""), valueToAdd)
        if valueFromNASA:
            workingTableExoplanets.at[index, valueToAdd] = valueFromNASA

print(f"Enriched(1) planets with data from NASA: {countEnrichedPlanets}")

print("[3/3] Cross-checking with NASA, checking...")

def findAndFixErrors(starName, planetName, paramName, paramItself, paramErrorMin, paramErrorMax):
    valItself = None
    valErrorMin = None
    valErrorMax = None
    if (
        pandas.isnull(paramItself)
        or
        (
            (pandas.isnull(paramErrorMin) or paramErrorMin == 0.) or (pandas.isnull(paramErrorMax) or paramErrorMax == 0.)
        )
    ):
        valItself = lookForParameterInNASA(
            starName,
            planetName,
            dictForCheckingplanet[paramName]
        )
        valErrorMin = lookForParameterInNASA(
            starName,
            planetName,
            f"{dictForCheckingplanet[paramName]}err2"
        )
        valErrorMax = lookForParameterInNASA(
            starName,
            planetName,
            f"{dictForCheckingplanet[paramName]}err1"
        )
    return valItself, valErrorMin, valErrorMax

workingTableExoplanets.to_pickle(f"./merg/{queryStar}.pkl")
