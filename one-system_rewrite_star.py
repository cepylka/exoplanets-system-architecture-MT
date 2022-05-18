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

curves = pandas.read_fwf("../coaz/python/exoplanets/curves_all.txt", index_col=0)#, header=None)

# Exoplanets

print("[1/6] Querying Exoplanets...")

# For individual planetary systems: data from Exoplanet.eu and cross-cheking with NASA
tableNameExoplanets = "exoplanet.epn_core"

fieldsExoplanets = [
    "star_name",
    "star_spec_type",
    "star_teff",
    "star_radius",
    "star_mass",
    "star_metallicity",
    "star_age",
    "mass",
    "mass_error_min",
    "mass_error_max",
    "radius",
    "radius_error_min",
    "radius_error_max",
    "semi_major_axis",
    "semi_major_axis_error_min",
    "semi_major_axis_error_max",
    "period",
    "period_error_min",
    "period_error_max",
    "albedo",
    "albedo_error_min",
    "albedo_error_max",
    "mass_detection_type",
    "radius_detection_type",
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
                # f"SELECT TOP 1 {param}",
                f"SELECT {param}",
                f"FROM {tableNameNASA}",
                f"WHERE hostname = '{systemName}' AND {param} IS NOT NULL",
                "ORDER BY pl_pubdate DESC"
            ))
        )
        resultsNASA = serviceNASA.search(queryNASA)
        if len(resultsNASA) != 0:
            # print(f"All results from NASA for this parameter: {resultsNASA}")
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

print("[2/6] Cross-checking with NASA, adding...")

dictForCheckingplanet = {
    "mass": "pl_massj",
    "radius": "pl_radj",
    "semi_major_axis": "pl_orbsmax",
    "period": "pl_orbper"
}

dictForCheckingstar = {
    "star_teff": "st_teff",
    "star_radius": "st_rad",
    "star_mass": "st_mass",
    "star_metallicity": "st_met",
    "star_age": "st_age",
    "st_lum": "st_lum",
    "sy_snum": "sy_snum",
    "st_rotp": "st_rotp",
    "st_metratio": "st_metratio",
    "star_spec_type": "st_spectype"
}

dictForChecking = {}
dictForChecking.update(dictForCheckingplanet.items())
dictForChecking.update(dictForCheckingstar.items())

listForAdding = {
    "st_lum": 0,
    "sy_snum": 0,
    "st_rotp": 0,
    "st_metratio": 1,
    "sy_pnum": 0,
    "cb_flag": 0,
    "pl_massjlim": 0,
    "pl_radjlim": 0,
    "pl_orbsmaxlim": 0,
    "pl_orbperlim": 0,
}
for valueToAdd, valueType in listForAdding.items():
    # 1 - str
    if valueType == 1:
        workingTableExoplanets[valueToAdd] = numpy.array(numpy.NaN, dtype=str)
    # 0 - float
    else:
        workingTableExoplanets[valueToAdd] = numpy.array(numpy.NaN, dtype=float)

# fill [M/H] for those with values in metallicity
for index, row in workingTableExoplanets.iterrows():
    if not pandas.isnull(row["star_metallicity"]):
        workingTableExoplanets.at[index, "st_metratio"] = "[M/H]"

countEnrichedPlanets = 0
for index, row in workingTableExoplanets.iterrows():
    rowAlreadyCounted = False
    systemName = row["star_name"]
    print(f"Iterating {index}, star name: {systemName}")
    # automated checking
    for valueToCheck, valueToFind in dictForCheckingstar.items():
        if pandas.isnull(row[valueToCheck]):
            print(f"Value: {row[valueToCheck]}")
            # print(f"Is null: {pandas.isnull(row[valueToCheck])}")
            # print(f"Is NaN: {pandas.isna(row[valueToCheck])}")
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

print("[3/6] Cross-checking with NASA, checking...")


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

print("[4/6] Checking for NaNs...")

countEnrichedPlanets = 0
for index, row in workingTableExoplanets.iterrows():
    rowAlreadyCounted = False
    for paramNameExoplanets, paramNameNASA in dictForCheckingplanet.items():
        print(f"Checking parameter: {paramNameExoplanets}")
        valItself, valErrorMin, valErrorMax = findAndFixErrors(
            row["star_name"],
            index,
            paramNameExoplanets,
            row[paramNameExoplanets],
            row[f"{paramNameExoplanets}_error_min"],
            row[f"{paramNameExoplanets}_error_max"]
        )
        if (
            valItself is not None
        ):
            print(f"Got values from NASA: {valItself}, {valErrorMin}, {valErrorMax}")
            workingTableExoplanets.at[index, paramNameExoplanets] = valItself
            if valErrorMin is not None:
                workingTableExoplanets.at[index, f"{paramNameExoplanets}_error_min"] = abs(valErrorMin)
            if valErrorMax is not None:
                workingTableExoplanets.at[index, f"{paramNameExoplanets}_error_max"] = valErrorMax
            countEnrichedPlanets += 1

print(f"Enriched(2) planets with data from NASA: {countEnrichedPlanets}")


print("[5/6] Habitable zone calculation...")

# first we create 6 rows for HZ coefficients
workingTableExoplanets["st_recentVenus"] = numpy.array(numpy.NaN, dtype=float)
workingTableExoplanets["st_runawayGreenhouse"] = numpy.array(numpy.NaN, dtype=float)
workingTableExoplanets["st_maxGreenhouse"] = numpy.array(numpy.NaN, dtype=float)
workingTableExoplanets["st_earlyMars"] = numpy.array(numpy.NaN, dtype=float)
workingTableExoplanets["st_half_Earth"] = numpy.array(numpy.NaN, dtype=float)
workingTableExoplanets["st_five_Earth"] = numpy.array(numpy.NaN, dtype=float)

# Coeffcients to be used in the analytical expression to calculate habitable zone flux
# boundaries

seffsun  = [1.776, 1.107, 0.356, 0.320, 1.188, 0.99]
a = [2.136e-4, 1.332e-4, 6.171e-5, 5.547e-5, 1.433e-4, 1.209e-4]
b = [2.533e-8, 1.580e-8, 1.698e-9, 1.526e-9, 1.707e-8, 1.404e-8]
c = [-1.332e-11, -8.308e-12, -3.198e-12, -2.874e-12, -8.968e-12, -7.418e-12]
d = [-3.097e-15, -1.931e-15, -5.575e-16, -5.011e-16, -2.084e-15, -1.713e-15]

for index, row in workingTableExoplanets.iterrows():
    workingTableExoplanets.at[index, "st_recentVenus"]  = seffsun[0] + a[0] * (workingTableExoplanets.at[index, "star_teff"] - 5780.0) + b[0] * (workingTableExoplanets.at[index,"star_teff"] - 5780.0) ** 2 + c[0] * (workingTableExoplanets.at[index,"star_teff"] - 5780.0) ** 3 + d[0] * (workingTableExoplanets.at[index,"star_teff"] - 5780.0) ** 4
    workingTableExoplanets.at[index, "st_runawayGreenhouse"]  = seffsun[1] + a[1] * (workingTableExoplanets.at[index, "star_teff"] - 5780.0) + b[1] * (workingTableExoplanets.at[index,"star_teff"] - 5780.0) ** 2 + c[1] * (workingTableExoplanets.at[index,"star_teff"] - 5780.0) ** 3 + d[1] * (workingTableExoplanets.at[index,"star_teff"] - 5780.0) ** 4
    workingTableExoplanets.at[index, "st_maxGreenhouse"]  = seffsun[2] + a[2] * (workingTableExoplanets.at[index, "star_teff"] - 5780.0) + b[2] * (workingTableExoplanets.at[index,"star_teff"] - 5780.0) ** 2 + c[2] * (workingTableExoplanets.at[index,"star_teff"] - 5780.0) ** 3 + d[2] * (workingTableExoplanets.at[index,"star_teff"] - 5780.0) ** 4
    workingTableExoplanets.at[index, "st_earlyMars"]  = seffsun[3] + a[3] * (workingTableExoplanets.at[index, "star_teff"] - 5780.0) + b[3] * (workingTableExoplanets.at[index,"star_teff"] - 5780.0) ** 2 + c[3] * (workingTableExoplanets.at[index,"star_teff"] - 5780.0) ** 3 + d[3] * (workingTableExoplanets.at[index,"star_teff"] - 5780.0) ** 4
    workingTableExoplanets.at[index, "st_half_Earth"]  = seffsun[4] + a[4] * (workingTableExoplanets.at[index, "star_teff"] - 5780.0) + b[4] * (workingTableExoplanets.at[index,"star_teff"] - 5780.0) ** 2 + c[4] * (workingTableExoplanets.at[index,"star_teff"] - 5780.0) ** 3 + d[4] * (workingTableExoplanets.at[index,"star_teff"] - 5780.0) ** 4
    workingTableExoplanets.at[index, "st_five_Earth"]  = seffsun[5] + a[5] * (workingTableExoplanets.at[index, "star_teff"] - 5780.0) + b[5] * (workingTableExoplanets.at[index,"star_teff"] - 5780.0) ** 2 + c[5] * (workingTableExoplanets.at[index,"star_teff"] - 5780.0) ** 3 + d[5] * (workingTableExoplanets.at[index,"star_teff"] - 5780.0) ** 4

workingTableExoplanets.to_pickle(f"./data/star_pickles/{queryStar}.pkl")
