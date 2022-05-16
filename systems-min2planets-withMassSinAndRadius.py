import numpy
import pandas
import pyvo

serviceExoplanets = pyvo.dal.TAPService("http://voparis-tap-planeto.obspm.fr/tap")

tableNameExoplanets = "exoplanet.epn_core"

fieldsExoplanets = [
    "mass",
    "radius",
    "star_name"
]

queryExoplanets = (
    " ".join((
        "SELECT",
        ", ".join(fieldsExoplanets),
        f"FROM {tableNameExoplanets}",
        "WHERE star_name IN",
        "(SELECT star_name",
        f"FROM {tableNameExoplanets}",
        "GROUP BY star_name",
        "HAVING COUNT(*) >= 2)"
    ))
)
print("Query to execute: ", queryExoplanets, "\n")

resultsExoplanets = serviceExoplanets.search(queryExoplanets).to_table().to_pandas()


starsWithMinumumSatisfyingParameters = resultsExoplanets.query("mass.notna() & radius.notna()").groupby("star_name").filter(lambda x: len(x) > 1)["star_name"].unique()#.drop_duplicates()
print(len(starsWithMinumumSatisfyingParameters), type(starsWithMinumumSatisfyingParameters))

numpy.savetxt('./systems-min2planets-withMassSinAndRadius.txt', starsWithMinumumSatisfyingParameters, fmt="%s")
