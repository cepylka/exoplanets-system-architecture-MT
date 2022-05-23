import sys
import numpy
import pandas
import pyvo
import pathlib

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
from matplotlib.colors import LogNorm
from collections import OrderedDict
from matplotlib.lines import Line2D
import itertools
from functools import partial
import matplotlib.ticker as mticker
from cycler import cycler
from astropy import units
from astropy import constants
pandas.options.mode.use_inf_as_na = True


if len(sys.argv) < 2:
    raise SystemExit("No pickle file name provided")

if len(sys.argv) > 2:
    raise SystemExit("Too many arguments")

pickleFileName = pathlib.Path(sys.argv[1])
print(f"Processing pickle file: {pickleFileName}")

resultingName = pickleFileName.stem
pickleFileName = "./data/systems_for_chart" / pickleFileName
workingTableExoplanets = pandas.read_pickle(pickleFileName)

quantity_support()
set_matplotlib_formats('svg')

mass_coeff = constants.M_jup / constants.M_earth
radius_coeff = constants.R_jup / constants.R_earth
# interior cirves data
curves = pandas.read_fwf("./data/curves_all.txt", index_col=0)#, header=None)

print("[1/2] Making charts #1...")
# solar system data
solarsystemTable1 = pandas.read_pickle("./data/solarsystemE.pkl")
for index, row in workingTableExoplanets.iterrows():
    if not pandas.isnull(row["semi_major_axis"]):
        continue
    else:
        if not pandas.isnull(row["period"]):
            workingTableExoplanets.at[index, "semi_major_axis"] = numpy.power((workingTableExoplanets.at[index, "star_mass"] * constants.M_sun) * constants.G/(4 * numpy.pi **2) * ((workingTableExoplanets.at[index, "period"]* 86400 * units.s) ** 2 ) , 1/3) / constants.au


workingTableExoplanets["luminosity"] = numpy.array(numpy.NaN, dtype=float)
for index, row in workingTableExoplanets.iterrows():
    if not pandas.isnull(row["st_lum"]):
        workingTableExoplanets.at[index, "luminosity"] = 10 ** workingTableExoplanets.at[index, "st_lum"]
    else:
        workingTableExoplanets.at[index, "luminosity"] = (4 * numpy.pi * (workingTableExoplanets.at[index, "star_radius"] * constants.R_sun.value) ** 2 * constants.sigma_sb.value * (workingTableExoplanets.at[index, "star_teff"] ** 4)) / constants.L_sun.value


earth_albedo = 0.306
fig, ax = plt.subplots(figsize=(7, 7))

# conversion to Earth's radius
radius = workingTableExoplanets["radius"] * radius_coeff
# conversion to Earth's mass
mass = workingTableExoplanets["mass"] * mass_coeff

i_units = "kg m-3"
def label_point(x, y, val, ax):
    a = pandas.concat({'x': x, 'y': y, 'val': val}, axis=1)
    for i, point in a.iterrows():
        ax.text(point['x']*1.05, point['y']*0.98, str(point['val']), fontsize=8, zorder=3)

star_t = workingTableExoplanets["star_teff"]
star_radius = workingTableExoplanets["star_radius"]
planet_distance = workingTableExoplanets["semi_major_axis"]

radius = workingTableExoplanets["radius"] * radius_coeff

mass = workingTableExoplanets["mass"] * mass_coeff
orbital = workingTableExoplanets["period"]
distance = workingTableExoplanets["semi_major_axis"]
if not numpy.isnan(star_t[0]) and not numpy.isnan(star_radius[0]):
    print(type(star_t[0]),star_radius[0])
    seaborn.set_style("white")

    cmap = seaborn.cubehelix_palette(rot=-.2, as_cmap=True)
    g = seaborn.relplot(
        facet_kws=dict(despine=False),
        x=distance, y=mass,
        palette=cmap,
    )
    if numpy.nanmax(mass)/numpy.nanmin(mass) >= 50. or numpy.nanmax(mass) >= 500.:
        ax.set_yscale('log')

    recentVenus = numpy.power((workingTableExoplanets["luminosity"] / workingTableExoplanets["st_recentVenus"]), 0.5)
    runawayGreenhouse = numpy.power((workingTableExoplanets["luminosity"] / workingTableExoplanets["st_runawayGreenhouse"]), 0.5)
    maxGreenhouse = numpy.power((workingTableExoplanets["luminosity"] / workingTableExoplanets["st_maxGreenhouse"]), 0.5)
    earlyMars = numpy.power((workingTableExoplanets["luminosity"] / workingTableExoplanets["st_earlyMars"]), 0.5)
    half_Earth = numpy.power((workingTableExoplanets["luminosity"] / workingTableExoplanets["st_half_Earth"]), 0.5)
    five_Earth = numpy.power((workingTableExoplanets["luminosity"] / workingTableExoplanets["st_five_Earth"]), 0.5)
    plt.ylabel(r"Planet mass [M$_\oplus$]")
    plt.xlabel(r"Orbital distance [au] ")
    # plottong HZ edges
    plt.plot([numpy.nanmax(recentVenus), numpy.nanmax(recentVenus)], [0, numpy.nanmax(mass)+0.5], color='r', label="recent Venus")
    plt.plot([numpy.nanmax(runawayGreenhouse), numpy.nanmax(runawayGreenhouse)], [0, numpy.nanmax(mass)+0.5], color='g', label="runaway Greenhouse")
    plt.plot([numpy.nanmax(half_Earth), numpy.nanmax(half_Earth)], [0, numpy.nanmax(mass)+0.5], color='g', linestyle="dotted", label="0.5 Earth runaway Greenhouse")
    plt.plot([numpy.nanmax(five_Earth), numpy.nanmax(five_Earth)], [0, numpy.nanmax(mass)+0.5], color='g', linestyle="dashed", label="5 Earth runaway Greenhouse")
    plt.plot([numpy.nanmax(maxGreenhouse), numpy.nanmax(maxGreenhouse)], [0, numpy.nanmax(mass)+0.5], color='g', label="max Greenhouse")
    plt.plot([numpy.nanmax(earlyMars), numpy.nanmax(earlyMars)], [0, numpy.nanmax(mass)+0.5], color='r', label="early Mars")

    label_point(distance, mass, workingTableExoplanets.index.to_series(name='index'), plt.gca())
    plt.tick_params(axis='both', which='major', labelsize=12)
    g.set(xscale="log")

    # system data to table
    if numpy.nanmax(earlyMars) < numpy.nanmax(distance):
        g.set(xlim=(min(distance)*0.5, numpy.nanmax(distance)*2.5))
    elif min(recentVenus) < min(distance):
        g.set(xlim=(min(recentVenus)*0.5, numpy.nanmax(distance)*2.5))
    else:
        g.set(xlim=(min(distance)*0.5, numpy.nanmax(earlyMars)*2.5))
    star = pandas.DataFrame(
    )
    star["value"] = numpy.array(numpy.NaN, dtype=str)
    star.at["Star name", "value"] = workingTableExoplanets['star_name'][0]
    star.at["Planets", "value"] = workingTableExoplanets['sy_pnum'][0] if not pandas.isna(workingTableExoplanets['sy_pnum'][0]) else ""
    star.at["Star Teff", "value"] = workingTableExoplanets['star_teff'][0] if not pandas.isna(workingTableExoplanets['star_teff'][0]) else ""
    star.at["Star spectral type", "value"] = workingTableExoplanets['star_spec_type'][0] if not pandas.isna(workingTableExoplanets['star_spec_type'][0]) else ""
    star.at["Star mass", "value"] = workingTableExoplanets['star_mass'][0] if not pandas.isna(workingTableExoplanets['star_mass'][0]) else ""
    star.at["Star age", "value"] = workingTableExoplanets['star_age'][0] if not pandas.isna(workingTableExoplanets['star_age'][0]) else ""
    star.at["Star metallicity", "value"] = workingTableExoplanets['star_metallicity'][0] if not pandas.isna(workingTableExoplanets['star_metallicity'][0]) else ""
    if not pandas.isnull(workingTableExoplanets['star_metallicity'][0]) and not pandas.isna(workingTableExoplanets['st_metratio'][0]):
        sitsRow = pandas.DataFrame({"value": workingTableExoplanets['st_metratio'][0]}, index=["Met. ratio"])
        star = star.append(sitsRow)
    if workingTableExoplanets['sy_snum'][0] > 1:
        sitsRow1 = pandas.DataFrame({"value": workingTableExoplanets['sy_snum'][0]}, index=["Stars in the system"])
        star = star.append(sitsRow1)
    c=0
    for index, row in workingTableExoplanets.iterrows():
        if not pandas.isnull (row["pl_massjlim"]):
            c+=1
    if c !=0:
        sitsRow2 = pandas.DataFrame({"value": c}, index=["Mass uncertain in"])
        star = star.append(sitsRow2)
    d=0
    for index, row in workingTableExoplanets.iterrows():
        if not pandas.isnull(row["pl_radjlim"]):
            d+=1
    if d !=0:
        sitsRow3 = pandas.DataFrame({"value": d}, index=["Radius uncertain in"])
        star = star.append(sitsRow3)

    ax = plt.table(cellText=star.values, rowLabels=star.index, bbox=(1.4, 0.65, 0.35, 0.035*len(star.index)))
    ax.auto_set_font_size(False)
    ax.set_fontsize(10)
    plt.grid(False)

    plt.title("HZ graph in %s system" % workingTableExoplanets['star_name'][0])
    plt.legend(bbox_to_anchor=[1.05,.09,.3,.5], loc='center left', frameon=False)
    # saving plot in .png and .svg
    plt.savefig(f"./charts/{resultingName}_with_HZ_edges.png", bbox_inches="tight")
    plt.savefig(f"./charts/{resultingName}_with_HZ_edges.svg", bbox_inches="tight")

print("[2/2] Making charts #2...")
seaborn.set_style("white")
vmin = 150.
vmax = 10000.
fig, ax = plt.subplots(figsize=(7, 7))
plt.grid(False)
ax.set_xscale('log')

#interior curves
levelsZ = ["100%Fe", "30%Fe", "Rocky", "40%H2O","100%H2O","ColdH2/He"]
colors=['#E7340E','#C16228', '#DF9706', '#1EB70C', '#0398FF','#add8e6']

for i in range(len(levelsZ)):
    ax.plot(curves[levelsZ[i]], color=colors[i], zorder=0) # Mearth Rearth
legend_elements = ['%s' % value for value in levelsZ]

#solar system
massJ = solarsystemTable1["massJ"] * mass_coeff
radiusJ = solarsystemTable1["radiusJ"] * radius_coeff
words = solarsystemTable1.index.values
solarsys = ax.scatter(massJ[:], radiusJ[:], marker="o", facecolor='tab:green', zorder=2)
# Solar system terrestrial included
if numpy.nanmax(mass) <= massJ[2] * 2:
    label_point(massJ[:4]*1.01, radiusJ[:4], solarsystemTable1.index.to_series(name='index'), plt.gca())
    xL = 0.0055
    xR = 10.
    yB = 0.055
    yT = 2.5
    ax.legend(legend_elements,  fontsize=9, facecolor='white',edgecolor='white', framealpha=1, bbox_to_anchor=(0.85, 0.35))
# Solar system big rocky and gice giants included
elif massJ[2] * 2  < numpy.nanmax(mass) < massJ[7] * 2 and numpy.nanmax(radius) < 4.5:
    xL = 0.125
    xR = 70.
    yB = 0.055
    yT = 5.5
    plt.text(massJ[1]*1.08, radiusJ[1]*.93, "Venus", color="k", fontsize=8, zorder=4)
    plt.text(massJ[2]*0.7, radiusJ[2]*1.01, "Earth", color="k", fontsize=8, zorder=4)
    plt.text(massJ[6]*0.615, radiusJ[6]*.99, "Uranus", color="k", fontsize=8, zorder=4)
    plt.text(massJ[7]*1.065, radiusJ[7]*.99, "Neptune", color="k", fontsize=8, zorder=4)
    # label_point(massJ[1:3], radiusJ[1:3]*0.99, solarsystemTable1.index.to_series(name='index'), plt.gca())
    # label_point(massJ[6:], radiusJ[6:], solarsystemTable1.index.to_series(name='index'), plt.gca())

    ax.legend(legend_elements, fontsize=9, facecolor='white',edgecolor='white', framealpha=1, bbox_to_anchor=(0.35, 0.35))
# Solar system gas giants included
elif numpy.nanmax(mass) >= massJ[7] * 2 or numpy.nanmax(radius) >= 4.5:
    plt.text(massJ[4]*1.04, radiusJ[4]*0.97, "Jupiter", color="k", fontsize=8, zorder=4)
    plt.text(massJ[5]*1.09, radiusJ[5]*0.97, "Saturn", color="k", fontsize=8, zorder=4)
    plt.text(massJ[1]*0.95, radiusJ[1]*.56, "Venus", color="k", fontsize=8, zorder=4)
    plt.text(massJ[2]*0.9, radiusJ[2]*1.2, "Earth", color="k", fontsize=8, zorder=4)
    plt.text(massJ[6]*0.39, radiusJ[6]*0.99, "Uranus", color="k", fontsize=8, zorder=4)
    plt.text(massJ[7]*1.11, radiusJ[7]*0.99, "Neptune", color="k", fontsize=8, zorder=4)
    xL = 0.125
    if numpy.nanmax(mass) > massJ[4] * 2:
        xR = numpy.nanmax(mass)*4.5
    else:
        xR = 5e4
    yB = 0.055
    if numpy.nanmax(radius) >= 15.:
        yT = numpy.nanmax(radius) + 0.5
    else:
        yT = 17.5
    ax.legend(legend_elements, fontsize=9, facecolor='white',edgecolor='white', framealpha=1, loc="best")

ax.set_xlim((xL,xR))
ax.set_ylim((yB,yT))
# fig, ax = plt.subplots()
# density curves plot
nx, ny = (100, 100)
massgrid = numpy.logspace(numpy.log10(xL), numpy.log10(xR), nx)
raduisgrid = numpy.linspace(yB, yT, ny)
massgridX, raduisgridY = numpy.meshgrid(massgrid, raduisgrid)

rho = [100,300,1000,3000,10000,30000]
#[0.1, 0.3, 1., 3., 10., 30.]
density = (3 * massgridX * constants.M_earth)/ (4 * numpy.pi * (raduisgridY * constants.R_earth) ** 3)  / units.Quantity(1,unit=i_units)
contourplotD = ax.contour(massgridX, raduisgridY, density, levels=rho, colors=['#D3D3D3'], extend='both',linestyles='dashed',zorder=-1)
# manual_locations = [(0.4, 0.4), (0.4, 0.4), (0.4, 0.4), (0.4, 0.4), (0.4, 0.4), (0.4, 0.4)]
ax.clabel(contourplotD, contourplotD.levels, fontsize=9, inline=1, inline_spacing=7, fmt='%1.0f kg m-3', use_clabeltext=True,zorder=-1)
# interior curves plot


#     plt.text(0.32990, curves.at[0.32990, levelsZ[i]] - 0.1 + 0.06 * i, levelsZ[i], color=colors[i], fontsize=9, bbox=dict(facecolor='white', edgecolor='none', zorder=0))
# plt.text(0.32990, 0.65, levelsZ[1], color=colors[1], fontsize=9, bbox=dict(facecolor='white', edgecolor='none', zorder=0))
    # plt.text(0.32990, 0.5 + 0.6 * i, levelsZ[i], color=colors[i], fontsize=9, bbox=dict(facecolor='white', edgecolor='none'))
xerr = numpy.nan_to_num(workingTableExoplanets[["mass_error_min", "mass_error_max"]].to_numpy().T, posinf=0.) * mass_coeff
yerr = numpy.nan_to_num(workingTableExoplanets[["radius_error_min", "radius_error_max"]].to_numpy().T, posinf=0.) * radius_coeff
# planet's names

ax.errorbar(mass, radius, xerr=numpy.abs(xerr), yerr=numpy.abs(yerr), ls='none', fmt='0.8', ecolor='tab:gray', elinewidth=0.8, capsize=None, barsabove=True, zorder=1)
scatter = ax.scatter(mass, radius, marker="o", facecolor='tab:gray', zorder=1)
# for cifer in range(len(workingTableExoplanets.index)):
#     plt.text(mass[cifer]-.2, radius[cifer]-.15, workingTableExoplanets.index[cifer], fontsize=9, zorder=2)

# flux temperature calculations
earth_albedo = 0.306
planet_temperature = star_t  * numpy.power((star_radius * constants.R_sun.value /(2 * planet_distance * constants.au.value)), (1/2)) * numpy.power((1 - earth_albedo),(1/4))
# im = ax.scatter(mass, radius, c=planet_temperature, norm=matplotlib.colors.LogNorm(), marker="o", cmap='jet', zorder=1) # for logarithmic temperatures
im = ax.scatter(mass, radius, c=planet_temperature, norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax), marker="o", cmap='jet', zorder=2)
# for cifer in range(len(workingTableExoplanets.index)):
#     ax.annotate(workingTableExoplanets.index[cifer], (mass[cifer], radius[cifer]), fontsize=9)
label_point(mass, radius, workingTableExoplanets.index.to_series(name='index'), plt.gca())

legend_elements = ['%s' % value for value in levelsZ]
ax.legend(legend_elements, fontsize=9, facecolor='white',edgecolor='white', framealpha=1, bbox_to_anchor=(0.25, 0.75))

plt.xlabel(f"Mass [M$_\oplus$]")
plt.ylabel(f"Radius [R$_\oplus$]")



plt.title(f"MMR diagram for {resultingName}")

ax4 = fig.add_axes([0.65, 0.1, 0.19, 0.2])
ax4.get_xaxis().set_visible(False)
ax4.get_yaxis().set_visible(False)
ax4.spines["right"].set_visible(False)
ax4.spines["top"].set_visible(False)
ax4.spines["bottom"].set_visible(False)
ax4.spines["left"].set_visible(False)
ax4.patch.set_facecolor('white')
ax4.patch.set_alpha(0.)
ax4.tick_params(top=False)
ax4.tick_params(labeltop=False)
ax4.tick_params(right=False)
ax4.tick_params(labelright=False)
ax4.tick_params(bottom=False)
ax4.tick_params(labelbottom=False)
ax4.tick_params(left=False)
ax4.tick_params(labelleft=False)
cbar = plt.colorbar(im, ax=ax4, aspect=3.8, shrink=0.5)
cbar.ax.tick_params(labelsize=9)
cbar.ax.tick_params(which='minor', color='white')
cbar.ax.set_title('T,K')
plt.grid(False)
# saving plot in .png and .svg
plt.savefig(f"./charts/{resultingName}.png",bbox_inches="tight")
plt.savefig(f"./charts/{resultingName}.svg",bbox_inches="tight")
