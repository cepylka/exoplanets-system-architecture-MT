# Exoplanets system architecture

## Preparing the sample

For retrieving information from working databases initially we need a list of potentially interesting systems, which can be formed by following script: `systems-min2planets-withMassSinAndRadius.py`.

This script allows to retrieve star names for the systems, where at least two planets have both mass and radius defined, by fetching parameter `star_name` from Exoplanet.eu database.

``` sh
$ python ./systems-min2planets-withMassSinAndRadius.py
```

You could run this script in the terminal. It will write in the `systems-min2planets-withMassSinAndRadius.txt` file, which could be edited manually, if you already knew what system you want to exclude. It will save some processing time.

For retrieving the whole information about the system in interest, you can run `task-manager.sh` which will work through the planetary system list in `systems-min2planets-withMassSinAndRadius.txt` , using script `one-system_rewrite_star.py` and write one `.pkl` file for each system in the designated folder, in this case `./data`

``` sh
$ ./task-manager.sh
```

At this stage you could visualise systems, and for that you would need data of solar system planets. The shell script you can run, `task-manager-charts.sh`, will load it from the folder `./data/solarsystem.pkl` and the interior curves data file `curves.txt` from [Zeng et al, 2019](#Zeng).

``` sh
$ ./task-manager-charts.sh
```

This script will work through the planetary systems files  `"system name".pkl` in the folder `./data/systems_for_chart/` using the other script `one-system-charts_hz.py`, and making figures for Habitable zone edges [1] and mass-radius diagrams [2] for every individual planetary system. There are the examples for Kepler-107 system:

![](./img/Kepler-107_with_HZ_edges.png)

![](./img/Kepler-107.png)

The HZ edges are calculated according to [Kopparapu et al., 2014](#Kopparapu).

In order to further investigate the data you need to concatenate all chosen system files to one, with which the certain jupyter notebooks will work, also this combined file in `.pkl` form will be needed to create working files for internal systems structure analyses. You can do the concatenation by running the script `merg-edges.py`, but you need to put all needed system's `.pkl` in the folder `./merg` first. The script creates file `all_my_systems.pkl`.

``` sh
$ python ./merg-edges.py
```

In aforementioned folders: `./data`, `./data/systems_for_chart/`, `./merg` I put example planetary system .pkl files for checking how the scripts work. After doing the routine from the start, you will obtain such files of the new systems that meet your selection criteria. You can aso make your own file tree for your convenience, but do not forget to change paths if you do.

## Analyses of the sample

### System architecture: masses, radii similarities

For analysing "peas in the pod" tendency among masses and radii of adjacent planets in the system, you can create a new `.pkl` file from existing sample file `all_my_systems.pkl` by using `pairs_MRP_stellar.py`, which you can run in the terminal.

``` sh
$ python ./pairs_MRP_stellar.py
```

It will write a new file `MRP_data_sets.pkl` in the `data` folder. The file consists now paired data for two adjacent planets in the system as one row in the pandas data frame.

For analyses of similarity in parameters of adjacent planet pairs and its depending on stellar parameters you can use the following script: `MRD_adjacent_planets.py`.

``` sh
$ python ./MRD_adjacent_planets.py
```

First, this script will create general plots for parameters: mass, radius, and density, and for two populations: P<sub>i</sub> for the inner planet in the pair of adjacent planets, P<sub>i+1</sub> for the outer planet in this pair, where P is one of the parameters and i=1,2 for first and second planet in the system by ascending of their semi-major axis. Note, that one system can produce multiple pairs, if it has more that 2 planets with well-defined parameters. The Pearson correlation test will be conducted and R- and P-values will be calculated, and the result along with the ordering in pairs calculations will be output to stdout.

Then, the script calculates the Pearson coefficient distributions for the Monte Carlo random uniform sampling populations from error intervals (10<sup>5</sup> attempts), and resulting figures will represent R- and P-values distributions, the normal distribution function fitting. And the mean and median coefficient values and the standard deviation $\sigma$ value will be printed in every figure.

In the third part, the script will make "moving window" test for different stellar parameters for sub-samples of x=40 data points with n=2 step. It will produce graphs of R- and P-value dynamics in stellar mass, radius, metallicity, T<sub>eff</sub> and age ranges.


For analysing periods ratio, you can create a new `.pkl` file from existing sample file `all_my_systems.pkl` by using `triples_MRPD.py`, which you can run in the terminal.

``` sh
$ python ./triples_MRPD.py
```

It will write a new file `triples_MR.pkl` in the `data` folder. The file consists now data for three adjacent planets in the system as one row in the pandas data frame.

For analyses of similarity in parameters of adjacent planet triples, and for assessing how these similarity trends depend on stellar parameters, you can use the following script: `P_adjacent_planets.py`.

``` sh
$ python ./P_adjacent_planets.py
```

Analogously with the script for masses, radii and densities (`MRD_adjacent_planets.py`), the `P_adjacent_planets.py` script calculates R- and P-values and output to stdout, creates general plots for period ratios, and for two populations: P<sub>i+1</sub>/P<sub>i</sub> for the inner planet pair of adjacent planets, P<sub>i+2</sub>/P<sub>i+1</sub>  for the outer planet pair, where P is the orbital period and i=1,2,3 for first, second and third planet in the system by ascending of their semi-major axis. Note, that one system can produce multiple triples, if it has more that 3 planets with well-defined parameters.

Second, the script calculates the Pearson coefficient distributions from random uniform simulations from error intervals (10<sup>5</sup> attempts), and plots figures of R- and P-values distributions, the normal distribution function fitting along with the mean and median coefficient values and the standard deviation $\sigma$ value.

After that, the script will make "moving window" test for different stellar parameters for sub-samples of x=40 data points with n=2 step. It will produce graphs of R- and P-value dynamics in stellar mass, radius, metallicity, T<sub>eff</sub> and age ranges.



References:

<a name="Kopparapu"></a> Kopparapu, R. K., Ramirez, R. M., SchottelKotte, J., Kasting, J. F., Domagal-Goldman, S. and Eymet, V. (2014). "Habitable zones around main-sequence stars: dependence
on planetary mass". In: The Astrophysical Journal Letters vol. 787, no. 2, p. L29.

<a name="Zeng"></a> Zeng, L., Jacobsen, S. B., Sasselov, D. D., Petaev, M. I., Vanderburg, A., Lopez-Morales, M., Perez-Mercader, J., Mattsson, T. R., Li, G., Heising, M. Z. et al. (2019). "Growth model interpretation of planet size distribution". In: Proceedings of the National Academy of Sciences vol. 116, no. 20, pp. 9723–9728.
