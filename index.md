# The Difference between Composite and Regressed OLR on RMM index 
 
## Kai Huang

## Introduction

My CLIM 680 project investigates the difference between the composited OLR and regressed OLR on the RMM1 and RMM2 index of Madden-Julian Oscillation (MJO). I choose this topic because the composite and linear regression have been widely used in climate research, but I wanted to see how they are different in representing the OLR structure of MJO. Here I consider the linearly regressed OLR as the linear mode of MJO because it only capture the OLR with same periodicity, while the composited OLR contains variations of all time scales thus both linear and nonlinear structures are included.

## Data

The datasets used in my project are:

__Era-Interim zonal wind at 850hPa and 200hPa__
The Era-Interim reanalysis (stored in `/shared/ccsm4/khuang/obs/era-interim/ on COLA servers`) dataset is daily global data from Jan, 1979 to Dec, 2017. It has horizontal resoluiton of 2.5 deg longitude by 2.5 deg latitude.
 
__NOAA Outgoing Longwave Radiation__
The NOAA Outgoing Longwave Radiation (stored in `/shared/ccsm4/khuang/obs/NOAA-OLR/ on COLA servers`) is daily global Outgoing Longwave Radiation at the top of atmosphere from Jun 1979 to Dec 2019. It is on a 2.5 deg longitude by 2.5 deg latitude gride.

## Results and Codes

### Functions
I created a set of functions in `clim_utils.py` for doing common tasks used throughout my analysis, including:
* `MonthlyClimatology(ds)` and `DailyClimatology(ds)` for calculating monthly and daily climatology
* `label_latlon(ax,lons,lats)` for labelling plots
* `DailyAnomaly(ds)` for calculating daily anomalies
* `LinearRegression_Map(v,p,ndimlat,ndimlon)` for calculating the linearly regressed 2-dimensional v on predictor p
* `rmm(ds_u850,ds_u200,ds_olr)` for calculating the two leading EOFs and corresponding PCs of RMM index


### Conda Environment

The environment.yml file is provided to define the environment needed to run all codes.

### Figures

The figures from my notebooks are saved in the `Figs` subdirectory and shown in my notebooks.

### Analyses and Notebooks
 
#### [Annual Mean](https://github.com/KaiHuang94/CLIM680/blob/master/project.annual.mean.ipynb)  

The annual mean of OLR. High value in tropics and decreases poleward. 

#### [Climatology](https://github.com/KaiHuang94/CLIM680/blob/master/project.monthly.climatology.ipynb)

The monthly climatology of OLR. Similar patterns as the annual mean, but with relatively higher values in the summer hemisphere. 

#### [Calculate the RMM index](https://github.com/KaiHuang94/CLIM680/blob/master/project.rmm.ipynb)

Two leading EOFs and corresponding PCs based on the meridional means (15S-15N) of u850, u200, and OLR are constructed. 
EOF1 illustrates the mode of MJO with its convection active center over the Indian center, and a baroclinic circulation with convergence at low troposphere and divergence at high troposphere. 
EOF2 illustrates the mode of MJO with its convection avtive center over the west Pacific while its depressed convection center over Indian Ocean, also accompained with a baroclinic circulation.

#### [Calculate OLR Composites of EOF1 and EOF2](https://github.com/KaiHuang94/CLIM680/blob/master/project.rmm.composite.ipynb)

Link to notebook and describe results

####  [Calculate the Mean Difference between Composites with Significance]()

Link to notebook and describe results

#### Calculate and compare Composites with NAO index  

Link to notebook and describe results

#### Correlation of ENSO Composites for two time periods

Link to notebook and describe results

#### Regression of Nino3.4 with Temperature and Precipitation

Link to notebook and describe results

#### EOFs of monsthly SSTs

Link to notebook and describe results

## Summary

Provide short summary of what you learned from your analysis of your data (both scientific and technical), what you plan to do next, and any challenges or issues you encountered/overcame. 
