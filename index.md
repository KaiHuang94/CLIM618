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

The OLR composites of two EOFs are done by compositing the daily OLR anomalies when the corresponding PC is greater than 1 standard deviation. 
Composited OLR anomalies for RMM1 shows a strong organized large-scale convection over the Indian Ocean, also shows some active convections in ITCZ and SPCZ regions. 
Composited OLR anomalies for RMM2 shows a dipole mode with depressed convections over west Indian Ocean, and active convections over west Pacific. The active convections in ITCZ and SPCZ are also evident. 

####  [Calculate the Mean Difference between Composites with Significance](https://github.com/KaiHuang94/CLIM680/blob/master/project.rmm.composite.ipynb)

The Difference between strong RMM1, RMM2, and no-MJO OLR composites are shown. 
Most negative/positive OLR values with great magnitude reach the 95% signigicance level. 

#### [Regression of OLR against PC1 and PC2](https://github.com/KaiHuang94/CLIM680/blob/master/project.rmm.regression.ipynb)

The regressed OLR patterns against PC1 and PC2 are similar to the composited OLR patterns, but with a smaller magnitude. 

#### [Differences between composite and regressed OLR](https://github.com/KaiHuang94/CLIM680/blob/master/project.rmm.difference.ipynb)

The differences between composite and regressed OLR is calculated by taking their difference directly. Results show that the difference is very like the linearly regressed pattern, indicating the main contributer to it is the magnitude difference. It makes sense becaues the convections of MJO contains both linear mode and nonlinear mode. Therefore, directly taking the difference does not telling much information. 
To better demostrate the structure difference of composited and regressed OLR, before taking the difference, each pattern is devided by its pattern standard deviation. 
After being divided by its pattern standard deviation, the OLR fields becomes unitless and they represents the structures with same magnitude.
Both the direct taking difference and taking difference after normalization are done for OLR fields only with values passing 95% significance level and for OLR fields with all values. 
The differences after normalization show when the MJO convection center is over the Indian Ocean, the nonlinear mode depress the convection over Maritime Continent, enhancing the convections in ITCZ. When the MJO convection center is over the Maritime Continent, the nonlinear mode is relatively weak, but still with a depressed convection over MC and central Pacific. 

## Summary
### Findings
The composited OLR field on RMM index is considered as the total mode of MJO convection. It contains both the linear and nonlinear mode. The linearly regressed OLR field on RMM index is considered as the linear mode of MJO convection. 
The direct difference between comsited and regressed OLR field is very similar to that of the original composite field. Therefore, the nonlinear mode amplifies the MJO convections. 
By normalizing composited and regressed OLR, the unitless fields reflects the structures of linear and nonlinear mode with the same magnitude. Therefore, by taking the differences between normalized fields, the structural differences are clearly shown. 
Over the Maritime Continent, the nonlinear mode, compared with linear mode, is always with a weaker convection magnitude. 
Over the ITCZ and west Indian, the nonlinear mode is with a relatively stronger convection magnitude. 
It is interesting that over the MC, the nonlinear mode is always with a relatively weaker convection. It may be realted to the advection terms, because the surface temperature and moisture is always high there due to the warm pool. Considering that, the advections in or out of that region always lead to negative temperature/moisture tendencies. But this is only a hypothesis. Over MC, there are many nonlinear processes, such as the sub-grid convection and cloud processes, radiations, and the atmosphere/ocean mixed layer processes. 
Therefore, more dynamic processes should be diagnosed to investigate the above hypothesis. 

### One technical problem
I tried to read in multiple data files using `xarray.open_mfdataset`. 
The dataset is very large with daily 4D (37 vertial layers, 2.5deg horizontal resolution) outputs from 1979 to 2017 (50+ GB). It exceeds the memory size and therefore I could not read in the whole dataset.
Is there any way in Python that I can read in the whole dataset, mabe by chunking the data? 