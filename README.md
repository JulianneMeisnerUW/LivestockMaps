# LivestockMaps
A time series of cattle and pig density (ratio of livestock/humans) maps (0.017 decimal degrees) for Uganda (2000-2020), Malawi (2000-2020), and DRC (2008-2015), median and width of posterior 95% credible interval. County-level 2008 density estimates for South Sudan, median and sd. 

## Data sources for livestock maps:

### Livestock data
* [IHSN survey catalog](https://catalog.ihsn.org/catalog)
* [IPUMS](https://ipums.org)
* [IPUMS DHS](https://www.idhsdata.org/idhs/)

### Predictors
* [NOAA nighttime lights](https://www.worldpop.org/geodata/listing?id=75)
* [WorldPop population density](https://ngdc.noaa.gov/eog/dmsp/downloadV4composites.html)
* [Elevation (GMTED2010)](https://www.usgs.gov/centers/eros/science/usgs-eros-archive-digital-elevation-global-multi-resolution-terrain-elevation?qt-science_center_objects=0#qt-science_center_objects)
* [Water bodies (Global Lakes and Wetlands)](https://www.worldwildlife.org/pages/global-lakes-and-wetlands-database)
* [Protected areas (World Database of Protected Areas)](https://www.protectedplanet.net/en)
* [WorldPop under US$/day](https://www.worldpop.org/project/categories?id=9)

### Shapefiles
* [GADM](https://gadm.org/download_country_v3.html)
* [Settlements (DRC)](https://cod-data.forest-atlas.org/datasets/eaa138a4908b4d1588e6ba3d21ea5698_0/data)
- (may be missing some)

## Data sources for g-formula and regression implementation

### Outcome (HAT) data source:
* [WHO Atlas of HAT](https://www.who.int/trypanosomiasis_african/country/foci_AFRO/en/)

### Denominator (offset; population)
* [WorldPop population counts](https://www.worldpop.org/project/categories?id=3)

### Confounders and mediators
* [NDVI and LST](https://ladsweb.modaps.eosdis.nasa.gov)
* [Armed conflict (UCDP)](https://ucdp.uu.se)
* [Disasters (EM-DAT)](https://www.emdat.be)

## Recommended file structure:
### Codes:
```bash
codes/
```
### Data:

## Order for running codes: 

### SPDE countries (Malawi, Uganda, DRC):

- (1) data_processing (country-specific)
- (2) dataset_prep (country-specific)
- (3) urbanicity_prep
- (4) urbanicity_validation
- (5) urbanicity_regression_surface
- (6) workingfile, which calls models_cattle, then models_pigs, then prediction_cattle, then prediction_pigs, then validation
- (7) livestock_model selection
- (8) data_processing_wealth
- (9) workingfile_wealth, which calls models_wealth, and validation_wealth
- (10) wealth_model_selection
- (11) data_processing_regression, choose the two models
- (12) GLW regression section of livestock_model_selection, and WorldPop regression section of wealth_model_selection
- (13) regression_prep
- (14) index_events (country-specific)
- (15) remote_sense (country-specific)
- (16) gform_prep 
- (17) plots
- (18) gmethod_nonsp (country-specific)
- (19) gmethod_nonsp_med (country-specific)
- (20) gform_plot
- (21) gformula_results

### Functions code:
- (1) my_functions_mw (used for data processing)
- (2) functions_prediction (used for working with INLA models)
- (3) gformula_functions (used for gformula stuff)

### ICAR (South Sudan) 
- (1) data_processing_ss
- (2) models_ss
- (3) data_processing_wealth_ss
- (4) outcome_data_processing_SS
- (5) remote_sense_ss
- (6) index_events_ss
- (7) regression_ss
- (8) plots (same one as for other countries, above)
