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

### Shapefiles
* [GADM](https://gadm.org/download_country_v3.html)
* [Settlements (DRC)] (https://cod-data.forest-atlas.org/datasets/eaa138a4908b4d1588e6ba3d21ea5698_0/data)

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
