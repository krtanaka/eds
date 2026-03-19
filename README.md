
# eds

<!-- badges: start -->
<!-- badges: end -->

The eds package provides an automated framework for the spatiotemporal co-location and temporal summarization of gridded satellite and oceanographic model datasets with time-stamped and georeferenced field observations. It is designed to bridge the gap between biological surveys (e.g., in situ fish or coral counts) and the environmental conditions that drive ecosystem dynamics.

Rather than relying on broad regional averages, eds automates the extraction of localized environmental data for specific survey coordinates and timestamps.

Core Capabilities:

Data Integration: Seamlessly merges field survey data with hundreds of external environmental datasets. Supported variables include Sea Surface Temperature (SST), Chlorophyll-a (Chl-a), Photosynthetically Active Radiation (PAR), Turbidity (Kd490), Wave Height, Precipitation, Population Density, and Effluent.

Temporal Summarization: Ecosystems are shaped by past conditions, not just a single-day snapshot. eds calculates the exact environmental conditions leading up to your survey date using customizable time windows (e.g., the previous 7 days, 30 days, or 10 years).

Automated Metrics: Computes standard historical summaries (Mean, Standard Deviation, Maximum, Minimum, Median, and Sum) to format the extracted environmental data into an analysis-ready structure for ecological and statistical modeling.


## Installation

You can install the development version of eds from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("krtanaka/eds")
```

## Example

``` r
library(eds)

# Define the parameters required for the EDS. This data frame includes information about each dataset to be downloaded and processed.
eds_parameter <- data.frame(
  Dataset = c("Bathymetry_ETOPO_2022_v1_15s", "Sea_Surface_Temperature_OISST_Monthly"),
  Download = c("YES", "YES"),
  Frequency = c("Climatology", "Monthly"),
  URL = c("https://coastwatch.pfeg.noaa.gov/erddap/", "https://upwell.pfeg.noaa.gov/erddap/"),
  Dataset_ID = c("ETOPO_2022_v1_15s", "noaa_psl_4af9_4ab0_ab10"),
  Fields = c("z", "sst"),
  Summaries = c(NA, "mean;q05;q95;sd"),
  Mask = c(FALSE, FALSE)
)

# Save the data frame to a CSV file on the desktop
write.csv(eds_parameter, file = file.path("/Users/", Sys.info()[7], "Desktop", "eds_parameters.csv"), row.names = FALSE)

# Load example dataset (NCRMP 2010-2022)
df <- subset(df, region == "MHI")

run_eds(lon = df$lon,
        lat = df$lat,
        unit = df$island,
        time = df$date)
```

