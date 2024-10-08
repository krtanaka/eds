
# eds

<!-- badges: start -->
<!-- badges: end -->

This package lets you run Environmental Data Summary (EDS) that merges survey data with external environmental data, such as chlorophyll A, photosynthetically active radiation (PAR), turbidity (Kd490), wave height, and sea surface temperature (SST). 

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

