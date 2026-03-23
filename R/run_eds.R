#' Run Environmental Data Synthesis (EDS)
#'
#' Downloads and processes environmental netCDF files based on provided parameters.
#' This function merges multiple netCDF files into a single file and extracts relevant
#' climatological and time-series data for specified geographic points.
#'
#' @param lon Numeric vector of longitudes for the data points.
#' @param lat Numeric vector of latitudes for the data points.
#' @param unit Character vector specifying the unit (e.g., island name) for each point.
#' @param time Vector of dates corresponding to each data point.
#' @param buffer Numeric value specifying the buffer distance for spatial operations. Default is 0.25.
#'
#' @return The function writes merged netCDF files and saves extracted climatological
#' and time-series data to specified directories. It returns `NULL` invisibly.
#'
#' @importFrom raster raster crs<- merge writeRaster as.data.frame shift crs extent rotate
#' @importFrom sp coordinates
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom methods slot
#' @importFrom stats complete.cases quantile var
#' @importFrom utils head install.packages read.csv tail
#'
#' @export
run_eds = function(lon, lat, unit, time, buffer = 0.25) {

  length_no_na = function(x){return(length(x[!is.na(x)]))}

  expanding_extract = function(r, SpDF, Dists = c(0, 500, 1000)){

    require(raster)
    require(spatial)

    OutDF = data.frame(values = rep(NA, nrow(SpDF)),
                       Dist = rep(NA, nrow(SpDF)),
                       N = rep(NA, nrow(SpDF)))

    nDists = length(Dists)
    cnt = 1

    NAi = which(is.na(OutDF$values))
    NAsLeft = length(NAi) > 0

    while(cnt <= nDists & NAsLeft){

      NAi = which(is.na(OutDF$values))

      pull = raster::extract(x = r, y = SpDF[NAi, ],
                             buffer = Dists[cnt],
                             small = TRUE,
                             na.rm = TRUE)

      Nper = unlist(lapply(pull, length_no_na))

      OutDF$values[NAi] = unlist(lapply(pull, mean, na.rm = TRUE))
      OutDF$Dist[NAi] = Dists[cnt]
      OutDF$N[NAi] = Nper

      NAi = which(is.na(OutDF$values))
      NAsLeft = length(NAi) > 0

      cnt = cnt + 1
    }

    return(OutDF)

  }

  maskfun = function(x, na.rm = F, depth_threshold = -30, percent_threshold = 5){

    denom = length(x)
    numer = length(which(x > depth_threshold))
    outp = numer/denom

    if(outp > (percent_threshold/100)){
      return(NA)
    } else {
      return(1)
    }

  }

  merge_times_nc = function(filename1, filename2, variable_name, outfilename){

    file1  <- nc_open(filename1)
    file2  <- nc_open(filename2)

    #Get Dimensions for output
    f1Dims = NcReadDims(file1)

    # Just for one variable for now
    mat1 = ncvar_get(file1, variable_name)
    mat2 = ncvar_get(file2, variable_name)

    #Check for Temporal Overlap, remove any overlapping times from mat2
    t1 = ncvar_get(file1, "time")
    t2 = ncvar_get(file2, "time")

    tINT = intersect(t1, t2)

    if(length(tINT) > 0) {
      tINTi = match(tINT, t2)
      mat2 = mat2[, , -tINTi]
    }

    dat_new = abind(mat1, mat2, along = 3)

    eval(parse(text = paste0("var = file1$var[variable_name]$", variable_name)))

    #Get dimensions - Keep in Same Orientation as first input file
    DEC_LON=file1$dim$longitude$vals[1]>file1$dim$longitude$vals[2]
    Dx = ncdim_def("longitude", "degrees_east", sort(unique(c(file1$dim$longitude$vals, file2$dim$longitude$vals)), decreasing = DEC_LON))

    DEC_LAT=file1$dim$latitude$vals[1]>file1$dim$latitude$vals[2]
    Dy = ncdim_def("latitude", "degrees_north", sort(unique(c(file1$dim$latitude$vals, file2$dim$latitude$vals)), decreasing = DEC_LAT))

    DEC_T=file1$dim$time$vals[1]>file1$dim$time$vals[2]
    if (is.na(DEC_T) == T) DEC_T = FALSE
    Dt = ncdim_def("time", "seconds since 1970-01-01T00:00:00Z", sort(unique(c(file1$dim$time$vals, file2$dim$time$vals)), decreasing = DEC_T))

    file_new  <- nc_create(
      filename = outfilename,
      vars = ncvar_def(
        name = variable_name,
        units = var$units,
        dim = list(Dx, Dy, Dt)))

    ncvar_put(
      nc = file_new,
      varid = variable_name,
      vals = dat_new)

    nc_close(file1)
    nc_close(file2)
    nc_close(file_new)
  }

  merge_times_nc_list = function(infilenamelist, variable_name, outfilename){

    Nfiles = length(infilenamelist)

    if(Nfiles == 0){
      print("No files to join")
      return(1)
    } else if (Nfiles == 1){
      cat("Single filename provided, copying/renaming file.\n")
      file.copy(from = infilenamelist[[1]], to = outfilename, overwrite = TRUE)
      return(1)
    } else if (Nfiles == 2){
      cat("Two filenames provided, calling 'merge_times_nc'.\n")
      status = merge_times_nc(filename1 = infilenamelist[[1]],
                              filename2 = infilenamelist[[2]],
                              variable_name = variable_name,
                              outfilename = outfilename)
      return(1)
    } else if (Nfiles > 2){
      cat("More than two filenames provided, calling 'merge_times_nc' successively.\n")
      status = merge_times_nc(filename1 = infilenamelist[[1]],
                              filename2 = infilenamelist[[2]],
                              variable_name = variable_name,
                              outfilename = outfilename)

      for(i in 3:length(infilenamelist)){
        print(paste0("First ", i-1, " files merged.."))
        status = merge_times_nc(filename1 = outfilename,
                                filename2 = infilenamelist[[i]],
                                variable_name = variable_name,
                                outfilename = outfilename)
      }
      print(paste0("All ", Nfiles, " complete. Written to ", outfilename))
    }
  }

  points_in_polys = function(pts.x, pts.y, bb){

    nbox = nrow(bb)
    npts = length(pts.x)
    InMat = matrix(0, ncol = nbox, nrow = npts)
    colnames(InMat) = bb$unit
    DataISL = rep("NONE_ASSIGNED", npts)

    for(i in 1:nbox){
      InMat[,i] = point.in.polygon(pts.x,
                                   pts.y,
                                   bb[i,c("x_min","x_max","x_max","x_min")],
                                   bb[i,c("y_min","y_min","y_max","y_max")])
      DataISL[which(InMat[,i] >= 1)] = as.vector(bb$unit[i])
    }

    out = list(DataISL,InMat)
    names(out) = c("DATA_UNIT","IN_MATRIX")

    return(out)
  }

  xyt_to_ijk = function(xyt_df, x_grid, y_grid, t_grid, lon, lat, t){

    xyt_df = data.frame(x = as.numeric(xyt_df[, 1]),
                        y = as.numeric(xyt_df[, 2]),
                        t = as.Date(xyt_df[, 3]))

    xyt_df$x[xyt_df$x > 180] = xyt_df$x[xyt_df$x > 180] - 360

    x_grid = as.numeric(x_grid)
    x_grid[x_grid>180] = x_grid[x_grid > 180] - 360
    y_grid = as.numeric(y_grid)
    t_grid = as.Date(t_grid)

    n_pts = nrow(xyt_df)
    n_xg = length(x_grid)
    n_yg = length(y_grid)
    n_tg = length(t_grid)

    x_buffer = mean(diff(lon))/2
    y_buffer = mean(diff(lat))/2
    t_buffer = mean(diff(t))/2

    #X Match
    x_gMat = outer(rep(1, n_pts), x_grid)
    x_pMat = outer(xyt_df$x, rep(1, n_xg))
    x_i = apply(abs(x_pMat-x_gMat), 1, which.min)
    oob_x = which(xyt_df$x < (min(x_grid, na.rm = T) - x_buffer) | xyt_df$x > (max(x_grid, na.rm = T) + x_buffer))
    x_i[oob_x] = NA

    #Y Match
    y_gMat = outer(rep(1, n_pts), y_grid)
    y_pMat = outer(xyt_df$y, rep(1, n_yg))
    y_j = apply(abs(y_pMat-y_gMat), 1, which.min)
    oob_y = which(xyt_df$y < (min(y_grid, na.rm = T) - y_buffer) | xyt_df$y > (max(y_grid, na.rm = T) + y_buffer))
    y_j[oob_y] = NA

    #T Match
    t_gMat = outer(rep(1, n_pts), t_grid)
    t_pMat = outer(xyt_df$t, rep(1, n_tg))
    t_k = apply(abs(t_pMat-t_gMat), 1, which.min)
    oob_t = which(xyt_df$t < (min(t_grid, na.rm = T) - t_buffer) | xyt_df$t > (max(t_grid, na.rm = T) + t_buffer))
    t_k[oob_t] = NA

    ijk = data.frame(x_i, y_j, t_k)

    return(ijk)
  }

  na_stack_count = function(x){
    return(length(which(is.na(x))))
  }

  q05 = function(x, na.rm = T){
    return(quantile(x, .05, na.rm = T))
  }

  q95 = function(x, na.rm = T){
    return(quantile(x, .95, na.rm = T))
  }

  options(warn = -1)

  if (!is.numeric(lon) || !is.numeric(lat) || !is.character(unit)) stop("Invalid input parameters. lon and lat should be numeric, unit should be character.")

  ### download ERDDAP ###
  # PLYR BEFORE DPLYR to prevent masking errors
  packages_to_load <- c("rerddap", "readr", "zoo", "ncdf4", "RNetCDF", "easyNCDF", "raster", "lubridate", "abind", "acss", "plyr", "dplyr", "doParallel", "foreach", "spatial", "data.table", "splitstackshape")

  install_and_load_package <- function(package_name) {
    if (!requireNamespace(package_name, quietly = TRUE)) {
      install.packages(package_name, dependencies = TRUE)
    }
    library(package_name, character.only = TRUE)
  }

  invisible(sapply(packages_to_load, install_and_load_package))

  if ("package:reshape2" %in% search()) {
    unloadNamespace("reshape2")
    install_and_load_package("reshape2")
  }

  select = dplyr::select

  cache_setup(temp_dir = T)
  cache_delete_all()
  closeAllConnections()

  if(.Platform$OS.type == "unix") {
    eds_path <- file.path("~/Desktop/eds/")
    source_path <- file.path("~/Desktop", "eds_parameters.csv")
    destination_path <- file.path("~/Desktop/eds", "eds_parameters.csv")
  } else {
    eds_path <- paste0("C:/Users/", Sys.info()[7], "/Desktop/eds/")
    source_path <- paste0("C:/Users/", Sys.info()[7], "/Desktop/", "eds_parameters.csv")
    destination_path <- paste0("C:/Users/", Sys.info()[7], "/Desktop/eds/", "eds_parameters.csv")
  }

  if (!dir.exists(eds_path)) dir.create(eds_path, recursive = T)
  if (!file.exists(source_path)) stop("The file eds_parameters.csv does not exist on your desktop.")

  file.copy(from = source_path, to = destination_path, overwrite = TRUE)

  ParamDF <- read_csv(file.path(paste0(eds_path, "/eds_parameters.csv")), show_col_types = FALSE) %>% dplyr::filter(Download == "YES")
  uP <- ParamDF %>% dplyr::pull(Dataset) %>% unique()

  df <- as.data.frame(list(lon = lon, lat = lat, unit = unit, date_r = time))
  df <- df %>% dplyr::filter(!is.na(lon) & !is.na(lat) & !is.na(date_r))

  df$unit = gsub(" ", "_", unit)

  # Explicit DPLYR calls to prevent PLYR overlap
  boxes <- df %>%
    dplyr::group_by(unit) %>%
    dplyr::summarise(x_min = min(lon) - buffer,
                     x_max = max(lon) + buffer,
                     y_min = min(lat) - buffer,
                     y_max = max(lat) + buffer, .groups = 'drop') %>%
    dplyr::mutate(dplyr::across(c(x_max, x_min, y_max, y_min), round, digits = 2))

  df_frame <- df %>%
    dplyr::select(unit) %>%
    dplyr::distinct()

  bbox <- dplyr::left_join(boxes, df_frame, by = "unit") %>% dplyr::mutate(unit = gsub(" ", "_", unit))

  uI <- dplyr::distinct(bbox, unit)$unit

  ### download ERDDAP netCDF files ###
  cat("\nEDS is currently downloading netCDF files...\n")

  for (iP in 1:length(uP)){

    thisp = subset(ParamDF, Dataset == uP[iP])
    cat(paste0("Downloading ", thisp$Dataset, "\n"))

    thisinfo <- tryCatch({
      info(datasetid = thisp$Dataset_ID, url = thisp$URL)
    }, error = function(e) {
      cat("An error occurred: ", conditionMessage(e), "\n")
      cat("Unable to find the ERDDAP data by its ID...\n")
      return(NULL)
    })

    if (is.null(thisinfo)) {
      cat("Dataset not found on ERDDAP. Skipping to the next dataset...\n")
      next
    }

    longinfo = thisinfo$alldata$longitude
    longrange = longinfo[which(longinfo$attribute_name == "actual_range"), "value"]
    longrange_num = as.numeric(unlist(strsplit(longrange, ",")))

    if (min(longrange_num) < 0){ long180or360 = 180 } else { long180or360 = 360 }

    if (thisp$Frequency == "Climatology"){

      paramoutpath = file.path(eds_path, uP[iP])
      if (!dir.exists(paramoutpath)) dir.create(paramoutpath, recursive = T)

      pib_path = file.path(paramoutpath, "Block_Level_Data")
      if (!dir.exists(pib_path)) dir.create(pib_path)

      if (length(list.files(paramoutpath, pattern = "all_units")) > 0) {
        cat(paste0("NetCDF files for ", thisp$Dataset, " already exists.\n"))
        next
      }

      num_cores <- max(1, detectCores()/2)
      registerDoParallel(cores = num_cores)

      indices <- 1:length(uI)

      foreach(ii = indices, .packages = c("rerddap")) %dopar% {

        this_unit <- subset(bbox, unit == uI[ii])

        lon_vals <- as.numeric(c(this_unit$x_min, this_unit$x_max))
        if (long180or360 == 360) {
          thislong <- ifelse(lon_vals < 0, lon_vals + 360, lon_vals)
        } else {
          thislong <- lon_vals
        }

        # Clean Numeric Vectors explicitly sorted
        thislong_vec <- sort(as.numeric(thislong))
        thislat_vec <- sort(as.numeric(c(this_unit$y_min, this_unit$y_max)))

        targetfilename <- file.path(pib_path, paste0(this_unit$unit, "_", thisp$Dataset, ".nc"))

        if (!file.exists(targetfilename)) {
          tryCatch({
            thisIP <- griddap(
              datasetx  = thisp$Dataset_ID,
              url       = thisp$URL,
              fields    = thisp$Fields,
              longitude = thislong_vec,
              latitude  = thislat_vec,  # FIXED
              fmt       = "nc",
              store     = disk(path = pib_path)
            )
            file.rename(thisIP$summary$filename, targetfilename)
            cat(paste0("Spatial unit ", this_unit$unit, " written to disk. ", ii, " of ", length(uI), "\n"))
          }, error = function(e) {
            cat(paste0("Error in unit ", this_unit$unit, ": ", e$message, ". Skipping...\n"))
          })
        }
      }

      stopImplicitCluster()

      cat(paste0("Completed ", thisp$Dataset, ". Check: ", length(uI), " units data present.\n"))

      MergedRasterFileName = file.path(paramoutpath, paste0(uP[iP], "_all_units.nc"))

      if (!file.exists(MergedRasterFileName)){
        cat(paste0("Completed each unit. Merging netCDF files now...\n"))
        ILnc = list.files(pib_path, pattern = "*.nc", full.names = T)

        if (length(ILnc) > 0) {
          raster_list <- lapply(ILnc, function(nc_file) {
            r <- raster::raster(nc_file)
            crs(r) <- "+proj=longlat +datum=WGS84"
            return(r)
          })

          r <- do.call(merge, c(raster_list, fun = mean))

          raster::writeRaster(x = r,
                              filename = MergedRasterFileName,
                              format = "CDF",
                              overwrite = T)
        }
      }
      cat(paste0("Completed ", thisp$Dataset, ". Check: all units merged netCDF file present.\n"))
    }

    if (thisp$Frequency != "Climatology"){

      paramoutpath = file.path(eds_path, uP[iP])
      if (!dir.exists(paramoutpath)) dir.create(paramoutpath, recursive = T)

      pib_path = file.path(paramoutpath, "Block_Level_Data")
      if (!dir.exists(pib_path)) dir.create(pib_path)

      for (ii in 1:length(uI)){

        this_unit <- subset(bbox, unit == uI[ii])

        targetfilename <- file.path(pib_path, paste0(this_unit$unit, "_", thisp$Dataset_ID, ".nc"))

        if (file.exists(targetfilename)) {
          cat(paste0("NetCDF file for ", this_unit$unit, " already exists. Skipping...\n"))
          next
        }

        lon_vals <- as.numeric(c(this_unit$x_min, this_unit$x_max))
        if (long180or360 == 360) {
          thislong <- ifelse(lon_vals < 0, lon_vals + 360, lon_vals)
        } else {
          thislong <- lon_vals
        }

        # Vector format for standard download
        thislong_vec <- sort(as.numeric(thislong))
        thislat_vec <- sort(as.numeric(c(this_unit$y_min, this_unit$y_max)))

        cat(paste0("Processing unit ", this_unit$unit, " (", ii, " of ", length(uI), ")... "))

        testIP <- tryCatch({
          griddap(
            datasetx  = thisp$Dataset_ID,
            url       = thisp$URL,
            fields    = trimws(thisp$Fields),
            time      = c('last', 'last'),
            longitude = thislong_vec,
            latitude  = thislat_vec, # FIXED
            store     = memory()
          )
        }, error = function(e) {
          cat(paste0("\nGRIDDAP ERROR: ", e$message, "\n"))
          return(NULL)
        })

        if (is.null(testIP)) {
          cat("One or both longitude values outside data range. Skipping this unit...\n")
          next
        }

        NCG = thisinfo$alldata$NC_GLOBAL

        if (thisp$Frequency == "Monthly") timestep = 30.42
        if (thisp$Frequency == "14day") timestep = 14
        if (thisp$Frequency == "8day") timestep = 8
        if (thisp$Frequency == "5day") timestep = 5
        if (thisp$Frequency == "Weekly") timestep = 7
        if (thisp$Frequency == "Daily") timestep = 1
        if (thisp$Frequency == "1hour") timestep = round(1/6, 2)

        ts_start = NCG %>%
          dplyr::filter(attribute_name == "time_coverage_start") %>%
          dplyr::select(value) %>%
          dplyr::mutate(value = substring(value, 1, 10)) %>%
          lubridate::parse_date_time(orders = c("ymd", "mdy", "dmy"), tz = "UTC")

        ts_start <- ts_start + lubridate::days(1)

        ts_end = NCG %>%
          dplyr::filter(attribute_name == "time_coverage_end") %>%
          dplyr::select(value) %>%
          dplyr::mutate(value = substring(value, 1, 10)) %>%
          lubridate::parse_date_time(orders = c("ymd", "mdy", "dmy"), tz = "UTC")

        ts_end <- format(as.Date(ts_end) - lubridate::days(as.numeric(format(ts_end, "%d"))), "%Y-%m-%d UTC")
        ts_end = ts_end %>% lubridate::parse_date_time(orders = c("ymd", "mdy", "dmy"), tz = "UTC")

        total_timesteps = round((ts_end - ts_start)/timestep, 0) %>% as.numeric()

        block_step <- timestep * 1000

        Nblocks <- ceiling(as.numeric((ts_end - ts_start) / block_step))

        if (Nblocks %% 1 != 0) {
          Nblocks <- ceiling(Nblocks)
        }

        block_step <- (ts_end - ts_start) / Nblocks

        num_cores <- min(Nblocks, max(1, detectCores()/2))
        registerDoParallel(cores = num_cores)

        indices <- 1:Nblocks

        # Parallel block execution
        foreach(blockI = indices, .packages = c("lubridate", "rerddap")) %dopar% {

          this_start = floor_date(ts_start + (blockI-1) * block_step, unit = "day")
          if (blockI > 1) this_start = this_start + days(1)

          this_end = floor_date(ts_start + ((blockI) * block_step) - timestep, unit = "day")
          if (this_end > ts_end | blockI == Nblocks) this_end = ts_end

          targetfilename = file.path(pib_path, paste0(this_unit$unit, "_", thisp$Dataset, "_", this_start, "_", this_end, ".nc"))

          this_start_val = as.character(floor_date(ts_start + (blockI-1) * block_step, unit = "day"))
          this_end_val   = as.character(floor_date(ts_start + ((blockI) * block_step) - timestep, unit = "day"))

          if (!file.exists(targetfilename)) {

            continue_loop = TRUE
            retry_count = 0

            while(continue_loop && retry_count < 3) {
              tryCatch({
                thisIP = griddap(
                  datasetx  = thisp$Dataset_ID,
                  url       = thisp$URL,
                  fields    = c(thisp$Fields),
                  time      = c(this_start_val, this_end_val),
                  longitude = thislong_vec, # FIXED
                  latitude  = thislat_vec,  # FIXED
                  fmt       = "nc",
                  store     = disk(path = pib_path),
                  read      = TRUE
                )
                continue_loop = FALSE
                file.rename(thisIP$summary$filename, targetfilename)

              }, error = function(e) {
                retry_count <- retry_count + 1
                Sys.sleep(3)
              })
            }
          }
        }

        stopImplicitCluster()

        cat(paste0("Completed ", thisp$Dataset, ". Check: ", Nblocks, " blocks present for ", this_unit$unit, ". ", ii, " of ", length(uI), " units.\n"))

        pi_path = file.path(paramoutpath, "Unit_Level_Data")
        if (!dir.exists(pi_path)) dir.create(pi_path)

        outfile = file.path(pi_path, paste0(this_unit$unit, "_", thisp$Dataset, "_", floor_date(ts_start, unit = "day"), "_", floor_date(ts_end, unit = "day"), ".nc"))

        if(!file.exists(outfile)){
          AllBlock = list.files(pib_path, full.names = T, pattern = this_unit$unit)

          if(length(AllBlock) > 0) {
            out = merge_times_nc_list(infilenamelist = as.list(AllBlock),
                                      variable_name = thisp$Fields,
                                      outfilename = outfile)
          }
        }
        cat(paste0("Completed ", thisp$Dataset, ". Merged .nc time-series present for ", this_unit$unit, ".\n"))
      }
    }
  }

  ### extract climatology ###
  cat("\nEDS is currently extracting climatological data...\n")

  df_sp = df[,c("lon", "lat", "unit")]
  df_sp$lon = ifelse(df_sp$lon < 0, df_sp$lon + 360, df_sp$lon)
  df_sp = df_sp[complete.cases(df_sp[,c("lon", "lat")]), ]
  coordinates(df_sp) = ~lon + lat

  rasterlist = list.files(eds_path, recursive = T, pattern = "_all_units.nc", full.names = T)

  ParamDF <- read_csv(file.path(eds_path, "eds_parameters.csv"), show_col_types = FALSE) %>% dplyr::filter(Download == "YES")
  uP <- ParamDF %>% dplyr::filter(Frequency == "Climatology") %>% dplyr::pull(Dataset) %>% unique()

  rasterlist <- rasterlist[grep(paste(uP, collapse = "|"), rasterlist)]
  shift = raster::shift

  for (raster_i in 1:length(rasterlist)){

    rasname_full = rasterlist[raster_i]
    rasname_sp = strsplit(rasname_full, "/")[[1]]
    rasname = rasname_sp[length(rasname_sp)]
    rasname = gsub(rasname, pattern = "-", replacement = ".")
    rasname = gsub(rasname, pattern = "_all_units.nc", replacement = "")

    this_r = raster(rasterlist[raster_i])

    if (grepl("Bathymetry", rasname)) {
      this_r[this_r > 0] <- NA
    }

    if (this_r@extent@xmin < 0) this_r = shift(rotate(shift(this_r, 180)), 180)

    crs(df_sp) = crs(this_r)
    cat(paste0("Step ", raster_i, " of ", length(rasterlist), ": ", rasname, "\n"))

    this_Ex = expanding_extract(this_r, df_sp, Dists = seq(0, 100, 5))
    eval(parse(text = paste0("df_sp$", rasname, " = this_Ex$values")))
    cat(paste0("Step ", raster_i, " of ", length(rasterlist), ": Extraction Complete.\n"))

  }

  df_clim = as.data.frame(df_sp)
  df_clim$lon = ifelse(df_clim$lon > 180, df_clim$lon - 360, df_clim$lon)
  df_clim = df_clim %>% dplyr::mutate_all(~ifelse(is.nan(.), NA, .))

  ### extract timeseries ###
  cat("\nEDS is currently extracting timeseries data...\n")

  df$date_r = lubridate::mdy(time)

  paramdir = eds_path
  parameters = list.files(path = paramdir, full.names = F)

  Parameter_Table = read.csv(file.path(eds_path, "eds_parameters.csv")) %>% dplyr::filter(Download == "YES" & Frequency != "Climatology")
  parameters = unique(Parameter_Table$Dataset)

  PT = points_in_polys(df$lon, df$lat, bbox)
  df$DATA_UNIT = df$unit

  cat(paste("Dropping", length(which(df$DATA_UNIT == "NONE_ASSIGNED")),
            "points out of", nrow(df),
            "entered points due to being outside of the geographic scope.\n"))

  df = subset(df, DATA_UNIT != "NONE_ASSIGNED")
  df = df[!(is.na(df$lat) & is.na(df$lon)) & !is.na(df$date_r), ]
  df = df[!duplicated(df[c("lat", "lon", "date_r")]), ]

  unique_units = sort(unique(df$DATA_UNIT))

  for(parameter_i in 1:length(parameters)){

    param.name = parameters[parameter_i]
    this_param_i = which(Parameter_Table$Dataset == param.name)
    godir = file.path(paramdir, param.name, "Unit_Level_Data")
    paramsum = unlist(strsplit(as.vector(Parameter_Table[this_param_i, "Summaries"]), ";"))

    for(unit_i in 1:length(unique_units)){

      ncfile = list.files(godir, pattern = paste0(unique_units[unit_i], "_"), full.names = T)

      if(length(ncfile) == 0){
        cat(paste0("Skipping ", unique_units[unit_i], " as there is no available data.\n"))
        next
      }

      nc_p = nc_open(ncfile)
      rawvar = ncvar_get(nc = nc_p, varid = as.vector(Parameter_Table$Fields[this_param_i]))
      lon = ncvar_get(nc = nc_p, varid = "longitude")
      lat = ncvar_get(nc = nc_p, varid = "latitude")
      rawt = ncvar_get(nc = nc_p, varid = "time")
      nc_close(nc_p)

      t = lubridate::as_date(lubridate::as_datetime(as.numeric(rawt), origin = lubridate::ymd("1970/1/1")))

      df_i = which(df$DATA_UNIT == unique_units[unit_i])

      if (max(lon) > 180) {
        xyt_df_i = df %>% dplyr::mutate(lon = ifelse(df$lon < 0, df$lon + 360, df$lon))
      } else{
        xyt_df_i = df
      }

      ijk = xyt_to_ijk(xyt_df = as.data.frame(xyt_df_i[df_i,c("lon","lat","date_r")]),
                       x_grid = lon, y_grid = lat, t_grid = t,
                       lon = lon, lat = lat, t = t)

      droprows = which(is.na(ijk), arr.ind = T)[,1]

      if (length(droprows) > 0){
        ijk = ijk[-droprows,]
        df_i = df_i[-droprows]
      }

      if (length(df_i) == 0) {
        cat(paste0("No timeseries data found for ", unique_units[unit_i], ". Consider trying higher-resolution data. Skipping to the next...\n"))
        next()
      } else {

        if (length(dim(rawvar)) == 1) {
          dim1 <- 1
          dim2 <- 1
          dim3 <- length(rawvar) / (dim1 * dim2)
          rawvar <- array(rawvar, dim = c(dim1, dim2, dim3))
        }

        naP_xy = plyr::aaply(rawvar, c(1,2), na_stack_count)/dim(rawvar)[3]
        i_masked = which(naP_xy[cbind(ijk$x_i, ijk$y_j)] > 0.8)
        cnt = 1

        for(i_infill in i_masked){
          pNA = naP_xy[cbind(ijk$x_i[i_infill], ijk$y_j[i_infill])]
          ij_ex = 1

          while(pNA > 0.8 & ij_ex < 3){
            ij_vec = -ij_ex:ij_ex
            max_x = dim(rawvar)[1]
            max_y = dim(rawvar)[2]

            ts_x = ijk$x_i[i_infill]+ij_vec
            ts_y = ijk$y_j[i_infill]+ij_vec

            ts_x = subset(ts_x, ts_x <= max_x & ts_x > 0)
            ts_y = subset(ts_y, ts_y <= max_y & ts_y > 0)

            ts = plyr::aaply(rawvar[ts_x, ts_y,], c(3), mean, na.rm = T)
            pNA = length(which(is.na(ts)))/length(ts)

            if(pNA < 0.9){
              rawvar[ijk$x_i[i_infill], ijk$y_j[i_infill],] = ts
              naP_xy = plyr::aaply(rawvar,c(1,2), na_stack_count)/dim(rawvar)[3]
              cat(paste("In-fill complete", cnt, "of", length(i_masked), ". Pixel +/-", ij_ex), "\n")
            }
            ij_ex = ij_ex + 1
          }
          cnt = cnt + 1
        }

        if (Parameter_Table$Frequency[this_param_i] == "Monthly") tstep = 30.42
        if (Parameter_Table$Frequency[this_param_i] == "14day") tstep = 14
        if (Parameter_Table$Frequency[this_param_i] == "8day") tstep = 8
        if (Parameter_Table$Frequency[this_param_i] == "5day") tstep = 5
        if (Parameter_Table$Frequency[this_param_i] == "Weekly") tstep = 7
        if (Parameter_Table$Frequency[this_param_i] == "Daily") tstep = 1

        ijk$t_01dy = ijk$t_k - (1/tstep-1)
        ijk$t_01wk = ijk$t_k - (7/tstep-1)
        ijk$t_01mo = round(ijk$t_k - (30.42/tstep-1))
        ijk$t_03mo = round(ijk$t_k - (91.25/tstep-1))
        ijk$t_06mo = round(ijk$t_k - (182.5/tstep-1))
        ijk$t_01yr = round(ijk$t_k - (1*365.25/tstep-1))
        ijk$t_03yr = round(ijk$t_k - (3*365.25/tstep-1))
        ijk$t_05yr = round(ijk$t_k - (5*365.25/tstep-1))
        ijk$t_10yr = round(ijk$t_k - (10*365.25/tstep-1))

        ijk[,c("t_k","t_01dy","t_01wk","t_01mo","t_03mo","t_06mo","t_01yr","t_03yr","t_05yr","t_10yr")][which(ijk[,c("t_k","t_01dy","t_01wk","t_01mo","t_03mo","t_06mo","t_01yr","t_03yr","t_05yr","t_10yr")] < 1, arr.ind = T)] = 1

        for(sum_i in 1:length(paramsum)){

          paramsum.name = paste0(paramsum[sum_i], "_", param.name)

          if(!paramsum.name %in% substr(names(df), 1, nchar(paramsum.name))){
            eval(parse(text = paste0("df$",paramsum.name,"_DY01=-9991")))
            eval(parse(text = paste0("df$",paramsum.name,"_WK01=-9991")))
            eval(parse(text = paste0("df$",paramsum.name,"_MO01=-9991")))
            eval(parse(text = paste0("df$",paramsum.name,"_MO03=-9991")))
            eval(parse(text = paste0("df$",paramsum.name,"_MO06=-9991")))
            eval(parse(text = paste0("df$",paramsum.name,"_YR01=-9991")))
            eval(parse(text = paste0("df$",paramsum.name,"_YR03=-9991")))
            eval(parse(text = paste0("df$",paramsum.name,"_YR05=-9991")))
            eval(parse(text = paste0("df$",paramsum.name,"_YR10=-9991")))
            eval(parse(text = paste0("df$",paramsum.name,"_YR10YR01=-9991")))
            eval(parse(text = paste0("df$",paramsum.name,"_ALLB4=-9991")))
          }

          for(sumpt_i in 1:length(df_i)){
            ts_01dy = rawvar[ijk$x_i[sumpt_i], ijk$y_j[sumpt_i], ijk$t_01dy[sumpt_i]:ijk$t_k[sumpt_i]]
            ts_01wk = rawvar[ijk$x_i[sumpt_i], ijk$y_j[sumpt_i], ijk$t_01wk[sumpt_i]:ijk$t_k[sumpt_i]]
            ts_01mo = rawvar[ijk$x_i[sumpt_i], ijk$y_j[sumpt_i], ijk$t_01mo[sumpt_i]:ijk$t_k[sumpt_i]]
            ts_03mo = rawvar[ijk$x_i[sumpt_i], ijk$y_j[sumpt_i], ijk$t_03mo[sumpt_i]:ijk$t_k[sumpt_i]]
            ts_06mo = rawvar[ijk$x_i[sumpt_i], ijk$y_j[sumpt_i], ijk$t_06mo[sumpt_i]:ijk$t_k[sumpt_i]]
            ts_01yr = rawvar[ijk$x_i[sumpt_i], ijk$y_j[sumpt_i], ijk$t_01yr[sumpt_i]:ijk$t_k[sumpt_i]]
            ts_03yr = rawvar[ijk$x_i[sumpt_i], ijk$y_j[sumpt_i], ijk$t_03yr[sumpt_i]:ijk$t_k[sumpt_i]]
            ts_05yr = rawvar[ijk$x_i[sumpt_i], ijk$y_j[sumpt_i], ijk$t_05yr[sumpt_i]:ijk$t_k[sumpt_i]]
            ts_10yr = rawvar[ijk$x_i[sumpt_i], ijk$y_j[sumpt_i], ijk$t_10yr[sumpt_i]:ijk$t_k[sumpt_i]]
            ts_10yr01yr = rawvar[ijk$x_i[sumpt_i], ijk$y_j[sumpt_i], ijk$t_10yr[sumpt_i]:ijk$t_01yr[sumpt_i]]
            ts_ALLB4 = rawvar[ijk$x_i[sumpt_i], ijk$y_j[sumpt_i], 1:ijk$t_k[sumpt_i]]

            t_01dy = t[ijk$t_01dy[sumpt_i]:ijk$t_k[sumpt_i]]
            t_01wk = t[ijk$t_01wk[sumpt_i]:ijk$t_k[sumpt_i]]
            t_01mo = t[ijk$t_01mo[sumpt_i]:ijk$t_k[sumpt_i]]
            t_03mo = t[ijk$t_03mo[sumpt_i]:ijk$t_k[sumpt_i]]
            t_06mo = t[ijk$t_06mo[sumpt_i]:ijk$t_k[sumpt_i]]
            t_01yr = t[ijk$t_01yr[sumpt_i]:ijk$t_k[sumpt_i]]
            t_03yr = t[ijk$t_03yr[sumpt_i]:ijk$t_k[sumpt_i]]
            t_05yr = t[ijk$t_05yr[sumpt_i]:ijk$t_k[sumpt_i]]
            t_10yr = t[ijk$t_10yr[sumpt_i]:ijk$t_k[sumpt_i]]
            t_10yr01yr = t[ijk$t_10yr[sumpt_i]:ijk$t_01yr[sumpt_i]]
            t_ALLB4 = t[1:ijk$t_k[sumpt_i]]

            if(paramsum[sum_i] %in% c("mean", "q05", "q95","sd")){
              eval(parse(text = paste0("df$", paramsum.name, "_DY01[df_i[sumpt_i]] = ", paramsum[sum_i], "(x = ts_01dy, na.rm = T)")))
              eval(parse(text = paste0("df$", paramsum.name, "_WK01[df_i[sumpt_i]] = ", paramsum[sum_i], "(x = ts_01wk, na.rm = T)")))
              eval(parse(text = paste0("df$", paramsum.name, "_MO01[df_i[sumpt_i]] = ", paramsum[sum_i], "(x = ts_01mo, na.rm = T)")))
              eval(parse(text = paste0("df$", paramsum.name, "_MO03[df_i[sumpt_i]] = ", paramsum[sum_i], "(x = ts_03mo, na.rm = T)")))
              eval(parse(text = paste0("df$", paramsum.name, "_MO06[df_i[sumpt_i]] = ", paramsum[sum_i], "(x = ts_06mo, na.rm = T)")))
              eval(parse(text = paste0("df$", paramsum.name, "_YR01[df_i[sumpt_i]] = ", paramsum[sum_i], "(x = ts_01yr, na.rm = T)")))
              eval(parse(text = paste0("df$", paramsum.name, "_YR03[df_i[sumpt_i]] = ", paramsum[sum_i], "(x = ts_03yr, na.rm = T)")))
              eval(parse(text = paste0("df$", paramsum.name, "_YR05[df_i[sumpt_i]] = ", paramsum[sum_i], "(x = ts_05yr, na.rm = T)")))
              eval(parse(text = paste0("df$", paramsum.name, "_YR10[df_i[sumpt_i]] = ", paramsum[sum_i], "(x = ts_10yr, na.rm = T)")))
              eval(parse(text = paste0("df$", paramsum.name, "_YR10YR01[df_i[sumpt_i]] = ", paramsum[sum_i], "(x = ts_10yr01yr, na.rm = T)")))
              eval(parse(text = paste0("df$", paramsum.name, "_ALLB4[df_i[sumpt_i]] = ", paramsum[sum_i], "(x = ts_ALLB4, na.rm = T)")))
            }else{
              eval(parse(text = paste0("df$", paramsum.name, "_DY01[df_i[sumpt_i]] = ", paramsum[sum_i], "(x = ts_01dy, na.rm = T, t = t_01dy)")))
              eval(parse(text = paste0("df$", paramsum.name, "_WK01[df_i[sumpt_i]] = ", paramsum[sum_i], "(x = ts_01wk, na.rm = T, t = t_01wk)")))
              eval(parse(text = paste0("df$", paramsum.name, "_MO01[df_i[sumpt_i]] = ", paramsum[sum_i], "(x = ts_01mo, na.rm = T, t = t_01mo)")))
              eval(parse(text = paste0("df$", paramsum.name, "_MO03[df_i[sumpt_i]] = ", paramsum[sum_i], "(x = ts_03mo, na.rm = T, t = t_03mo)")))
              eval(parse(text = paste0("df$", paramsum.name, "_MO06[df_i[sumpt_i]] = ", paramsum[sum_i], "(x = ts_06mo, na.rm = T, t = t_06mo)")))
              eval(parse(text = paste0("df$", paramsum.name, "_YR01[df_i[sumpt_i]] = ", paramsum[sum_i], "(x = ts_01yr, na.rm = T, t = t_01yr)")))
              eval(parse(text = paste0("df$", paramsum.name, "_YR03[df_i[sumpt_i]] = ", paramsum[sum_i], "(x = ts_03yr, na.rm = T, t = t_03yr)")))
              eval(parse(text = paste0("df$", paramsum.name, "_YR05[df_i[sumpt_i]] = ", paramsum[sum_i], "(x = ts_05yr, na.rm = T, t = t_05yr)")))
              eval(parse(text = paste0("df$", paramsum.name, "_YR10[df_i[sumpt_i]] = ", paramsum[sum_i], "(x = ts_10yr, na.rm = T, t = t_10yr)")))
              eval(parse(text = paste0("df$", paramsum.name, "_YR10YR01[df_i[sumpt_i]] = ", paramsum[sum_i], "(x = ts_10yr01yr, na.rm = T, t = t_10yr01yr)")))
              eval(parse(text = paste0("df$", paramsum.name, "_ALLB4[df_i[sumpt_i]] = ", paramsum[sum_i], "(x = ts_ALLB4, na.rm = T, t = t_ALLB4)")))
            }
          }

          cat(paste("Processing", sumpt_i, "points for unit", unique_units[unit_i], "out of", length(unique_units), "units.", paramsum.name, "completed.\n"))
          save(df, file = file.path(eds_path, "eds_time.Rdata"))
        }
      }
    }
  }

  df[df == -9991] <- NA

  names(df) <- gsub("Daily", "daily", names(df))
  names(df) <- gsub("Weekly", "weekly", names(df))
  names(df) <- gsub("_8Day_", "_8days_", names(df))
  names(df) <- gsub("Monthly", "monthly", names(df))
  names(df) <- gsub("DHW.", "", names(df))
  names(df) <- gsub("Major_", "major_", names(df))
  names(df) <- gsub("Np10y_", "", names(df))
  names(df) <- gsub("MeanMax_", "mean_maximum_", names(df))
  names(df) <- gsub("CI95Max_", "ci95_maximum_", names(df))
  names(df) <- gsub("MeanDur_", "mean_durnal_", names(df))
  names(df) <- gsub("MaxMax_", "maximum_maximum_", names(df))
  names(df) <- gsub("_DY01", "_01dy", names(df))
  names(df) <- gsub("_WK01", "_01wk", names(df))
  names(df) <- gsub("_MO01", "_01mo", names(df))
  names(df) <- gsub("_MO03", "_03mo", names(df))
  names(df) <- gsub("_MO06", "_06mo", names(df))
  names(df) <- gsub("_YR10YR01", "_10yr_01yr", names(df))
  names(df) <- gsub("_YR01", "_01yr", names(df))
  names(df) <- gsub("_YR03", "_03yr", names(df))
  names(df) <- gsub("_YR05", "_05yr", names(df))
  names(df) <- gsub("_YR10", "_10yr", names(df))
  names(df) <- gsub("_ALLB4", "_all_before", names(df))

  df = df %>% dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ round(., 3)))

  df_time = df

  save(df_clim, file = file.path(eds_path, "eds_clim.RData"))
  save(df_time, file = file.path(eds_path, "eds_time.RData"))

  options(warn = 0)

  cat("\nEDS run has been successfully completed.\n")
  cat("climatology output saved to:", file.path(eds_path, "eds_clim.RData"), "\n")
  cat("temporally summarized output saved to:", file.path(eds_path, "eds_time.RData"), "\n")

}
