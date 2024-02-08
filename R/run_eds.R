#' importFrom(stats, complete.cases, filter, quantile, var)
#' importFrom(utils, head, install.packages, read.csv, tail)
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
    dim(dat_new)
    eval(parse(text = paste0("var = file1$var[variable_name]$", variable_name)))

    #Get dimensions - Keep in Same Orientation as first input file
    #Lon
    DEC_LON=file1$dim$longitude$vals[1]>file1$dim$longitude$vals[2]
    Dx = ncdim_def("longitude", "degrees_east", sort(unique(c(file1$dim$longitude$vals, file2$dim$longitude$vals)), decreasing = DEC_LON))

    #Lat
    DEC_LAT=file1$dim$latitude$vals[1]>file1$dim$latitude$vals[2]
    Dy = ncdim_def("latitude", "degrees_north", sort(unique(c(file1$dim$latitude$vals, file2$dim$latitude$vals)), decreasing = DEC_LAT))

    #Time
    DEC_T=file1$dim$time$vals[1]>file1$dim$time$vals[2]
    if (is.na(DEC_T) == T) DEC_T = FALSE
    Dt = ncdim_def("time", "seconds since 1970-01-01T00:00:00Z", sort(unique(c(file1$dim$time$vals, file2$dim$time$vals)), decreasing = DEC_T))

    # Create a new file
    file_new  <- nc_create(
      filename = outfilename,
      # We need to define the variables here
      vars = ncvar_def(
        name = variable_name,
        units = var$units,
        dim = list(Dx, Dy, Dt)))

    #copy over all the metadata, but not the dimensions...
    #  modifyNcdfCopyMetadata(file.con.orig = file1, file.con.copy = file_new, glob.atts = T, dimensions = F)

    # And write to it
    ncvar_put(
      nc = file_new,
      varid = variable_name,
      vals = dat_new)

    # Finally, close the files
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
      file.copy(from = infilenamelist[[1]], to = outfilename, )
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

    #Make sure data are apples to apples class-wise
    xyt_df = data.frame(x = as.numeric(xyt_df[, 1]),
                        y = as.numeric(xyt_df[, 2]),
                        t = as.Date(xyt_df[, 3]))

    xyt_df$x[xyt_df$x > 180] = xyt_df$x[xyt_df$x > 180] - 360

    x_grid = as.numeric(x_grid)
    x_grid[x_grid>180] = x_grid[x_grid > 180] - 360
    y_grid = as.numeric(y_grid)
    t_grid = as.Date(t_grid)

    #Sizes
    n_pts = nrow(xyt_df)
    n_xg = length(x_grid)
    n_yg = length(y_grid)
    n_tg = length(t_grid)

    # set x, y, t buffers
    x_buffer = mean(diff(lon))/2; x_buffer
    y_buffer = mean(diff(lat))/2; y_buffer
    t_buffer = mean(diff(t))/2; t_buffer

    #X Match
    x_gMat = outer(rep(1, n_pts), x_grid)
    x_pMat = outer(xyt_df$x, rep(1, n_xg))
    x_i = apply(abs(x_pMat-x_gMat), 1, which.min)
    #check for out of bound points
    oob_x = which(xyt_df$x < (min(x_grid, na.rm = T) - x_buffer) | xyt_df$x > (max(x_grid, na.rm = T) + x_buffer))
    x_i[oob_x] = NA

    #Y Match
    y_gMat = outer(rep(1, n_pts), y_grid)
    y_pMat = outer(xyt_df$y, rep(1, n_yg))
    y_j = apply(abs(y_pMat-y_gMat), 1, which.min)
    #check for out of bound points
    oob_y = which(xyt_df$y < (min(y_grid, na.rm = T) - y_buffer) | xyt_df$y > (max(y_grid, na.rm = T) + y_buffer))
    y_j[oob_y] = NA

    #T Match
    t_gMat = outer(rep(1, n_pts), t_grid)
    t_pMat = outer(xyt_df$t, rep(1, n_tg))
    t_k = apply(abs(t_pMat-t_gMat), 1, which.min)
    #check for out of bound points
    oob_t = which(xyt_df$t < (min(t_grid, na.rm = T) - t_buffer) | xyt_df$t > (max(t_grid, na.rm = T) + t_buffer))
    t_k[oob_t] = NA

    #Build output df
    ijk = data.frame(x_i, y_j, t_k)

    #return...
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

  # load("data/SURVEY_MASTER.RData")
  # df <- subset(df, region == "MHI")
  # # df = df[sample(1:nrow(df), 300, replace = F),]
  # lon = df$lon
  # lat = df$lat
  # unit = df$island
  # time = df$date
  # buffer = 0.25

  if (!is.numeric(lon) || !is.numeric(lat) || !is.character(unit)) stop("Invalid input parameters. lon and lat should be numeric, unit should be character.")

  ### download ERDDAP ###
  # List of packages to load
  packages_to_load <- c("rerddap", "readr", "zoo", "ncdf4", "RNetCDF", "easyNCDF", "raster", "lubridate", "abind", "acss", "dplyr", "plyr", "doParallel", "foreach", "spatial", "data.table", "splitstackshape")

  # Function to install and load a package
  install_and_load_package <- function(package_name) {

    if (!requireNamespace(package_name, quietly = TRUE)) {
      install.packages(package_name, dependencies = TRUE)
    }

    library(package_name, character.only = TRUE)

  }

  invisible(sapply(packages_to_load, install_and_load_package))

  # if ("package:plyr" %in% search()) {
  #   unloadNamespace("plyr")
  # }

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

  ParamDF <- read_csv(file.path(paste0(eds_path, "/eds_parameters.csv")), show_col_types = FALSE) %>% filter(Download == "YES")
  uP <- ParamDF %>% pull(Dataset) %>% unique()

  df <- as.data.frame(list(lon = lon, lat = lat, unit = unit, date_r = time))
  df <- df %>% filter(!is.na(lon) & !is.na(lat) & !is.na(date_r))

  df$unit = gsub(" ", "_", unit)

  boxes <- df %>%
    group_by(unit) %>%
    dplyr::summarise(x_min = min(lon) - buffer,
                     x_max = max(lon) + buffer,
                     y_min = min(lat) - buffer,
                     y_max = max(lat) + buffer) %>%
    mutate(across(c(x_max, x_min, y_max, y_min), round, digits = 2))

  df_frame <- df %>%
    select(unit) %>%
    distinct()

  bbox <- left_join(boxes, df_frame, by = "unit") %>% mutate(unit = gsub(" ", "_", unit))

  uI <- distinct(bbox, unit)$unit

  ### download ERDDAP netCDF files ###
  cat("\nEDS is currently downloading netCDF files...\n")

  for (iP in 1:length(uP)){

    # iP = 2

    # Select dataset
    thisp = subset(ParamDF, Dataset == uP[iP])

    cat(paste0("Downloading ", thisp$Dataset, "\n"))

    # Fetch dataset information from ERDDAP
    thisinfo <- tryCatch({

      info(datasetid = thisp$Dataset_ID, url = thisp$URL)

    }, error = function(e) {

      cat("An error occurred: ", conditionMessage(e), "\n")
      cat("Unable to find the ERDDAP data by its ID...\n")
      return(NULL)  # Return NULL instead of printing and using next

    })

    # Check if dataset exists in ERDDAP
    if (is.null(thisinfo)) {

      cat("Dataset not found on ERDDAP. Skipping to the next dataset...\n")
      next

    }

    # Check how longitude is stored, determine if 180 or 360 style
    longinfo = thisinfo$alldata$longitude
    longrange = longinfo[which(longinfo$attribute_name == "actual_range"), "value"]
    longrange_num = as.numeric(unlist(strsplit(longrange, ",")))

    if (min(longrange_num) < 0){ long180or360 = 180 } else { long180or360 = 360 }

    if (thisp$Frequency == "Climatology"){

      # Find or create output directory
      paramoutpath = paste0(eds_path, uP[iP])
      if (!dir.exists(paramoutpath)) dir.create(paramoutpath, recursive = T)

      pib_path = paste0(paramoutpath,"/Block_Level_Data")
      if (!dir.exists(pib_path)) dir.create(pib_path)

      # Check if any files match the pattern
      if (length(list.files(paste0(paramoutpath, "/"), pattern = "all_units")) > 0) {

        cat(paste0("NetCDF files for ", thisp$Dataset, " already exists.\n"))
        next

      }

      # Set the number of cores to use
      num_cores <- detectCores()/2  # Change this to the desired number of cores

      # Initialize parallel backend
      registerDoParallel(cores = num_cores)

      # Create a list of indices for parallel processing
      indices <- 1:length(uI)

      # Loop through each unit
      foreach(ii = indices, .packages = c("rerddap")) %dopar% {

        # ii = 1

        # Select unit's data
        this_unit = subset(bbox, unit == uI[ii])

        # Get appropriate Longitude span
        if (long180or360 == 360){

          this_unit[, c("x_min", "x_max")] = ifelse(this_unit[, c("x_min", "x_max")] < 0,
                                                    this_unit[, c("x_min", "x_max")] + 360,
                                                    this_unit[, c("x_min", "x_max")])

          thislong = this_unit[, c("x_min","x_max")]

        } else {

          thislong = this_unit[, c("x_min","x_max")]

        }

        # Make a griddap call to pull data from server
        targetfilename = paste0(pib_path, "/",
                                this_unit$unit, "_",
                                thisp$Dataset, ".nc")

        if (!file.exists(targetfilename)){

          tryCatch({

            thisIP = griddap(datasetx = thisp$Dataset_ID,
                             url = thisp$URL,
                             fields = c(thisp$Fields),
                             longitude = thislong,
                             latitude = this_unit[, c("y_min","y_max")],
                             fmt = "nc",
                             store = disk(path = pib_path))

            ncstatus = file.rename(thisIP$summary$filename,
                                   targetfilename)

            cat(paste0("Spatial unit", this_unit$unit,
                       " written to disk. ",
                       ii,
                       ' of ',
                       length(uI), "\n"))

          }, error = function(e) {

            cat(paste0(e$message, ". Skipping ", this_unit$unit, "...\n"))

          })
        }

      }

      # Stop the parallel backend
      stopImplicitCluster()

      cat(paste0("Completed ", thisp$Dataset, ". Check: ", length(uI), " units data present.\n"))

      MergedRasterFileName = paste0(paramoutpath, "/", uP[iP], "_all_units.nc")

      if (!file.exists(MergedRasterFileName)){

        # Re-read each file, merge into single ncdf, output...
        cat(paste0("Completed each unit. Merging netCDF files now...\n"))

        # Get list of island-level ncdfs
        ILnc = list.files(pib_path,pattern = "*.nc", full.names = T)

        # skip if previous step fails to produce summary .nc files
        if (length(ILnc) == 0) next

        # Set number of parallel workers
        num_workers <- detectCores()/2  # Adjust the number of workers based on your system's capacity

        # Register parallel backend
        cl <- makeCluster(num_workers)
        registerDoParallel(cl)

        # Load and merge raster files using parallel processing
        r <- foreach(i = 1:length(ILnc), .combine = merge, .packages = c("raster")) %dopar% {

          r <- raster(ILnc[i])
          crs(r) <- "+proj=longlat +datum=WGS84"
          r

        }

        # Stop parallel processing and clean up
        stopCluster(cl)
        registerDoSEQ()

        r = mean(r)

        # Write Raster as nc file
        raster::writeRaster(x = readAll(r),
                            filename = MergedRasterFileName,
                            format = "CDF",
                            overwrite = T)

      } # End of Merged Raster export

      cat(paste0("Completed ", thisp$Dataset, ". Check: all units merged netCDF file present.\n"))

    }

    if (thisp$Frequency != "Climatology"){

      # Find or create output directory
      paramoutpath = paste0(eds_path,"/", uP[iP])
      if (!dir.exists(paramoutpath)) dir.create(paramoutpath, recursive = T)

      pib_path = paste0(paramoutpath, "/Block_Level_Data")
      if (!dir.exists(pib_path)) dir.create(pib_path)

      # Loop through each unit
      for (ii in 1:length(uI)){

        # ii = 1

        # Select unit's bounding box data
        this_unit = subset(bbox, unit == uI[ii])

        # Check if any files match the pattern
        if (length(list.files(paste0(paramoutpath, "/Unit_Level_Data/"), pattern = this_unit$unit)) > 0) {

          cat(paste0("NetCDF files for ", this_unit$unit, " already exists.\n"))
          next

        }

        # Get appropriate Longitude span
        if (long180or360 == 360){

          this_unit[, c("x_min", "x_max")] = ifelse(this_unit[, c("x_min", "x_max")] < 0,
                                                    this_unit[, c("x_min", "x_max")] + 360,
                                                    this_unit[, c("x_min", "x_max")])

          thislong = this_unit[, c("x_min","x_max")]

        } else {

          thislong = this_unit[, c("x_min","x_max")]

        }

        # call griddap() to pull test data from server
        # skip if an unit is outside of data range
        testIP <- tryCatch({

          griddap(datasetx = thisp$Dataset_ID,
                  url = thisp$URL,
                  fields = c(trim(thisp$Fields)),
                  time = c('last', 'last'),
                  longitude = thislong,
                  latitude = this_unit[, c("y_min", "y_max")],
                  store = memory())

        }, error = function(e) {

          cat("GRIDDAP ERROR...\n")
          return(NULL)

        })

        if (is.null(testIP)) {

          cat("One or both longitude values outside data range. Skipping this unit...\n")
          next

        }

        # Get Metadata
        NCG = thisinfo$alldata$NC_GLOBAL

        if (thisp$Frequency == "Monthly") timestep = 30.42 # (60*60*24*30.42)
        if (thisp$Frequency == "14day") timestep = 14 # (60*60*24*14)
        if (thisp$Frequency == "8day") timestep = 8 # (60*60*24*8)
        if (thisp$Frequency == "5day") timestep = 5 # (60*60*24*5)
        if (thisp$Frequency == "Weekly") timestep = 7 # (60*60*24*7)
        if (thisp$Frequency == "Daily") timestep = 1 # (60*60*24*1)
        if (thisp$Frequency == "1hour") timestep = round(1/6, 2) # (60*60*24*1)

        # Set start and end dates for each block
        ts_start = NCG %>%
          filter(attribute_name == "time_coverage_start") %>%
          dplyr::select(value) %>%
          mutate(value = substring(value, 1, 10)) %>%  # extract only date part
          parse_date_time(orders = c("ymd", "mdy", "dmy"), tz = "UTC"); ts_start

        # Avoid start dates that start at midnight)
        ts_start <- ts_start + days(1)

        ts_end = NCG %>%
          filter(attribute_name == "time_coverage_end") %>%
          dplyr::select(value) %>%
          mutate(value = substring(value, 1, 10)) %>%  # extract only date part
          parse_date_time(orders = c("ymd", "mdy", "dmy"), tz = "UTC"); ts_end

        # Point to the last day of the previous month
        ts_end <- format(as.Date(ts_end) - days(as.numeric(format(ts_end, "%d"))), "%Y-%m-%d UTC")
        ts_end = ts_end %>% parse_date_time(orders = c("ymd", "mdy", "dmy"), tz = "UTC"); ts_end

        singlestep = nrow(testIP$data); singlestep

        total_timesteps = round((ts_end - ts_start)/timestep, 0) %>% as.numeric(); total_timesteps

        # Create Nblocks as an integer where block_step is always 100 times longer than timestep
        block_step <- timestep * 1000

        # Calculate Nblocks based on the updated block_step and ensure it's an integer
        Nblocks <- ceiling(as.numeric((ts_end - ts_start) / block_step))

        if (Nblocks %% 1 != 0) {

          Nblocks <- ceiling(Nblocks)

        }

        # Recalculate block_step based on the updated Nblocks
        block_step <- (ts_end - ts_start) / Nblocks

        block_step
        Nblocks

        # Set the number of cores to use
        num_cores <- min(Nblocks, detectCores()/2)

        # Initialize parallel backend
        registerDoParallel(cores = num_cores)
        # cl <- makeCluster(num_cores)
        # registerDoParallel(cl)

        # Create a list of indices for parallel processing
        indices <- 1:Nblocks

        # Parallel loop
        foreach(blockI = indices, .packages = c("lubridate", "rerddap")) %dopar% {

          # blockI = 1

          this_start = floor_date(ts_start + (blockI-1) * block_step, unit = "day")
          if (blockI > 1) this_start = this_start + days(1)

          this_end = floor_date(ts_start + ((blockI) * block_step) - timestep, unit = "day")
          if (this_end > ts_end | blockI == Nblocks) this_end = ts_end

          targetfilename = paste0(pib_path, "/",
                                  this_unit$unit, "_",
                                  thisp$Dataset, "_",
                                  this_start, "_",
                                  this_end, ".nc")

          this_start = as.character(this_start)
          this_end = as.character(this_end)

          # If the targetfile doesn't already exist, call griddap
          if (!file.exists(targetfilename)) {

            # Define a variable to keep track of whether to continue the loop
            continue_loop = TRUE

            while(continue_loop) {

              tryCatch({

                thisIP = griddap(datasetx = thisp$Dataset_ID,
                                 url = thisp$URL,
                                 fields = c(thisp$Fields),
                                 time = c(this_start, this_end),
                                 longitude = thislong,
                                 latitude = this_unit[, c("y_min", "y_max")],
                                 fmt = "nc",
                                 store = disk(path = pib_path),
                                 read = TRUE)

                continue_loop = FALSE # If no error occurs, set continue_loop to FALSE to exit the loop

              }, error = function(e) {

                cat("GRIDDAP ERROR...\n") # The loop will continue to run until no error occurs

              })
            }

            # Once the griddap call works, rename the file
            ncstatus = file.rename(thisIP$summary$filename, targetfilename)

            cat(paste0(this_unit$unit, ", block #",
                       blockI, " of ",
                       Nblocks, " written to disk. Unit #",
                       ii, ' of ',
                       length(uI), "\n"))

          }

        }

        # Stop the parallel backend
        stopImplicitCluster()
        # stopCluster(cl)
        # registerDoSEQ()

        cat(paste0("Completed ", thisp$Dataset,
                   ". Check: ", Nblocks,
                   " blocks present for ", this_unit$unit,
                   ". ", ii,
                   " of ", length(uI),
                   " units.\n"))

        #For each unit - set up single unit folder
        pi_path = paste0(paramoutpath, "/Unit_Level_Data")
        if (!dir.exists(pi_path)) dir.create(pi_path)

        outfile = paste0(pi_path, "/",
                         this_unit$unit, "_",
                         thisp$Dataset, "_",
                         floor_date(ts_start, unit = "day"), "_",
                         floor_date(ts_end, unit = "day"), ".nc")

        if(!file.exists(outfile)){

          AllBlock = list.files(pib_path,
                                full.names = T,
                                pattern = this_unit$unit)

          out = merge_times_nc_list(infilenamelist = as.list(AllBlock),
                                    variable_name = thisp$Fields,
                                    outfilename = outfile)
        }

        cat(paste0("Completed ",
                   thisp$Dataset,
                   ". Merged .nc time-series present for ",
                   this_unit$unit, ".\n"))

      }

    }

  }

  ### extract climatology ###
  cat("\nEDS is currently extracting climatological data...\n")

  df_sp = df[,c("lon", "lat", "unit")]
  df_sp$lon = ifelse(df_sp$lon < 0, df_sp$lon + 360, df_sp$lon)
  df_sp = df_sp[complete.cases(df_sp[,c("lon", "lat")]), ]
  coordinates(df_sp) = ~lon + lat

  # Create list of climatology raster files
  rasterlist = list.files(eds_path,
                          recursive = T,
                          pattern = "_all_units.nc",
                          full.names = T)

  # Read Parameter and Time Series Summary Definitions
  ParamDF <- read_csv(paste0(eds_path, "eds_parameters.csv"), show_col_types = FALSE) %>% filter(Download == "YES")

  # Extract unique datasets marked as "Climatology" from ParamDF
  uP <- ParamDF %>% filter(Frequency == "Climatology") %>% pull(Dataset) %>% unique()

  # Filter rasterlist based on the uP keywords
  rasterlist <- rasterlist[grep(paste(uP, collapse = "|"), rasterlist)]

  shift = raster::shift

  for (raster_i in 1:length(rasterlist)){

    # raster_i = 1

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
  df_clim = df_clim %>% mutate_all(~ifelse(is.nan(.), NA, .))

  ### extract timeseries ###
  cat("\nEDS is currently extracting timeseries data...\n")

  df$date_r = mdy(time)

  # Build list of target environmental variables
  paramdir = eds_path
  parameters = list.files(path = paramdir, full.names = F); parameters

  # Read EDS Parameter/Variable Table
  Parameter_Table = read.csv(paste0(eds_path, "/eds_parameters.csv")) %>% filter(Download == "YES" & Frequency != "Climatology")
  names(Parameter_Table)
  parameters = unique(Parameter_Table$Dataset)

  # locate each point in an unit bounding box
  PT = points_in_polys(df$lon,
                       df$lat,
                       bbox)

  df$DATA_UNIT = df$unit

  # Drop points outside target area and report
  cat(paste("Dropping", length(which(df$DATA_UNIT == "NONE_ASSIGNED")),
            "points out of", nrow(df),
            "entered points due to being outside of the geographic scope.\n"))

df = subset(df, DATA_UNIT != "NONE_ASSIGNED")
df = df[!(is.na(df$lat) & is.na(df$lon)) & !is.na(df$date_r), ]
df = df[!duplicated(df[c("lat", "lon", "date_r")]), ]

  # List of spatial units
  unique_units = sort(unique(df$DATA_UNIT)); unique_units

  # Extract Time Series Variables for every survey point and loop for each parameter
  for(parameter_i in 1:length(parameters)){

    # parameter_i = 1

    # Get Unit_Level_Data Directory for this Param
    param.name = parameters[parameter_i]; param.name

    this_param_i = which(Parameter_Table$Dataset == param.name); this_param_i

    godir = paste(paramdir, param.name, "/Unit_Level_Data", sep = ""); godir

    paramsum = unlist(strsplit(as.vector(Parameter_Table[this_param_i, "Summaries"]), ";")); paramsum

    # For each unit
    for(unit_i in 1:length(unique_units)){

      # unit_i = 1

      unique_units[unit_i]

      # retrieve raw netcdf data for each unit
      ncfile = list.files(godir, pattern = paste0(unique_units[unit_i], "_"), full.names = T)

      # if there are no data then move to next unit
      if(length(ncfile) == 0){

        cat(paste0("Skipping ", unique_units[unit_i], " as there is no available data.\n"))
        next

      }

      nc_p = nc_open(ncfile)

      # Pull var array
      rawvar = ncvar_get(nc = nc_p, varid = as.vector(Parameter_Table$Fields[this_param_i]))

      # Pull dim vectors
      lon = ncvar_get(nc = nc_p, varid = "longitude"); lon
      lat = ncvar_get(nc = nc_p, varid = "latitude"); lat
      rawt = ncvar_get(nc = nc_p, varid = "time"); rawt

      # Close nc
      nc_close(nc_p)

      t = as_date(as_datetime(as.numeric(rawt), origin = ymd("1970/1/1")))
      head(t); tail(t)

      # Subset to Unit
      df_i = which(df$DATA_UNIT == unique_units[unit_i])

      if (max(lon) > 180) {

        xyt_df_i = df %>% mutate(lon = ifelse(df$lon < 0, df$lon + 360, df$lon))

      } else{

        xyt_df_i = df

      }

      # Locate all points in rawvar array, flagging any out of bound with NA in "ijk"
      ijk = xyt_to_ijk(xyt_df = as.data.frame(xyt_df_i[df_i,c("lon","lat","date_r")]),
                    x_grid = lon,
                    y_grid = lat,
                    t_grid = t,
                    lon = lon,
                    lat = lat,
                    t = t)
      ijk

      # Check for any out of bound points (stored as NA), we want to drop them from both ijk and df_i
      droprows = which(is.na(ijk), arr.ind = T)[,1] #Finds points as rows in ijk

      # If there are any NA, drops rows from ijk and indices from df_i
      if (length(droprows) > 0){

        ijk = ijk[-droprows,]
        df_i = df_i[-droprows]

      }

      if (length(df_i) == 0) {

        cat(paste0("No timeseries data found for ", unique_units[unit_i], ". Consider trying higher-resolution data. Skipping to the next...\n"))
        next()

      } else {

        # make sure rawvar is 3d array
        if (length(dim(rawvar)) == 1) {

          dim1 <- 1
          dim2 <- 1
          dim3 <- length(rawvar) / (dim1 * dim2)
          rawvar <- array(rawvar, dim = c(dim1, dim2, dim3))

        }

        # Count NA in var array (will use to solve NA issues)
        naP_xy = plyr::aaply(rawvar,
                             c(1,2),
                             na_stack_count)/dim(rawvar)[3]

        naP_xy

        # Id points sitting on NA-heavy timeseries
        i_masked = which(naP_xy[cbind(ijk$x_i, ijk$y_j)] > 0.8); i_masked

        # Infill selected points with spatial interpolation
        cnt = 1

        for(i_infill in i_masked){

          # Update NA blocks
          pNA = naP_xy[cbind(ijk$x_i[i_infill],
                             ijk$y_j[i_infill])]

          # Selected NA timeseries +/- x pixel steps
          ij_ex = 1

          while(pNA > 0.8 & ij_ex < 3){

            ij_vec = -ij_ex:ij_ex

            # Make sure selected NA timeseries +/- x pixel steps fits within the size of rawvar

            max_x = dim(rawvar)[1]
            max_y = dim(rawvar)[2]

            ts_x = ijk$x_i[i_infill]+ij_vec
            ts_y = ijk$y_j[i_infill]+ij_vec

            ts_x = subset(ts_x, ts_x <= max_x)
            ts_x = subset(ts_x, ts_x > 0)
            ts_y = subset(ts_y, ts_y <= max_y)
            ts_y = subset(ts_y, ts_y > 0)

            # Generates "infill" time series
            ts = aaply(rawvar[ts_x,
                              ts_y,],
                       c(3),
                       mean, na.rm = T)

            pNA = length(which(is.na(ts)))/length(ts)

            if(pNA < 0.9){

              # Update rawvar
              rawvar[ijk$x_i[i_infill],
                     ijk$y_j[i_infill],] = ts

              # Update naP
              naP_xy = aaply(rawvar,c(1,2),
                             na_stack_count)/dim(rawvar)[3]

              cat(paste("In-fill complete", cnt, "of", length(i_masked), ". Pixel +/-", ij_ex), "\n")

            }

            ij_ex = ij_ex + 1

          } # Close While

          # Fill in interpolated data
          cnt = cnt + 1

        } # Close infill for

        # Set Time Step
        if (Parameter_Table$Frequency[this_param_i] == "Monthly") tstep = 30.42 # (60*60*24*30.42)
        if (Parameter_Table$Frequency[this_param_i] == "14day") tstep = 14 # (60*60*24*14)
        if (Parameter_Table$Frequency[this_param_i] == "8day") tstep = 8 # (60*60*24*8)
        if (Parameter_Table$Frequency[this_param_i] == "5day") tstep = 5 # (60*60*24*5)
        if (Parameter_Table$Frequency[this_param_i] == "Weekly") tstep = 7 # (60*60*24*7)
        if (Parameter_Table$Frequency[this_param_i] == "Daily") tstep = 1 # (60*60*24*1)

        # TimeSeries Pull Indices
        ijk$t_01dy = ijk$t_k - (1/tstep-1)
        ijk$t_01wk = ijk$t_k - (7/tstep-1)
        ijk$t_01mo = round(ijk$t_k - (30.42/tstep-1))
        ijk$t_03mo = round(ijk$t_k - (91.25/tstep-1))
        ijk$t_06mo = round(ijk$t_k - (182.5/tstep-1))
        ijk$t_01yr = round(ijk$t_k - (1*365.25/tstep-1))
        ijk$t_03yr = round(ijk$t_k - (3*365.25/tstep-1))
        ijk$t_05yr = round(ijk$t_k - (5*365.25/tstep-1))
        ijk$t_10yr = round(ijk$t_k - (10*365.25/tstep-1))

        ijk[,c("t_k",
               "t_01dy",
               "t_01wk",
               "t_01mo",
               "t_03mo",
               "t_06mo",
               "t_01yr",
               "t_03yr",
               "t_05yr",
               "t_10yr")][which(ijk[,c("t_k",
                                       "t_01dy",
                                       "t_01wk",
                                       "t_01mo",
                                       "t_03mo",
                                       "t_06mo",
                                       "t_01yr",
                                       "t_03yr",
                                       "t_05yr",
                                       "t_10yr")] < 1, arr.ind = T)] = 1

        # Apply Summaries to Timeseries
        for(sum_i in 1:length(paramsum)){

          # sum_i = 1

          paramsum.name = paste0(paramsum[sum_i], "_", param.name); paramsum.name

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

          # For each point in df_i
          for(sumpt_i in 1:length(df_i)){

            # sumpt_i = 1

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

            } # END if

          } # END Loop over this unit's points (for 1:length(df_i))

          cat(paste("Processing", sumpt_i,
                    "points for unit", unique_units[unit_i],
                    "out of", length(unique_units),
                    "units.", paramsum.name, "completed.\n"))

          save(df, file = paste0(eds_path, "/eds_time.Rdata"))


        } # END Loop over each summary function  (for 1:length(paramsum))

      }

    } # END Loop over each unit

  }

  df[df == -9991] <- NA

  names(df) <- gsub("Daily", "daily", names(df)); names(df)
  names(df) <- gsub("Weekly", "weekly", names(df)); names(df)
  names(df) <- gsub("_8Day_", "_8days_", names(df)); names(df)
  names(df) <- gsub("Monthly", "monthly", names(df)); names(df)

  names(df) <- gsub("DHW.", "", names(df)); names(df)
  names(df) <- gsub("Major_", "major_", names(df)); names(df)
  names(df) <- gsub("Np10y_", "", names(df)); names(df)
  names(df) <- gsub("MeanMax_", "mean_maximum_", names(df)); names(df)
  names(df) <- gsub("CI95Max_", "ci95_maximum_", names(df)); names(df)
  names(df) <- gsub("MeanDur_", "mean_durnal_", names(df)); names(df)
  names(df) <- gsub("MaxMax_", "maximum_maximum_", names(df)); names(df)

  names(df) <- gsub("_DY01", "_01dy", names(df)); names(df)
  names(df) <- gsub("_WK01", "_01wk", names(df)); names(df)
  names(df) <- gsub("_MO01", "_01mo", names(df)); names(df)
  names(df) <- gsub("_MO03", "_03mo", names(df)); names(df)
  names(df) <- gsub("_MO06", "_06mo", names(df)); names(df)
  names(df) <- gsub("_YR10YR01", "_10yr_01yr", names(df)); names(df)
  names(df) <- gsub("_YR01", "_01yr", names(df)); names(df)
  names(df) <- gsub("_YR03", "_03yr", names(df)); names(df)
  names(df) <- gsub("_YR05", "_05yr", names(df)); names(df)
  names(df) <- gsub("_YR10", "_10yr", names(df)); names(df)
  names(df) <- gsub("_ALLB4", "_all_before", names(df)); names(df)

  df = df %>% mutate(across(where(is.numeric), ~ round(., 3)))

  df_time = df

  save(df_clim, file = file.path(eds_path, "eds_clim.RData"))
  save(df_time, file = file.path(eds_path, "eds_time.RData"))

  options(warn = 0)

  cat("\nEDS run has been successfully completed.\n")

}
