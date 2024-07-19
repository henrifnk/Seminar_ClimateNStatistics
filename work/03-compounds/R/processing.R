library(sp)
library(dplyr)

processing <- function(string1, string2) {
  data_1881 <- readLines(paste0("work/03-compounds/data/", string1, string2, "1881_01.asc"))
  
  xll_value_1881 <- NULL
  yll_value_1881 <- NULL
  corner_type_1881 <- NULL
  
  for (line in data_1881) {
    line <- trimws(line)
    
    if (startsWith(line, "XLL")) {
      parts <- strsplit(line, "\\s+")[[1]]
      if (length(parts) > 1) {
        xll_value_1881 <- as.numeric(parts[2])
      }
      if (startsWith(line, "XLLCENTER")) {
        corner_type_1881 <- "CENTER"
      } else if (startsWith(line, "XLLCORNER")) {
        corner_type_1881 <- "CORNER"
      }
    }
    
    if (startsWith(line, "YLL")) {
      parts <- strsplit(line, "\\s+")[[1]]
      if (length(parts) > 1) {
        yll_value_1881 <- as.numeric(parts[2])
      }
    }
    
    if (startsWith(line, "NODATA")) {
      break
    }
  }
  print(paste("XLL Value:", xll_value_1881))
  print(paste("YLL Value:", yll_value_1881))
  print(paste("Corner Type:", corner_type_1881))
  
  data_lines <- grep("^\\s*-?\\d+", data_1881) # Find the first line that starts with a number
  
  if (length(data_lines) > 0) {
    grid_data <- matrix(nrow = length(data_lines), ncol = length(strsplit(data_1881[[data_lines[1]]], "\\s+")[[1]]), 
                        byrow = TRUE)
    
    for (i in seq_along(data_lines)) {
      grid_data[i, ] <- as.numeric(strsplit(data_1881[[data_lines[i]]], "\\s+")[[1]])
    }
    
    
    grid_data[grid_data == -999] <- NA
    
    row_coords <- seq(yll_value_1881, length.out = nrow(grid_data), by = 1000)
    col_coords <- seq(xll_value_1881, length.out = ncol(grid_data), by = 1000)  
    
    df <- as.data.frame(grid_data)
    rownames(df) <- rev(row_coords)
    colnames(df) <- col_coords
  } else {
    print("No valid data lines found.")
  }
  
  map_data <- data.frame(
    Longitude = as.numeric(rep(colnames(df), times = nrow(df))),
    Latitude = as.numeric(rep(rownames(df), each = ncol(df))),
    Value = as.vector(t(df))  # 将矩阵转换为向量
  )
  
  # Creating data
  GK <- map_data[, 1:2]
  
  #Spatial Information, Coordinates Transform
  coordinates(GK) <- c("Longitude", "Latitude")
  proj4string(GK) <- CRS("+init=epsg:31467") # Defining Gauss Krüger
  GK.WGS84 <- spTransform(GK, CRS("+init=epsg:4326")) # tranforming to WGS84 longlat
  
  data <- data.frame(
    Longitude.GK = map_data[, 1],
    Latitude.GK = map_data[, 2],
    Longitude = coordinates(GK.WGS84)[,1],
    Latitude = coordinates(GK.WGS84)[,2],
    "1881" = map_data$Value)
  
  for (j in 1882:2023) {
    path <- paste0(string2, j, "_01.asc")
    te_i <- readLines(paste0("work/03-compounds/data/", string1, path))
    
    xll_value <- NULL
    yll_value <- NULL
    corner_type <- NULL
    
    for (line in te_i) {
      line <- trimws(line)
      
      if (startsWith(line, "XLL")) {
        parts <- strsplit(line, "\\s+")[[1]]
        if (length(parts) > 1) {
          xll_value <- as.numeric(parts[2])
        }
        if (startsWith(line, "XLLCENTER")) {
          corner_type <- "CENTER"
        } else if (startsWith(line, "XLLCORNER")) {
          corner_type <- "CORNER"
        }
      }
      
      if (startsWith(line, "YLL")) {
        parts <- strsplit(line, "\\s+")[[1]]
        if (length(parts) > 1) {
          yll_value <- as.numeric(parts[2])
        }
      }
      
      if (startsWith(line, "NODATA")) {
        break
      }
    }
    
    if (!(xll_value == xll_value_1881 & yll_value == yll_value_1881 & corner_type == corner_type_1881)) {
      break
    }
    
    data_lines <- grep("^\\s*-?\\d+", te_i)
    
    if (length(data_lines) > 0) {
      grid_data <- matrix(nrow = length(data_lines), ncol = length(strsplit(te_i[[data_lines[1]]], "\\s+")[[1]]), 
                          byrow = TRUE)
      
      for (i in seq_along(data_lines)) {
        grid_data[i, ] <- as.numeric(strsplit(te_i[[data_lines[i]]], "\\s+")[[1]])
      }
      
      grid_data[grid_data == -999] <- NA
      
      row_coords <- seq(yll_value, length.out = nrow(grid_data), by = 1000)  # 行坐标
      col_coords <- seq(xll_value, length.out = ncol(grid_data), by = 1000)  # 列坐标
      
      df <- as.data.frame(grid_data)
      rownames(df) <- rev(row_coords)
      colnames(df) <- col_coords
    } else {
      print("No valid data lines found.")
    }
    
    map_data <- data.frame(
      Longitude.GK = as.numeric(rep(colnames(df), times = nrow(df))),
      Latitude.GK = as.numeric(rep(rownames(df), each = ncol(df))),
      j = as.vector(t(df))) 
      
      colnames(map_data)[3] <- paste0("X", j)
      
      data <- left_join(data, map_data, by = c("Longitude.GK", "Latitude.GK"))
  }
  data[, 5:ncol(data)] <- data[, 5:ncol(data)]
  
  return(data)
}

data_te_MAM <- processing("Temperature/MAM/", "TAMM_13_")
data_te_JJA <- processing("Temperature/JJA/", "TAMM_14_")
data_te_SON <- processing("Temperature/SON/", "TAMM_15_")

all(data_te_JJA[, 1:4] == data_te_MAM[, 1:4]) & all(data_te_MAM[, 1:4] == data_te_SON[, 1:4])

data_te <- data_te_JJA
data_te[, 5:ncol(data_te)] <- (92*data_te_MAM[, 5:ncol(data_te)]
                               + 92*data_te_JJA[, 5:ncol(data_te)]
                               + 91*data_te_SON[, 5:ncol(data_te)])/2750

data_pr_MAM <- processing("Precipitation/MAM/", "RSMS_13_")
data_pr_JJA <- processing("Precipitation/JJA/", "RSMS_14_")
data_pr_SON <- processing("Precipitation/SON/", "RSMS_15_")

all(data_pr_JJA[, 1:4] == data_pr_MAM[, 1:4]) & all(data_pr_MAM[, 1:4] == data_pr_SON[, 1:4])
data_pr <- data_pr_JJA
data_pr[, 5:ncol(data_pr)] <- (data_pr_MAM[, 5:ncol(data_pr)]
                               + data_pr_JJA[, 5:ncol(data_pr)]
                               + data_pr_SON[, 5:ncol(data_pr)])
