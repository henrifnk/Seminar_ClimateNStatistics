library(ggplot2)
library(sp)
library(dplyr)
library(magick)
library(gganimate)

te_1881 <- readLines("work/03-compounds/data/Temperature/MAM/TAMM_13_1881_01.asc")

xll_value_1881_t <- NULL
yll_value_1881_t <- NULL
corner_type_1881_t <- NULL

for (line in te_1881) {
  line <- trimws(line)
  
  if (startsWith(line, "XLL")) {
    parts <- strsplit(line, "\\s+")[[1]]
    if (length(parts) > 1) {
      xll_value_1881_t <- as.numeric(parts[2])
    }
    if (startsWith(line, "XLLCENTER")) {
      corner_type_1881_t <- "CENTER"
    } else if (startsWith(line, "XLLCORNER")) {
      corner_type_1881_t <- "CORNER"
    }
  }
  
  if (startsWith(line, "YLL")) {
    parts <- strsplit(line, "\\s+")[[1]]
    if (length(parts) > 1) {
      yll_value_1881_t <- as.numeric(parts[2])
    }
  }
  
  if (startsWith(line, "NODATA")) {
    break
  }
}

print(paste("XLL Value:", xll_value_1881_t))
print(paste("YLL Value:", yll_value_1881_t))
print(paste("Corner Type:", corner_type_1881_t))

data_lines <- grep("^\\s*-?\\d+", te_1881)  # 找到第一个以数字开头的行

if (length(data_lines) > 0) {
  grid_data <- matrix(nrow = length(data_lines), ncol = length(strsplit(te_1881[[data_lines[1]]], "\\s+")[[1]]), 
                      byrow = TRUE)
  
  for (i in seq_along(data_lines)) {
    grid_data[i, ] <- as.numeric(strsplit(te_1881[[data_lines[i]]], "\\s+")[[1]])
  }
  
  # 将未被占用的点（-999）替换为NA
  grid_data[grid_data == -999] <- NA
  
  row_coords <- seq(yll_value_1881_t, length.out = nrow(grid_data), by = 1000)  # 行坐标
  col_coords <- seq(xll_value_1881_t, length.out = ncol(grid_data), by = 1000)  # 列坐标
  
  # 创建带有坐标的数据框
  df_t <- as.data.frame(grid_data)
  rownames(df_t) <- rev(row_coords)
  colnames(df_t) <- col_coords
} else {
  print("No valid data lines found.")
}

map_data <- data.frame(
  Longitude = as.numeric(rep(colnames(df_t), times = nrow(df_t))),
  Latitude = as.numeric(rep(rownames(df_t), each = ncol(df_t))),
  Value = as.vector(t(df_t))  # 将矩阵转换为向量
)

# Creating data
GK <- map_data[, 1:2]

#Spatial Information, Coordinates Transform
coordinates(GK) <- c("Longitude", "Latitude")
proj4string(GK) <- CRS("+init=epsg:31467") # Defining Gauss Krüger
GK.WGS84 <- spTransform(GK, CRS("+init=epsg:4326")) # tranforming to WGS84 longlat

data_te_MAM <- data.frame(
  Longitude.GK = map_data[, 1],
  Latitude.GK = map_data[, 2],
  Longitude = coordinates(GK.WGS84)[,1],
  Latitude = coordinates(GK.WGS84)[,2],
  "1881" = map_data$Value)  # 假设这是点的大小数据


for (j in 1882:2023) {
  path <- paste0("TAMM_13_", j, "_01.asc")
  te_i <- readLines(paste0("work/03-compounds/data/Temperature/MAM/", path))
  
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
  
  if (!(xll_value == xll_value_1881_t & yll_value == yll_value_1881_t & corner_type == corner_type_1881_t)) {
    break
  }
  
  data_lines <- grep("^\\s*-?\\d+", te_i)  # 找到第一个以数字开头的行
  
  if (length(data_lines) > 0) {
    grid_data <- matrix(nrow = length(data_lines), ncol = length(strsplit(te_i[[data_lines[1]]], "\\s+")[[1]]), 
                        byrow = TRUE)
    
    for (i in seq_along(data_lines)) {
      grid_data[i, ] <- as.numeric(strsplit(te_i[[data_lines[i]]], "\\s+")[[1]])
    }
    
    # 将未被占用的点（-999）替换为NA
    grid_data[grid_data == -999] <- NA
    
    row_coords <- seq(yll_value, length.out = nrow(grid_data), by = 1000)  # 行坐标
    col_coords <- seq(xll_value, length.out = ncol(grid_data), by = 1000)  # 列坐标
    
    # 创建带有坐标的数据框
    df_t <- as.data.frame(grid_data)
    rownames(df_t) <- rev(row_coords)
    colnames(df_t) <- col_coords
  } else {
    print("No valid data lines found.")
  }
  
  map_data <- data.frame(
    Longitude.GK = as.numeric(rep(colnames(df_t), times = nrow(df_t))),
    Latitude.GK = as.numeric(rep(rownames(df_t), each = ncol(df_t))),
    j = as.vector(t(df_t))  # 将矩阵转换为向量
  )
  colnames(map_data)[3] <- paste0("X", j)
  
  data_te_MAM <- left_join(data_te_MAM, map_data, by = c("Longitude.GK", "Latitude.GK"))
}


te_1881 <- readLines("work/03-compounds/data/Temperature/SON/TAMM_15_1881_01.asc")

xll_value_1881_t <- NULL
yll_value_1881_t <- NULL
corner_type_1881_t <- NULL

for (line in te_1881) {
  line <- trimws(line)
  
  if (startsWith(line, "XLL")) {
    parts <- strsplit(line, "\\s+")[[1]]
    if (length(parts) > 1) {
      xll_value_1881_t <- as.numeric(parts[2])
    }
    if (startsWith(line, "XLLCENTER")) {
      corner_type_1881_t <- "CENTER"
    } else if (startsWith(line, "XLLCORNER")) {
      corner_type_1881_t <- "CORNER"
    }
  }
  
  if (startsWith(line, "YLL")) {
    parts <- strsplit(line, "\\s+")[[1]]
    if (length(parts) > 1) {
      yll_value_1881_t <- as.numeric(parts[2])
    }
  }
  
  if (startsWith(line, "NODATA")) {
    break
  }
}

print(paste("XLL Value:", xll_value_1881_t))
print(paste("YLL Value:", yll_value_1881_t))
print(paste("Corner Type:", corner_type_1881_t))

data_lines <- grep("^\\s*-?\\d+", te_1881)  # 找到第一个以数字开头的行

if (length(data_lines) > 0) {
  grid_data <- matrix(nrow = length(data_lines), ncol = length(strsplit(te_1881[[data_lines[1]]], "\\s+")[[1]]), 
                      byrow = TRUE)
  
  for (i in seq_along(data_lines)) {
    grid_data[i, ] <- as.numeric(strsplit(te_1881[[data_lines[i]]], "\\s+")[[1]])
  }
  
  # 将未被占用的点（-999）替换为NA
  grid_data[grid_data == -999] <- NA
  
  row_coords <- seq(yll_value, length.out = nrow(grid_data), by = 1000)  # 行坐标
  col_coords <- seq(xll_value, length.out = ncol(grid_data), by = 1000)  # 列坐标
  
  # 创建带有坐标的数据框
  df_t <- as.data.frame(grid_data)
  rownames(df_t) <- rev(row_coords)
  colnames(df_t) <- col_coords
} else {
  print("No valid data lines found.")
}

map_data <- data.frame(
  Longitude = as.numeric(rep(colnames(df_t), times = nrow(df_t))),
  Latitude = as.numeric(rep(rownames(df_t), each = ncol(df_t))),
  Value = as.vector(t(df_t))  # 将矩阵转换为向量
)

# Creating data
GK <- map_data[, 1:2]

#Spatial Information, Coordinates Transform
coordinates(GK) <- c("Longitude", "Latitude")
proj4string(GK) <- CRS("+init=epsg:31467") # Defining Gauss Krüger
GK.WGS84 <- spTransform(GK, CRS("+init=epsg:4326")) # tranforming to WGS84 longlat

data_te_SON <- data.frame(
  Longitude.GK = map_data[, 1],
  Latitude.GK = map_data[, 2],
  Longitude = coordinates(GK.WGS84)[,1],
  Latitude = coordinates(GK.WGS84)[,2],
  "1881" = map_data$Value)  # 假设这是点的大小数据


for (j in 1882:2023) {
  path <- paste0("TAMM_15_", j, "_01.asc")
  te_i <- readLines(paste0("work/03-compounds/data/Temperature/SON/", path))
  
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
  
  if (!(xll_value == xll_value_1881_t & yll_value == yll_value_1881_t & corner_type == corner_type_1881_t)) {
    break
  }
  
  data_lines <- grep("^\\s*-?\\d+", te_i)  # 找到第一个以数字开头的行
  
  if (length(data_lines) > 0) {
    grid_data <- matrix(nrow = length(data_lines), ncol = length(strsplit(te_i[[data_lines[1]]], "\\s+")[[1]]), 
                        byrow = TRUE)
    
    for (i in seq_along(data_lines)) {
      grid_data[i, ] <- as.numeric(strsplit(te_i[[data_lines[i]]], "\\s+")[[1]])
    }
    
    # 将未被占用的点（-999）替换为NA
    grid_data[grid_data == -999] <- NA
    
    row_coords <- seq(yll_value, length.out = nrow(grid_data), by = 1000)  # 行坐标
    col_coords <- seq(xll_value, length.out = ncol(grid_data), by = 1000)  # 列坐标
    
    # 创建带有坐标的数据框
    df_t <- as.data.frame(grid_data)
    rownames(df_t) <- rev(row_coords)
    colnames(df_t) <- col_coords
  } else {
    print("No valid data lines found.")
  }
  
  map_data <- data.frame(
    Longitude.GK = as.numeric(rep(colnames(df_t), times = nrow(df_t))),
    Latitude.GK = as.numeric(rep(rownames(df_t), each = ncol(df_t))),
    j = as.vector(t(df_t))  # 将矩阵转换为向量
  )
  colnames(map_data)[3] <- paste0("X", j)
  
  data_te_SON <- left_join(data_te_SON, map_data, by = c("Longitude.GK", "Latitude.GK"))
}

all(data_te_JJA[, 1:4] == data_te_MAM[, 1:4]) & all(data_te_MAM[, 1:4] == data_te_SON[, 1:4])
data_te <- data_te_JJA
data_te[, 5:ncol(data_te)] <- (92*data_te_MAM[, 5:ncol(data_te)] + 92*data_te_JJA[, 5:ncol(data_te)]*10 + 91*data_te_SON[, 5:ncol(data_te)])/2750

year <- 1881:2023

for (i in 5:ncol(data_te)) {
  mapsubset <- data_te[, c(3, 4, i)]
  colnames(mapsubset)[3] <- "Value"
  mapsubset <- na.omit(mapsubset)
  mapsubset <- mutate(mapsubset, level = NA)
  mapsubset$level <-
    cut(
      mapsubset$Value,
      breaks = c(-5, 0, 5, 8, 10, 12, 15, 20),
      labels = c("<= 0", "0-5", "5-8", "8-10", "10-12", "12-15",">=15"),
      order = TRUE,
      include.lowest = TRUE,
      right = TRUE
    )
  p <- ggplot() +
    coord_fixed(ratio = 1.3) +
    labs(title = "Temperature in Germany", subtitle = year[i - 4], caption = "Data source: DWD Climate Data Center (CDC)") + 
    geom_point(data = mapsubset, aes(x = Longitude, y = Latitude, color = level), size = 0.1) +
    scale_color_manual(
      breaks = c("<= 0", "0-5", "5-8", "8-10", "10-12", "12-15",">=15"),
      values = c("blue3", "lightblue","#FFFFCC","#FED976","#FEB24C","#FD8D3C","#E31A1C"),
      limits = c("<= 0", "0-5", "5-8", "8-10", "10-12", "12-15",">=15"))+
    theme_minimal()+
    guides(color = guide_legend(override.aes = list(size = 5), title = "Temperature in °C", reverse = TRUE))
  
  ggsave(
    plot = p,
    filename = paste0("work/03-compounds/figures/Temperature/", year[i - 4], ".png"),
    width = 20,
    height = 20,
    units = "cm",
    bg = "white"
  )
}
animate_p2 <-
  image_animate(image = image_read(path = paste0("work/03-compounds/figures/Temperature/",
                                                 year, ".png")), delay = 60)
anim_save(filename = "work/03-compounds/figures/Temperature/pre.gif", animation = animate_p2)
