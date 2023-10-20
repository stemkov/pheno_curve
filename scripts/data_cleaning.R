
library(lubridate)
library(mgsub)
library(data.table)

`%!in%` <- Negate(`%in%`)

setwd("/home/michael/Documents/Grad School/Research Projects/abundance_curve/for_submission")

### loading in raw data

bees <- read.csv("raw_data/bee_data.csv", stringsAsFactors = FALSE)
flowers <- read.csv("raw_data/flower_data.csv")
effort <- read.csv("raw_data/sampling_effort.csv")
locations <- read.csv("raw_data/plot_locations.csv")

bees_2021 <- read.csv("raw_data/bee_data_2021.csv", stringsAsFactors = FALSE)
flowers_2021 <- read.csv("raw_data/flower_data_2021.csv")
effort_2021 <- read.csv("raw_data/sampling_effort_2021.csv")

### combining years

bees <- rbind(bees, bees_2021)

effort <- rbind(effort, effort_2021[,-c(18:35)]) # removing coords from 2021 data to combine

# flowers needs a bit of work because I switched to counting all flowers in 2021 rather than number of plants and # flowers on 10 plants

flowers_2021$Number.of.plants <- 10
flowers_2021$Number.of.flowers.1 <- flowers_2021$Number.of.flowers/10
flowers_2021$Number.of.flowers.2 <- flowers_2021$Number.of.flowers/10
flowers_2021$Number.of.flowers.3 <- flowers_2021$Number.of.flowers/10
flowers_2021$Number.of.flowers.4 <- flowers_2021$Number.of.flowers/10
flowers_2021$Number.of.flowers.5 <- flowers_2021$Number.of.flowers/10
flowers_2021$Number.of.flowers.6 <- flowers_2021$Number.of.flowers/10
flowers_2021$Number.of.flowers.7 <- flowers_2021$Number.of.flowers/10
flowers_2021$Number.of.flowers.8 <- flowers_2021$Number.of.flowers/10
flowers_2021$Number.of.flowers.9 <- flowers_2021$Number.of.flowers/10
flowers_2021$Number.of.flowers.10 <- flowers_2021$Number.of.flowers/10
# unneeded columns
flowers$X <- NULL; flowers$X.1 <- NULL
flowers_2021$Number.of.flowers <- NULL
flowers <- rbind(flowers, flowers_2021)

### checking data

# checks to run to make sure I entered data correctly:

# bee data
# spelling of bee names
# whether date and site exist in sampling effort
# sex correct for solitaries, and caste for bumbles

sort(unique(bees$Field.ID))
table(bees$Field.ID)

sort(unique(bees$Caught.on))

# flower data
# spelling of flower names
# whether date and site exist in sampling effort
# NAs for "out of quadrat" and "out of transect" species
# plots are 1-6
# quadrats are A-D
# number of plants and number flr counts matches
# there's plant data for every sampled transect

sort(unique(flowers$Species))
table(flowers$Species)

# sampling effort/location
# effort and plot location data matches
# AM is before PM

### cleaning data

# put dates in ymd format

clean_bees <- bees
clean_bees$Date <- mdy(bees$Date) # problem here
clean_bees$Field.sex <- toupper(bees$Field.sex)

get.collector <- function(x){
  if(substr(x,1,2) == "AF") return("A. Fife")
  if(substr(x,1,2) == "MS") return("M. Stemkovski")
  else return("")
}

clean_bees$collector <- sapply(bees$ID.number, get.collector)

pinned <- clean_bees[which(clean_bees$Collected. == "Y"),]
unique_id <- paste("2019",1:nrow(pinned),sep="_")

clean_bees$unique_id <- ""
clean_bees$unique_id[which(clean_bees$Collected. == "Y")] <- unique_id

ordered_dets <- dets[order(dets$id),]
clean_bees[which(clean_bees$unique_id %in% dets$unique_id),"Lab.ID"] <- ordered_dets$determination
clean_bees[which(clean_bees$unique_id %in% dets$unique_id),"Lab.sex"] <- toupper(ordered_dets$sex)

clean_bees$year <- year(clean_bees$Date)
head(clean_bees)

# relabelling two late-season (August) points because they are reemergent females
clean_bees[which(clean_bees$Field.ID == "Halictus rubicundus" & clean_bees$year == 2019), ]
clean_bees[which(clean_bees$Field.ID == "Halictus rubicundus" & clean_bees$year == 2019 & clean_bees$Field.sex == "F" & (clean_bees$Date == "2019-08-15" | clean_bees$Date == "2019-08-12") ), "Field.ID"] <- c("Halictus rubicundus - reemergent")

write.csv(clean_bees, "clean_data/clean_bee_data.csv", row.names = FALSE)

clean_flowers <- flowers
clean_flowers$Date <- mdy(flowers$Date)
clean_flowers$year <- year(clean_flowers$Date)
clean_flowers$Species <- mgsub(clean_flowers$Species, c("Potentilla pulcherimma"), c("Potentilla pulcherrima"))

write.csv(clean_flowers, "clean_data/clean_flower_data.csv", row.names = FALSE)

clean_effort <- effort
clean_effort$Date <- mdy(effort$Date)
clean_effort$year <- year(clean_effort$Date)

write.csv(clean_effort, "clean_data/clean_sampling_effort.csv", row.names = FALSE)

clean_locations <- locations
clean_locations$Date <- mdy(locations$Date)

write.csv(clean_locations, "clean_data/clean_plot_locations.csv", row.names = FALSE)


### pulling out coords for map
# can skip this: not used in analysis

bees_loc <- fread("clean_data/clean_bee_data.csv", stringsAsFactors = FALSE)
locs_2019 <- fread("clean_data/clean_plot_locations.csv", stringsAsFactors = FALSE)
locs_2021 <- fread("raw_data/sampling_effort_2021.csv")

get.coords <- function(date_var, site_var, plot_var){
  
  if(site_var %!in% c("401 Trail", "Waterfall")) return(c(lat = NA_real_, lon = NA_real_))
  if(is.na(plot_var)) return(c(lat = NA_real_, lon = NA_real_))
  if(is.na(date_var)) return(c(lat = NA_real_, lon = NA_real_))
  if(is.na(site_var)) return(c(lat = NA_real_, lon = NA_real_))
  #if(any(is.na(c(date_var, site_var, plot_var)))) return(c(lat = NA_real_, lon = NA_real_))
  
  year <- year(ymd(date_var))
  
  if(year == 2019){
    coords <- locs_2019[Date == date_var & Plot.number == plot_var & Site. == site_var, Plot.Coordinates]
    # doing this because lat and lon are sometimes out of order
    coords_vec <- strsplit(coords, "\\s")
    coords_vec_clean <- c(sapply(coords_vec, function(x) gsub("[^0-9.-]","", x)))
    where_lat <- grep("^3", coords_vec_clean)
    lat <- coords_vec_clean[where_lat]
    lon <- coords_vec_clean[-where_lat]
  }
  if(year == 2021){
    lat_col <- paste("Plot", plot_var, "lat")
    lon_col <- paste("Plot", plot_var, "lon")
    lat <- as.numeric(subset(locs_2021, date == date_var & Site == site_var, select=lat_col)[1])
    lon <- as.numeric(subset(locs_2021, date == date_var & Site == site_var, select=lon_col)[1])
  }
  

  if(any(is.na(c(lat, lon)))) return(c(lat = NA_real_, lon = NA_real_))
  if(length(c(as.numeric(lat), as.numeric(lon))) != 2) return(c(lat = NA_real_, lon = NA_real_))
  
  return(c(lat = as.numeric(lat), lon = as.numeric(lon)))
}

coords <- mapply(get.coords, bees_loc$Date, bees_loc$Site, bees_loc$Plot)

bees_loc$lat <- coords[1,]
bees_loc$lon <- -abs(coords[2,]) # making sure lon is negative

# getting rid of a few typos
bees_loc[round(lon) == -108, lon := NA]
bees_loc[lat < 38.97, lat := NA]

plot(lat ~ lon, data = bees_loc) # yay

# n hrub by day/plot/site
hrub_loc <- bees_loc[, .(n_hrub = sum(Field.ID %in% c("Halictus rubicundus", "Halictus rubicundus ")),
                         lat = unique(lat),
                         lon = unique(lon)),
                           by=.(Date, Site, Plot)]


### map figure

library(sf)
site_map <- st_read("raw_data/rmbl_sites.kmz")

col_2019 <- "#dd5d2e"
col_2021 <- "#073f78"

hrub_loc[, year := year(Date), by=.(Date)]
hrub_loc[, col := adjustcolor("black", alpha.f=0.25)]
hrub_loc[year == 2019 & n_hrub >= 1, col := adjustcolor(col_2019, alpha.f=0.85)]
hrub_loc[year == 2021 & n_hrub >= 1, col := adjustcolor(col_2021, alpha.f=0.85)]


svg("figures/map.svg", width=9, height=5)
  par(mfrow=c(1,2))
  plot(site_map$geometry[1], lwd=2)
  points(lat ~ lon, data = hrub_loc, pch=19, col=col)
  points(lat ~ lon, data = hrub_loc, pch=1, col=col)
  
  plot(site_map$geometry[2], lwd=2)
  points(lat ~ lon, data = hrub_loc, pch=19, col=col)
  points(lat ~ lon, data = hrub_loc, pch=1, col=col)
dev.off()



