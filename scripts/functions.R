`%!in%` <- Negate(`%in%`)

### simulation functions

verhulst.func <- function(t_given, B0, e, f, m, t_onset, get="P"){
  
  if(FALSE){
    B0 <- 40; e <- 0.2; f <- 0.004; m <- 0.1
    t_given <- 111; t_onset <- 110.1
  }
  
  t_length <- 100
  t_vec <- c(1:t_length) - 1
  P <- rep(NA, t_length)
  B <- rep(NA, t_length)
  
  if(t_onset >= t_given & get == "P") return(0)
  if(t_onset >= t_given & get == "B") return(B0)
  
  P[1] <- 0
  B[1] <- B0
  
  for(t in 2:t_length){
    B[t] <- B[t-1] - e*B[t-1] + f*B[t-1]^2;
    P[t] <- P[t-1] + e*B[t-1] - f*B[t-1]^2 - m*P[t-1];
  }
  
  t_out <- t_given-t_onset
  P_out <- approx(t_vec, P, t_out)
  B_out <- approx(t_vec, B, t_out)
  
  if(get == "P") return(P_out$y)
  if(get == "B") return(B_out$y)
}

add.covars <- function(a_flr, flower_data){
  out <- a_flr*flower_data
  return(out)
}

verhulst.wrapper <- function(t, B0, e, f, m, t_onset,
                             a_flr=NA, flower_data=NA,
                             plot=F, t_data=NA, P_data=NA,
                             plot_lines=F){
  P_covars <- NA # just so that return has something if no covars included
  P <- sapply(t, verhulst.func, B0, e, f, m, t_onset, "P")
  B <- sapply(t, verhulst.func, B0, e, f, m, t_onset, "B")
  if(!is.na(a_flr)) P_covars <- P[which(t %in% t_data)] + add.covars(a_flr, flower_data)
  if(plot & !plot_lines){
    plot(P ~ t, type="l", lwd=2, col="red")
    if(!is.na(a_flr)){
      arrows(x0 = t_data, y0 = P[which(t %in% t_data)], y1 = P_covars, col="red", lwd=2, code=0)
      points(P_covars ~ t_data, col="red", pch=19)
    } 
    if(!any(is.na(c(t_data, P_data)))) points(P_data ~ t_data)
    
  }
  if(plot_lines) lines(P ~ t, lwd=1, col=adjustcolor("red", 0.25)) # just plot curve for ppc
  return(list(P=P, B=B, t=t, B0=B0, e=e, f=f, m=m, t_onset=t_onset, a_flr=a_flr, flower_data=flower_data, t_data=t_data, P_data=P_data, P_covars=P_covars))
}


### data collation functions
# apologies, there are many

get.sampling.dates <- function(site, year=2019){ ymd(effort[which(effort$Site == site & effort$year == year),"Date"]) }

get.effort <- function(site, date, year=2019, plots="all"){
  this_effort <- effort[which(effort$Date == date & effort$Site == site & effort$year == year),]
  
  if("all" %in% plots){plots <- 1:6}
  
  which_plots <- paste("Plot",rep(plots,2),"start",rep(c("AM","PM"),each=length(plots)),sep=".")
  total_effort <- sum(this_effort[,which_plots] != "")*10
  names(total_effort) <- date
  return(total_effort)
}

# species=="all" is a special case that returns a table of all species caught
get.abundance <- function(species, site, date, year=2019, plots="all", sex="all", effort_adjust=FALSE){
  if("all" %in% plots){plots <- 1:6}
  this_period <- bees[which(bees$Date == date & bees$Site == site & bees$Plot %in% plots & bees$year == year),]
  if(species == "all"){species <- unique(this_period$Field.ID)}
  if(sex == "all"){sex <- c("","F","M","Q","W")}
  this_sp <- this_period[which(this_period$Field.ID %in% species & this_period$Field.sex %in% sex),]
  ab_table <- table(this_sp$Field.ID, this_sp$Field.sex)
  if(effort_adjust){ab_table <- ab_table/(get.effort(site, date, year, plots)/60)} # returns bees caught per hour
  if(nrow(ab_table) == 0){return(0)}
  return(ab_table)
}

# species=="all" returns sum off all flowers regardless of species
get.flowers <- function(species, site, year=2019, date, plots="all", quads="all"){
  if("all" %in% plots){plots <- 1:6}
  if(species == "all"){
    this_period <- flowers[which(flowers$Date == date & flowers$Site == site & flowers$Plot %in% plots & flowers$In.out.quadrats == "in" & flowers$In.out.transects == "in" & flowers$year == year),]
  } else{
    this_period <- flowers[which(flowers$Date == date & flowers$Site == site & flowers$Plot %in% plots & flowers$Species == species & flowers$In.out.quadrats == "in" & flowers$In.out.transects == "in" & flowers$year == year),]
  }
  if(nrow(this_period) == 0){return(0)}
  flr_cols <- paste("Number.of.flowers.",1:10,sep="")
  total_flrs <- sapply(1:nrow(this_period), function(x){this_period[x,"Number.of.plants"]*mean(as.numeric(this_period[x,flr_cols]), na.rm = TRUE)})
  total_flrs[which(is.na(total_flrs))] <- 0
  all_flrs <- sum(total_flrs)
  return(all_flrs)
}

plot.timeseries <- function(species, site, year=2019, legend_position = "topright", effort_adjust=FALSE){
  days <- get.sampling.dates(site, year)
  doys <- yday(days)
  this_effort <- sapply(days, get.effort, site=site, year=year)
  
  ab_f <- mapply(get.abundance, species=species, site=site, date=days, year=year, sex="F", effort_adjust=effort_adjust, USE.NAMES = FALSE)
  ab_m <- mapply(get.abundance, species=species, site=site, date=days, year=year, sex="M", effort_adjust=effort_adjust, USE.NAMES = FALSE)
  
  plot_range <- c(min(ab_f, ab_m),max(ab_f, ab_m))
  
  plot(NA, xlim=range(doys), ylim=plot_range, axes=F, xlab="Day of year", ylab="Abundance")
  axis(1, seq(range(doys)[1], range(doys)[2], 10), seq(range(doys)[1], range(doys)[2], 10), las=2)
  axis(2, seq(plot_range[1], plot_range[2]), seq(plot_range[1], plot_range[2]), las=2)
  points(ab_f ~ doys, pch=20, col="blue")
  lines(ab_f ~ doys, lty=2, col="blue")
  
  points(ab_m ~ doys, pch=20, col="red")
  lines(ab_m ~ doys, lty=2, col="red")
  
  legend(legend_position, inset=0.05, c("females", "males"), col=c("blue", "red"), pch=20, lty=2, cex=0.8,
         title=paste(species," \n ",site," site",sep=""), bty="n")
  box(bty="L",lwd=1.5)
}

get.timeseries <- function(species, site, year=2019, effort_adjust=FALSE){
  days <- get.sampling.dates(site, year)
  doys <- yday(days)
  this_effort <- sapply(days, get.effort, site=site, year=year)
  
  ab_f <- mapply(get.abundance, species=species, site=site, date=days, year=year, sex="F", effort_adjust=effort_adjust, USE.NAMES = FALSE)
  ab_m <- mapply(get.abundance, species=species, site=site, date=days, year=year, sex="M", effort_adjust=effort_adjust, USE.NAMES = FALSE)
  
  output <- data.frame(date = days,
                       doy = doys,
                       effort = this_effort,
                       females = ab_f,
                       males = ab_m,
                       row.names = 1:length(doys))
  return(output)
}

plot.flowers <- function(species, site, year=2019, legend_position="topright"){
  days <- get.sampling.dates(site, year)
  doys <- yday(days)
  this_effort <- sapply(days, get.effort, site=site, year=year)
  
  flr_frame <- sapply(species, function(sp){mapply(get.flowers, species=sp, site=site, date=days, year=year, USE.NAMES = FALSE)})
  
  plot_range <- range(flr_frame)
  
  plot(NA, xlim=range(doys), ylim=plot_range, axes=F, xlab="Day of year", ylab="Abundance")
  axis(1, seq(range(doys)[1], range(doys)[2], 10), seq(range(doys)[1], range(doys)[2], 10), las=2)
  axis(2, seq(plot_range[1], plot_range[2], 300), seq(plot_range[1], plot_range[2], 300), las=2)
  points(flr_frame[,2] ~ doys, pch=20, col="blue")
  lines(flr_frame[,2] ~ doys, lty=2, col="blue")
  
  par(new=TRUE)
  plot(NA, xlim=range(doys), ylim=range(flr_frame[,1]), axes=F, xlab="Day of year", ylab="Abundance")
  points(flr_frame[,1] ~ doys, pch=20, col="red")
  lines(flr_frame[,1] ~ doys, lty=2, col="red")
  axis(4, seq(range(flr_frame[,1])[1], range(flr_frame[,1])[2], 4), seq(range(flr_frame[,1])[1], range(flr_frame[,1])[2], 4), las=2)
  
  legend(legend_position, inset=0.05, c("Potentilla p.", "Taraxacum o."), col=c("blue", "red"), pch=20, lty=2, cex=0.8,
         title=paste(site," site",sep=""), bty="n")
  box(bty="U",lwd=1.5)
}

get.flower.timeseries <- function(species, site, year=2019){
  days <- get.sampling.dates(site, year)
  doys <- yday(days)
  this_effort <- sapply(days, get.effort, site=site, year=year)
  
  flr_frame <- as.data.frame(sapply(species, function(sp){mapply(get.flowers, species=sp, site=site, date=days, year=year, USE.NAMES = FALSE)}))
  colnames(flr_frame) <- c("P")
  flr_frame$date <- as.character(days)
  flr_frame$doy <- doys
  
  return(flr_frame)
}

get.temp <- function(site, location, year){
  if(FALSE){
    site <- "401 Trail"
    location <- "air"
    year <- 2019
  }
  site_var <- site; location_var <- location; year_var <- year
  if(location %!in% c("soil", "air")) stop("Invalid location")
  dates <- get.sampling.dates(site, year)
  temps <- rep(NA, length(dates))
  for(i in seq_along(dates)){
    d <- dates[i]
    # between 6am and 6pm
    day_data <- temp[site == site_var & location == location_var & date(date) == d & hour(date) >= 6 & hour(date) < 18, ]
    day_avg <- median(day_data$temp_c, na.rm=T)
    temps[i] <- day_avg
  }
  return(temps)
}

# wrapper to get all relevant data for timeseries
assemble.data <- function(site, year, time_shift=NA, plot=F){
  species <- "Halictus rubicundus"
  if(FALSE){
    site <- "Waterfall"
    year <- 2019
    time_shift <- 1
  }
  
  raw_data <- get.timeseries(species, site, year=year)
  if(is.na(time_shift)){
    data <- data.frame(time = raw_data$doy,
                       P = raw_data$females)
  } else{
    data <- data.frame(time = raw_data$doy - (min(raw_data$doy)-time_shift),
                       P = raw_data$females)
  }
  
  data$taraxacum <- get.flower.timeseries("Taraxacum officinale", site, year)$P
  # if soil and air temperatures are desired:
  # data$air <- get.temp(site, "air", year)
  # data$soil <- get.temp(site, "soil", year)
  data$site <- site; data$year <- year; data$time_shift <- time_shift
  
  if(plot) plot(P ~ time, data=data)
  
  return(data)
}



