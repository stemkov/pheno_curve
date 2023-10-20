library(rstan)
library(bayesplot)
library(plotrix)

setwd("/home/michael/Documents/Grad School/Research Projects/abundance_curve/for_submission")

bees <- read.csv("clean_data/clean_bee_data.csv", stringsAsFactors = FALSE)
flowers <- read.csv("clean_data/clean_flower_data.csv", stringsAsFactors = FALSE)
effort <- read.csv("clean_data/clean_sampling_effort.csv", stringsAsFactors = FALSE)
locations <- read.csv("clean_data/clean_plot_locations.csv", stringsAsFactors = FALSE)

### assembling full dataset to analyze
sites <- c("401 Trail", "Waterfall", "401 Trail", "Waterfall")
years <- c(2019, 2019, 2021, 2021)
# sites <- c("401 Trail", "401 Trail")
# years <- c(2019, 2021)
all_data_list <- mapply(assemble.data, sites, years, SIMPLIFY = F, USE.NAMES = F)
all_data <- rbindlist(all_data_list)
all_data[, site_year_label := as.factor(paste(site, year, sep="_"))]
all_data[, site_year := as.numeric(site_year_label)]

Ps <- all_data$P
times <- all_data$time
ts <- all_data$site_year # for B0 to vary by timeseries
year_id <- as.numeric(as.factor(all_data$year)) # for onset to vary by year
flowers <- all_data$taraxacum

data <- list(P = Ps, time = times, flowers = flowers, n = length(Ps), ts = ts, n_ts = max(ts), year = year_id, n_year = max(year_id))

### fitting model

model <- stan_model(file='/home/michael/Documents/Grad School/Research Projects/abundance_curve/scripts/model.stan')
fit <- sampling(model, data, iter=10000, chains=4, cores=4, control=list(adapt_delta=0.95)) 
saveRDS(fit, "scripts/model_fit.RDS")
# fit <- readRDS("scripts/model_fit.RDS") # read in previous fit if you want to skip fitting, which takes a while

### Just ran the task above.
# maybe try centering flower data... but I guess that would give me negative predictions unless I also have some intercept... not sure how to interpret the intercept

summary(fit)$summary
posterior <- extract(fit, permute=F)

### diagnostics

png("figures/trace_B0.png", width=700, height=700)
  par(mfrow=c(2,2))
  matplot(posterior[,,"B0[1]"], type="l", lty=1)
  matplot(posterior[,,"B0[2]"], type="l", lty=1)
  matplot(posterior[,,"B0[3]"], type="l", lty=1)
  matplot(posterior[,,"B0[4]"], type="l", lty=1)
dev.off()

png("figures/trace_e_f_m.png", width=600, height=900)
 par(mfrow=c(3,2))
  matplot(posterior[,,"e[1]"], type="l", lty=1)
  matplot(posterior[,,"e[2]"], type="l", lty=1)
  matplot(posterior[,,"f[1]"], type="l", lty=1)
  matplot(posterior[,,"f[2]"], type="l", lty=1)
  matplot(posterior[,,"m[1]"], type="l", lty=1)
  matplot(posterior[,,"m[2]"], type="l", lty=1)
dev.off()

png("figures/trace_onset.png", width=400, height=700)
  par(mfrow=c(2,1))
  matplot(posterior[,,"onset[1]"], type="l", lty=1)
  matplot(posterior[,,"onset[2]"], type="l", lty=1)
dev.off()

png("figures/trace_flower_sigma.png", width=400, height=700)
  par(mfrow=c(2,1))
  matplot(posterior[,,"a_flr"], type="l", lty=1)
  matplot(posterior[,,"sigma"], type="l", lty=1)
dev.off()

matplot(posterior[,,"lp__"], type="l", lty=1)
matplot(posterior[,,"onset_mu"], type="l", lty=1)
matplot(posterior[,,"onset_sigma"], type="l", lty=1)


ts1_data <- all_data[site_year==1,]
ts2_data <- all_data[site_year==2,]
ts3_data <- all_data[site_year==3,]
ts4_data <- all_data[site_year==4,]

# posterior predictive check
png("figures/ppc.png", width=700, height=700)
  par(mfrow=c(2,2))
  ppc_samples <- sample(1:(3000*4), 200)
  
  plot(NA, xlab="time", ylab="P", xlim=c(130,220), ylim=c(0,8))
  sapply(ppc_samples, function(x) verhulst.wrapper(c(130:220),B0 = posterior[,,"B0[1]"][x],e = posterior[,,"e[1]"][x],f = posterior[,,"f[1]"][x],m = posterior[,,"m[1]"][x],t_onset = posterior[,,"onset[1]"][x], a_flr = posterior[,,"a_flr"][x], plot_lines=T))
  points(P ~ time, data = ts1_data)
  
  plot(NA, xlab="time", ylab="P", xlim=c(130,220), ylim=c(0,8))
  sapply(ppc_samples, function(x) verhulst.wrapper(c(130:220),B0 = posterior[,,"B0[2]"][x],e = posterior[,,"e[2]"][x],f = posterior[,,"f[2]"][x],m = posterior[,,"m[2]"][x],t_onset = posterior[,,"onset[2]"][x], a_flr = posterior[,,"a_flr"][x], plot_lines=T))
  points(P ~ time, data = ts2_data)
  
  plot(NA, xlab="time", ylab="P", xlim=c(130,220), ylim=c(0,8))
  sapply(ppc_samples, function(x) verhulst.wrapper(c(130:220),B0 = posterior[,,"B0[3]"][x],e = posterior[,,"e[1]"][x],f = posterior[,,"f[1]"][x],m = posterior[,,"m[1]"][x],t_onset = posterior[,,"onset[1]"][x], a_flr = posterior[,,"a_flr"][x], plot_lines=T))
  points(P ~ time, data = ts3_data)
  
  plot(NA, xlab="time", ylab="P", xlim=c(130,220), ylim=c(0,8))
  sapply(ppc_samples, function(x) verhulst.wrapper(c(130:220),B0 = posterior[,,"B0[4]"][x],e = posterior[,,"e[2]"][x],f = posterior[,,"f[2]"][x],m = posterior[,,"m[2]"][x],t_onset = posterior[,,"onset[2]"][x], a_flr = posterior[,,"a_flr"][x], plot_lines=T))
  points(P ~ time, data = ts4_data)
dev.off()

# priors vs posteriors
par(mfrow=c(2,2))
hist(rnorm(20000, 40, 10), 100, border="white", col=adjustcolor("red", alpha.f = 0.5), main="B0_mu")
hist(posterior[,,"B0_mu"], 50, add=T, border="white", col=adjustcolor("blue", alpha.f = 0.5))
hist(rlnorm(20000, 0, 1), 100, border="white", col=adjustcolor("red", alpha.f = 0.5), main="B0_sigma")
hist(posterior[,,"B0_sigma"], 100, add=T, border="white", col=adjustcolor("blue", alpha.f = 0.5))

hist(rnorm(20000, 160, 10), 100, border="white", col=adjustcolor("red", alpha.f = 0.5), main="onset_mu")
hist(posterior[,,"onset_mu"], 50, add=T, border="white", col=adjustcolor("blue", alpha.f = 0.5))
hist(rlnorm(20000, 1, 1), 100, border="white", col=adjustcolor("red", alpha.f = 0.5), main="onset_sigma")
hist(posterior[,,"onset_sigma"], 75, add=T, border="white", col=adjustcolor("blue", alpha.f = 0.5))

hist(rnorm(20000, 0, 1), 100, border="white", col=adjustcolor("red", alpha.f = 0.5), main="e_tilde_mu")
hist(posterior[,,"e_tilde_mu"], 50, add=T, border="white", col=adjustcolor("blue", alpha.f = 0.5))
hist(rnorm(20000, 0, 1), 100, border="white", col=adjustcolor("red", alpha.f = 0.5), main="f_tilde_mu")
hist(posterior[,,"f_tilde_mu"], 100, add=T, border="white", col=adjustcolor("blue", alpha.f = 0.5))
hist(rnorm(20000, 0, 1), 100, border="white", col=adjustcolor("red", alpha.f = 0.5), main="m_tilde_mu")
hist(posterior[,,"m_tilde_mu"], 100, add=T, border="white", col=adjustcolor("blue", alpha.f = 0.5))

hist(rnorm(20000, 0, 5), 100, border="white", col=adjustcolor("red", alpha.f = 0.5), main="a_flr")
hist(posterior[,,"a_flr"], 1, add=T, border="white", col=adjustcolor("blue", alpha.f = 0.5))
hist(rlnorm(20000, 0, 1), 100, border="white", col=adjustcolor("red", alpha.f = 0.5), main="sigma")
hist(posterior[,,"sigma"], 2, add=T, border="white", col=adjustcolor("blue", alpha.f = 0.5))



par(mfrow=c(2,2))
hist(posterior[,,"B0[1]"], 100)
hist(posterior[,,"B0[2]"], 100)
hist(posterior[,,"B0[3]"], 100)
hist(posterior[,,"B0[4]"], 100)

hist(posterior[,,"e[1]"], 100, col=adjustcolor("red", 0.25))
hist(posterior[,,"e[2]"], 100, col=adjustcolor("blue", 0.25), add=T)
hist(posterior[,,"f[1]"], 100, col=adjustcolor("red", 0.25))
hist(posterior[,,"f[2]"], 100, col=adjustcolor("blue", 0.25), add=T)
hist(posterior[,,"m[1]"], 100, col=adjustcolor("red", 0.25))
hist(posterior[,,"m[2]"], 100, col=adjustcolor("blue", 0.25), add=T)

hist(posterior[,,"onset[1]"], 100, col=adjustcolor("red", 0.25), xlim=c(130,185), border=adjustcolor("white", 0))
hist(posterior[,,"onset[2]"], 100, col=adjustcolor("blue", 0.25), add=T, border=adjustcolor("white", 0))

par(mfrow=c(1,1))
hist(posterior[,,"e"], 100)
hist(posterior[,,"f"], 100)
hist(posterior[,,"m"], 100)
hist(posterior[,,"onset[1]"], 100)
abline(v=mean(posterior[,,"onset[1]"]), col="red", lwd=2)
mean(posterior[,,"onset[1]"])
as.Date(mean(posterior[,,"onset[1]"]), origin = "2019-01-01")
hist(posterior[,,"onset[2]"], 100)
abline(v=mean(posterior[,,"onset[2]"]), col="red", lwd=2)
mean(posterior[,,"onset[2]"])
as.Date(mean(posterior[,,"onset[2]"]), origin = "2021-01-01")

hist(posterior[,,"sigma"], 100)
hist(posterior[,,"a_flr"], 100)

hist(posterior[,,"onset_mu"], 100)
hist(posterior[,,"onset_sigma"], 100)
hist(posterior[,,"B0_mu"], 100)
hist(posterior[,,"B0_sigma"], 100) # some truncation because I use a normal dist

hist(posterior[,,"e_tilde_mu"], 100)
hist(posterior[,,"e_tilde_sigma"], 100) # some truncation because I use a normal dist
hist(posterior[,,"a_flr"], 100)
hist(posterior[,,"sigma"], 100)

par(mfrow=c(2,2))

# chain 1 fits
verhulst.wrapper(c(120:250),B0 = mean(posterior[,1,"B0[1]"]),e = mean(posterior[,1,"e[1]"]),f = mean(posterior[,1,"f[1]"]),m = mean(posterior[,1,"m[1]"]),t_onset = mean(posterior[,1,"onset[1]"]),a_flr = mean(posterior[,1,"a_flr"]),plot=T, t_data = ts1_data$time, P_data = ts1_data$P, flower_data = ts1_data$taraxacum)
verhulst.wrapper(c(120:250),B0 = mean(posterior[,1,"B0[2]"]),e = mean(posterior[,1,"e[2]"]),f = mean(posterior[,1,"f[2]"]),m = mean(posterior[,1,"m[2]"]),t_onset = mean(posterior[,1,"onset[2]"]), a_flr=mean(posterior[,1,"a_flr"]), plot=T, t_data = ts2_data$time, P_data = ts2_data$P, flower_data = ts2_data$taraxacum)
verhulst.wrapper(c(120:250),B0 = mean(posterior[,1,"B0[3]"]),e = mean(posterior[,1,"e[1]"]),f = mean(posterior[,1,"f[1]"]),m = mean(posterior[,1,"m[1]"]),t_onset = mean(posterior[,1,"onset[1]"]),a_flr = mean(posterior[,1,"a_flr"]),plot=T, t_data = ts3_data$time, P_data = ts3_data$P, flower_data = ts3_data$taraxacum)
verhulst.wrapper(c(120:250),B0 = mean(posterior[,1,"B0[4]"]),e = mean(posterior[,1,"e[2]"]),f = mean(posterior[,1,"f[2]"]),m = mean(posterior[,1,"m[2]"]),t_onset = mean(posterior[,1,"onset[2]"]), a_flr=mean(posterior[,1,"a_flr"]), plot=T, t_data = ts4_data$time, P_data = ts4_data$P, flower_data = ts4_data$taraxacum)

# chain 2
verhulst.wrapper(c(120:250),B0 = mean(posterior[,2,"B0[1]"]),e = mean(posterior[,2,"e[1]"]),f = mean(posterior[,2,"f[1]"]),m = mean(posterior[,2,"m[1]"]),t_onset = mean(posterior[,2,"onset[1]"]),a_flr = mean(posterior[,2,"a_flr"]),plot=T, t_data = ts1_data$time, P_data = ts1_data$P, flower_data = ts1_data$taraxacum)
verhulst.wrapper(c(120:250),B0 = mean(posterior[,2,"B0[2]"]),e = mean(posterior[,2,"e[2]"]),f = mean(posterior[,2,"f[2]"]),m = mean(posterior[,2,"m[2]"]),t_onset = mean(posterior[,2,"onset[2]"]), a_flr=mean(posterior[,2,"a_flr"]), plot=T, t_data = ts2_data$time, P_data = ts2_data$P, flower_data = ts2_data$taraxacum)
verhulst.wrapper(c(120:250),B0 = mean(posterior[,2,"B0[3]"]),e = mean(posterior[,2,"e[1]"]),f = mean(posterior[,2,"f[1]"]),m = mean(posterior[,2,"m[1]"]),t_onset = mean(posterior[,2,"onset[1]"]),a_flr = mean(posterior[,2,"a_flr"]),plot=T, t_data = ts3_data$time, P_data = ts3_data$P, flower_data = ts3_data$taraxacum)
verhulst.wrapper(c(120:250),B0 = mean(posterior[,2,"B0[4]"]),e = mean(posterior[,2,"e[2]"]),f = mean(posterior[,2,"f[2]"]),m = mean(posterior[,2,"m[2]"]),t_onset = mean(posterior[,2,"onset[2]"]), a_flr=mean(posterior[,2,"a_flr"]), plot=T, t_data = ts4_data$time, P_data = ts4_data$P, flower_data = ts4_data$taraxacum)

# chain 3
verhulst.wrapper(c(120:250),B0 = mean(posterior[,3,"B0[1]"]),e = mean(posterior[,3,"e[1]"]),f = mean(posterior[,3,"f[1]"]),m = mean(posterior[,3,"m[1]"]),t_onset = mean(posterior[,3,"onset[1]"]),a_flr = mean(posterior[,3,"a_flr"]),plot=T, t_data = ts1_data$time, P_data = ts1_data$P, flower_data = ts1_data$taraxacum)
verhulst.wrapper(c(120:250),B0 = mean(posterior[,3,"B0[2]"]),e = mean(posterior[,3,"e[2]"]),f = mean(posterior[,3,"f[2]"]),m = mean(posterior[,3,"m[2]"]),t_onset = mean(posterior[,3,"onset[2]"]), a_flr=mean(posterior[,3,"a_flr"]), plot=T, t_data = ts2_data$time, P_data = ts2_data$P, flower_data = ts2_data$taraxacum)
verhulst.wrapper(c(120:250),B0 = mean(posterior[,3,"B0[3]"]),e = mean(posterior[,3,"e[1]"]),f = mean(posterior[,3,"f[1]"]),m = mean(posterior[,3,"m[1]"]),t_onset = mean(posterior[,3,"onset[1]"]),a_flr = mean(posterior[,3,"a_flr"]),plot=T, t_data = ts3_data$time, P_data = ts3_data$P, flower_data = ts3_data$taraxacum)
verhulst.wrapper(c(120:250),B0 = mean(posterior[,3,"B0[4]"]),e = mean(posterior[,3,"e[2]"]),f = mean(posterior[,3,"f[2]"]),m = mean(posterior[,3,"m[2]"]),t_onset = mean(posterior[,3,"onset[2]"]), a_flr=mean(posterior[,3,"a_flr"]), plot=T, t_data = ts4_data$time, P_data = ts4_data$P, flower_data = ts4_data$taraxacum)

# chain 4
verhulst.wrapper(c(120:250),B0 = mean(posterior[,4,"B0[1]"]),e = mean(posterior[,4,"e[1]"]),f = mean(posterior[,4,"f[1]"]),m = mean(posterior[,4,"m[1]"]),t_onset = mean(posterior[,4,"onset[1]"]),a_flr = mean(posterior[,4,"a_flr"]),plot=T, t_data = ts1_data$time, P_data = ts1_data$P, flower_data = ts1_data$taraxacum)
verhulst.wrapper(c(120:250),B0 = mean(posterior[,4,"B0[2]"]),e = mean(posterior[,4,"e[2]"]),f = mean(posterior[,4,"f[2]"]),m = mean(posterior[,4,"m[2]"]),t_onset = mean(posterior[,4,"onset[2]"]), a_flr=mean(posterior[,4,"a_flr"]), plot=T, t_data = ts2_data$time, P_data = ts2_data$P, flower_data = ts2_data$taraxacum)
verhulst.wrapper(c(120:250),B0 = mean(posterior[,4,"B0[3]"]),e = mean(posterior[,4,"e[1]"]),f = mean(posterior[,4,"f[1]"]),m = mean(posterior[,4,"m[1]"]),t_onset = mean(posterior[,4,"onset[1]"]),a_flr = mean(posterior[,4,"a_flr"]),plot=T, t_data = ts3_data$time, P_data = ts3_data$P, flower_data = ts3_data$taraxacum)
verhulst.wrapper(c(120:250),B0 = mean(posterior[,4,"B0[4]"]),e = mean(posterior[,4,"e[2]"]),f = mean(posterior[,4,"f[2]"]),m = mean(posterior[,4,"m[2]"]),t_onset = mean(posterior[,4,"onset[2]"]), a_flr=mean(posterior[,4,"a_flr"]), plot=T, t_data = ts4_data$time, P_data = ts4_data$P, flower_data = ts4_data$taraxacum)



###############
### figures ###
###############

### pretty timeseries figure
trail_2019 <- verhulst.wrapper(c(120:250),B0 = mean(c(posterior[,,"B0[1]"])),e = mean(c(posterior[,,"e[1]"])),f = mean(c(posterior[,,"f[1]"])),m = mean(c(posterior[,,"m[1]"])),t_onset = mean(c(posterior[,,"onset[1]"])), a_flr = mean(c(posterior[,,"a_flr"])),plot=T, t_data = ts1_data$time, P_data = ts1_data$P, flower_data = ts1_data$taraxacum)
trail_2021 <- verhulst.wrapper(c(120:250),B0 = mean(c(posterior[,,"B0[2]"])),e = mean(c(posterior[,,"e[2]"])),f = mean(c(posterior[,,"f[2]"])),m = mean(c(posterior[,,"m[2]"])),t_onset = mean(c(posterior[,,"onset[2]"])), a_flr = mean(c(posterior[,,"a_flr"])), plot=T, t_data = ts2_data$time, P_data = ts2_data$P, flower_data = ts2_data$taraxacum)
waterfall_2019 <- verhulst.wrapper(c(120:250),B0 = mean(c(posterior[,,"B0[3]"])),e = mean(c(posterior[,,"e[1]"])),f = mean(c(posterior[,,"f[1]"])),m = mean(c(posterior[,,"m[1]"])),t_onset = mean(c(posterior[,,"onset[1]"])), a_flr = mean(c(posterior[,,"a_flr"])),plot=T, t_data = ts3_data$time, P_data = ts3_data$P, flower_data = ts3_data$taraxacum)
waterfall_2021 <- verhulst.wrapper(c(120:250),B0 = mean(c(posterior[,,"B0[4]"])),e = mean(c(posterior[,,"e[2]"])),f = mean(c(posterior[,,"f[2]"])),m = mean(c(posterior[,,"m[2]"])),t_onset = mean(c(posterior[,,"onset[2]"])), a_flr = mean(c(posterior[,,"a_flr"])), plot=T, t_data = ts4_data$time, P_data = ts4_data$P, flower_data = ts4_data$taraxacum)


col_2019 <- "#dd5d2e"
col_2021 <- "#073f78"

# takes output of verhulst.wrapper, outputs pretty text for figure
print.pars <- function(x){
  e <- x$e; m <- x$m; f <- x$f; B0 <- x$B0; onset <- x$t_onset
  text <- paste0("emergence: ", round(e,4), "   \n",
                 "friction: ", round(f,4), "   \n",
                 "mortality: ", round(m,4), "   \n",
                 "population: ", round(B0), "   \n",
                 "onset DOY: ", round(onset), "   \n")
  return(text)
}

png("figures/fig1.png", width=8, height=7, units="in", res=300)
# svg("fig1.svg", width=8, height=7)

  B_max_trail <- max(c(trail_2019$B, trail_2021$B), na.rm=T) # na.rm=T because there are some at the end of the sim
  B_max_waterfall <- max(c(waterfall_2019$B, waterfall_2021$B), na.rm=T)

  par(mfrow=c(2,1))
  par(mar=c(0,4,3,4), oma=c(0,0,0,0))
  
  # Trail P curves
  plot(NA, xlim=c(120,250), ylim=c(0,8), xlab="", ylab="", xaxt="n")
  with(trail_2019, {
    lines(P ~ t, col=col_2019, lwd=3)
    polygon(x = c(t,rev(t)), y = c(P,rep(0,length(P))),
            border=NA, col=adjustcolor(col_2019, 0.25))
    arrows(x0 = t_data, y0 = P[which(t %in% t_data)], y1 = P_covars, col=col_2019, lwd=2, code=0)
    points(P_data ~ t_data, pch=19, col=col_2019)
  })
  with(trail_2021, {

    lines(P ~ t, col=col_2021, lwd=3)
    polygon(x = c(t,rev(t)), y = c(P,rep(0,length(P))),
            border=NA, col=adjustcolor(col_2021, 0.25))
    arrows(x0 = t_data, y0 = P[which(t %in% t_data)], y1 = P_covars, col=col_2021, lwd=2, code=0)
    points(P_data ~ t_data, pch=19, col=col_2021)
  })
  
  # Trail B overlays
  par(new=T)
  with(trail_2019, {
    plot(B ~ t, type="l", col=col_2019, lty=2, lwd=2,
         xlim=c(120,250), ylim=c(0,B_max_trail), xlab="", ylab="", axes=F, bty="n")
  })
  par(new=T)
  with(trail_2021, {
    plot(B ~ t, type="l", col=col_2021, lty=2, lwd=2,
         xlim=c(120,250), ylim=c(0,B_max_trail), xlab="", ylab="", axes=F, bty="n")
  })
  axis(4, pretty(c(0,B_max_trail)))
  legend("left", inset=0.025, lwd=c(3,3,2,2), lty=c(1,1,1,2), col=c(col_2019, col_2021, "black", "black"), legend=c("2019", "2021", "Emerged", "Unemerged"), cex=0.75, bty="n")
  mtext(print.pars(trail_2019), line=-5, adj=1, cex=0.75, col=col_2019)
  mtext(print.pars(trail_2021), line=-10, adj=1, cex=0.75, col=col_2021)
  
  mtext("Trail site", side=3, line=-1.25)
  
  par(mar=c(4,4,0,4))
  
  # Waterfall P curves
  plot(NA, xlim=c(120,250), ylim=c(0,6), xlab="Day of year", ylab="", xaxt="n")
  doy_date <- data.frame(doy = c(120:250), date = as.Date(c(120:250), origin=c("2019-01-01")))
  doy_date$month <- month.abb[month(doy_date$date)]; doy_date$day <- day(doy_date$date)
  doy_date$day_month <- paste(doy_date$day, doy_date$month)
  doy_date_axis <- doy_date[which(doy_date$day %in% c(1,11,21)),]
  axis(1, at = doy_date_axis$doy, labels = doy_date_axis$day_month, las=2)
  
  with(waterfall_2019, {
    lines(P ~ t, col=col_2019, lwd=3)
    polygon(x = c(t,rev(t)), y = c(P,rep(0,length(P))),
            border=NA, col=adjustcolor(col_2019, 0.25))
    arrows(x0 = t_data, y0 = P[which(t %in% t_data)], y1 = P_covars, col=col_2019, lwd=2, code=0)
    points(P_data ~ t_data, pch=19, col=col_2019)
  })
  with(waterfall_2021, {
    lines(P ~ t, col=col_2021, lwd=3)
    polygon(x = c(t,rev(t)), y = c(P,rep(0,length(P))),
            border=NA, col=adjustcolor(col_2021, 0.25))
    arrows(x0 = t_data, y0 = P[which(t %in% t_data)], y1 = P_covars, col=col_2021, lwd=2, code=0)
    points(P_data ~ t_data, pch=19, col=col_2021)
  })

  # Waterfall B overlays
  par(new=T)
  with(waterfall_2019, {
    plot(B ~ t, type="l", col=col_2019, lty=2, lwd=2,
         xlim=c(120,250), ylim=c(0,B_max_waterfall), xlab="", ylab="", axes=F, bty="n")
  })
  par(new=T)
  with(waterfall_2021, {
    plot(B ~ t, type="l", col=col_2021, lty=2, lwd=2,
         xlim=c(120,250), ylim=c(0,B_max_waterfall), xlab="", ylab="", axes=F, bty="n")
  })
  axis(4, pretty(c(0,B_max_waterfall)))
  mtext(print.pars(waterfall_2019), line=-5, adj=1, cex=0.75, col=col_2019)
  mtext(print.pars(waterfall_2021), line=-10, adj=1, cex=0.75, col=col_2021)
  
  mtext("Waterfall site", side=3, line=-1.25)
  
  # axis labels
  par(new=T)
  #par(fig=c(0,1,0,0.75))
  par(mfrow=c(1,1), mar=c(4.5,4.1,1.5,1.5))
  #par(mar=c(0,4,6,4))
  mtext("Unemerged abundance", 4, line=2.5) 
  mtext("Emerged abundance", 2, line=2.5) 
  
dev.off()


### m vs. B0 figure - posterior correlations
# e/m/f: odds are 2019, evens are 2021
# B0s: 1 = trail 2019, 2 = trail 2021, 3 = waterfall 2019, 4 = waterfall 2021

e_pal <- colorRampPalette(c("#ffe369", "#9e1000"))
B0s <- c(posterior[,,"B0[1]"])
ms <- c(posterior[,,"m[1]"])
es <- c(posterior[,,"e[1]"])
fs <- c(posterior[,,"f[1]"])
onsets <- c(posterior[,,"onset[1]"])
es_order <- findInterval(es, sort(es))


png("figures/fig2.png", width=8, height=5, units="in", res=175)
# svg("figures/fig2_less_points.svg", width=8, height=5)

  xlim <- c(0,1); ylim <- c(0,80)
  
  layout(matrix(c(1,1,1,2,3,4), nrow=3, ncol=2), widths=c(3,1), heights=c(1,1,1))
  
  par(mar=c(4.5,4.1,1.5,0))
  par(mar=c(0,0,0,0), oma=c(4.5,4.5,1.5,4.5))
  
  post_sample <- seq(1, dim(posterior)[1]) # for all posterior samples
  post_sample <- sample(seq(1, dim(posterior)[1]), 1000) # for subsample of posterior
  
  # main panel
  plot(posterior[,,"B0[1]"] ~ posterior[,,"m[1]"],
       pch=19, cex=0.3, col=adjustcolor(col_2019, 0.25),
       xlim=xlim, ylim=ylim)
  viz_model <- lm(B0s ~ ms + c(ms^2) + c(ms^3)) # for visualization only
  ms_pred <- seq(min(ms),max(ms),length.out=1000)
  viz_pred <- predict(viz_model, data.frame(ms=ms_pred))
  lines(viz_pred ~ ms_pred, lwd=3, col=adjustcolor(col_2019, offset = c(-0.2,-0.2,-0.2,0)))
  mtext("   2019 Trail site", side=3, adj=0, line=-2)
  
  # panel 2
  plot(posterior[post_sample,,"B0[2]"] ~ posterior[post_sample,,"m[2]"], cex=0.2, col=adjustcolor(col_2021, 0.25), xlim=xlim, ylim=ylim, xaxt="n", yaxt="n", ylab="", xlab="")
  axis(4, at=seq(0,75, length.out=4), labels=seq(0,75, length.out=4))
  mtext("   2019 Waterfall site", side=3, adj=0, line=-1.75, cex=0.7)
  
  # panel 3
  plot(posterior[post_sample,,"B0[3]"] ~ posterior[post_sample,,"m[1]"], cex=0.2, col=adjustcolor(col_2019, 0.25), xlim=xlim, ylim=ylim, xaxt="n", yaxt="n", ylab="", xlab="")
  axis(4, at=seq(0,75, length.out=4), labels=seq(0,75, length.out=4))
  mtext("   2021 Trail site", side=3, adj=0, line=-1.75, cex=0.7)
  
  # panel 4
  plot(posterior[post_sample,,"B0[4]"] ~ posterior[post_sample,,"m[2]"], cex=0.2, col=adjustcolor(col_2021, 0.25), xlim=xlim, ylim=ylim, xaxt="n", yaxt="n", ylab="", xlab="")
  axis(4, at=seq(0,75, length.out=4), labels=seq(0,75, length.out=4))
  axis(1, at=seq(0,1.2, length.out=4), labels=seq(0,1.2, length.out=4))
  mtext("   2021 Waterfall site", side=3, adj=0, line=-1.75, cex=0.7)
  
  # axes labels
  mtext("Senescence rate", 1, line=2.75, outer=T) 
  mtext("Population size", 2, line=2.75, outer=T) 
  mtext("Population size", 4, line=2.75, outer=T) 
  
  # new plotting window for bubbles
  par(fig = c(0,1,0,1), mar=c(0,0,0,0), oma=c(4.5,4.5,1.5,4.5), new=T)
  plot.new()
  plot.window(xlim=xlim*c(0,1.40), ylim=ylim)
  
  # parameter points to demonstrate
  m_demo <- c(0.15, 0.3, 0.6)
  B0_demo <- predict(viz_model, data.frame(ms=m_demo))
  
  # lines to bubbles
  bubble_coords_x <- c(0.165, 0.58, .81)
  bubble_coords_y <- c(45, 63, 11)
  segments(x0 = m_demo, x1 = bubble_coords_x, y0 = B0_demo, y1 = bubble_coords_y)
  
  # drawing bubbles
  draw.circle(bubble_coords_x[1], bubble_coords_y[1], .125, col="#f0f3f7")
  draw.circle(bubble_coords_x[2], bubble_coords_y[2], .125, col="#f0f3f7")
  draw.circle(bubble_coords_x[3], bubble_coords_y[3], .125, col="#f0f3f7")
  
  
  # bubble subplots with curve
  curve.plot <- function(x1, x2, y1, y2,
                         t = c(120:250), B0, e, f, m, t_onset,
                         col = "red", lwd=3){
    
    if(FALSE){
      t = c(120:250); B0 = 25; e = 0.2; f = 0.004; m = 0.3; t_onset = 150
    }
    
    par(fig = c(0,1,0,1), mar=c(0,0,0,0), oma=c(4.5,4.5,1.5,4.5))
    par(fig = c(x1,x2,y1,y2), mar=c(1.5,0,1.5,0), new = T, xpd=TRUE) 
    
    curve <- verhulst.wrapper(t, B0, e, f, m, t_onset, plot=F)
    plot(P ~ t, data = curve, type="l", col=col, lwd=lwd, axes=F, xlab="", ylab="")
    polygon(x = c(curve$t,rev(curve$t)), y = c(curve$P,rep(0,length(curve$P))),
            border=NA, col=adjustcolor(col, 0.25))
    #print(max(curve$P))
    text(curve$t[which(curve$P == max(curve$P))], max(curve$P)*1, paste("           —", round(max(curve$P), 1)), cex=1.25)
    t_onset_date <- doy_date[which(doy_date$doy == round(t_onset)), "day_month"]
    text(round(t_onset), max(curve$P)*-0.15, paste("|\n          ", t_onset_date), cex=1.25)
    #points(median(curve$t), mean(range(curve$P, na.rm=T)), cex=20) # one big circle
    
  }
  
  # demo parameters
  demo_model_e <- lm(es ~ ms + B0s)
  demo_model_f <- lm(fs ~ ms + B0s)
  demo_model_onset <- lm(onsets ~ ms + B0s)
  e_demo <- predict(demo_model_e, data.frame(ms = m_demo, B0s = B0_demo))
  f_demo <- predict(demo_model_f, data.frame(ms = m_demo, B0s = B0_demo))
  onset_demo <- predict(demo_model_onset, data.frame(ms = m_demo, B0s = B0_demo))
  
  0.35; curve2_y <- 0.66
  curve3_x <- 0.5; curve3_y <- 0.06
  curve_width <- 0.15
  curve_height <- 0.225
  
  # drawing curves
  curve.plot(x1=curve1_x, x2=curve1_x+curve_width, y1=curve1_y, y2=curve1_y+curve_height,
             t = c(170:225), B0=B0_demo[1], e=e_demo[1], f=f_demo[1], m=m_demo[1], t_onset=onset_demo[1],
             col = col_2019, lwd=3)
  
  curve.plot(x1=curve2_x, x2=curve2_x+curve_width, y1=curve2_y, y2=curve2_y+curve_height,
             t =  c(170:225), B0=B0_demo[2], e=e_demo[2], f=f_demo[2], m=m_demo[2], t_onset=onset_demo[2],
             col = col_2019, lwd=3)
  
  curve.plot(x1=curve3_x, x2=curve3_x+curve_width, y1=curve3_y, y2=curve3_y+curve_height,
             t =  c(170:225), B0=B0_demo[3], e=e_demo[3], f=f_demo[3], m=m_demo[3], t_onset=onset_demo[3],
             col = col_2019, lwd=3)

dev.off()



### parameter scans to show that the curve is flexible

cols <- colorRampPalette(c("#ffd000","#02d9ce"))(10)

es <- seq(0.1,0.5,length.out=10)
fs <- seq(0,0.0078,length.out=10)
ms <- seq(0,0.5, length.out=10)
B0s <- seq(2,45, length.out=10)
onsets <- seq(105,123, by=2)

svg("figures/par_scan.svg", width=8, height=8)
#png("figures/par_scan.png", width=600, height=600)

  layout(matrix(c(1,2,5,3,4,5),nrow=3,ncol=2), widths=c(1,1,0.5))
  par(mar=c(2,4,1,2))
  
  plot(NA, xlim=c(100,150), ylim=c(0,15), xlab="", ylab="", yaxt="n", xaxt="n")
  for(i in seq_along(es)){
    e <- es[i]
    col <- cols[i]
    sim <-  verhulst.wrapper(c(100:150), B0=25, e=e, f=0.004, m=0.3, t_onset=108)
    lines(P ~ t, data=sim, col=col, lwd=2)
  }
  legend("topright", inset=0.05, legend=round(es, 2), col=cols, lwd=2, bty="n", cex=0.75)
  mtext("Emergence rate (e)", line=-2, cex=0.9)
  
  plot(NA, xlim=c(100,150), ylim=c(0,9), xlab="", ylab="", yaxt="n", xaxt="n")
  for(i in seq_along(fs)){
    f <- fs[i]
    col <- cols[i]
    sim <-  verhulst.wrapper(c(100:150), B0=25, e=0.2, f=f, m=0.3, t_onset=108)
    lines(P ~ t, data=sim, col=col, lwd=2)
  }
  legend("topright", inset=0.05, legend=round(fs, 5), col=cols, lwd=2, bty="n", cex=0.75)
  mtext("Emergence friction (f)", line=-2, cex=0.9)
  mtext("Abundance", side=2, line=1.8, cex=1)
  
  plot(NA, xlim=c(100,150), ylim=c(0,28), xlab="", ylab="", yaxt="n", xaxt="n")
  for(i in seq_along(ms)){
    m <- ms[i]
    col <- cols[i]
    sim <-  verhulst.wrapper(c(100:150), B0=25, e=0.2, f=0.004, m=m, t_onset=108)
    lines(P ~ t, data=sim, col=col, lwd=2)
  }
  legend("topleft", inset=0.05, legend=round(ms, 2), col=cols, lwd=2, bty="n", cex=0.75)
  mtext("Mortality rate (m)", line=-2, cex=0.9)
  
  plot(NA, xlim=c(100,150), ylim=c(0,9), xlab="", ylab="", yaxt="n", xaxt="n")
  for(i in seq_along(B0s)){
    B0 <- B0s[i]
    col <- cols[i]
    sim <-  verhulst.wrapper(c(100:150), B0=B0, e=0.2, f=0.004, m=0.3, t_onset=108)
    lines(P ~ t, data=sim, col=col, lwd=2)
  }
  legend("topright", inset=0.05, legend=round(B0s, 0), col=cols, lwd=2, bty="n", cex=0.75)
  mtext("Population size (u)", line=-2, cex=0.9)
  
  par(mar=c(3,17,0,17))
  
  plot(NA, xlim=c(100,150), ylim=c(0,8), xlab="", ylab="", yaxt="n", xaxt="n")
  for(i in seq_along(onsets)){
    onset <- onsets[i]
    col <- cols[i]
    sim <-  verhulst.wrapper(c(100:150), B0=25, e=0.2, f=0.004, m=0.3, t_onset=onset)
    lines(P ~ t, data=sim, col=col, lwd=2)
  }
  legend("topright", inset=0.05, legend=round(onsets, 0), col=cols, lwd=2, bty="n", cex=0.75)
  mtext("Onset phenology (h)", line=-2, cex=0.9)
  mtext("Time", side=1, line=1.8, cex=1)

dev.off()
