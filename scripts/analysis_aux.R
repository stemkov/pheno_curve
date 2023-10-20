### running auxiliary analyses

setwd("/home/michael/Documents/Grad School/Research Projects/abundance_curve")

# bees <- fread("clean_data/clean_bee_data.csv", stringsAsFactors = FALSE)
# flowers <- fread("clean_data/clean_flower_data.csv", stringsAsFactors = FALSE)
# effort <- fread("clean_data/clean_sampling_effort.csv", stringsAsFactors = FALSE)
# locations <- fread("clean_data/clean_plot_locations.csv", stringsAsFactors = FALSE)
# temp <- fread("clean_data/temp_data.csv")

sites <- c("401 Trail", "Waterfall", "401 Trail", "Waterfall")
years <- c(2019, 2019, 2021, 2021)
all_data_list <- mapply(assemble.data, sites, years, SIMPLIFY = F, USE.NAMES = F)
all_data <- rbindlist(all_data_list)
all_data[, date := as.Date(time, origin = paste0(year, "-01-01"))]
all_data[, site_year_label := as.factor(paste(site, year, sep="_"))]
all_data[, site_year := as.numeric(site_year_label)]
trail_2019 <- all_data[site=="401 Trail" & year == 2019,]
waterfall_2019 <- all_data[site=="Waterfall" & year == 2019,]
trail_2021 <- all_data[site=="401 Trail" & year == 2021,]
waterfall_2021 <- all_data[site=="Waterfall" & year == 2021,]

### summary statistics

# start & end dates
effort[year == 2019 & Site == "401 Trail", Date]
effort[year == 2019 & Site == "Waterfall", Date]
effort[year == 2021 & Site == "401 Trail", Date]
effort[year == 2021 & Site == "Waterfall", Date]

# total effort
sum(sapply(get.sampling.dates("401 Trail", 2019), function(x) get.effort(site="401 Trail", date=x, year=2019)))
sum(sapply(get.sampling.dates("Waterfall", 2019), function(x) get.effort(site="Waterfall", date=x, year=2019)))
sum(sapply(get.sampling.dates("401 Trail", 2021), function(x) get.effort(site="401 Trail", date=x, year=2021)))
sum(sapply(get.sampling.dates("Waterfall", 2021), function(x) get.effort(site="Waterfall", date=x, year=2021)))

# number of transects walked per site per year
sum(sapply(get.sampling.dates("401 Trail", 2019), function(x) get.effort(site="401 Trail", date=x, year=2019))/10)
sum(sapply(get.sampling.dates("Waterfall", 2019), function(x) get.effort(site="Waterfall", date=x, year=2019))/10)
sum(sapply(get.sampling.dates("401 Trail", 2021), function(x) get.effort(site="401 Trail", date=x, year=2021))/10)
sum(sapply(get.sampling.dates("Waterfall", 2021), function(x) get.effort(site="Waterfall", date=x, year=2021))/10)


# total H rub caught
paste("Total hrub at Trail 2019:",sum(all_data_list[[1]]$P))
paste("Total hrub at Waterfall 2019:",sum(all_data_list[[2]]$P))
paste("Total hrub at Trail 2021:",sum(all_data_list[[3]]$P))
paste("Total hrub at Waterfall 2021:",sum(all_data_list[[4]]$P))

paste("Total hrub overall:" , sum(sapply(1:4, function(x) sum(all_data_list[[x]]$P))),
      ", across unique days:", sum(sapply(1:4, function(x) length(which(all_data_list[[x]]$P != 0)))))

# total T officinale counted
paste("Total Toffic at Trail 2019:",sum(all_data_list[[1]]$taraxacum))
paste("Total hrub at Waterfall 2019:",sum(all_data_list[[2]]$taraxacum))
paste("Total hrub at Trail 2021:",sum(all_data_list[[3]]$taraxacum))
paste("Total hrub at Waterfall 2021:",sum(all_data_list[[4]]$taraxacum))


# coarse look at phenology onset
trail_2019[P != 0,][1, date]
waterfall_2019[P != 0,][1, date]
trail_2021[P != 0,][1, date]
waterfall_2021[P != 0,][1, date]

trail_2019[taraxacum != 0,][1, date]
waterfall_2019[taraxacum != 0,][1, date]
trail_2021[taraxacum != 0,][1, date]
waterfall_2021[taraxacum != 0,][1, date]

# coarse look at phenology end
trail_2019[P != 0,][,date]
waterfall_2019[P != 0,][,date]
trail_2021[P != 0,][,date]
waterfall_2021[P != 0,][,date]

trail_2019[taraxacum != 0,][,date]
waterfall_2019[taraxacum != 0,][,date]
trail_2021[taraxacum != 0,][,date]
waterfall_2021[taraxacum != 0,][,date]

par(mfrow=c(2,1))

# onset parameter estimates
hist(posterior[,,"onset[1]"], 100)
abline(v=mean(posterior[,,"onset[1]"]), col="red", lwd=2)
mean(posterior[,,"onset[1]"])
quantile(posterior[,,"onset[1]"], c(0.025, 0.975))
as.Date(mean(posterior[,,"onset[1]"]), origin = "2019-01-01")
hist(posterior[,,"onset[2]"], 100)
abline(v=mean(posterior[,,"onset[2]"]), col="red", lwd=2)
mean(posterior[,,"onset[2]"])
quantile(posterior[,,"onset[2]"], c(0.025, 0.975))
as.Date(mean(posterior[,,"onset[2]"]), origin = "2021-01-01")

# B0 parameter estimates
hist(posterior[,,"B0[1]"], 100)
abline(v=mean(posterior[,,"B0[1]"]), col="red", lwd=2)
mean(posterior[,,"B0[1]"])
quantile(posterior[,,"B0[1]"], c(0.025, 0.975))
hist(posterior[,,"B0[2]"], 100)
abline(v=mean(posterior[,,"B0[2]"]), col="red", lwd=2)
mean(posterior[,,"B0[2]"])
quantile(posterior[,,"B0[2]"], c(0.025, 0.975))
hist(posterior[,,"B0[3]"], 100)
abline(v=mean(posterior[,,"B0[3]"]), col="red", lwd=2)
mean(posterior[,,"B0[3]"])
quantile(posterior[,,"B0[3]"], c(0.025, 0.975))
hist(posterior[,,"B0[4]"], 100)
abline(v=mean(posterior[,,"B0[4]"]), col="red", lwd=2)
mean(posterior[,,"B0[4]"])
quantile(posterior[,,"B0[4]"], c(0.025, 0.975))



# emergence parameter estimates
hist(posterior[,,"e[1]"], 100)
abline(v=mean(posterior[,,"e[1]"]), col="red", lwd=2)
mean(posterior[,,"e[1]"])
quantile(posterior[,,"e[1]"], c(0.025, 0.975))
hist(posterior[,,"e[2]"], 100)
abline(v=mean(posterior[,,"e[2]"]), col="red", lwd=2)
mean(posterior[,,"e[2]"])
quantile(posterior[,,"e[2]"], c(0.025, 0.975))

# mortality parameter estimates
hist(posterior[,,"m[1]"], 100)
abline(v=mean(posterior[,,"m[1]"]), col="red", lwd=2)
mean(posterior[,,"m[1]"])
quantile(posterior[,,"m[1]"], c(0.025, 0.975))
hist(posterior[,,"m[2]"], 100)
abline(v=mean(posterior[,,"m[2]"]), col="red", lwd=2)
mean(posterior[,,"m[2]"])
quantile(posterior[,,"m[2]"], c(0.025, 0.975))

# mortality parameter estimates
hist(posterior[,,"f[1]"], 100)
abline(v=mean(posterior[,,"f[1]"]), col="red", lwd=2)
mean(posterior[,,"f[1]"])
quantile(posterior[,,"f[1]"], c(0.025, 0.975))
hist(posterior[,,"f[2]"], 100)
abline(v=mean(posterior[,,"f[2]"]), col="red", lwd=2)
mean(posterior[,,"f[2]"])
quantile(posterior[,,"f[1]"], c(0.025, 0.975))

# flower covariate
hist(posterior[,,"a_flr"], 100)
abline(v=mean(posterior[,,"a_flr"]), col="red", lwd=2)
mean(posterior[,,"a_flr"])
quantile(posterior[,,"a_flr"], c(0.025, 0.975))
sum(posterior[,,"a_flr"] > 0) / length(posterior[,,"a_flr"])


# hypothesis tests
hist(posterior[,,"e[1]"], 100, col=adjustcolor("red", 0.25))
hist(posterior[,,"e[2]"], 100, col=adjustcolor("blue", 0.25), add=T)
hist(posterior[,,"f[1]"], 100, col=adjustcolor("red", 0.25))
hist(posterior[,,"f[2]"], 100, col=adjustcolor("blue", 0.25), add=T)
hist(posterior[,,"m[1]"], 100, col=adjustcolor("red", 0.25))
hist(posterior[,,"m[2]"], 100, col=adjustcolor("blue", 0.25), add=T)

e_dif <- posterior[,,"e[1]"] - posterior[,,"e[2]"]
hist(e_dif, 100)
sum(e_dif > 0) / length(e_dif)

f_dif <- posterior[,,"f[1]"] - posterior[,,"f[2]"]
hist(f_dif, 100)
sum(f_dif > 0) / length(f_dif)

m_dif <- posterior[,,"m[1]"] - posterior[,,"m[2]"]
hist(m_dif, 100)
sum(m_dif > 0) / length(m_dif)

h_dif <- posterior[,,"onset[1]"] - posterior[,,"onset[2]"]
hist(h_dif, 100)
sum(h_dif > 0) / length(h_dif) # proportion of onset differences greater than 0


### predictions from model
ts1_data <- all_data[site_year==1,]
ts2_data <- all_data[site_year==2,]
ts3_data <- all_data[site_year==3,]
ts4_data <- all_data[site_year==4,]

trail_2019 <- verhulst.wrapper(c(120:250),B0 = mean(c(posterior[,,"B0[1]"])),e = mean(c(posterior[,,"e[1]"])),f = mean(c(posterior[,,"f[1]"])),m = mean(c(posterior[,,"m[1]"])),t_onset = mean(c(posterior[,,"onset[1]"])), a_flr = mean(c(posterior[,,"a_flr"])),plot=T, t_data = ts1_data$time, P_data = ts1_data$P, flower_data = ts1_data$taraxacum)
trail_2021 <- verhulst.wrapper(c(120:250),B0 = mean(c(posterior[,,"B0[2]"])),e = mean(c(posterior[,,"e[2]"])),f = mean(c(posterior[,,"f[2]"])),m = mean(c(posterior[,,"m[2]"])),t_onset = mean(c(posterior[,,"onset[2]"])), a_flr = mean(c(posterior[,,"a_flr"])), plot=T, t_data = ts2_data$time, P_data = ts2_data$P, flower_data = ts2_data$taraxacum)
waterfall_2019 <- verhulst.wrapper(c(120:250),B0 = mean(c(posterior[,,"B0[3]"])),e = mean(c(posterior[,,"e[1]"])),f = mean(c(posterior[,,"f[1]"])),m = mean(c(posterior[,,"m[1]"])),t_onset = mean(c(posterior[,,"onset[1]"])), a_flr = mean(c(posterior[,,"a_flr"])),plot=T, t_data = ts3_data$time, P_data = ts3_data$P, flower_data = ts3_data$taraxacum)
waterfall_2021 <- verhulst.wrapper(c(120:250),B0 = mean(c(posterior[,,"B0[4]"])),e = mean(c(posterior[,,"e[2]"])),f = mean(c(posterior[,,"f[2]"])),m = mean(c(posterior[,,"m[2]"])),t_onset = mean(c(posterior[,,"onset[2]"])), a_flr = mean(c(posterior[,,"a_flr"])), plot=T, t_data = ts4_data$time, P_data = ts4_data$P, flower_data = ts4_data$taraxacum)

# peak dates
doy_max_2019t <- trail_2019$t[which.max(trail_2019$P)]
doy_max_2019w <- trail_2019$t[which.max(waterfall_2019$P)]
doy_max_2021t <- trail_2019$t[which.max(trail_2021$P)]
doy_max_2021w <- trail_2019$t[which.max(waterfall_2021$P)]

as.Date(doy_max_2019t, origin="2019-01-01")
as.Date(doy_max_2019w, origin="2019-01-01")
as.Date(doy_max_2021t, origin="2021-01-01")
as.Date(doy_max_2021w, origin="2021-01-01")

# MAEs
with(trail_2019, mean(abs(P_data - P[which(t %in% t_data)])))
with(waterfall_2019, mean(abs(P_data - P[which(t %in% t_data)])))
with(trail_2021, mean(abs(P_data - P[which(t %in% t_data)])))
with(waterfall_2021, mean(abs(P_data - P[which(t %in% t_data)])))




