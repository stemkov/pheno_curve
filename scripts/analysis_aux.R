### running auxiliary analyses

setwd(wd_path)

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


# flowering plant species richness
flowers <- fread("clean_data/clean_flower_data.csv", stringsAsFactors = FALSE)
length(unique(flowers[Site == "401 Trail" & year == 2019 , Species]))
length(unique(flowers[Site == "401 Trail" & year == 2021 , Species]))
length(unique(flowers[Site == "Waterfall" & year == 2019 , Species]))
length(unique(flowers[Site == "Waterfall" & year == 2021 , Species]))

# june 15 - July 15 temperatures
mean(temp[site == "401 Trail" & date >= ymd("2019-06-15") & date <= ymd("2019-07-15"), temp_c], na.rm=T)
mean(temp[site == "401 Trail" & date >= ymd("2021-06-15") & date <= ymd("2021-07-15"), temp_c], na.rm=T)
mean(temp[site == "Waterfall" & date >= ymd("2019-06-15") & date <= ymd("2019-07-15"), temp_c], na.rm=T)
mean(temp[site == "Waterfall" & date >= ymd("2021-06-15") & date <= ymd("2021-07-15"), temp_c], na.rm=T)

# total bees - not just Hrub
bees <- fread("clean_data/clean_bee_data.csv", stringsAsFactors = FALSE)
bees[Site == "401 Trail" & year == 2019 & Field.ID %!in% c("Halictus rubicundus", "none"), .N]
bees[Site == "401 Trail" & year == 2021 & Field.ID %!in% c("Halictus rubicundus", "none"), .N]
bees[Site == "Waterfall" & year == 2019 & Field.ID %!in% c("Halictus rubicundus", "none"), .N]
bees[Site == "Waterfall" & year == 2021 & Field.ID %!in% c("Halictus rubicundus", "none"), .N]

unique(bees[Site == "401 Trail" & year == 2019 & Field.ID %!in% c("Halictus rubicundus", "none"), gsub(" .*", "", Field.ID)]) # 12
unique(bees[Site == "401 Trail" & year == 2021 & Field.ID %!in% c("Halictus rubicundus", "none"), gsub(" .*", "", Field.ID)]) # 9
unique(bees[Site == "Waterfall" & year == 2019 & Field.ID %!in% c("Halictus rubicundus", "none"), gsub(" .*", "", Field.ID)]) # 17
unique(bees[Site == "Waterfall" & year == 2021 & Field.ID %!in% c("Halictus rubicundus", "none"), gsub(" .*", "", Field.ID)]) # 8


### parameter summary table
names(posterior[1,1,])

posterior$parameters
str(posterior)

parameter_names <- names(posterior[1,1,])

# Calculate summary statistics (e.g., mean, median, quantiles)
summary_stats <- sapply(parameter_names, function(param) {
  quantiles <- quantile(posterior[,,param], c(0.025, 0.5, 0.975))
  names(quantiles) <- c("2.5%", "Median", "97.5%")
  return(quantiles)
})

# Convert summary statistics to data frame
summary_df <- as.data.frame(t(summary_stats))

# Add parameter names as a column
summary_df$Parameter <- parameter_names

summary_to_output <- summary_df[c("onset[1]", "onset[2]",
                                  "B0[1]", "B0[2]", "B0[3]", "B0[4]",
                                  "e[1]", "e[2]",
                                  "f[1]", "f[2]",
                                  "m[1]", "m[2]",
                                  "a_flr", "sigma"),
                                c(4,2,1,3)]
summary_to_output$Parameter <- c("h_2019", "h_2021",
                                 "u_Trail,2019", "u_Trail,2021", "u_Waterfall,2019", "u_Waterfall,2021",
                                 "a_2019", "a_2021",
                                 "b_2019", "b_2021",
                                 "m_2019", "m_2021",
                                 "beta_F", "sigma")


library(flextable)
ft <- flextable(summary_to_output)
save_as_html(ft, path = "summary_table.html")


### showing how traditional phenophase estimation methods can give misleading results, Supplement 2

library(mgcv)

set.seed(0)

t <- c(100:200)
B0 <- 100
e <- .15
f <- 0.0001
m <- 0.1
t_onset <- 110

test.method <- function(t, B0, e, f, m, t_onset, plot=F){
  
  sample_i <- round(seq(runif(1,1,20), 100, length.out=10))
  sample_t <- t[sample_i]
  
  sim_data <- verhulst.wrapper(t=t, B0=B0, e=e, f=f, m=m, t_onset=t_onset, plot=plot)
  P_obs <- round(sim_data$P[sample_i] + rnorm(10, 0, 2))
  
  first_est <- t[sample_i[which(P_obs > 0)[1]]]
  quant_est <- quantile(rep(sample_t, time = pmax(P_obs,0)),.05)
  gam_model <- gam(pmax(P_obs,0) ~ s(sample_t), family = nb())
  gam_pred <- predict(gam_model, newdata = data.frame(sample_t=c(50:250)))
  gam_est_0 <- c(50:250)[which(gam_pred >= 1)[1]]
  gam_est_05 <- c(50:250)[which(gam_pred >= max(gam_pred)*0.05)[1]]
  
  return(c(t_onset = t_onset, first_est = first_est, quant_est = as.numeric(quant_est), gam_est_0 = gam_est_0, gam_est_05 = gam_est_05))
}

if(FALSE){
  test.method(t, B0, e=1, f, m, t_onset, plot=T) # demonstrate single run
}

t_onsets <- seq(100,150)
es <- runif(length(t_onsets), .03, 1)
tests <- mapply(test.method, t_onset = t_onsets, e = es, MoreArgs = list(t = t, B0 = B0, f = f, m = m))
tests_df <- data.table(t(tests))

library(RColorBrewer)
pal <- brewer.pal(5,"Dark2")

#points(gam_est_05 ~ t_onset, data=tests_df, pch=20, col = adjustcolor(pal[4], 0.5))
#abline(lm(gam_est_0 ~ t_onset, data=tests_df), lwd=1, col=pal[4])

model_for_test <- stan_model(file='scripts/model_for_testing.stan')

test.verhulst <- function(t, B0, e, f, m, t_onset){
  sample_i <- round(seq(runif(1,1,20), 100, length.out=10))
  sample_t <- t[sample_i]
  sim_data <- verhulst.wrapper(t=t, B0=B0, e=e, f=f, m=m, t_onset=t_onset)
  P_obs <- round(sim_data$P[sample_i] + rnorm(10, 0, 2))
  test_data <- list(P = P_obs, time = sample_t, n = length(P_obs))
  
  test_fit <- sampling(model_for_test, test_data, iter=4000, chains=4, cores=4, control=list(adapt_delta=0.95)) 
  verhulst_est <- summary(test_fit)$summary[1,1]
  return(verhulst_est)
}

#subsample_datasets <- round(seq(1,50, length.out=10)) # running verhult test on just 10 simulated datasets rather than 50
subsample_datasets <- c(1:51)
verhulst_tests <- mapply(test.verhulst, t_onset = t_onsets[subsample_datasets], e = es[subsample_datasets], MoreArgs = list(t = t, B0 = B0, f = f, m = m))
verhulst_tests_df <- data.table(t(verhulst_tests))


png("methods_test_50.png", width=700, height=600)
#svg("methods_test_50.svg", width=7, height=6) # if a vector graphic is desired

  plot(NA, xlim=range(t_onsets), ylim=range(t_onsets), xlab="Actual onset (DOY)", ylab="Estimated onset (DOY)")
  abline(0, 1, lwd=3, col="black", lty=2)
  
  points(first_est ~ t_onset, data=tests_df, pch=20, col = adjustcolor(pal[1], 0.5))
  abline(lm(first_est ~ t_onset, data=tests_df), lwd=2, col=pal[1])
  
  points(quant_est ~ t_onset, data=tests_df, pch=20, col = adjustcolor(pal[2], 0.5))
  abline(lm(quant_est ~ t_onset, data=tests_df), lwd=2, col=pal[2])
  
  points(gam_est_05 ~ t_onset, data=tests_df, pch=20, col = adjustcolor(pal[3], 0.5))
  abline(lm(gam_est_05 ~ t_onset, data=tests_df), lwd=2, col=pal[3])
  
  points(verhulst_tests ~ t_onsets[subsample_datasets], data=tests_df, pch=20, col = adjustcolor(pal[4], 0.5))
  abline(lm(verhulst_tests ~ t_onsets[subsample_datasets], data=tests_df), lwd=2, col=pal[4])

  legend("topleft", inset=0.05,
         legend=c("First observation", "Empirical quantile", "Generalized Additive Model", "Bayesian phenological curve model", "One-to-one line"),
         col=c(pal[1:4],"black"), lwd=2, lty=c(1,1,1,1,2), cex=0.9)
  
dev.off()

all_tests <- list(tests_df = tests_df, t_onsets = t_onsets, es = es, subsample_datasets = subsample_datasets, verhulst_tests = verhulst_tests)

saveRDS(all_tests, "all_tests.RDS")
# all_tests <- readRDS("all_tests.RDS") # run this and not the line above to save time
tests_df <- all_tests$tests_df
t_onsets <- all_tests$t_onsets
es <- all_tests$es
subsample_datasets <- all_tests$subsample_datasets
verhulst_tests <- all_tests$verhulst_tests

summary(lm((t_onsets - tests_df$first_est) ~ 1))
summary(lm((t_onsets - tests_df$quant_est) ~ 1))
summary(lm((t_onsets - tests_df$gam_est_05) ~ 1))
summary(lm((t_onsets - verhulst_tests) ~ 1))

library(car)
linearHypothesis(lm(first_est ~ t_onset, data=tests_df), "t_onset = 1")
linearHypothesis(lm(quant_est ~ t_onset, data=tests_df), "t_onset = 1")
linearHypothesis(lm(gam_est_05 ~ t_onset, data=tests_df), "t_onset = 1")
linearHypothesis(lm(verhulst_tests ~ t_onset, data=tests_df), "t_onset = 1")


