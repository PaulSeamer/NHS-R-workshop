

# required packages
library("season")  # monthglm()
library("zoo")  # as.yearmon()
library("tidyverse")
# also requires "Epi", "lmtest" and "tsModel" - these just need to be installed not loaded




# read Sicily dataset ---------------------------------------------------------
# read Sicily dataset from csv file
# YOU NEED TO POINT R TO THE "sicilyDataset.csv" FILE ON YOUR MACHINE
sicily <- read_csv("C:/.../sicilyDataset.csv")
sicily
# View(sicily)

# variables
# time = elapsed time since the start of the study
# aces = monthly count of acute coronary episodes (the outcome)
# smokban = smoking ban (the intervention) coded 0 before intervention, 1 after
# pop = total population
# stdpop = age standardised population




# explore the Sicily dataste --------------------------------------------------
# create date variable - year-month format using zoo::as.yearmon()
sicily <- sicily %>%
  mutate(dt = as.yearmon(paste(sicily$year, sicily$month, sep = "-")))

# plot admissions prior to intervention
ggplot(data = sicily %>% filter(smokban == 0)) +
  geom_point(aes(x = dt, y = aces), size = 1, shape = 16, color = "red") +
  scale_x_yearmon(name = "time", format = "%b-%y") +
  scale_y_continuous(limits = c(0, 1200))

# plot full series and intervention
ggplot(data = sicily) +
  geom_point(aes(x = time, y = aces), size = 1, shape = 16, color = "red") +
  geom_vline(xintercept = 36.5, color = "grey") +
  scale_x_continuous(name = "time") +
  scale_y_continuous(limits = c(0, 1200))

# plot population
ggplot(data = sicily) +
  geom_point(aes(x = dt, y = pop), size = 1, shape = 16, color = "red") +
  scale_x_yearmon(name = "time", format = "%b-%y") +
  scale_y_continuous(limits = c(0, max(sicily$pop)))




# investiagte seasonality -----------------------------------------------------
# plot admissions series by month prior to intervention -
# you need to set the admission series as a time series variables using ts()
monthplot(ts(sicily %>% filter(smokban == 0) %>% pull(aces), frequency = 12))

# use stl() to seasonally decompose the admissions series -
# you need to set the admission series as a time series variables using ts()
plot(stl(ts(sicily %>% filter(smokban == 0) %>% pull(aces), frequency = 12), s.window = "periodic"))

# use a loess smoother to help see underlying patterns in admissions series
ggplot(sicily %>% filter(smokban == 0)) +
  geom_line(aes(x = dt, y = aces), size = 0.1) +
  geom_smooth(aes(x = dt, y = aces) 
              , method = "loess", span = 0.4, size = 1, color = "#5881c1") +
  scale_x_yearmon(name = "time", format = "%b-%y") +
  scale_y_continuous(limits = c(0, 1200))

# use a poisson glm model to estimate the seasonal effect of month of year
sicilyMoy <- monthglm(
  aces ~ year
  , data = sicily %>% filter(smokban == 0)
  , family = poisson()
  , refmonth = 1
  , monthvar = "month"
  , offsetmonth = TRUE)

# extract the coefficients from the above model
sicilyMoyCoef <- tibble(
  coef = names(sicilyMoy$glm$coefficients)
  , names2 = c("(Intercept)", "year", tail(month.abb, -1L))
  , value = sicilyMoy$glm$coefficients)

# give the coefficients nice names ready for plotting
sicilyMoyPlotData <- sicilyMoyCoef %>% 
  filter(!coef %in% c("year", "(Intercept)")) %>% 
  select(-coef) %>% 
  mutate(names2 = factor(names2, levels = month.abb))

# plot month of year coefficients - these show the seasonal pattern in admissions by month
# e.g. admissions in November & December are c.5% higher than in January, and 9% lower in August than in January
ggplot(sicilyMoyPlotData) +
  geom_point(aes(x = names2, y = value)) +
  geom_hline(aes(yintercept = .0), size = .5) 




# level change impact model ----------------------------------------------------
# use segemented regression to test for a level change in admissions associated with the smoking ban

# fit a simple linear trend model to the pre intervention period
modelPre <- glm(aces ~ time, family = gaussian, data = sicily %>% filter(smokban == 0))
summary(modelPre)

sicilyPre <- sicily %>% 
  filter(smokban == 0) %>% 
  mutate(preTrend = predict(modelPre, type = "response"))

# plot the linear trend model
ggplot(data = sicilyPre) +
  geom_point(aes(x = time, y = aces)) +
  geom_line(aes(x = time, y = preTrend), color = "red")
  # # if you want prettier x-axis labels add these 2 lines to the plot
  # scale_x_continuous(name = "time", breaks = seq(.5, 60.5, by = 12), labels = 2002:2007) +
  # scale_y_continuous(limits = c(0, 1200))

# use the above model to predict admissions as if the ban was never implemented i.e. to generate a counterfactual
sicilyCf <- sicily %>%
  # note you need to supply the entire dataset in newdata = 
  mutate(preCf = predict(modelPre, type = "response", newdata = sicily))

# plot the counterfactual
ggplot(data = sicilyCf) +
  geom_point(aes(x = time, y = aces)) +
  geom_line(aes(x = time, y = preCf), color = "red", linetype = "dashed") +
  geom_vline(xintercept = 36.5, color = "grey")
  # # add these extra lines if you want a prettier x-axis
  # scale_x_continuous(name = "time", breaks = seq(.5, 60.5, by = 6), labels = c("", rbind(2002:2006, ""))) +
  # scale_y_continuous(limits = c(0, 1200)) +
  # theme(
  #   panel.grid.minor.x = element_blank()
  #   , panel.grid.major.x = element_line(colour = c("white", rep(c(NA, "white"), t = 6)))
  #   , axis.ticks.x=element_line(colour=c("black", rep(c(NA, "black"), t = 6)))
  # )

# impact model for a change in level
modelLevelChg <- glm(aces ~ time + smokban, family = gaussian, data = sicily)
# note - smokan coefficient suggests a drop of 89 admissions associated with the smoking ban
summary(modelLevelChg)
# print confidence intervals
round(confint(modelLevelChg), 3)

# calculate predicted admissions from the level change model
sicilyLevelChg <- sicily %>% 
  mutate(levelChg = predict(modelLevelChg, type = "response"))

# plot the level change model
ggplot(data = sicilyLevelChg) +
  geom_point(aes(x = time, y = aces)) +
  geom_line(aes(x = time, y = levelChg), color = "red") +
  geom_vline(xintercept = 36.5, color = "grey") +
  scale_x_continuous(
    name = "time", breaks = seq(.5, 60.5, by = 6)
    , labels = c("", rbind(2002:2006, ""))) +
  scale_y_continuous(limits = c(0, 1200)) +
  theme(
    panel.grid.minor.x = element_blank()
    , panel.grid.major.x = element_line(
      colour = c("white", rep(c(NA, "white"), t = 6)))
    , axis.ticks.x=element_line(colour=c("black", rep(c(NA, "black"), t = 6)))
  )

# create a counterfactual using the level change model
sicilyLevelChg <- sicilyLevelChg %>% 
  mutate(levelChgCf = predict(
    modelLevelChg, type = "response", newdata = sicily %>% mutate(smokban = 0)))

# plot counterfactual and level change impact model
ggplot(data = sicilyLevelChg) +
  geom_point(aes(x = time, y = aces)) +
  geom_line(aes(x = time, y = levelChg), color = "red") +
  geom_line(aes(x = time, y = levelChgCf), color = "blue") +
  geom_vline(xintercept = 36.5, color = "grey") +
  scale_x_continuous(name = "time", breaks = seq(.5, 60.5, by = 6), 
                     labels = c("", rbind(2002:2006, ""))) +
  scale_y_continuous(limits = c(0, 1200)) +
  theme(
    panel.grid.minor.x = element_blank()
    , panel.grid.major.x = element_line(
      colour = c("white", rep(c(NA, "white"), t = 6)))
    , axis.ticks.x=element_line(colour=c("black", rep(c(NA, "black"), t = 6)))
  )



# model checks ----------------------------------------------------------------
# check level change model residuals by plotting against time
res <- residuals(modelLevelChg, type = "deviance")

ggplot() +
  geom_point(aes(x = sicilyLevelChg$time, y = res)) +
  scale_x_continuous() +
  scale_y_continuous(limits = c(min(res * 1.5), max(res * 1.5)))

# check for autocorrelation
acf(res)
pacf(res)
lmtest::dwtest(modelLevelChg)  # Durbin-Watson test suggests presence of autocorrelation

# more diagnostic plots
plot(modelLevelChg)



# level and slope change impact model -----------------------------------------
# sicily <- sicily %>%
#   mutate(time2 = c(rep(0, 36), seq(37, 59)))

modelSlopeLevel <- glm(aces ~ time + smokban*time, data = sicily)
summary(modelSlopeLevel)
round(confint(modelSlopeLevel), 3)
# neither level change or slope change coefficient are statistically significant

# generate predictions from this model for plotting
sicilySlopeLevel <- sicily %>% 
  mutate(slopeLevel = predict(modelSlopeLevel, type = "response"))

ggplot(data = sicilySlopeLevel) +
  geom_point(aes(x = time, y = aces)) +
  geom_line(aes(x = time, y = slopeLevel), color = "red") +
  geom_vline(xintercept = 36.5, color = "grey") +
  scale_x_continuous(name = "time", breaks = seq(.5, 60.5, by = 6), labels = c("", rbind(2002:2006, ""))) +
  scale_y_continuous(limits = c(0, 1200)) +
  theme(
    panel.grid.minor.x = element_blank()
    , panel.grid.major.x = element_line(colour = c("white", rep(c(NA, "white"), t = 6)))
    , axis.ticks.x=element_line(colour=c("black", rep(c(NA, "black"), t = 6)))
  )




# more models -----------------------------------------------------------------
# month as a categorical variable to account for seasonality
m01 <- glm(aces ~ time + smokban + as.factor(month), family = gaussian, data = sicily)
summary(m01)

# poisson model for count variables
m02 <- glm(aces ~ time + smokban, family = poisson, sicily)
# coefficents are changes in log of expected counts
summary(m02)
# exponentiate to return incidence rate ratios
# e.g. effect of smoking ban is c.10% drop in admissions 
round(Epi::ci.lin(m02, Exp = T), 3)

# poisson model with harmoinic for seasonality
m03 <- glm(aces ~ time + smokban + tsModel::harmonic(month, 2, 12), family = poisson, sicily)
summary(m03)
round(Epi::ci.lin(m03, Exp = T), 3)

# suggests presence of over dispersion
summary(m03)$dispersion

# try quasipoisson to account for overdispersion
m04 <- glm(aces ~ time + smokban + tsModel::harmonic(month, 2, 12), family = quasipoisson, sicily)
summary(m04)
round(Epi::ci.lin(m04, Exp = T), 3)

# better
summary(m04)$dispersion

# diagnostics
res <- residuals(m04, type = "deviance")

ggplot() +
  geom_point(aes(x = sicily$time, y = res)) +
  scale_x_continuous() +
  scale_y_continuous(limits = c(min(res * 1.5), max(res * 1.5)))

acf(res)
pacf(res)
lmtest::dwtest(m04)  # autocorrelation no longer an issue

# plot m04
ggplot(data = sicily %>% mutate(m04 = predict(m04, type = "response", data = sicily))) +
  geom_point(aes(x = time, y = aces)) +
  geom_line(aes(x = time, y = m04), color = "red") +
  geom_vline(xintercept = 36.5, color = "grey") +
  scale_x_continuous(name = "time", breaks = seq(.5, 60.5, by = 6), labels = c("", rbind(2002:2006, ""))) +
  scale_y_continuous(limits = c(0, 1200)) +
  theme(
    panel.grid.minor.x = element_blank()
    , panel.grid.major.x = element_line(colour = c("white", rep(c(NA, "white"), t = 6)))
    , axis.ticks.x=element_line(colour=c("black", rep(c(NA, "black"), t = 6)))
  )

