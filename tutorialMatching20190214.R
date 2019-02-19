
# required packages
library("tidyverse")
library("Hmisc")
library("zoo")
library("Matching")
library("tableone")
library("pROC")




# read rhc dataset ------------------------------------------------------------
rhc <- readRDS("rhc.RDS")




# exercise 1 ------------------------------------------------------------------
# explore the dataset
is.tibble(rhc)

# identify exposure and outcome variables
# exposure = swang1; outcome = dth30
table(rhc$swang1, useNA = "ifany")
table(rhc$dth30, useNA = "ifany")
table(rhc$swang1, rhc$dth30, useNA = "ifany")

addmargins(table(rhc$swang1, rhc$dth30, useNA = "ifany"))
addmargins(table(rhc$swang1, rhc$dth30, useNA = "ifany"), margin = 2)

prop.table(table(rhc$swang1, useNA = "ifany"))
prop.table(table(rhc$swang1, rhc$dth30, useNA = "ifany"), margin = 1)

# plot counts
ggplot(rhc) +
  geom_bar(aes(x = swang1, fill = dth30))

# plot proportions
ggplot(rhc %>%
         group_by(swang1, dth30) %>% 
         summarise (n = n()) %>%
         mutate(freq = n / sum(n))) +
  geom_bar(aes(x = swang1, y = freq, fill = dth30), stat = "identity") +
  coord_flip()




# exercise 2 ------------------------------------------------------------------
# why exact matching is problematic
# select a small number of key matching variables
matchVars <- c("age", "sex", "cat1", "meanbp1")

# filter on key matching variables
rhc <- rhc %>%
  dplyr::select_at(vars(swang1, dth30, matchVars)) %>% 
  mutate(age = as.integer(age))

# use tableone
unmatchedTbl1 <- CreateTableOne(vars = matchVars, strata = "swang1", data = rhc, test = FALSE)
print(unmatchedTbl1)
print(unmatchedTbl1, smd = TRUE)

# Match() requires numeric variables
rhc <- rhc %>% 
  mutate(
    swang1 = case_when(swang1 == "No RHC" ~ 0, TRUE ~ 1)
    , dth30 = case_when(dth30 == "No" ~ 0, TRUE ~ 1))

# Match() requires numeric variables
matchVarsDf <- rhc %>% 
  dplyr::select(-dth30, -swang1) %>%
  mutate_if(is.character, funs(as.numeric(as.factor(.))))
  
# why exact matching is not feasible
# 
exactMatch <- Match(
  Y = rhc$dth30
  , Tr = rhc$swang1
  , M = 1
  , X = matchVarsDf
  , ties = FALSE
  , replace = FALSE  # ties are randomly broken when replace = FALSE
  , estimand = "ATT"
  , exact = colnames(matchVarsDf) %in% c("age", "sex", "cat1", "meanbp1") )
  
# recover matched dataset
exactMatch$mdata[["Y"]]
exactMatch$mdata[["Tr"]]
exactMatch$mdata[["X"]]

# exactMatchDf <- bind_cols(
#   tibble(dth30 = exactMatch$mdata[["Y"]], swang1 = exactMatch$mdata[["Tr"]])
#   , as_tibble(exactMatch$mdata[["X"]]))

# more useful to have variables as character/factor 
exactMatchDf <- rhc[unlist(exactMatch[c("index.treated", "index.control")]), ]

# only 208 matches (from 2184 cases)
table(exactMatchDf$swang1)

# check balance
matchedTbl1 <- CreateTableOne(vars = matchVars, strata = "swang1", data = exactMatchDf , test = FALSE)
print(matchedTbl1, smd = TRUE)

# examine outcome
prop.table(table(exactMatchDf$swang1, exactMatchDf$dth30, useNA = "ifany"), margin = 1)

# test outcome
trtDth30 <- exactMatchDf %>% filter(swang1 == 1) %>% pull(dth30)
conDth30 <- exactMatchDf %>% filter(swang1 == 0) %>% pull(dth30)
t.test(trtDth30, conDth30, alternative = c("two.sided"), paired = TRUE)

# try running it again - what result do you get?
# https://stackoverflow.com/questions/13605271/reasons-for-using-the-set-seed-function




# exercise 3 ------------------------------------------------------------------
# build a propensity score model

# fit a propensity score model (logistic regression)
psmM1 <- glm(swang1 ~ age + sex + cat1 + meanbp1
               , family = binomial()
               , data = rhc)

# investigate coefficients and mdoel statistics
summary(psmM1)

# predictions from model - probability of death within 30 days
rhc$psmM1 <- predict.glm(psmM1, type = c("response"))

# calculate a roc curve to evaluate models predictive ability
# psmM1Roc <- roc(response = rhc$swang1, predictor = rhc$psmM1)
psmM1Roc <- roc(swang1 ~ psmM1, data = rhc)

# plot the roc curve
plot(psmM1Roc)

# area under the curve
auc(psmM1Roc)
ci.auc(psmM1Roc)

# use propensity score for matching
pScoreM1 <- predict.glm(psmM1, type = c("link"))

# run Match()
pScoreM1Match <- Match(
  Y = rhc$dth30
  , Tr = rhc$swang1
  , M = 1
  , X = pScoreM1
  , ties = FALSE
  , replace = FALSE  # ties are randomly broken when replace = FALSE
  , estimand = "ATT"
  , exact = NULL)

# more useful to have variables as character/factor 
pScoreM1MatchDf <- rhc[unlist(pScoreM1Match[c("index.treated", "index.control")]), ]

# only 2184 matches (from 2184 cases)
table(pScoreM1MatchDf$swang1)

# check balance
matchedTbl1 <- CreateTableOne(vars = matchVars, strata = "swang1", data = pScoreM1MatchDf , test = FALSE)
print(matchedTbl1, smd = TRUE)

# examine outcome
prop.table(table(pScoreM1MatchDf$swang1, pScoreM1MatchDf$dth30, useNA = "ifany"), margin = 1)

# test outcome
trtDth30 <- pScoreM1MatchDf %>% filter(swang1 == 1) %>% pull(dth30)
conDth30 <- pScoreM1MatchDf %>% filter(swang1 == 0) %>% pull(dth30)
t.test(trtDth30, conDth30, alternative = c("two.sided"), paired = TRUE)

# produce plot showing balance pre-post matching
# construct df with variable name and smd 
lovePlotData <- data.frame(variable = dimnames(ExtractSmd(unmatchedTbl1))[[1]]
                           , unmatched = unname(ExtractSmd(unmatchedTbl1))
                           , matched = unname(ExtractSmd(matchedTbl1)))

# wrangle long-form data for ggplot2
lovePlotData <- lovePlotData %>% 
  gather(key = "method", value = "smd", -variable)

# plot
ggplot(lovePlotData, aes(x = variable, y = smd, group = method, color = method)) +
  geom_point(shape = 16, size = 3) +
  geom_hline(yintercept = 0.1, color = "#2c2825", size = 0.4, linetype = "33") +
  coord_flip() +
  scale_color_manual(values = c("#ec6555", "#5881c1")) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.05, 1.05), name = "mean difference", breaks = c(seq(0, 1, 0.1))) +
  scale_x_discrete(expand = c(0, 1), name = element_blank())

# test outcome using regression model (clean up any residual unbalance)
# gaussian model - no covariates - should match t.test() result
glmM1 <- glm(formula = dth30 ~ swang1
             , family = gaussian()
             , data = pScoreM1MatchDf)
summary(glmM1)

# binomial model (appropriate for binary outcome) with covariates
glmM2 <- glm(formula = as.formula(paste0("dth30 ~ swang1 +", paste(matchVars, collapse = "+"))) 
             , family = binomial()
             , data = pScoreM1MatchDf)
summary(glmM2)

exp(cbind(OddsRatios = coef(glmM2), confint(glmM2)))





MatchBalance(
  as.formula(checkBalFm)
  , data = matchingPop
  , match.out = matchOut
  , ks = TRUE
  , nboots = 1000
  , weights = NULL
  , digits = 5
  , paired = TRUE
  , print.level = 1)














# F; T = 2078
# F; F = 1678  # ties are randomly broken when replace = FALSE
# T; F = 1678  # ties are randomly broken when replace = FALSE
# T; T = 16620 (2078)

matched <- matched %>%
  mutate(obsWeight = rep(exactMatch$weights, 2))

sum(exactMatch$weights)
matchVarsDf %>% group_by_all(.) %>% summarise(n = n()) %>% arrange(-n)
exactMatch$index.dropped
exactMatch$ndrops.matches
exactMatch$orig.nobs
exactMatch$orig.treated.nobs
exactMatch$nobs
exactMatch$wnobs

y_trt <- matched$dth30[matched$swang1 == 1] * matched$obsWeight[matched$swang1 == 1]
y_con <- matched$dth30[matched$swang1 == 0] * matched$obsWeight[matched$swang1 == 0]




# pairwise difference
diffy <- y_trt - y_con
# 651 / 1678 = 39%
# 486 / 1678 = 29%
# paired t-test
t.test(diffy)



table(y_trt, y_con)
mcnemar.test(table(y_trt, y_con))


library("survey")
matchedSvy <- svydesign(ids = ~ 1, weights = ~ obsWeight, data = matched)
# table counts of outcome - unweighted
table(matched[matched$swang1 == 1, ]$dth30)
table(matched[matched$swang1 == 0, ]$dth30)
# crude (unweighted) odds
6095 / (10525 + 6095) # 36.7
4676 / (11944 + 4676) # 28.1
# correct (weighted) odds
svytable(~ dth30 + swang1, design = matchedSvy)
606.675 / (1471.325 + 606.675) # 29.1
791 / (1287 + 791) # 38.1



glmM1 <- glm(formula = dth30 ~ swang1
             , family = gaussian()
             , data = matched)
summary(glmM1)

glmM1 <- glm(formula = dth30 ~ swang1
             , family = gaussian()
             , data = rhc)
summary(glmM1)

glmM2 <- glm(formula = dth30 ~ swang1
             , family = binomial()
             , data = matched)
summary(glmM2)
exp(cbind(OddsRatios = coef(glmM2), confint(glmM2)))


glmM2 <- glm(formula = dth30 ~ swang1
             , family = binomial()
             , data = matched
             , weights = )
summary(glmM2)
exp(cbind(OddsRatios = coef(glmM2), confint(glmM2)))


glmM2 <- glm(formula = dth30 ~ swang1
             , family = binomial()
             , data = rhc)
summary(glmM2)
exp(cbind(OddsRatios = coef(glmM2), confint(glmM2)))


# try Mahanalobis distance




# psm ---------------------------------------------------------------------



psmMatch <- Match(
  Tr = rhc$swang1
  , M=3
  , X= pscore2
  , replace = FALSE
  , caliper = .2)

matched <- rhc[unlist(psmMatch[c("index.treated", "index.control")]), ]
matchedTbl1 <- CreateTableOne(vars = matchVars, strata = "swang1", data = matched, test = FALSE)
print(matchedTbl1, smd = TRUE)
length(psmMatch$weights)

y_trt <- matched$dth30[matched$swang1 == 1]
y_con <- matched$dth30[matched$swang1 == 0]

# pairwise difference
diffy <- y_trt - y_con
# 760 / 2000 = 38%
# 585 / 2000 = 29%
# paired t-test
t.test(diffy)
table(y_trt, y_con)


a <- MatchBalance(
  swang1 ~ age + sex + cat1
  , data = rhc
  , match.out = matched
  , ks = TRUE
  , nboots = 1000
  , weights = NULL
  , digits = 5
  , paired = TRUE
  , print.level = 1)


















###################
#RHC Example

#install packages
install.packages("tableone")
install.packages("Matching")

#load packages
library(tableone)
library(Matching)

#read in data
load(url("http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/rhc.sav"))
#view data
View(rhc)

#treatment variables is swang1
#x variables that we will use
#cat1: primary disease category
#age
#sex
#meanbp1: mean blood pressure

#create a data set with just these variables, for simplicity
ARF<-as.numeric(rhc$cat1=='ARF')
CHF<-as.numeric(rhc$cat1=='CHF')
Cirr<-as.numeric(rhc$cat1=='Cirrhosis')
colcan<-as.numeric(rhc$cat1=='Colon Cancer')
Coma<-as.numeric(rhc$cat1=='Coma')
COPD<-as.numeric(rhc$cat1=='COPD')
lungcan<-as.numeric(rhc$cat1=='Lung Cancer')
MOSF<-as.numeric(rhc$cat1=='MOSF w/Malignancy')
sepsis<-as.numeric(rhc$cat1=='MOSF w/Sepsis')
female<-as.numeric(rhc$sex=='Female')
died<-as.numeric(rhc$death=='Yes')
age<-rhc$age
treatment<-as.numeric(rhc$swang1=='RHC')
meanbp1<-rhc$meanbp1

#new dataset
mydata<-cbind(ARF,CHF,Cirr,colcan,Coma,lungcan,MOSF,sepsis,
              age,female,meanbp1,treatment,died)
mydata<-data.frame(mydata)

#covariates we will use (shorter list than you would use in practice)
xvars<-c("ARF","CHF","Cirr","colcan","Coma","lungcan","MOSF","sepsis",
         "age","female","meanbp1")

#look at a table 1
table1<- CreateTableOne(vars=xvars,strata="treatment", data=mydata, test=FALSE)
## include standardized mean difference (SMD)
print(table1,smd=TRUE)


############################################
#do greedy matching on Mahalanobis distance
############################################

greedymatch<-Match(Tr=treatment,M=1,X=mydata[xvars],replace=FALSE)
matched<-mydata[unlist(greedymatch[c("index.treated","index.control")]), ]

#get table 1 for matched data with standardized differences
matchedtab1<-CreateTableOne(vars=xvars, strata ="treatment", 
                            data=matched, test = FALSE)
print(matchedtab1, smd = TRUE)

#outcome analysis
y_trt<-matched$died[matched$treatment==1]
y_con<-matched$died[matched$treatment==0]

#pairwise difference
diffy<-y_trt-y_con

#paired t-test
t.test(diffy)

#McNemar test
table(y_trt,y_con)

mcnemar.test(matrix(c(973,513,395,303),2,2))



##########################
#propensity score matching
#########################

#fit a propensity score model. logistic regression

psmodel<-glm(treatment~ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+
               sepsis+age+female+meanbp1+aps,
             family=binomial(),data=mydata)

#show coefficients etc
summary(psmodel)
#create propensity score
pscore<-psmodel$fitted.values


#do greedy matching on logit(PS) using Match with a caliper

logit <- function(p) {log(p)-log(1-p)}
psmatch<-Match(Tr=mydata$treatment,M=1,X=logit(pscore),replace=FALSE,caliper=.2)
matched<-mydata[unlist(psmatch[c("index.treated","index.control")]), ]
xvars<-c("ARF","CHF","Cirr","colcan","Coma","lungcan","MOSF","sepsis",
         "age","female","meanbp1")

#get standardized differences
matchedtab1<-CreateTableOne(vars=xvars, strata ="treatment", 
                            data=matched, test = FALSE)
print(matchedtab1, smd = TRUE)

#outcome analysis
y_trt<-matched$died[matched$treatment==1]
y_con<-matched$died[matched$treatment==0]

#pairwise difference
diffy<-y_trt-y_con

#paired t-test
t.test(diffy)


