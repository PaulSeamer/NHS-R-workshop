

# required packages
library("tidyverse")
library("Hmisc")
library("zoo")
library("Matching")
library("tableone")
library("pROC")




# read rhc dataset ------------------------------------------------------------
# YOU NEED TO POINT R TO THE "rhc.RDS" FILE ON YOUR MACHINE
rhc <- readRDS("rhcDataset.RDS")
rhc
# View(rhc)




# explore the rhc dataset -----------------------------------------------------
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
         summarise(n = n()) %>%
         mutate(freq = n / sum(n))) +
  geom_bar(aes(x = swang1, y = freq, fill = dth30), stat = "identity") +
  coord_flip()

# plot age
ggplot(rhc) +
  geom_density(aes(x = age, group = swang1, fill = swang1), alpha = .4)




# exact matching --------------------------------------------------------------
# why exact matching is problematic
# select a small number of key matching variables
matchVars <- c("age", "sex", "cat1", "meanbp1")

# filter on key matching variables
rhc <- rhc %>%
  dplyr::select_at(vars(swang1, dth30, matchVars)) %>% 
  mutate(age = as.integer(age))

# use tableone
unmatchedTbl1 <- CreateTableOne(
  vars = matchVars, strata = "swang1", data = rhc, test = FALSE)
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
  
# why exact matching is problematci
set.seed(12345)
exactMatch <- Match(
  Y = rhc$dth30
  , Tr = rhc$swang1
  , M = 1
  , X = matchVarsDf
  , ties = FALSE
  , replace = FALSE  # ties are randomly broken when replace = FALSE
  , estimand = "ATT"
  , exact = colnames(matchVarsDf) %in% c("age", "sex", "cat1", "meanbp1"))
  
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




# psm matching ----------------------------------------------------------------
# fit a propensity score model (logistic regression)
psmM1 <- glm(swang1 ~ age + sex + cat1 + meanbp1
               , family = binomial()
               , data = rhc)

# check psm model statistics
summary(psmM1)

# predictions from model - probability of assignment to treatement i.e. receiving rhc
rhc$psmM1 <- predict.glm(psmM1, type = c("response"))

# plot a roc curve to evaluate models predictive ability
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

# 2184 matches (from 2184 cases)
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




# greedy matching on Mahanalobis ----------------------------------------------
# add-in some more covariates
rhc <- readRDS("rhc.RDS")
matchVars <- c("age", "sex", "cat1", "meanbp1"
               , "surv2md1", "hrt1", "resp1", "renalhx", "liverhx")

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

greedyMatch1 <- Match(
  Y = rhc$dth30
  , Tr = rhc$swang1
  , M = 1
  , X = matchVarsDf
  , ties = FALSE
  , replace = FALSE  # ties are randomly broken when replace = FALSE
  , estimand = "ATT"
  , exact = colnames(matchVarsDf) %in% c("cat1"))


# more useful to have variables as character/factor 
greedyMatch1Df <- rhc[unlist(greedyMatch1[c("index.treated", "index.control")]), ]

# 2011 matches (from 2184 cases)
table(greedyMatch1Df$swang1)

# check balance
matchedTbl1 <- CreateTableOne(vars = matchVars, strata = "swang1", data = greedyMatch1Df , test = FALSE)
print(matchedTbl1, smd = TRUE)

# examine outcome
prop.table(table(greedyMatch1Df$swang1, greedyMatch1Df$dth30, useNA = "ifany"), margin = 1)

# test outcome
trtDth30 <- greedyMatch1Df %>% filter(swang1 == 1) %>% pull(dth30)
conDth30 <- greedyMatch1Df %>% filter(swang1 == 0) %>% pull(dth30)
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

