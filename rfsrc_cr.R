# Clean workspace
rm(list=ls())

library(dplyr)
library(ggplot2)
library(data.table)
library(tidyr)
library(tibble)
library(gridExtra)
library(iBST)
library(tidyverse)
library(reshape2)
library(randomForestSRC) # the main package
library(pec) # crps, and ibs are used for plots. selectFGR and selectCox are used for the CoxPH model.
library(readxl)
library(SurvMetrics)
library(dplyr)
library(scales)
library(tidyr)

# data preprocessing
#Car_Leases <- read_excel("/Users/seansexton/Desktop/Research/Survival_Analysis/Data/sample comp-risks (cols order changed).xlsx")
Car_Leases <- read_excel("~/Desktop/PhD/CR_IML/sample comp-risks (cols order changed).xlsx")

Car_Leases <- Car_Leases[Car_Leases$NUMBER_OF_INSTALMENTS == 36,] # only keep leases with 36 month terms (9,081 observations)

# Same list as paper
namelist <- c("CAR_TYPE","NEW_USED_CAR","PHONE","CUSTOMER_BEFORE", # Factors 
              "NET_ASSETS","DOWNPAYMENT_TO_CAR_PRICE","ANNUAL_TURNOVER_LAST_YEAR", # Continuous 
              "ANNUAL_INCOME_LAST_YEAR","CAR_PRICE", "RESIDUAL_VALUE", "RESIDUAL_VALUE_TO_CAR_PRICE","ANNUAL_COSTS",
              "COMPANY_AGE", "NUMBER_OF_ALL_ACTIVE_PREVIOUS_CO","NUMBER_OF_ALL_PREVIOUS_CONTRACTS", # Discrete
              "NUMBER_OF_ALL_BAD_PREVIOUS_CONTR", "NUMBER_OF_ALL_COMPLETED_PREVIOUS",
              "car_age","NUMBER_OF_INSTALMENTS","NUMBER_OF_EMPLOYEES","length_coop",
              "firm_LTD_PARTN","firm_joint_stock","firm_other", # legal form of the company, ref: freelancer
              "REGION1","REGION2","REGION3","REGION4", # Geographical Region; ref: Region E.
              "BRANCH_services","BRANCH_construction","BRANCH_Sale","BRANCH_Production", # Branch, ref: not specified
              "BRANCH_health","BRANCH_other","duration","status")

events <- Car_Leases[Car_Leases$event == "D" | Car_Leases$event == "E",] # dataset of events
non_events <- Car_Leases[Car_Leases$event != "D" & Car_Leases$event != "E",] # dataset of non-events

set.seed(1337)
sample <- non_events[sample(nrow(non_events), 2665), ] # randomly sample 2665 leases from non_events
CL <- rbind(events, sample) # 'more balanced' dataset. C/0=current, D/1=default,E/2=early prepayment, P/3=paid on time
CL <- CL[,namelist]

CL$NUMBER_OF_INSTALMENTS <- NULL # singular
CL$NUMBER_OF_ALL_ACTIVE_PREVIOUS_CO <- NULL # correlation of almost 1 with NUMBER_OF_ALL_ACTIVE_PREVIOUS_CO
CL$NEW_USED_CAR <- NULL # we have car age 
# CL$length_coop <- NULL # high correlation with CUSTOMER_BEFORE and showed poor p-values earlier.
# CL$RESIDUAL_VALUE <- NULL # We have residual value to car price
# CL$CAR_PRICE <- NULL # we have residual value to car price
CL$REGION4 <- NULL # dummy variable trap. Keep REGION1, 2, and 3.

CL <- CL[CL$duration <= 36,]
CL <- CL[-which(CL$status == 0 & CL$duration == 36),]
CL <- CL[-which(CL$status == 2 & CL$duration == 36),]


# fit rfsrc
temp <- CL

# When the goal is t-year prediction, then the CIF is the appropriate measure 
# and in this case the appropriate splitting rule is logrankCR
v.obj5 <- rfsrc(Surv(duration, status) ~., data = temp, 
                ntree = 1000,
                nodesize = 15, # Note that generally, competing risks requires a larger nodesize than survival settings
                mtry = 6,
                nsplit = 10,
                samptype = "swr", # sample with replacement
                splitrule = "logrankCR", # note that Gray's rule is used here. 
                importance = TRUE, # logrank
                statistics = TRUE)

print(v.obj5)

# variable importance 
vimp.results <- round(as.data.frame(v.obj5$importance),3)
vimp.results <- tibble::rownames_to_column(vimp.results, "Variables")
colnames(vimp.results)[2:3] <- c("default_importance","prepayment_importance")
vimp.results <- vimp.results[order(-vimp.results$default_importance),]
vimp.results

# Minimal Depth & find threshold.

maxsubtreeinfo2 <- max.subtree(v.obj5)

# strong variables have minimal depth less than or equal to the following threshold
# print(maxsubtreeinfo2$threshold) # different # than the paper: 8.943 vs 9.375

# this corresponds to the set of variables
topvars2 <- maxsubtreeinfo2$topvars

md2 <- maxsubtreeinfo2$order[, 1] # the minimal depth is the first order depth
md2 <- md2[topvars2]
md2 <- as.data.frame(md2)
md2 <- tibble::rownames_to_column(md2, "Variables")

top_variables <- merge(vimp.results, md2, by=c("Variables"))
top_variables <- top_variables[order(-top_variables$default_importance),]

print(top_variables)


## Create dataset with only the top variables
final_variables3 <- as.character(top_variables$Variables)
namelist <- c("status", "duration", final_variables3) # include 
temp <- CL[,namelist]

v.obj6 <- rfsrc(Surv(duration, status) ~., data = temp, 
                ntree = 1000, # "B" bootstraps with replacement. Alternative is to sample with replacement 0.63 × N = 0.63(1596) = 1008.672
                nodesize = 15, # min # of observations in terminal node
                mtry = 6, # of variables to possibly split at each node. sqrt(p)=sqrt(32) rounded up.
                nsplit = 10, # specifies number of random splits for splitting a variable
                samptype = "swr", # sample with replacement
                splitrule = "logrankCR",
                importance = TRUE, # logrank
                statistics = TRUE) # used for later statistic info

print(v.obj6)

# cumulative median default rate vs median prepayment rate
df <- as.data.frame(v.obj6$cif)
default.cif <- as.data.frame(df[,1:32]) # obtain the cif for default
prepay.cif <- as.data.frame(df[,33:64]) # obtain cif for prepayments
default.cif.median <- apply(default.cif, 2, median)
prepay.cif.median <- apply(prepay.cif, 2, median) 

# Create a dataframe with the vectors
df <- data.frame(Default = default.cif.median,
                 Prepay = prepay.cif.median)

# Add an 'Index' column representing the sequence of indices
df$Index <- seq_along(default.cif.median)

# Plot using ggplot2
prepay.vs.default <- ggplot(df, aes(x = Index)) +
  geom_line(aes(y = Default, color = "Default", linetype = "Default")) +
  geom_line(aes(y = Prepay, color = "Prepay", linetype = "Prepay")) +
  labs(x = "Index", y = "Value", title = "Comparison of default and prepay") +
  scale_color_manual(values = c(Default = "blue", Prepay = "red"),
                     labels = c("Default", "Prepay")) +
  scale_linetype_manual(values = c(Default = "solid", Prepay = "dashed"),
                        labels = c("Default", "Prepay")) +
  theme_minimal() +
  theme(legend.position = c(0.2, 0.90)) +  # Position legend in the top-left corner
  labs(color = "", linetype = "") +  # Remove legend titles
  scale_y_continuous(limits = c(0, 0.06), labels = percent_format())  # Set y-axis limits and format labels as percentages

prepay.vs.default

# pdp 
plot.variable(v.obj6, 
              xvar.names = c("NUMBER_OF_ALL_BAD_PREVIOUS_CONTR", "NET_ASSETS", "DOWNPAYMENT_TO_CAR_PRICE"), # Names of the x-variables to be used
              target = 1, # indicates the event of interest. Default is 1
              time = NULL, # do we need to put something here?
              surv.type = c("years.lost"), # cif and chf do not work.
              partial = TRUE, 
              oob = TRUE, # OOB (TRUE) or in-bag (FALSE) predicted values.
              show.plots = TRUE,
              plots.per.page = NULL, # Integer value controlling page layout.
              sorted = FALSE, # sort variables by importance
              smooth.lines = TRUE)

CL <- CL[-which(CL$status == 3),]