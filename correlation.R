library(corrplot)
library(RColorBrewer)
library(ggplot2)
feature_data <- temp[, (names(temp) %in% 
                          c("ANNUAL_TURNOVER_LAST_YEAR", 
                            "ANNUAL_INCOME_LAST_YEAR", 
                            "BRANCH_other",
                            "BRANCH_services",
                            "PHONE",
                            "car_age",
                            "length_coop",
                            "CAR_PRICE",
                            "COMPANY_AGE",
                            "CUSTOMER_BEFORE",
                            "DOWNPAYMENT_TO_CAR_PRICE",
                            "NET_ASSETS",
                            "NUMBER_OF_ALL_BAD_PREVIOUS_CONTR",
                            "NUMBER_OF_EMPLOYEES"))]
M <- cor(feature_data)
head(round(M,2))
corrplot(M, order = 'AOE', col = COL2('RdBu', 10), tl.cex = 0.4)

