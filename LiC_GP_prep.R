###################################################################
######### Life in Conservation: Goal progress preparation #########
###################################################################

# The following describes the primary steps used to prepare the data used in the analysis presented in the article 'Balancing making a difference with making a living in the conservation sector'. 

######### Steps #########
# 1) Set up the environment 
# 2) Multiple imputation 
# 3) Post imputation manipulation 
# 4) Estimate plausible values of dispositional optimism 

###### 1) Set up the environment ######
### Load packages ### 
library(mice)
library(fastDummies)
library(plyr)
library(psych)
library(kableExtra)
library(ltm)
library(lavaan)
library(semTools)
library(semPlot)
library(ggplot2)
library(ggpubr)
library(lavaan)

### Set seed ###
set.seed(123)

### Load data ###
# DF.GP.1.Rdata should be in the same folder as the code  
load("DF.GP.1.Rdata") 

### Load functions - save this ### 
# DF.GP.1.Rdata should be in the same folder as the code  
source("LiC_GP_functions.R")

###


###### 2) Multiple imputation ######
### Subset to the imputed variables ###
# Imputed variables 
imptuted_variables <- c("LOTR_1","LOTR_2","LOTR_3","LOTR_4", "LOTR_5","LOTR_6","age_year_scaled", "years_cons_scaled", 
"gender" , "CV", "position_simple", "education_simple", "SO_Region" , "WH_scaled", "Environment")

# Subset the DF to include only the variables to be imputed
DF.GP.1_imp <- DF.GP.1[imptuted_variables]


### Function that ensures the data is in the correct format, drops "other" as a category from gender (because there are too few observations for stable estimation), determines the imputation method from the data type, and performs the imputation using MICE.
MI_function_1 <-function(DF1) {
  
  ### Correct data type ###
  # As ordinal
  DF1[,c("LOTR_1","LOTR_3","LOTR_6","LOTR_2","LOTR_4" , "LOTR_5")] <- lapply(DF1[,c("LOTR_1","LOTR_3","LOTR_6","LOTR_2","LOTR_4" , "LOTR_5")], as.ordered)
  
  # As numeric
  DF1[,c("age_year_scaled","years_cons_scaled", "WH_scaled")] <- lapply(DF1[,c("age_year_scaled","years_cons_scaled", "WH_scaled")], as.numeric)
  
  # As factor 
  DF1[,c("education_simple", "gender",  "position_simple", "SO_Region" , "Environment")] <- lapply(DF1[,c("education_simple", "gender","position_simple", "SO_Region",  "Environment")], as.factor)
  
  ### Drop "other" as a category in gender (there were too few observations to model)
  DF1<-DF1[!(DF1$gender=="Other"),]
  DF1<- droplevels(DF1)
  
  ### Number of imputed DF ###
  # 10 imputed DF
  N.Imp = 10
  
  ### determine the imputation method from the data type ##
  # The data type for each variable 
  str_out <- data.frame(capture.output(str(DF1)))
  
  # Delete the first and last row
  str_out <- data.frame(str_output = matrix(str_out[2:(nrow(str_out)-1),]))
  
  # Create a column that contain the appropriate model for each variable - this only works if the variable is of the correct data type in the first place 
  str_out$type <- ifelse(grepl("Ord.factor", str_out$str_output, fixed = TRUE)==T, "polr", 
                         ifelse(grepl("num", str_out$str_output, fixed = TRUE)==T, "pmm", 
                                ifelse(grepl("Factor", str_out$str_output, fixed = TRUE)==T, "polyreg",
                                       ifelse(grepl("int", str_out$str_output, fixed = TRUE)==T, "logreg", "ERROR"))))
  
  # Conduct the MI - with the number of datasets specified by N.Imp, and the estimation type specified by str_out$type (derived from the above)
  DF1_imp <- mice(DF1, m = N.Imp, method = str_out$type )
  
  # Print the first 50 logged events, if they occur 
  print(head(DF1_imp$loggedEvents, 50))
  
  # Return the imputed data
  return(DF1_imp)
}

### Conduct the imputation ### 
# Conduct the imputation using DF.GP.1 (excluding country code)
mice.imp.GP_int <- MI_function_1(DF = DF.GP.1_imp)

# Save the imputed data
save(mice.imp.GP_int, file = "mice.imp.GP_int.Rdata")


###


###### 3) Post imputation manipulation ######
### Function to add goal endorsement and satisfaction variable, SO_CountryCode & respondent ID, change the reference level for region
MI_function_2 <-function(DF1, DF_imp) {
  
  ### Number of imputed DF ###
  # 10 imputed DF
  N.Imp = 10
  
  ### Drop "other" as a category in gender (there were too few observations to model)
  DF1<-DF1[!(DF1$gender=="Other"),]
  DF1<- droplevels(DF1)
  
  ### Extract each imputed dataset and perform additional manipulation ###
  # Create a list to store the imputed datasets 
  mice.imp <- list()
  
  # For i in each dataset
  for(i in 1:N.Imp) {
    
    ### Extract the imputed data
    mice.imp[[i]] <- mice::complete(DF_imp, action= i, inc=FALSE)
    
    ### Add data from DF1 
    
    # Add goal endorsement and progress
    mice.imp[[i]] <- cbind(mice.imp[[i]], DF1[c("GP_1_a", "GP_2_a", "GP_3_a", "GP_4_a","GP_5_a", "GP_6_a", "GP_7_a", "GP_8_a", "GP_9_a","GP_10_a",
                                                "GP_1_b", "GP_2_b", "GP_3_b", "GP_4_b","GP_5_b", "GP_6_b", "GP_7_b", "GP_8_b", "GP_9_b","GP_10_b")])
    
    # Add country code 
    mice.imp[[i]]$SO_CountryCode <- DF1$SO_CountryCode
    
    # Add respondent ID
    mice.imp[[i]]$ID <- DF1$ID
    
    # Change the reference level for region
    mice.imp[[i]]$SO_Region <- as.factor(mice.imp[[i]]$SO_Region)
    mice.imp[[i]]$SO_Region <- relevel(mice.imp[[i]]$SO_Region, ref="Europe and Northern America")
    
    # Change the reference level for Environment
    mice.imp[[i]]$Environment <- as.factor(mice.imp[[i]]$Environment)
    mice.imp[[i]]$Environment <- relevel(mice.imp[[i]]$Environment, ref="Terrestrial")
    
    # Add the 'Conservation' variable - which indicates if individuals currently consider themselves to be working or studying in conservation 
    mice.imp[[i]]$Conservation <- DF1$Conservation
  }
  
  # Return the manipulated DF 
  return(mice.imp)
}

# Implement the function 
mice.imp.GP <- MI_function_2(DF1 = DF.GP.1, DF_imp = mice.imp.GP_int)

# Save the data 
save(mice.imp.GP, file = "mice.imp.GP.Rdata")


###


###### 4) Estimate plausible values of dispositional optimism ###### 

# The structure of the DO model 
model_OP_method <- '
###### Dispositional optimism 
# Dispositional optimism 
OP =~ LOTR_1 + LOTR_2 + LOTR_3 +  LOTR_4 + LOTR_5 + LOTR_6

# The method effect 
method =~  LOTR_1 + LOTR_3 + LOTR_6

# OP and the method effect are orthognal 
OP ~~0*method
'

# LOTR row names
LOTR_row_name <- c("LOTR_1","LOTR_2","LOTR_3","LOTR_4","LOTR_5","LOTR_6")

# Convert to numeric, across all imputed DF 
for (i in seq_along(1:length(mice.imp.GP))){
  mice.imp.GP[[i]][,LOTR_row_name] <- apply(mice.imp.GP[[i]][,LOTR_row_name], 2, as.numeric)
}

# Run the SEM on all imputed DF 
model_GP_1_MI <- runMI(model_OP_method, data = mice.imp.GP, fun="sem", estimator = "WLSMVS", FUN = fitMeasures)

# Extract the plausible values from each imputed DF - we're doing just one draw from each of the ten imputed DF given the computational time 
pv_DO <- plausibleValues(model_GP_1_MI, nDraws = 1, seed = 123)

# Match each PV to each DF - technically, it should be ten PV for each imputed datasets (totalling 100), but this is likely massive overkill given we have complete cases for the LOTR questions, and means the resulting computation is too much. 
mice.imp.GP.sim <- mice.imp.GP
for (i in seq_along(1:length(mice.imp.GP))){
  mice.imp.GP.sim[[i]]$DO <- pv_DO[[i]]$OP
}

# Save this DF
save(mice.imp.GP.sim, file ="mice.imp.GP.sim.RData")
