################################################################
######### Life in Conservation: Goal progress analysis #########
################################################################

# The following describes the primary analysis steps used to generate the results presented in the article 'Balancing making a difference with making a living in the conservation sector'.  

# Users can substitute 'XXXX' (as shown below) with their own AWS S3 bucket details, in order to run this code, where appropriate.
# This code was run using an R Studio Server Amazon Machine Image (see https://www.louisaslett.com/RStudio_AMI/), using an AWS EC2 instance type with 61 GiB of memory, 8 vCPUs, and a 64-bit platform (instance type name 'r4.2xlarge').


######### Steps #########
# 1) Set up the environment 
# 2) Preparing the data and specifying the model 
# 3) Run test models using the first of the ten imputed datasets
# 4) Perform model diagnostics following the WAMBS check list 
# 5) Run the model on all 10 datasets
# 6) Prepare the results including making predictions and creating figures 
# 7) Repeat the analysis with a conservative definition of conservationists 
# 8) Repeat the analysis with the inclusion of a COVID-19 dummy variable 
# 9) Repeat the analysis with the inclusion of a work hours dummy variable
# 10) Repeat the analysis disaggregating the 'collective' goals into those relating to altruistic and biospheric values
# 11) Characteristics associated with goal endorsement 
# 12) Association between satisfaction and goal type

######### 1) Set up the environment #########

# Set seed
set.seed(123)

# Libraries 
library(dplyr)
library(plyr)
library(parallel)
library(brms)
library(purrr)
library(aws.s3)
library(psych)
library(ggmcmc)
library(coda)
library(reshape2)
library(bayesplot)
library(ggridges)
library(bayestestR)
library(tidybayes)
library(tidyverse)
library(ggpubr)

# Set AWS settings 
s3BucketName <- "XXXX"
Sys.setenv("AWS_ACCESS_KEY_ID" = "XXXX",
           "AWS_SECRET_ACCESS_KEY" = "XXXX",
           "AWS_DEFAULT_REGION" = "XXXX")

# How many cores
detectCores() 
options(mc.cores = parallel::detectCores())

# Load data  (these files should be in the working directory)
load("mice.imp.GP.sim.RData")
load("DF.GP.1.Rdata") 

# My pallet 
mypal <- c("#fe9929", "#006837")
mypal.bw <- c("#bdbdbd", "#636363")

# ggplot themes 
theme_set(theme_pubr()) 

######### 2) Preparing the data and specifying the model #########

# A function that returns complete cases for each goal, progress towards each goal, the name of that goal, and converting back from numeric to ordered factor 
complete.goals <- function(DF, goal){
  DF_goal <- DF[complete.cases(DF[,goal]),]
  DF_goal$goal_num <- as.factor(goal) 
  DF_goal$goal_prog <- DF_goal[,goal]
  DF_goal$goal_prog <- as.ordered(DF_goal$goal_prog)
  return(DF_goal)
}

# Lists of the stacked DF (stacked by complete case of each goal) for each MI dataset - individual goals
mice.imp_individual.b <- list_along(1:length(mice.imp.GP.sim))

for (i in seq_along(1:10)){
  # Personal 
  mice.imp_individual.b[[i]] <- rbind(complete.goals(mice.imp.GP.sim[[i]], "GP_1_b" ),
                                      complete.goals(mice.imp.GP.sim[[i]], "GP_2_b" ), 
                                      complete.goals(mice.imp.GP.sim[[i]], "GP_3_b" ), 
                                      complete.goals(mice.imp.GP.sim[[i]], "GP_4_b" ))
  mice.imp_individual.b[[i]]$goal_prog <- as.ordered(mice.imp_individual.b[[i]]$goal_prog)
  mice.imp_individual.b[[i]]$goal_prog <-  plyr::mapvalues(mice.imp_individual.b[[i]]$goal_prog, from = c("1", "2","3", "4", "5"), to = c("Very dissatisfied" , "Dissatisfied" , "Neutral" , "Satisfied" , "Very satisfied"))
  mice.imp_individual.b[[i]]$SO_Region <- relevel(mice.imp_individual.b[[i]]$SO_Region, ref = "Europe and Northern America") 
  mice.imp_individual.b[[i]]$DO <- as.numeric(scale(mice.imp_individual.b[[i]]$DO, scale = T, center = T))
}

# Lists of the stacked DF (stacked by complete case of each goal) for each MI dataset - collective goals
mice.imp_collective.b <- list_along(1:length(mice.imp.GP.sim))

for (i in seq_along(1:10)){
  # Personal 
  mice.imp_collective.b[[i]] <- rbind(complete.goals(mice.imp.GP.sim[[i]], "GP_5_b" ),
                                      complete.goals(mice.imp.GP.sim[[i]], "GP_6_b" ),
                                      complete.goals(mice.imp.GP.sim[[i]], "GP_7_b" ),
                                      complete.goals(mice.imp.GP.sim[[i]], "GP_8_b" ),
                                      complete.goals(mice.imp.GP.sim[[i]], "GP_9_b" ),
                                      complete.goals(mice.imp.GP.sim[[i]], "GP_10_b"))
  mice.imp_collective.b[[i]]$goal_prog <- as.ordered(mice.imp_collective.b[[i]]$goal_prog)
  mice.imp_collective.b[[i]]$goal_prog <-  plyr::mapvalues(mice.imp_collective.b[[i]]$goal_prog, from = c("1", "2","3", "4", "5"), to = c("Very dissatisfied" , "Dissatisfied" , "Neutral" , "Satisfied" , "Very satisfied"))
  mice.imp_collective.b[[i]]$SO_Region <- relevel(mice.imp_collective.b[[i]]$SO_Region, ref = "Europe and Northern America") 
  mice.imp_collective.b[[i]]$DO <- as.numeric(scale(mice.imp_collective.b[[i]]$DO, scale = T, center = T))
}


### Checking positive covariance between all variables ###

# Select the first of the imputed datasets 
mice.imp.GP.one <- mice.imp.GP.sim[[1]]

# Convert from ordered factor to numeric
GP_row_name <- c("GP_1_b", "GP_2_b","GP_3_b","GP_4_b", "GP_5_b","GP_6_b" , "GP_7_b","GP_8_b", "GP_9_b" , "GP_10_b" )
mice.imp.GP.one[GP_row_name] <- apply(mice.imp.GP.one[GP_row_name],2, as.numeric)

### Polychoric correlation between variables ### 
poly_cor = psych::polychoric(mice.imp.GP.one[GP_row_name], na.rm=TRUE)
cor.plot(poly_cor$rho, numbers=F, upper=FALSE, main = "Polychoric correlation", show.legend = T)


### Model set-up ###
# Burn in and iterations for main model 
burn_in <- 4000
iterations <- 8000

# Burn in and iterations for double model
burn_in.d <- 8000
iterations.d <- 16000

# Set prior 
prior_all <- c(
  set_prior(prior = "normal(0,10)", class = "Intercept"),
  set_prior(prior = "normal(0,10)", class = "b"),
  set_prior(prior = "student_t(3, 0, 2.5)", class = "sd"))

# The model formula
brm_form.b <- formula("goal_prog ~ goal_num + years_cons_scaled + WH_scaled + DO  + gender + education_simple + SO_Region + position_simple + Environment + (1|ID)")

######### 3) Run test models using the first of the ten imputed datasets #########

### The test models used in the diagnostics 
# Single  - individual 
start_time_individual <- Sys.time()
individual_model.b <- brm(brm_form.b, family = cumulative(link = "logit", threshold = "flexible"), data = mice.imp_individual.b[[1]], warmup = burn_in,  iter  = iterations,  chains = 4, seed = 123, prior = prior_all, control = list(adapt_delta = 0.95), refresh = 0) 
save(individual_model.b, file ="individual_model.b.RData")
put_object(file = "individual_model.b.RData", object = "individual_model.b.RData", bucket = "XXXX", multipart=TRUE)
try(remove(list=c("individual_model.b")))
try(unlink("individual_model.b.RData"))
end_time_individual <- Sys.time()
end_time_individual - start_time_individual

# Single  - collective 
start_time_collective <- Sys.time()
collective_model.b <- brm(brm_form.b, family = cumulative(link = "logit", threshold = "flexible"), data = mice.imp_collective.b[[1]], warmup = burn_in,  iter  = iterations,  chains = 4, seed = 123,  prior = prior_all, control = list(adapt_delta = 0.95), refresh = 0)
save(collective_model.b, file ="collective_model.b.RData")
put_object(file = "collective_model.b.RData", object = "collective_model.b.RData", bucket = "XXXX", multipart=TRUE)
try(remove(list=c("collective_model.b")))
try(unlink("collective_model.b.RData"))
end_time_collective <- Sys.time()
end_time_collective - start_time_collective

### The test models with double the iterations, used in the diagnostics 

# Double - individual 
start_time_individual.d <- Sys.time()
individual_model.b.d <- brm(brm_form.b, family = cumulative(link = "logit", threshold = "flexible"), data = mice.imp_individual.b[[1]], warmup = burn_in.d,  iter  = iterations.d,  chains = 4, seed = 123,  prior = prior_all, control = list(adapt_delta = 0.95), refresh = 0)
save(individual_model.b.d, file ="individual_model.b.d.RData")
put_object(file = "individual_model.b.d.RData", object = "individual_model.b.d.RData", bucket = "XXXX", multipart=TRUE)
try(remove(list=c("individual_model.b.d")))
try(unlink("individual_model.b.d.RData"))
end_time_individual.d <- Sys.time()
end_time_individual.d - start_time_individual.d

# Double - collective 
start_time_collective.d <- Sys.time()
collective_model.b.d <- brm(brm_form.b, family = cumulative(link = "logit", threshold = "flexible"), data = mice.imp_collective.b[[1]], warmup = burn_in.d,  iter  = iterations.d,  chains = 4, seed = 123,  prior = prior_all, control = list(adapt_delta = 0.95), refresh = 0)
save(collective_model.b.d, file ="collective_model.b.d.RData")
put_object(file = "collective_model.b.d.RData", object = "collective_model.b.d.RData", bucket = "XXXX", multipart=TRUE)
try(remove(list=c("collective_model.b.d")))
try(unlink("individual_model.b.d.RData"))
end_time_collective.d <- Sys.time()
end_time_collective.d - start_time_collective.d

######### 4) Perform model diagnostics following the WAMBS check list #########

### Further details describing the diagnostics can be found by following the below links.
# - Individual model:
# - Collective model: 

###### Diagnostics for the individual model ######

# Load the personal data 
s3load("individual_model.b.RData", bucket = "XXXX")
# load("individual_model.b.RData")

### Examining convergence ### 
# How many parameters are we looking at?
param_num <- individual_model.b$fit@par_dims$b

# How many burn ins and iterations
burn_in <- individual_model.b$fit@sim$warmup
iterations <- individual_model.b$fit@sim$iter

# Checking the trace plots for each model 
individual_model_plot.b <- ggmcmc::ggs(individual_model.b)

# Plot 
ggplot(filter(individual_model_plot.b, Parameter %in% c(names(individual_model.b$fit)[1:param_num])),
       aes(x   = Iteration,
           y   = value, 
           col = as.factor(Chain)))+
  geom_line() +
  geom_vline(xintercept = burn_in)+
  facet_grid(Parameter ~ . ,
             scale  = 'free_y',
             switch = 'y')+
  labs(title = "Trace-plots: Individual goal", 
       col   = "Chains")

# Convert to MCMC object 
individual_modelposterior.b <- as.mcmc(individual_model.b)

# Gelman-Rubin Diagnostic
gelman.plot(individual_modelposterior.b[,1:param_num])

#  Geweke diagnostic
geweke.plot(individual_modelposterior.b[, 1:param_num])

### 3) Examination of convergence after doubling the number of iterations ###

# Load double model
s3load("individual_model.b.d.RData", bucket = "XXXX")
# load("individual_model.b.d.RData")

# Burn in and iterations for double model 
burn_in.d <- individual_model.b.d$fit@sim$warmup
iterations.d <- individual_model.b.d$fit@sim$iter

# Checking the trace plots for each model 
individual_model_plot_d.b <- ggmcmc::ggs(individual_model.b.d)

# Plot 
ggplot(filter(individual_model_plot_d.b, Parameter %in% c(names(individual_model.b.d$fit)[1:param_num])),
       aes(x   = Iteration,
           y   = value, 
           col = as.factor(Chain)))+
  geom_line() +
  geom_vline(xintercept = burn_in)+
  facet_grid(Parameter ~ . ,
             scale  = 'free_y',
             switch = 'y')+
  labs(title = "Trace-plots (double): Individual centered goal", 
       col   = "Chains")

# Convert to MCMC object 
individual_model_posterior_d <- as.mcmc(individual_model.b.d)

# Gelman-Rubin Diagnostic
gelman.plot(individual_model_posterior_d[,1:param_num])

# Geweke diagnostic
geweke.plot(individual_model_posterior_d[, 1:param_num])

# Remove used objects 
try(remove(list=c("individual_model.b.d")))
try(unlink("individual_model.b.d.RData"))

### 4) Examine the histogram of estimates in each iteration to determine if there were enough iterations to approximate the posterior ### 

# Histogram 
mcmc_plot(individual_model.b, type = "hist")

### 5) Examine the degree of autocorrelation between the chains ###

# Calculate lag
auto_lag <- data.frame(autocorr.diag(individual_model.b[,1:param_num], lags = c(0, 1,2,3,4, 5, 10, 50)), lag = factor(c("0", "1","2","3","4", "5", "10", "50"), levels = c("0", "1","2","3","4", "5", "10", "50")))

# Melt
auto_lag <- melt(auto_lag)

# Autocorrelation plot 
ggplot() + geom_pointrange(data=auto_lag, mapping=aes(x=lag, y=value, ymin=0, ymax=value)) + geom_hline(yintercept = 0)  +  facet_wrap(~ variable, ncol = 3) 

### 6) Evaluate if the posterior distribution is sensible ### 

# Personal 
mcmc_plot(individual_model.b, type = "dens", pars = "^b_")

# Remove used objects 
try(remove(list=c("individual_model.b")))
try(unlink("individual_model.b.RData"))

###### Diagnostics for the collective model ######

# Load the personal data 
s3load("collective_model.b.RData", bucket = "XXXX")
# load("collective_model.b.RData")

### Examining convergence ### 
# How many parameters are we looking at?
param_num <- collective_model.b$fit@par_dims$b

# How many burn ins and iterations
burn_in <- collective_model.b$fit@sim$warmup
iterations <- collective_model.b$fit@sim$iter

# Checking the trace plots for each model 
collective_model_plot.b <- ggmcmc::ggs(collective_model.b)

# Plot 
ggplot(filter(collective_model_plot.b, Parameter %in% c(names(collective_model.b$fit)[1:param_num])),
       aes(x   = Iteration,
           y   = value, 
           col = as.factor(Chain)))+
  geom_line() +
  geom_vline(xintercept = burn_in)+
  facet_grid(Parameter ~ . ,
             scale  = 'free_y',
             switch = 'y')+
  labs(title = "Trace-plots: collective goal", 
       col   = "Chains")

# Convert to MCMC object 
collective_modelposterior.b <- as.mcmc(collective_model.b)

# Gelman-Rubin Diagnostic
gelman.plot(collective_modelposterior.b[,1:param_num])

#  Geweke diagnostic
geweke.plot(collective_modelposterior.b[, 1:param_num])

### 3) Examination of convergence after doubling the number of iterations ###

# Load double model
s3load("collective_model.b.d.RData", bucket = "XXXX")
# load("collective_model.b.d.RData")

# Burn in and iterations for double model 
burn_in.d <- collective_model.b.d$fit@sim$warmup
iterations.d <- collective_model.b.d$fit@sim$iter

# Checking the trace plots for each model 
collective_model_plot_d.b <- ggmcmc::ggs(collective_model.b.d)

# Plot 
ggplot(filter(collective_model_plot_d.b, Parameter %in% c(names(collective_model.b.d$fit)[1:param_num])),
       aes(x   = Iteration,
           y   = value, 
           col = as.factor(Chain)))+
  geom_line() +
  geom_vline(xintercept = burn_in)+
  facet_grid(Parameter ~ . ,
             scale  = 'free_y',
             switch = 'y')+
  labs(title = "Trace-plots (double): collective centered goal", 
       col   = "Chains")

# Convert to MCMC object 
collective_model_posterior_d <- as.mcmc(collective_model.b.d)

# Gelman-Rubin Diagnostic
gelman.plot(collective_model_posterior_d[,1:param_num])

# Geweke diagnostic
geweke.plot(collective_model_posterior_d[, 1:param_num])

# Remove used objects 
try(remove(list=c("collective_model.b.d")))
try(unlink("collective_model.b.d.RData"))

### 4) Examine the histogram of estimates in each iteration to determine if there were enough iterations to approximate the posterior ### 

# Histogram 
mcmc_plot(collective_model.b, type = "hist")

### 5) Examine the degree of autocorrelation between the chains ###

# Calculate lag
auto_lag <- data.frame(autocorr.diag(collective_model.b[,1:param_num], lags = c(0, 1,2,3,4, 5, 10, 50)), lag = factor(c("0", "1","2","3","4", "5", "10", "50"), levels = c("0", "1","2","3","4", "5", "10", "50")))

# Melt
auto_lag <- melt(auto_lag)

# Autocorrelation plot 
ggplot() + geom_pointrange(data=auto_lag, mapping=aes(x=lag, y=value, ymin=0, ymax=value)) + geom_hline(yintercept = 0)  +  facet_wrap(~ variable, ncol = 3) 

### 6) Evaluate if the posterior distribution is sensible ### 

# Personal 
mcmc_plot(collective_model.b, type = "dens", pars = "^b_")

# Remove used objects 
try(remove(list=c("collective_model.b")))
try(unlink("collective_model.b.RData"))


######### 5) Run the model on all 10 datasets #########

###### Individual ######

### Run the model ###
individual_model_full.b <- brm_multiple(brm_form.b, family = cumulative(link = "logit", threshold = "flexible"), data = mice.imp_individual.b, warmup = burn_in,  iter   = iterations,  chains = 4, seed = 123,  prior = prior_all, control = list(adapt_delta = 0.95))
save(individual_model_full.b, file ="individual_model_full.b.RData")
put_object(file = "individual_model_full.b.RData", object = "individual_model_full.b.RData", bucket = "XXXX", multipart=TRUE)

# Extract the summary information
individual_model_full.b.summary <- summary(individual_model_full.b)

## Check for false convergence warning - individual should be >.99 & < 1.01
table(apply(data.frame(individual_model_full.b$rhats[1:18]), 2, min)<.99)
table(apply(data.frame(individual_model_full.b$rhats[1:18]), 2, max)>1.01)

# Save the summary data and put it in the S3 bucket
save(individual_model_full.b.summary, file ="individual_model_full.b.summary.RData")
put_object(file = "individual_model_full.b.summary.RData", object = "individual_model_full.b.summary.RData", bucket = "XXXX", multipart=TRUE)

# Remove the file after use
try(file.remove("individual_model_full.b.summary.RData"))
try(file.remove("individual_model_full.b.RData"))

# Remove used objects
try(remove(list=c("individual_model_full.b")))


###### Collective ######

### Run the model ###
collective_model_full.b <- brm_multiple(brm_form.b, family = cumulative(link = "logit", threshold = "flexible"), data = mice.imp_collective.b, warmup = burn_in,  iter   = iterations,  chains = 4, seed = 123,  prior = prior_all, control = list(adapt_delta = 0.95))
save(collective_model_full.b, file ="collective_model_full.b.RData")
put_object(file = "collective_model_full.b.RData", object = "collective_model_full.b.RData", bucket = "XXXX", multipart=TRUE)

# Extract the summary information
collective_model_full.b.summary <- summary(collective_model_full.b)

## Check for false convergence warning - collective should be >.99 & < 1.01
table(apply(data.frame(collective_model_full.b$rhats[1:18]), 2, min)<.99)
table(apply(data.frame(collective_model_full.b$rhats[1:18]), 2, max)>1.01)

# Save the summary data and put it in the S3 bucket
save(collective_model_full.b.summary, file ="collective_model_full.b.summary.RData")
put_object(file = "collective_model_full.b.summary.RData", object = "collective_model_full.b.summary.RData", bucket = "XXXX", multipart=TRUE)

# Remove the file after use
try(file.remove("collective_model_full.b.summary.RData"))
try(file.remove("collective_model_full.b.RData"))

# Remove used objects
try(remove(list=c("collective_model_full.b")))


######### 6) Prepare the results including making predictions and creating figures ######### 

###### Prediction ######

### Years in conservation ###
years_5   <- mean(DF.GP.1[DF.GP.1$years_cons == 5,]$years_cons_scaled, na.rm = T)
years_30  <- mean(DF.GP.1[DF.GP.1$years_cons == 30,]$years_cons_scaled, na.rm = T)

### Dispositional optimism ###
DO_1Sbelow <- mean(mice.imp_individual.b[[1]]$DO) - sd(mice.imp_individual.b[[1]]$DO)
DO_1SDabov <- mean(mice.imp_collective.b[[1]]$DO) + sd(mice.imp_collective.b[[1]]$DO)

### Work hours ###
WH_20   <- mean(DF.GP.1[DF.GP.1$WH == 20,]$WH_scaled, na.rm = T)
WH_40  <- mean(DF.GP.1[DF.GP.1$WH == 40,]$WH_scaled, na.rm = T)

# function to predict data
data_predict<- function(DF){
  return(data.frame(goal_num = c(rep(levels(DF$goal_num)[1],9)),
                    years_cons_scaled =  c(rep(years_5,1), rep(years_30,1), rep(years_5,7)),  # Years in conservation = 5 & 30 in SD
                    DO = c(rep(DO_1Sbelow, 2), rep(DO_1SDabov, 1) ,rep(DO_1Sbelow, 6)), # DO = mean and SD 
                    WH_scaled = c(rep(WH_20, 3), rep(WH_40,1), rep(WH_20, 5)), # Work hours = 20 & 40  in SD
                    gender = c(rep("Female",4), rep("Male",1), rep("Female",4)),
                    education_simple = c(rep("Non-university",5), rep("University",1), rep("Non-university",3)),
                    SO_Region = c(rep("Europe and Northern America",6), rep("Oceania",1), rep("Europe and Northern America",2)),
                    position_simple = c(rep("Academic",7), rep("Practice",1), rep("Academic",1)), 
                    Environment = c(rep("Terrestrial", 8), rep("Marine", 1))
  ))
}


# Function to extract fitted draws for satisfied (4) and very satisfied (5)
predict_funct <-  function(model, dummy_DF, prob_int) {
  
  # Create the conditional probability name 
  prob_name <- paste0("P(goal_prog|",  prob_int, ")")
  
  # Extract the fitted draws 
  out <- tibble(dummy_DF) %>%
    add_fitted_draws(model, re_formula = NA, scale = "response", value = prob_name, n=1000)
  
  # Subset to cat 4 and 5 
  out_4 <- out[ which(out$.category == 4), ]
  out_5 <- out[ which(out$.category == 5), ]
  
  # Create the predicted DF 
  predicted_df <- data.frame("prob_int" = c(out_4[prob_int]), "pred" = c(out_4[prob_name] + out_5[prob_name]))
  
  # Return
  return(predicted_df)
}

# How many increments in the predictor variable 
n_rows = 100 
n_rows_bin = 2
n_rows_cat = length(levels(mice.imp_individual.b[[1]]$SO_Region))

# Create dummy DF
dummy_df_individual.b <- data_predict(mice.imp_individual.b[[1]])
dummy_df_collective.b <- data_predict(mice.imp_collective.b[[1]])


######  Predicted outcomes ######
# Load models 
# # Use s3load after models are developed
s3load("individual_model_full.b.RData", bucket = "XXXX")
s3load("collective_model_full.b.RData", bucket = "XXXX")

### Individual predicted outcomes
# Predicted outcomes
individual.b.pred <- predict(individual_model_full.b, newdata = dummy_df_individual.b, re_formula = NA,
                             scale = "response",  probs = c(0.025, 0.975))

# Combine rows
dummy_df_individual.b <- cbind(dummy_df_individual.b, individual.b.pred )

# Save and remove 
save(dummy_df_individual.b, file ="dummy_df_individual.b.RData")
put_object(file = "dummy_df_individual.b.RData", object = "dummy_df_individual.b.RData", bucket = "XXXX", multipart=TRUE)
unlink("dummy_df_individual.b.RData")

### Collective predicted outcomes
# Predicted outcomes
collective.b.pred <- predict(collective_model_full.b, newdata = dummy_df_collective.b, re_formula = NA,
                             scale = "response",  probs = c(0.025, 0.975))

# Combine rows
dummy_df_collective.b <- cbind(dummy_df_collective.b, collective.b.pred )

# Save and remove 
save(dummy_df_collective.b, file ="dummy_df_collective.b.RData")
put_object(file = "dummy_df_collective.b.RData", object = "dummy_df_collective.b.RData", bucket = "XXXX", multipart=TRUE)
unlink("dummy_df_collective.b.RData")

###### Figures ######

# Extract 
individual_model_full.b.summary <- summary(individual_model_full.b)
collective_model_full.b.summary <- summary(collective_model_full.b)

# Extracting coefficients
coef_list_b <- list(
  
  # Personal 
  data.frame(Goal = "Individual",Variable=as.factor(rownames(individual_model_full.b.summary$fixed)), (individual_model_full.b.summary$fixed[,c(1,3,4)])),
  
  # People 
  data.frame(Goal = "Collective", Variable= as.factor(rownames(collective_model_full.b.summary$fixed)),  (collective_model_full.b.summary$fixed[,c(1,3,4)]))
)

# Create single DF
coef_df_b <- do.call(rbind, coef_list_b)

# Variables to retain
keep <- c("years_cons_scaled","DO","WH_scaled","genderMale", "EnvironmentCC", "EnvironmentMarine" , "education_simpleUniversity",
          "SO_RegionCentralandSouthernAsia" ,
          "SO_RegionEasternandSouthMEasternAsia" ,
          "SO_RegionSubMSaharanAfrica" ,
          "SO_RegionLatinAmericaandtheCaribbean" ,
          "SO_RegionNorthernAfricaandWesternAsia" ,
          "SO_RegionOceania",
          "position_simplePractice" )

# Keep variables and change to charecters
coef_df_b.s <- coef_df_b[ which(coef_df_b$Variable %in% keep), ]
coef_df_b.s$Variable <- as.character(coef_df_b.s$Variable)

# Recode variable 
coef_df_b.s$Variable <- recode(coef_df_b.s$Variable, "years_cons_scaled" = "Years in conservation",
                               "DO" = "Dispositional optimism",
                               "WH_scaled" = "Work hours",
                               "genderMale" = "Male",
                               "education_simpleUniversity" = "University",
                               "SO_RegionCentralandSouthernAsia" = "C. & S. Asia",
                               "SO_RegionEasternandSouthMEasternAsia" = "E. & S.E. Asia",
                               "SO_RegionSubMSaharanAfrica" = "Sub-Saharan Africa",
                               "SO_RegionLatinAmericaandtheCaribbean" = "Lat. America & Carib.",
                               "SO_RegionNorthernAfricaandWesternAsia" = "N. Africa & W. Asia",
                               "SO_RegionOceania"= "Oceania",
                               "position_simplePractice" = "Practice",
                               "EnvironmentCC" = "Cross-cutting", 
                               "EnvironmentMarine" = "Marine")



# Create reference level DF
Ref_level <- data.frame(Goal = c(rep("Reference", 5)), Variable = c("Female", "Non-university", "Europe & N. America", "Terrestrial", "Academia") , Estimate = c(rep(0, 5)), l.95..CI= c(rep(0, 5)),  u.95..CI= c(rep(0, 5)))

# Combine 
coef_df_b.s <- rbind(coef_df_b.s,Ref_level )

# Re-level 
coef_df_b.s$Variable <- factor(coef_df_b.s$Variable, levels = (c("Dispositional optimism","Years in conservation", "Work hours", "Female", "Male" ,
                                                                 "Non-university", "University", "Academia", "Practice", "Terrestrial", "Cross-cutting" ,
                                                                 "Marine", "Europe & N. America", "C. & S. Asia", "E. & S.E. Asia", "Sub-Saharan Africa",
                                                                 "Lat. America & Carib.", "N. Africa & W. Asia", "Oceania" )))

# Plot
Figure_2 <- ggplot(coef_df_b.s, aes(x= Variable, y= Estimate)) + ylab("Goal progress satisfaction (log-odds)")  + 
  geom_pointrange(aes(color = Goal, shape = Goal, ymax = u.95..CI, ymin= l.95..CI),
                  position=position_dodge(width=c(0.6,0.6)), fatten = 2) +
  geom_hline(yintercept = 0, linetype="dashed", colour="grey73")  +  coord_flip()    + scale_color_manual(values = c(mypal[1], mypal[2], "black" ) ) +  scale_shape_manual(values = c(16, 16,  1)  ) + 
  scale_x_discrete(breaks = rev(c(levels(coef_df_b.s$Variable))), 
                   limits = rev(c(levels(coef_df_b.s$Variable)[1:3], "skip", 
                                  levels(coef_df_b.s$Variable)[4:5], "skip", 
                                  levels(coef_df_b.s$Variable)[6:7], "skip",
                                  levels(coef_df_b.s$Variable)[8:9], "skip", 
                                  levels(coef_df_b.s$Variable)[10:12], "skip",
                                  levels(coef_df_b.s$Variable)[13:19]))) + 
  theme(plot.margin=unit(c(0,0,1.5,0),"cm"))+
  theme(legend.position=c(.72,-.15), legend.direction="horizontal")

# Save 
ggsave("Figure_2.jpeg", width = 120, height = 140, units = "mm",   dpi = 800)
put_object(file = "Figure_2.jpeg", object = "Figure_2.jpeg", bucket = "XXXX", multipart=TRUE)
unlink("Figure_2.jpeg")



# Plot
Figure_2.BW <- ggplot(coef_df_b.s, aes(x= Variable, y= Estimate)) + ylab("Goal progress satisfaction (log-odds)")  + 
  geom_pointrange(aes(color = Goal, shape = Goal, ymax = u.95..CI, ymin= l.95..CI),
                  position=position_dodge(width=c(0.6,0.6)), fatten = 2) +
  geom_hline(yintercept = 0, linetype="dashed", colour="grey73")  +  coord_flip()    + scale_color_manual(values = c(mypal.bw[1], mypal.bw[2], "black" ) ) +  scale_shape_manual(values = c(16, 16,  1)  ) + 
  scale_x_discrete(breaks = rev(c(levels(coef_df_b.s$Variable))), 
                   limits = rev(c(levels(coef_df_b.s$Variable)[1], "skip", 
                                  levels(coef_df_b.s$Variable)[2], "skip", 
                                  levels(coef_df_b.s$Variable)[3], "skip", 
                                  levels(coef_df_b.s$Variable)[4:5], "skip", 
                                  levels(coef_df_b.s$Variable)[6:7], "skip",
                                  levels(coef_df_b.s$Variable)[8:9], "skip", 
                                  levels(coef_df_b.s$Variable)[10:12], "skip",
                                  levels(coef_df_b.s$Variable)[13:19]))) + 
  theme(plot.margin=unit(c(0,0,1.5,0),"cm"))+
  theme(legend.position=c(.72,-.15), legend.direction="horizontal")

# Save 
ggsave("Figure_2.BW.jpeg", plot = Figure_2,  width = 120, height = 150, units = "mm",   dpi = 800)
put_object(file = "Figure_2.BW.jpeg", object = "Figure_2.BW.jpeg", bucket = "XXXX", multipart=TRUE)
unlink("Figure_2.BW.jpeg")

# Save 
ggsave("Figure_2.BW.eps", plot = Figure_2,  width = 120, height = 150, units = "mm",   dpi = 800)
put_object(file = "Figure_2.BW.eps", object = "Figure_2.BW.eps", bucket = "XXXX", multipart=TRUE)
unlink("Figure_2.BW.eps")


######### 7) Repeat the analysis with a conservative definition of conservationists ######### 

# Return ID of those who conformed to the stricter definition of conservationist 
DF.GP.1_def <- DF.GP.1[which(DF.GP.1$Conservation=='Oui' | 
                               DF.GP.1$Conservation=='S<ed>' |
                               DF.GP.1$Conservation=='Sí' |
                               DF.GP.1$Conservation=='Sim' |
                               DF.GP.1$Conservation=='Yes'), ]

# ID's corresponding to the strict definition 
ID_def <- DF.GP.1_def$ID

# Subset the imputed dataset - include only those individuals with the more conservative definition of conservationits 
mice.imp_individual.def <- list()
mice.imp_collective.def <- list()
for (i in seq_along(1:length(mice.imp.GP.sim))){
  mice.imp_individual.def[[i]] <- mice.imp_individual.b[[i]][which(mice.imp_individual.b[[i]]$ID %in% ID_def ), ] 
  mice.imp_collective.def[[i]] <- mice.imp_collective.b[[i]][which(mice.imp_collective.b[[i]]$ID %in% ID_def ), ] 
}

### Repeat the analysis with the first of the ten imputed datasets, using the conservative definition of conservationists 
# Single  - individual.def 
start_time_individual.def <- Sys.time()
individual_model.b.def <- brm(brm_form.b, family = cumulative(link = "logit", threshold = "flexible"), data = mice.imp_individual.def[[1]], warmup = burn_in,  iter  = iterations,  chains = 4, seed = 123, prior = prior_all, control = list(adapt_delta = 0.95), refresh = 0) 
save(individual_model.b.def, file ="individual_model.b.def.RData")
put_object(file = "individual_model.b.def.RData", object = "individual_model.b.def.RData", bucket = "XXXX", multipart=TRUE)
end_time_individual.def <- Sys.time()
end_time_individual.def - start_time_individual.def

# Single  - collective.def 
start_time_collective.def <- Sys.time()
collective_model.b.def <- brm(brm_form.b, family = cumulative(link = "logit", threshold = "flexible"), data = mice.imp_collective.def[[1]], warmup = burn_in,  iter  = iterations,  chains = 4, seed = 123,  prior = prior_all, control = list(adapt_delta = 0.95), refresh = 0)
save(collective_model.b.def, file ="collective_model.b.def.RData")
put_object(file = "collective_model.b.def.RData", object = "collective_model.b.def.RData", bucket = "XXXX", multipart=TRUE)
end_time_collective.def <- Sys.time()
end_time_collective.def - start_time_collective.def

### Create the figures ### 
# Extract 
individual_model_full.b.def.summary <- summary(individual_model.b.def)
collective_model_full.b.def.summary <- summary(collective_model.b.def)

# Extracting coefficients
coef_list_b <- list(
  
  # Personal 
  data.frame(Goal = "Individual",Variable=as.factor(rownames(individual_model_full.b.def.summary$fixed)), (individual_model_full.b.def.summary$fixed[,c(1,3,4)])),
  
  # People 
  data.frame(Goal = "Collective", Variable= as.factor(rownames(collective_model_full.b.def.summary$fixed)),  (collective_model_full.b.def.summary$fixed[,c(1,3,4)]))
)

# Create single DF
coef_df_b <- do.call(rbind, coef_list_b)

# Variables to retain
keep <- c("years_cons_scaled","DO","WH_scaled","genderMale", "EnvironmentCC", "EnvironmentMarine" , "education_simpleUniversity",
          "SO_RegionCentralandSouthernAsia" ,
          "SO_RegionEasternandSouthMEasternAsia" ,
          "SO_RegionSubMSaharanAfrica" ,
          "SO_RegionLatinAmericaandtheCaribbean" ,
          "SO_RegionNorthernAfricaandWesternAsia" ,
          "SO_RegionOceania",
          "position_simplePractice" )

# Keep variables and change to characters
coef_df_b.s <- coef_df_b[ which(coef_df_b$Variable %in% keep), ]
coef_df_b.s$Variable <- as.character(coef_df_b.s$Variable)

# Recode variable 
coef_df_b.s$Variable <- recode(coef_df_b.s$Variable, 
                               "years_cons_scaled" = "Years in conservation",
                               "DO" = "Dispositional optimism",
                               "WH_scaled" = "Work hours",
                               "genderMale" = "Male",
                               "education_simpleUniversity" = "University",
                               "SO_RegionCentralandSouthernAsia" = "C. & S. Asia",
                               "SO_RegionEasternandSouthMEasternAsia" = "E. & S.E. Asia",
                               "SO_RegionSubMSaharanAfrica" = "Sub-Saharan Africa",
                               "SO_RegionLatinAmericaandtheCaribbean" = "Lat. America & Carib.",
                               "SO_RegionNorthernAfricaandWesternAsia" = "N. Africa & W. Asia",
                               "SO_RegionOceania"= "Oceania",
                               "position_simplePractice" = "Practice",
                               "EnvironmentCC" = "Cross-cutting", 
                               "EnvironmentMarine" = "Marine")



# Create reference level DF
Ref_level <- data.frame(Goal = c(rep("Reference", 5)), Variable = c("Female", "Non-university", "Europe & N. America", "Terrestrial", "Academia") , Estimate = c(rep(0, 5)), l.95..CI= c(rep(0, 5)),  u.95..CI= c(rep(0, 5)))

# Combine 
coef_df_b.s <- rbind(coef_df_b.s,Ref_level )

# Re-level 
coef_df_b.s$Variable <- factor(coef_df_b.s$Variable, levels = (c("Dispositional optimism","Years in conservation", "Work hours", "Female", "Male" ,
                                                                 "Non-university", "University", "Academia", "Practice", "Terrestrial", "Cross-cutting" ,
                                                                 "Marine", "Europe & N. America", "C. & S. Asia", "E. & S.E. Asia", "Sub-Saharan Africa",
                                                                 "Lat. America & Carib.", "N. Africa & W. Asia", "Oceania" )))

# Plot
Figure_def <- ggplot(coef_df_b.s, aes(x= Variable, y= Estimate)) + ylab("Goal progress satisfaction (log-odds)")  + 
  geom_pointrange(aes(color = Goal, shape = Goal, ymax = u.95..CI, ymin= l.95..CI),
                  position=position_dodge(width=c(0.6,0.6)), fatten = 2) +
  geom_hline(yintercept = 0, linetype="dashed", colour="grey73")  +  coord_flip()    + scale_color_manual(values = c(mypal[1], mypal[2], "black" ) ) +  scale_shape_manual(values = c(16, 16,  1)  ) + 
  scale_x_discrete(breaks = rev(c(levels(coef_df_b.s$Variable))), 
                   limits = rev(c(levels(coef_df_b.s$Variable)[1:3], "skip", 
                                  levels(coef_df_b.s$Variable)[4:5], "skip", 
                                  levels(coef_df_b.s$Variable)[6:7], "skip",
                                  levels(coef_df_b.s$Variable)[8:9], "skip", 
                                  levels(coef_df_b.s$Variable)[10:12], "skip",
                                  levels(coef_df_b.s$Variable)[13:19]))) + 
  theme(plot.margin=unit(c(0,0,1.5,0),"cm"))+
  theme(legend.position=c(.72,-.15), legend.direction="horizontal")

# Save 
ggsave("Figure_def.jpeg", width = 120, height = 140, units = "mm",   dpi = 800)
put_object(file = "Figure_def.jpeg", object = "Figure_def.jpeg", bucket = "XXXX", multipart=TRUE)
unlink("Figure_def.jpeg")
try(remove(list=c("individual_model.b.def")))
try(unlink("individual_model.b.def.RData"))
try(remove(list=c("collective_model.b.def")))
try(unlink("collective_model.b.def.RData"))

#########  8) Repeat the analysis with the inclusion of a COVID-19 dummy variable ######### 

# The model formula with the COVID-19 dummy variable 
brm_form.b.CV <- formula("goal_prog ~ goal_num + CV + years_cons_scaled + WH_scaled + DO  + gender + education_simple + SO_Region + position_simple + Environment + (1|ID)")

# Single  - individual.CV 
start_time_individual.CV <- Sys.time()
individual_model.b.CV <- brm(brm_form.b.CV, family = cumulative(link = "logit", threshold = "flexible"), data = mice.imp_individual.b[[1]], warmup = burn_in,  iter  = iterations,  chains = 4, seed = 123, prior = prior_all, control = list(adapt_delta = 0.95), refresh = 0) 
save(individual_model.b.CV, file ="individual_model.b.CV.RData")
put_object(file = "individual_model.b.CV.RData", object = "individual_model.b.CV.RData", bucket = "XXXX", multipart=TRUE)
end_time_individual.CV <- Sys.time()
end_time_individual.CV - start_time_individual.CV

# Single  - collective.CV 
start_time_collective.CV <- Sys.time()
collective_model.b.CV <- brm(brm_form.b.CV, family = cumulative(link = "logit", threshold = "flexible"), data = mice.imp_collective.b[[1]], warmup = burn_in,  iter  = iterations,  chains = 4, seed = 123,  prior = prior_all, control = list(adapt_delta = 0.95), refresh = 0)
save(collective_model.b.CV, file ="collective_model.b.CV.RData")
put_object(file = "collective_model.b.CV.RData", object = "collective_model.b.CV.RData", bucket = "XXXX", multipart=TRUE)
end_time_collective.CV <- Sys.time()
end_time_collective.CV - start_time_collective.CV

### Create the figures ### 
# Extract 
individual_model_full.b.CV.summary <- summary(individual_model.b.CV)
collective_model_full.b.CV.summary <- summary(collective_model.b.CV)

# Extracting coefficients
coef_list_b <- list(
  
  # Personal 
  data.frame(Goal = "Individual",Variable=as.factor(rownames(individual_model_full.b.CV.summary$fixed)), (individual_model_full.b.CV.summary$fixed[,c(1,3,4)])),
  
  # People 
  data.frame(Goal = "Collective", Variable= as.factor(rownames(collective_model_full.b.CV.summary$fixed)),  (collective_model_full.b.CV.summary$fixed[,c(1,3,4)]))
)

# Create single DF
coef_df_b <- do.call(rbind, coef_list_b)

# Variables to retain
keep <- c("CV", "years_cons_scaled","DO","WH_scaled","genderMale", "EnvironmentCC", "EnvironmentMarine" , "education_simpleUniversity",
          "SO_RegionCentralandSouthernAsia" ,
          "SO_RegionEasternandSouthMEasternAsia" ,
          "SO_RegionSubMSaharanAfrica" ,
          "SO_RegionLatinAmericaandtheCaribbean" ,
          "SO_RegionNorthernAfricaandWesternAsia" ,
          "SO_RegionOceania",
          "position_simplePractice" )

# Keep variables and change to characters
coef_df_b.s <- coef_df_b[ which(coef_df_b$Variable %in% keep), ]
coef_df_b.s$Variable <- as.character(coef_df_b.s$Variable)

# Recode variable 
coef_df_b.s$Variable <- recode(coef_df_b.s$Variable,
                               "CV" = "COVID-19",
                               "years_cons_scaled" = "Years in conservation",
                               "DO" = "Dispositional optimism",
                               "WH_scaled" = "Work hours",
                               "genderMale" = "Male",
                               "education_simpleUniversity" = "University",
                               "SO_RegionCentralandSouthernAsia" = "C. & S. Asia",
                               "SO_RegionEasternandSouthMEasternAsia" = "E. & S.E. Asia",
                               "SO_RegionSubMSaharanAfrica" = "Sub-Saharan Africa",
                               "SO_RegionLatinAmericaandtheCaribbean" = "Lat. America & Carib.",
                               "SO_RegionNorthernAfricaandWesternAsia" = "N. Africa & W. Asia",
                               "SO_RegionOceania"= "Oceania",
                               "position_simplePractice" = "Practice",
                               "EnvironmentCC" = "Cross-cutting", 
                               "EnvironmentMarine" = "Marine")



# Create reference level DF
Ref_level <- data.frame(Goal = c(rep("Reference", 5)), Variable = c("Female", "Non-university", "Europe & N. America", "Terrestrial", "Academia") , Estimate = c(rep(0, 5)), l.95..CI= c(rep(0, 5)),  u.95..CI= c(rep(0, 5)))

# Combine 
coef_df_b.s <- rbind(coef_df_b.s,Ref_level )

# Re-level 
coef_df_b.s$Variable <- factor(coef_df_b.s$Variable, levels = (c("COVID-19", "Dispositional optimism","Years in conservation", "Work hours", "Female", "Male" ,
                                                                 "Non-university", "University", "Academia", "Practice", "Terrestrial", "Cross-cutting" ,
                                                                 "Marine", "Europe & N. America", "C. & S. Asia", "E. & S.E. Asia", "Sub-Saharan Africa",
                                                                 "Lat. America & Carib.", "N. Africa & W. Asia", "Oceania" )))

# Plot
Figure_CV <- ggplot(coef_df_b.s, aes(x= Variable, y= Estimate)) + ylab("Goal progress satisfaction (log-odds)")  + 
  geom_pointrange(aes(color = Goal, shape = Goal, ymax = u.95..CI, ymin= l.95..CI),
                  position=position_dodge(width=c(0.6,0.6)), fatten = 2) +
  geom_hline(yintercept = 0, linetype="dashed", colour="grey73")  +  coord_flip()    + scale_color_manual(values = c(mypal[1], mypal[2], "black" ) ) +  scale_shape_manual(values = c(16, 16,  1)  ) + 
  scale_x_discrete(breaks = rev(c(levels(coef_df_b.s$Variable))), 
                   limits = rev(c(levels(coef_df_b.s$Variable)[1:4], "skip", 
                                  levels(coef_df_b.s$Variable)[5:6], "skip", 
                                  levels(coef_df_b.s$Variable)[7:8], "skip",
                                  levels(coef_df_b.s$Variable)[9:10], "skip", 
                                  levels(coef_df_b.s$Variable)[11:13], "skip",
                                  levels(coef_df_b.s$Variable)[14:20]))) + 
  theme(plot.margin=unit(c(0,0,1.5,0),"cm"))+
  theme(legend.position=c(.72,-.15), legend.direction="horizontal")

# Save 
ggsave("Figure_CV.jpeg", width = 120, height = 140, units = "mm",   dpi = 800)
put_object(file = "Figure_CV.jpeg", object = "Figure_CV.jpeg", bucket = "XXXX", multipart=TRUE)
unlink("Figure_CV.jpeg")
try(remove(list=c("individual_model.b.CV")))
try(unlink("individual_model.b.CV.RData"))
try(remove(list=c("collective_model.b.CV")))
try(unlink("collective_model.b.CV.RData"))


#########  9) Repeat the analysis with the inclusion of a work hours dummy variable ######### 

# Create a binary variable indicating if WH is less than or greater than 40 hours a week 
mice.imp_individual.b[[1]]$WH_B <- ifelse(mice.imp_individual.b[[1]]$WH_scaled <= WH_40, 0,
                                          ifelse(mice.imp_individual.b[[1]]$WH_scaled > WH_40, 1, NA))


mice.imp_collective.b[[1]]$WH_B <- ifelse(mice.imp_collective.b[[1]]$WH_scaled <= WH_40, 0,
                                          ifelse(mice.imp_collective.b[[1]]$WH_scaled > WH_40, 1, NA))


# The model formula with the COVID-19 dummy variable 
brm_form.b.WHB <- formula("goal_prog ~ goal_num + years_cons_scaled + WH_B + DO  + gender + education_simple + SO_Region + position_simple + Environment + (1|ID)")

# Single  - individual.WHB 
start_time_individual.WHB <- Sys.time()
individual_model.b.WHB <- brm(brm_form.b.WHB, family = cumulative(link = "logit", threshold = "flexible"), data = mice.imp_individual.b[[1]], warmup = burn_in,  iter  = iterations,  chains = 4, seed = 123, prior = prior_all, control = list(adapt_delta = 0.95), refresh = 0) 
save(individual_model.b.WHB, file ="individual_model.b.WHB.RData")
put_object(file = "individual_model.b.WHB.RData", object = "individual_model.b.WHB.RData", bucket = "cor-analysis", multipart=TRUE)
end_time_individual.WHB <- Sys.time()
end_time_individual.WHB - start_time_individual.WHB

# Single  - collective.WHB 
start_time_collective.WHB <- Sys.time()
collective_model.b.WHB <- brm(brm_form.b.WHB, family = cumulative(link = "logit", threshold = "flexible"), data = mice.imp_collective.b[[1]], warmup = burn_in,  iter  = iterations,  chains = 4, seed = 123,  prior = prior_all, control = list(adapt_delta = 0.95), refresh = 0)
save(collective_model.b.WHB, file ="collective_model.b.WHB.RData")
put_object(file = "collective_model.b.WHB.RData", object = "collective_model.b.WHB.RData", bucket = "cor-analysis", multipart=TRUE)
end_time_collective.WHB <- Sys.time()
end_time_collective.WHB - start_time_collective.WHB

### Create the figures ### 
# Extract 
individual_model_full.b.WHB.summary <- summary(individual_model.b.WHB)
collective_model_full.b.WHB.summary <- summary(collective_model.b.WHB)

# Extracting coefficients
coef_list_b <- list(
  
  # Personal 
  data.frame(Goal = "Individual",Variable=as.factor(rownames(individual_model_full.b.WHB.summary$fixed)), (individual_model_full.b.WHB.summary$fixed[,c(1,3,4)])),
  
  # People 
  data.frame(Goal = "Collective", Variable= as.factor(rownames(collective_model_full.b.WHB.summary$fixed)),  (collective_model_full.b.WHB.summary$fixed[,c(1,3,4)]))
)

# Create single DF
coef_df_b <- do.call(rbind, coef_list_b)

# Variables to retain
keep <- c("years_cons_scaled","DO","WH_B","genderMale", "EnvironmentCC", "EnvironmentMarine" , "education_simpleUniversity",
          "SO_RegionCentralandSouthernAsia" ,
          "SO_RegionEasternandSouthMEasternAsia" ,
          "SO_RegionSubMSaharanAfrica" ,
          "SO_RegionLatinAmericaandtheCaribbean" ,
          "SO_RegionNorthernAfricaandWesternAsia" ,
          "SO_RegionOceania",
          "position_simplePractice" )

# Keep variables and change to characters
coef_df_b.s <- coef_df_b[ which(coef_df_b$Variable %in% keep), ]
coef_df_b.s$Variable <- as.character(coef_df_b.s$Variable)

# Recode variable 
coef_df_b.s$Variable <- recode(coef_df_b.s$Variable,
                               "WH_B" = "Work hours (binary)",
                               "years_cons_scaled" = "Years in conservation",
                               "DO" = "Dispositional optimism",
                               "WH_scaled" = "Work hours",
                               "genderMale" = "Male",
                               "education_simpleUniversity" = "University",
                               "SO_RegionCentralandSouthernAsia" = "C. & S. Asia",
                               "SO_RegionEasternandSouthMEasternAsia" = "E. & S.E. Asia",
                               "SO_RegionSubMSaharanAfrica" = "Sub-Saharan Africa",
                               "SO_RegionLatinAmericaandtheCaribbean" = "Lat. America & Carib.",
                               "SO_RegionNorthernAfricaandWesternAsia" = "N. Africa & W. Asia",
                               "SO_RegionOceania"= "Oceania",
                               "position_simplePractice" = "Practice",
                               "EnvironmentCC" = "Cross-cutting", 
                               "EnvironmentMarine" = "Marine")



# Create reference level DF
Ref_level <- data.frame(Goal = c(rep("Reference", 5)), Variable = c("Female", "Non-university", "Europe & N. America", "Terrestrial", "Academia") , Estimate = c(rep(0, 5)), l.95..CI= c(rep(0, 5)),  u.95..CI= c(rep(0, 5)))

# Combine 
coef_df_b.s <- rbind(coef_df_b.s,Ref_level )

# Re-level 
coef_df_b.s$Variable <- factor(coef_df_b.s$Variable, levels = (c("Dispositional optimism","Years in conservation", "Work hours (binary)", "Female", "Male" ,
                                                                 "Non-university", "University", "Academia", "Practice", "Terrestrial", "Cross-cutting" ,
                                                                 "Marine", "Europe & N. America", "C. & S. Asia", "E. & S.E. Asia", "Sub-Saharan Africa",
                                                                 "Lat. America & Carib.", "N. Africa & W. Asia", "Oceania" )))

# Plot
Figure_WHB <- ggplot(coef_df_b.s, aes(x= Variable, y= Estimate)) + ylab("Goal progress satisfaction (log-odds)")  + 
  geom_pointrange(aes(color = Goal, shape = Goal, ymax = u.95..CI, ymin= l.95..CI),
                  position=position_dodge(width=c(0.6,0.6)), fatten = 2) +
  geom_hline(yintercept = 0, linetype="dashed", colour="grey73")  +  coord_flip()    + scale_color_manual(values = c(mypal[1], mypal[2], "black" ) ) +  scale_shape_manual(values = c(16, 16,  1)  ) + 
  scale_x_discrete(breaks = rev(c(levels(coef_df_b.s$Variable))), 
                   limits = rev(c(levels(coef_df_b.s$Variable)[1:3], "skip", 
                                  levels(coef_df_b.s$Variable)[4:5], "skip", 
                                  levels(coef_df_b.s$Variable)[6:7], "skip",
                                  levels(coef_df_b.s$Variable)[8:9], "skip", 
                                  levels(coef_df_b.s$Variable)[10:12], "skip",
                                  levels(coef_df_b.s$Variable)[13:19]))) + 
  theme(plot.margin=unit(c(0,0,1.5,0),"cm"))+
  theme(legend.position=c(.72,-.15), legend.direction="horizontal")

# Save 
ggsave("Figure_WHB.jpeg", width = 120, height = 140, units = "mm",   dpi = 800)
put_object(file = "Figure_WHB.jpeg", object = "Figure_WHB.jpeg", bucket = "cor-analysis", multipart=TRUE)
unlink("Figure_WHB.jpeg")
try(remove(list=c("individual_model.b.WHB")))
try(unlink("individual_model.b.WHB.RData"))
try(remove(list=c("collective_model.b.WHB")))
try(unlink("collective_model.b.WHB.RData"))


######### 10) Repeat the analysis disaggregating the 'collective' goals into those relating to altruistic and biospheric values ######### 

# Load the previous model
s3load("collective_model_full.b.summary.RData", bucket = "XXXX")

# Lists of the stacked DF (stacked by complete case of each goal) for each MI dataset - biospheric goals
mice.imp_biospheric.b <- list_along(1:length(mice.imp.GP.sim))

for (i in seq_along(1:10)){
  # Personal 
  mice.imp_biospheric.b[[i]] <- rbind(complete.goals(mice.imp.GP.sim[[i]], "GP_5_b" ),
                                      complete.goals(mice.imp.GP.sim[[i]], "GP_10_b" ), 
                                      complete.goals(mice.imp.GP.sim[[i]], "GP_9_b" ))
  mice.imp_biospheric.b[[i]]$goal_prog <- as.ordered(mice.imp_biospheric.b[[i]]$goal_prog)
  mice.imp_biospheric.b[[i]]$goal_prog <-  plyr::mapvalues(mice.imp_biospheric.b[[i]]$goal_prog, from = c("1", "2","3", "4", "5"), to = c("Very dissatisfied" , "Dissatisfied" , "Neutral" , "Satisfied" , "Very satisfied"))
  mice.imp_biospheric.b[[i]]$SO_Region <- relevel(mice.imp_biospheric.b[[i]]$SO_Region, ref = "Europe and Northern America") 
  mice.imp_biospheric.b[[i]]$DO <- as.numeric(scale(mice.imp_biospheric.b[[i]]$DO, scale = T, center = T))
}

# Lists of the stacked DF (stacked by complete case of each goal) for each MI dataset - altruistic goals
mice.imp_altruistic.b <- list_along(1:length(mice.imp.GP.sim))

for (i in seq_along(1:10)){
  # Personal 
  mice.imp_altruistic.b[[i]] <- rbind(complete.goals(mice.imp.GP.sim[[i]], "GP_6_b" ),
                                      complete.goals(mice.imp.GP.sim[[i]], "GP_7_b" ),
                                      complete.goals(mice.imp.GP.sim[[i]], "GP_8_b" ))
  mice.imp_altruistic.b[[i]]$goal_prog <- as.ordered(mice.imp_altruistic.b[[i]]$goal_prog)
  mice.imp_altruistic.b[[i]]$goal_prog <-  plyr::mapvalues(mice.imp_altruistic.b[[i]]$goal_prog, from = c("1", "2","3", "4", "5"), to = c("Very dissatisfied" , "Dissatisfied" , "Neutral" , "Satisfied" , "Very satisfied"))
  mice.imp_altruistic.b[[i]]$SO_Region <- relevel(mice.imp_altruistic.b[[i]]$SO_Region, ref = "Europe and Northern America") 
  mice.imp_altruistic.b[[i]]$DO <- as.numeric(scale(mice.imp_altruistic.b[[i]]$DO, scale = T, center = T))
}

# The model formula
brm_form.b <- formula("goal_prog ~ goal_num + years_cons_scaled + WH_scaled + DO  + gender + education_simple + SO_Region + position_simple + Environment + (1|ID)")

### Logistic regression, using the first of the ten imputed datasets 
# Biospheric goals 
start_time_biospheric.b <- Sys.time()
biospheric_model.b <- brm(brm_form.b, family = cumulative(link = "logit", threshold = "flexible"), data = mice.imp_biospheric.b[[1]], warmup = burn_in,  iter  = iterations,  chains = 4, seed = 123,  prior = prior_all, control = list(adapt_delta = 0.95), refresh = 0)
save(biospheric_model.b, file ="biospheric_model.b.RData")
put_object(file = "biospheric_model.b.RData", object = "biospheric_model.b.RData", bucket = "XXXX", multipart=TRUE)
end_time_biospheric.b <- Sys.time()
end_time_biospheric.b - start_time_biospheric.b

# Altruistic goals   altruistic 
start_time_altruistic.b <- Sys.time()
altruistic_model.b <- brm(brm_form.b, family = cumulative(link = "logit", threshold = "flexible"), data = mice.imp_altruistic.b[[1]], warmup = burn_in,  iter  = iterations,  chains = 4, seed = 123,  prior = prior_all, control = list(adapt_delta = 0.95), refresh = 0)
save(altruistic_model.b, file ="altruistic_model.b.RData")
put_object(file = "altruistic_model.b.RData", object = "altruistic_model.b.RData", bucket = "XXXX", multipart=TRUE)
end_time_altruistic.b <- Sys.time()
end_time_altruistic.b - start_time_altruistic.b

### Figures ### 
# Extract 
biospheric_model_full.b.summary <- summary(biospheric_model.b)
altruistic_model_full.b.summary <- summary(altruistic_model.b)

# Extracting coefficients
coef_list_b <- list(
  # Personal 
  data.frame(Goal = "Collective", Variable=as.factor(rownames(collective_model_full.b.summary$fixed)), (collective_model_full.b.summary$fixed[,c(1,3,4)])),
  
  # Personal 
  data.frame(Goal = "Biospheric", Variable=as.factor(rownames(biospheric_model_full.b.summary$fixed)), (biospheric_model_full.b.summary$fixed[,c(1,3,4)])),
  
  # People 
  data.frame(Goal = "Altruistic", Variable= as.factor(rownames(altruistic_model_full.b.summary$fixed)),  (altruistic_model_full.b.summary$fixed[,c(1,3,4)]))
)

# Create single DF
coef_df_b <- do.call(rbind, coef_list_b)

# Variables to retain
keep <- c("years_cons_scaled","DO","WH_scaled","genderMale", "EnvironmentCC", "EnvironmentMarine" , "education_simpleUniversity",
          "SO_RegionCentralandSouthernAsia" ,
          "SO_RegionEasternandSouthMEasternAsia" ,
          "SO_RegionSubMSaharanAfrica" ,
          "SO_RegionLatinAmericaandtheCaribbean" ,
          "SO_RegionNorthernAfricaandWesternAsia" ,
          "SO_RegionOceania",
          "position_simplePractice" )

# Keep variables and change to characters
coef_df_b.s <- coef_df_b[ which(coef_df_b$Variable %in% keep), ]
coef_df_b.s$Variable <- as.character(coef_df_b.s$Variable)

# Recode variable 
coef_df_b.s$Variable <- recode(coef_df_b.s$Variable, "years_cons_scaled" = "Years in conservation",
                               "DO" = "Dispositional optimism",
                               "WH_scaled" = "Work hours",
                               "genderMale" = "Male",
                               "education_simpleUniversity" = "University",
                               "SO_RegionCentralandSouthernAsia" = "C. & S. Asia",
                               "SO_RegionEasternandSouthMEasternAsia" = "E. & S.E. Asia",
                               "SO_RegionSubMSaharanAfrica" = "Sub-Saharan Africa",
                               "SO_RegionLatinAmericaandtheCaribbean" = "Lat. America & Carib.",
                               "SO_RegionNorthernAfricaandWesternAsia" = "N. Africa & W. Asia",
                               "SO_RegionOceania"= "Oceania",
                               "position_simplePractice" = "Practice",
                               "EnvironmentCC" = "Cross-cutting", 
                               "EnvironmentMarine" = "Marine")



# Create reference level DF
Ref_level <- data.frame(Goal = c(rep("Reference", 5)), Variable = c("Female", "Non-university", "Europe & N. America", "Terrestrial", "Academia") , Estimate = c(rep(0, 5)), l.95..CI= c(rep(0, 5)),  u.95..CI= c(rep(0, 5)))

# Combine 
coef_df_b.s <- rbind(coef_df_b.s,Ref_level )

# Re-level 
coef_df_b.s$Variable <- factor(coef_df_b.s$Variable, levels = (c("Dispositional optimism","Years in conservation", "Work hours", "Female", "Male" ,
                                                                 "Non-university", "University", "Academia", "Practice", "Terrestrial", "Cross-cutting" ,
                                                                 "Marine", "Europe & N. America", "C. & S. Asia", "E. & S.E. Asia", "Sub-Saharan Africa",
                                                                 "Lat. America & Carib.", "N. Africa & W. Asia", "Oceania" )))

# Plot
Figure_goal.t_b <- ggplot(coef_df_b.s, aes(x= Variable, y= Estimate)) + ylab("Goal endorsement")  + 
  geom_pointrange(aes(color = Goal, shape = Goal, ymax = u.95..CI, ymin= l.95..CI),
                  position=position_dodge(width=c(0.6,0.6)), fatten = 2) +
  geom_hline(yintercept = 0, lnetype="dashed", colour="grey73")  +  coord_flip()    + scale_color_manual(values = c('#74c476', '#006d2c'  , mypal[1],  "black" ) ) +  scale_shape_manual(values = c(16, 16, 16,  1)  ) + 
  scale_x_discrete(breaks = rev(c(levels(coef_df_b.s$Variable))), 
                   limits = rev(c(levels(coef_df_b.s$Variable)[1:3], "skip", 
                                  levels(coef_df_b.s$Variable)[4:5], "skip", 
                                  levels(coef_df_b.s$Variable)[6:7], "skip",
                                  levels(coef_df_b.s$Variable)[8:9], "skip", 
                                  levels(coef_df_b.s$Variable)[10:12], "skip",
                                  levels(coef_df_b.s$Variable)[13:19]))) + 
  theme(plot.margin=unit(c(0,0,1.5,0),"cm"))+
  theme(legend.position=c(.55,-.15), legend.direction="horizontal")

# Save 
ggsave("Figure_goal.t_b.jpeg", width = 120, height = 140, units = "mm",   dpi = 800)
put_object(file = "Figure_goal.t_b.jpeg", object = "Figure_goal.t_b.jpeg", bucket = "XXXX", multipart=TRUE)
unlink("Figure_goal.t_b.jpeg")
try(remove(list=c("biospheric_model.b")))
try(unlink("biospheric_model.b.RData"))
try(remove(list=c("altruistic_model.b")))
try(unlink("altruistic_model.b.RData"))


######### 11) Characteristics associated with goal endorsement #########  

# A function that returns complete cases for each goal (they should and are all complete cases) and the name of that goal
complete.goals.a <- function(DF, goal){
  DF_goal <- DF[complete.cases(DF[,goal]),]
  DF_goal$goal_end <- DF_goal[,goal]
  DF_goal$goal_num <- as.factor(goal) 
  return(DF_goal)
}

# Lists of the stacked DF (stacked by complete case of each goal) for each MI dataset
mice.imp_individual.a <- list_along(1:length(mice.imp.GP.sim))
mice.imp_collective.a <- list_along(1:length(mice.imp.GP.sim))

for (i in seq_along(1:10)){
  
  # Individual  
  mice.imp_individual.a[[i]] <- rbind(complete.goals.a(mice.imp.GP.sim[[i]], "GP_1_a" ),
                                      complete.goals.a(mice.imp.GP.sim[[i]], "GP_2_a" ), 
                                      complete.goals.a(mice.imp.GP.sim[[i]], "GP_3_a" ), 
                                      complete.goals.a(mice.imp.GP.sim[[i]], "GP_4_a" ))
  
  # Collective
  mice.imp_collective.a[[i]] <- rbind(complete.goals.a(mice.imp.GP.sim[[i]], "GP_5_a" ),
                                      complete.goals.a(mice.imp.GP.sim[[i]], "GP_6_a" ),
                                      complete.goals.a(mice.imp.GP.sim[[i]], "GP_7_a" ),
                                      complete.goals.a(mice.imp.GP.sim[[i]], "GP_8_a" ),
                                      complete.goals.a(mice.imp.GP.sim[[i]], "GP_9_a" ),
                                      complete.goals.a(mice.imp.GP.sim[[i]], "GP_10_a"))
}

# The model formula
brm_form.a <- formula("goal_end ~ goal_num + years_cons_scaled + WH_scaled + DO  + gender + education_simple + SO_Region + position_simple + Environment + (1|ID)")

### Logistic regression, using the first of the ten imputed datasets 
# Individual goals 
start_time_individual.a <- Sys.time()
individual_model.a <- brm(brm_form.a, family = 'bernoulli', data = mice.imp_individual.a[[1]], warmup = burn_in,  iter  = iterations,  chains = 4, seed = 123,  prior = prior_all)
save(individual_model.a, file ="individual_model.a.RData")
put_object(file = "individual_model.a.RData", object = "individual_model.a.RData", bucket = "XXXX", multipart=TRUE)
end_time_individual.a <- Sys.time()
end_time_individual.a - start_time_individual.a

# Collective goals  
start_time_collective.a <- Sys.time()
collective_model.a <- brm(brm_form.a, family = 'bernoulli', data = mice.imp_collective.a[[1]], warmup = burn_in,  iter  = iterations,  chains = 4, seed = 123,  prior = prior_all)
save(collective_model.a, file ="collective_model.a.RData")
put_object(file = "collective_model.a.RData", object = "collective_model.a.RData", bucket = "XXXX", multipart=TRUE)
end_time_collective.a <- Sys.time()
end_time_collective.a - start_time_collective.a

### Figures ### 
# Extract 
individual_model_full.a.summary <- summary(individual_model.a)
collective_model_full.a.summary <- summary(collective_model.a)

# Extracting coefficients
coef_list_a <- list(
  
  # Personal 
  data.frame(Goal = "Individual", Variable=as.factor(rownames(individual_model_full.a.summary$fixed)), (individual_model_full.a.summary$fixed[,c(1,3,4)])),
  
  # People 
  data.frame(Goal = "Collective", Variable= as.factor(rownames(collective_model_full.a.summary$fixed)),  (collective_model_full.a.summary$fixed[,c(1,3,4)]))
)

# Create single DF
coef_df_a <- do.call(rbind, coef_list_a)

# Variables to retain
keep <- c("years_cons_scaled","DO","WH_scaled","genderMale", "EnvironmentCC", "EnvironmentMarine" , "education_simpleUniversity",
          "SO_RegionCentralandSouthernAsia" ,
          "SO_RegionEasternandSouthMEasternAsia" ,
          "SO_RegionSubMSaharanAfrica" ,
          "SO_RegionLatinAmericaandtheCaribbean" ,
          "SO_RegionNorthernAfricaandWesternAsia" ,
          "SO_RegionOceania",
          "position_simplePractice" )

# Keep variables and change to characters
coef_df_a.s <- coef_df_a[ which(coef_df_a$Variable %in% keep), ]
coef_df_a.s$Variable <- as.character(coef_df_a.s$Variable)

# Recode variable 
coef_df_a.s$Variable <- recode(coef_df_a.s$Variable, "years_cons_scaled" = "Years in conservation",
                               "DO" = "Dispositional optimism",
                               "WH_scaled" = "Work hours",
                               "genderMale" = "Male",
                               "education_simpleUniversity" = "University",
                               "SO_RegionCentralandSouthernAsia" = "C. & S. Asia",
                               "SO_RegionEasternandSouthMEasternAsia" = "E. & S.E. Asia",
                               "SO_RegionSubMSaharanAfrica" = "Sub-Saharan Africa",
                               "SO_RegionLatinAmericaandtheCaribbean" = "Lat. America & Carib.",
                               "SO_RegionNorthernAfricaandWesternAsia" = "N. Africa & W. Asia",
                               "SO_RegionOceania"= "Oceania",
                               "position_simplePractice" = "Practice",
                               "EnvironmentCC" = "Cross-cutting", 
                               "EnvironmentMarine" = "Marine")



# Create reference level DF
Ref_level <- data.frame(Goal = c(rep("Reference", 5)), Variable = c("Female", "Non-university", "Europe & N. America", "Terrestrial", "Academia") , Estimate = c(rep(0, 5)), l.95..CI= c(rep(0, 5)),  u.95..CI= c(rep(0, 5)))

# Combine 
coef_df_a.s <- rbind(coef_df_a.s,Ref_level )

# Re-level 
coef_df_a.s$Variable <- factor(coef_df_a.s$Variable, levels = (c("Dispositional optimism","Years in conservation", "Work hours", "Female", "Male" ,
                                                                 "Non-university", "University", "Academia", "Practice", "Terrestrial", "Cross-cutting" ,
                                                                 "Marine", "Europe & N. America", "C. & S. Asia", "E. & S.E. Asia", "Sub-Saharan Africa",
                                                                 "Lat. America & Carib.", "N. Africa & W. Asia", "Oceania" )))

# Plot
Figure_goal_a <- ggplot(coef_df_a.s, aes(x= Variable, y= Estimate)) + ylab("Goal endorsement")  + 
  geom_pointrange(aes(color = Goal, shape = Goal, ymax = u.95..CI, ymin= l.95..CI),
                  position=position_dodge(width=c(0.6,0.6)), fatten = 2) +
  geom_hline(yintercept = 0, linetype="dashed", colour="grey73")  +  coord_flip()    + scale_color_manual(values = c(mypal[1], mypal[2], "black" ) ) +  scale_shape_manual(values = c(16, 16,  1)  ) + 
  scale_x_discrete(breaks = rev(c(levels(coef_df_a.s$Variable))), 
                   limits = rev(c(levels(coef_df_a.s$Variable)[1:3], "skip", 
                                  levels(coef_df_a.s$Variable)[4:5], "skip", 
                                  levels(coef_df_a.s$Variable)[6:7], "skip",
                                  levels(coef_df_a.s$Variable)[8:9], "skip", 
                                  levels(coef_df_a.s$Variable)[10:12], "skip",
                                  levels(coef_df_a.s$Variable)[13:19]))) + 
  theme(plot.margin=unit(c(0,0,1.5,0),"cm"))+
  theme(legend.position=c(.72,-.15), legend.direction="horizontal")

# Save 
ggsave("Figure_goal_a.jpeg", width = 120, height = 140, units = "mm",   dpi = 800)
put_object(file = "Figure_goal_a.jpeg", object = "Figure_goal_a.jpeg", bucket = "XXXX", multipart=TRUE)
unlink("Figure_goal_a.jpeg")
try(remove(list=c("individual_model.a")))
try(unlink("individual_model.a.RData"))
try(remove(list=c("collective_model.a")))
try(unlink("collective_model.a.RData"))

######### 12) Association between satisfaction and goal type #########

# Lists of the stacked DF (stacked by complete case of each goal) for each MI dataset - for both the individual and collective goals
mice.imp_all.b <- list_along(1:length(mice.imp.GP.sim))

for (i in seq_along(1:10)){
  # Personal 
  mice.imp_all.b[[i]] <- rbind(complete.goals(mice.imp.GP.sim[[i]], "GP_1_b" ),
                               complete.goals(mice.imp.GP.sim[[i]], "GP_2_b" ), 
                               complete.goals(mice.imp.GP.sim[[i]], "GP_3_b" ), 
                               complete.goals(mice.imp.GP.sim[[i]], "GP_4_b" ),
                               complete.goals(mice.imp.GP.sim[[i]], "GP_5_b" ),
                               complete.goals(mice.imp.GP.sim[[i]], "GP_6_b" ),
                               complete.goals(mice.imp.GP.sim[[i]], "GP_7_b" ),
                               complete.goals(mice.imp.GP.sim[[i]], "GP_8_b" ),
                               complete.goals(mice.imp.GP.sim[[i]], "GP_9_b" ),
                               complete.goals(mice.imp.GP.sim[[i]], "GP_10_b"))
  
  mice.imp_all.b[[i]]$goal_prog <- as.ordered(mice.imp_all.b[[i]]$goal_prog)
  mice.imp_all.b[[i]]$goal_prog <-  plyr::mapvalues(mice.imp_all.b[[i]]$goal_prog, from = c("1", "2","3", "4", "5"), to = c("Very dissatisfied" , "Dissatisfied" , "Neutral" , "Satisfied" , "Very satisfied"))
  mice.imp_all.b[[i]]$SO_Region <- relevel(mice.imp_all.b[[i]]$SO_Region, ref = "Europe and Northern America") 
  mice.imp_all.b[[i]]$DO <- as.numeric(scale(mice.imp_all.b[[i]]$DO, scale = T, center = T))
  
  # Add goal type
  ind <- c("GP_1_b","GP_2_b","GP_3_b","GP_4_b")
  coll <- c("GP_5_b", "GP_6_b", "GP_7_b", "GP_8_b", "GP_9_b", "GP_10_b")
  mice.imp_all.b[[i]]$Goal_type <- ifelse( mice.imp_all.b[[i]]$goal_num %in% ind, "Individual", 
                                           ifelse(mice.imp_all.b[[i]]$goal_num %in% coll, "Collective", "Error")) 
}


# The model formula
brm_form.all <- formula("goal_prog ~ Goal_type + (1|ID)")

# Single  - all
start_time_all<- Sys.time()
all_model.b <- brm(brm_form.all, family = cumulative(link = "logit", threshold = "flexible"), data = mice.imp_all.b[[1]], warmup = burn_in,  iter  = iterations,  chains = 4, seed = 123, prior = prior_all, control = list(adapt_delta = 0.95), refresh = 0) 
save(all_model.b, file ="all_model.b.RData")
put_object(file = "all_model.b.RData", object = "all_model.b.RData", bucket = "XXXX", multipart=TRUE)
end_time_all<- Sys.time()
end_time_all- start_time_all

### Extract summary ### 
# Extract the summary information and save
all_model.b.summary <- summary(all_model.b)
save(all_model.b.summary, file ="all_model.b.summary.RData")
put_object(file = "all_model.b.summary.RData", object = "all_model.b.summary.RData", bucket = "XXXX", multipart=TRUE)


# Remove files
try(remove(list=c("all_model.b")))
try(unlink("all_model.b.RData"))
try(remove(list=c("all_model.b.summary")))
try(unlink("all_model.b.summary.RData"))
