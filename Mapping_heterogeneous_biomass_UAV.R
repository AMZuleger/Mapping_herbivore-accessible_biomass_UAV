######################################################################################################
################### Mapping Diverse Biomass Types Across a Heterogeneous Landscape ###################
######################## Using Multi-sensor High-Resolution UAV Data  ################################

################## Annika M. Zuleger, Martina M. Viti, Luise Quoss, Filipe S. Dias, ##################
##################### Luís Borda-de-Água, Miguel N. Bugalho & Henrique M. Pereira ####################

######################################## February, 2025 ##############################################

######## Load required packages ##########

library(dplyr)
library(psych)
library(mgcv)
library(visreg)
library(ggplot2)
library(gridExtra)

#### Load required data from GitHub ####

t <- tempdir()
setwd(t) # set temporal wd

# download data from github
url <- "https://github.com/AMZuleger/Mapping_heterogeneous_biomass_UAV/archive/refs/heads/main.zip"
download.file(url,destfile="Mapping_heterogeneous_biomass_UAV-main.zip")
unzip(zipfile="Mapping_heterogeneous_biomass_UAV-main.zip")

setwd(dir = file.path(t,"Mapping_heterogeneous_biomass_UAV-main"))
unzip(zipfile="Data.zip")

setwd(dir = file.path(t,"Mapping_heterogeneous_biomass_UAV-main/Data"))

biomass_data <- read.csv("biomass_data_2024.csv",header=T,sep=",")

biomass_data$Plot <- as.factor(biomass_data$Plot)
biomass_data$Transect <- as.factor(biomass_data$Transect)
biomass_data$Patch <- as.factor(biomass_data$Patch)
biomass_data$LULC_RF <- as.factor(biomass_data$LULC_RF)
biomass_data$Herb_biomass <- biomass_data$Forb_biomass + biomass_data$Grass_biomass

##### Summary statistics #### 

mean(biomass_data$Total_biomass,na.rm=T)
sd(biomass_data$Total_biomass, na.rm = TRUE) / sqrt(nrow(biomass_data))
biomass_data %>%
  group_by(LULC_RF) %>%
  summarise(mean_biomass = format(mean(Total_biomass, na.rm = TRUE), nsmall = 2),
            SE_biomass = format(sd(Total_biomass, na.rm = TRUE) / sqrt(n()), nsmall = 2),
            SD_biomass = format(sd(Total_biomass, na.rm = TRUE), nsmall = 2))

mean(biomass_data$Shrub_biomass,na.rm=T)
sd(biomass_data$Shrub_biomass, na.rm = TRUE) / sqrt(nrow(biomass_data))
biomass_data %>%
  group_by(LULC_RF) %>%
  summarise(mean_biomass = format(mean(Shrub_biomass, na.rm = TRUE), nsmall = 2),
            SE_biomass = format(sd(Shrub_biomass, na.rm = TRUE) / sqrt(n()), nsmall = 2),
            SD_biomass = format(sd(Shrub_biomass, na.rm = TRUE), nsmall = 2))

mean(biomass_data$Forb_biomass,na.rm=T)
sd(biomass_data$Forb_biomass, na.rm = TRUE) / sqrt(nrow(biomass_data))
biomass_data %>%
  group_by(LULC_RF) %>%
  summarise(mean_biomass = format(mean(Forb_biomass, na.rm = TRUE), nsmall = 2),
            SE_biomass = format(sd(Forb_biomass, na.rm = TRUE) / sqrt(n()), nsmall = 2),
            SD_biomass = format(sd(Forb_biomass, na.rm = TRUE), nsmall = 2))

mean(biomass_data$Grass_biomass,na.rm=T)
sd(biomass_data$Grass_biomass, na.rm = TRUE) / sqrt(nrow(biomass_data))
biomass_data %>%
  group_by(LULC_RF) %>%
  summarise(mean_biomass = format(mean(Grass_biomass, na.rm = TRUE), nsmall = 2),
            SE_biomass = format(sd(Grass_biomass, na.rm = TRUE) / sqrt(n()), nsmall = 2),
            SD_biomass = format(sd(Grass_biomass, na.rm = TRUE), nsmall = 2))

mean(biomass_data$Herb_biomass,na.rm=T)
sd(biomass_data$Herb_biomass, na.rm = TRUE) / sqrt(nrow(biomass_data))
biomass_data %>%
  group_by(LULC_RF) %>%
  summarise(mean_biomass = format(mean(Herb_biomass, na.rm = TRUE), nsmall = 2),
            SE_biomass = format(sd(Herb_biomass, na.rm = TRUE) / sqrt(n()), nsmall = 2),
            SD_biomass = format(sd(Herb_biomass, na.rm = TRUE), nsmall = 2))

#### Variable Correlation ####

corPlot(biomass_data[, c(4:7,29, 10:12,14:28)],method="spearman",use = "complete.obs", main="Correlation Plot All Covariates",xlas=3,n=20,cex=0.4,scale=F,upper=F,cex.axis=0.8)
corloads = cor(biomass_data[, c(4:7,29, 10:12,14:28)], use = "pairwise.complete.obs")
dissimilarity = 1 - corloads
distance = as.dist(dissimilarity) 
plot(hclust(distance), main="Dissimilarity = 1 - Correlation", xlab="") 

####### Fit biomass models to field data ######

##### Total Biomass #####

# Fit model
model_total <- gam(Total_biomass ~ s(Elevation_mean, k=5) + s(Slope_mean, k=5) +
                    s(Aspect_mean, bs = "cc", k = 5) +
                    s(NDVI_mean, by = LULC_RF, k = 5) + s(VHeight_min, by = LULC_RF, k = 5) + 
                    s(VHeight_max, by = LULC_RF, k = 5) + s(CC_min, k = 5) + s(vrd_0_mean, k = 5) + 
                    s(vrd_0.5_mean, k = 5) + s(vrd_1_mean, k = 5) + s(vrd_2_mean, k = 5) + s(vrd_2.5_mean, k = 5) +
                    s(Transect, bs = "re"),
                  data = biomass_data,
                  method = 'REML',
                  family = tw(),select=TRUE)

# Check model 
summary(model_total)
par(mfrow=c(2,2))
gam.check(model_total)

# Visualize variable effects
par(mfrow=c(3,2))
visreg(model_total, "NDVI_mean", by = "LULC_RF", type = "conditional",main="NDVI_mean")
visreg(model_total, "VHeight_min", by = "LULC_RF", type = "conditional",main="Vegetation Height (min)")
visreg(model_total, "VHeight_max", by = "LULC_RF", type = "conditional",main="Vegetation height (max)")
visreg(model_total, "Elevation_mean", type = "conditional",main="Elevation")
visreg(model_total, "Slope_mean", type = "conditional",main="Slope")
visreg(model_total, "vrd_1_mean", type = "conditional",main="VRD_1")
visreg(model_total, "vrd_2_mean", type = "conditional",main="VRD_2")
visreg(model_total, "vrd_2.5_mean", type = "conditional",main="VRD_2.5")


# Fit model including Land-use as a variable
model_total_lulc <- gam(Total_biomass ~ LULC_RF + s(Elevation_mean, k=5) + s(Slope_mean, k=5) +
                         s(Aspect_mean, bs = "cc", k = 5) +
                         s(NDVI_mean, by = LULC_RF, k = 5) + s(VHeight_min, by = LULC_RF, k = 5) + 
                         s(VHeight_max, by = LULC_RF, k = 5) + s(CC_min, k = 5) + s(vrd_0_mean, k = 5) + 
                         s(vrd_0.5_mean, k = 5) + s(vrd_1_mean, k = 5) + s(vrd_2_mean, k = 5) + s(vrd_2.5_mean, k = 5) +
                         s(Transect, bs = "re"),
                  data = biomass_data,
                  method = 'REML',
                  family = tw(),select=TRUE)

# Check model
summary(model_total_lulc)
par(mfrow=c(2,2))
gam.check(model_total_lulc)

# Compare models
AIC(model_total,model_total_lulc)
anova(model_total,model_total_lulc,test="Chisq")

# Visualize variable effects
par(mfrow=c(3,2))
visreg(model_total, "NDVI_mean", by = "LULC_RF", type = "conditional",main="NDVI_mean")
visreg(model_total, "VHeight_min", by = "LULC_RF", type = "conditional",main="Vegetation Height (min)")
visreg(model_total, "VHeight_max", by = "LULC_RF", type = "conditional",main="Vegetation height (max)")
visreg(model_total, "Elevation_mean", type = "conditional",main="Elevation")
visreg(model_total, "Slope_mean", type = "conditional",main="Slope")
visreg(model_total, "vrd_1_mean", type = "conditional",main="VRD_1")
visreg(model_total, "vrd_2_mean", type = "conditional",main="VRD_2")
visreg(model_total, "vrd_2.5_mean", type = "conditional",main="VRD_2.5")

### Cross-validation ####
cv_mod <- data.frame(stat=c("AIC","R²","Dev.Exp","REML","n","ß","mdf","rdf","MAE_80","MAE_20","RSME_80","RSME_20","Bias_80","Bias_20","R²_80","R²_20"))
n_iterations <- 10
lulc_classes <- unique(biomass_data$LULC_RF)

# Vector to store mean absolute errors (MAE)
mae_80_list <- numeric(n_iterations)
mae_20_list <- numeric(n_iterations)
rsme_80_list <- numeric(n_iterations)
rsme_20_list <- numeric(n_iterations)
bias_80_list <- numeric(n_iterations)
bias_20_list <- numeric(n_iterations)
rsq_80_list <- numeric(n_iterations)
rsq_20_list <- numeric(n_iterations)

combined_results_total <- data.frame(
  Iteration = integer(),
  Dataset = character(),
  Transect = integer(),
  Patch = factor(),
  Total_biomass = numeric(),
  Pred_biomass = numeric(),
  stringsAsFactors = FALSE
)

# Loop over the number of iterations
for (i in 1:n_iterations) {
  train_transects <- list()
  test_transects <- list()
  # Loop over each LULC_RF class to perform stratified sampling
  for (lulc in lulc_classes) {
    # Subset data for the current LULC_RF class
    class_data <- biomass_data[biomass_data$LULC_RF == lulc, ]
    # Get the unique transects for the current class
    unique_transects <- unique(class_data$Transect)
    # Perform stratified sampling: select 80% of the transects for training
    selected_train_transects <- sample(unique_transects, size = round(0.8 * length(unique_transects)), replace = FALSE)
    # Store the selected transects for training and testing
    train_transects[[lulc]] <- selected_train_transects
    test_transects[[lulc]] <- setdiff(unique_transects, selected_train_transects)
  }
  # Combine the selected transects across all LULC_RF classes for training and testing
  train_transects_combined <- unlist(train_transects)
  test_transects_combined <- unlist(test_transects)
  # Split the data into 80% training and 20% testing based on the selected transects
  data_80 <- biomass_data[biomass_data$Transect %in% train_transects_combined, ]
  data_20 <- biomass_data[biomass_data$Transect %in% test_transects_combined, ]
  # Ensure the levels of Patch are consistent in both training and testing data
  levels(data_80$Patch) <- levels(biomass_data$Patch)
  levels(data_20$Patch) <- levels(biomass_data$Patch)
  #   # Fit the model on the 80% training data
  model_rf_80 <- gam(model_total$formula,
                     data = data_80,
                     method = 'REML',
                     family = tw(),select=T)
  # Predict on the 80% training data
  data_80$Pred_biomass <- predict(model_rf_80, newdata = data_80, type = "response", exclude = "s(Transect)")
  # Calculate the prediction error (MAE and RSME)
  mae_80 <- mean(abs(data_80$Total_biomass - data_80$Pred_biomass), na.rm = TRUE)
  mae_80_list[i] <- mae_80
  rsme_80 <- sqrt(mean((data_80$Total_biomass - data_80$Pred_biomass)^2,na.rm = T))
  rsme_80_list[i] <- rsme_80
  # Calcualte bias
  bias_80 <- mean(data_80$Pred_biomass - data_80$Total_biomass,na.rm=T)
  bias_80_list[i] <- bias_80
  # Predict on the 20% training data
  data_20$Pred_biomass <- predict(model_rf_80, newdata = data_20, type = "response", exclude = "s(Transect)")
  # Calculate the prediction error (MAE and RSME)
  mae_20 <- mean(abs(data_20$Total_biomass - data_20$Pred_biomass), na.rm = TRUE)
  mae_20_list[i] <- mae_20
  rsme_20 <- sqrt(mean((data_20$Total_biomass - data_20$Pred_biomass)^2,na.rm = T))
  rsme_20_list[i] <- rsme_20
  # Calcualte bias
  bias_20 <- mean(data_20$Pred_biomass - data_20$Total_biomass,na.rm=T)
  bias_20_list[i] <- bias_20
  # # Plot fit
  # plot(data_20$Total_biomass~data_20$Pred_biomass,col="blue")
  # points(data_80$Total_biomass~data_80$Pred_biomass,col="red")
  # abline(lm(data_20$Total_biomass~data_20$Pred_biomass),col="blue")
  # abline(lm(data_80$Total_biomass~data_80$Pred_biomass),col="red")
  # abline(1,1,lty=2)
  rsq80 <- summary(lm(data_80$Total_biomass~data_80$Pred_biomass))$adj.r.squared
  rsq20 <- summary(lm(data_20$Total_biomass~data_20$Pred_biomass))$adj.r.squared
  # legend("bottomleft",legend = c(paste("R?(20):", round(rsq20, 2)),
  #                                paste("R?(80):", round(rsq80, 2))),
  #        col = c("blue", "red"),lty = 1)
  rsq_80_list[i] <- rsq80
  rsq_20_list[i] <- rsq20
  # Add a column for iteration and dataset type to data_80 and data_20
  data_80$Iteration <- i
  data_80$Dataset <- "Training"
  data_20$Iteration <- i
  data_20$Dataset <- "Testing"
  # Append the predicted and original values to combined_results
  combined_results_total <- rbind(
    combined_results_total,
    data_80[, c("Iteration", "Dataset", "Transect", "Patch", "Total_biomass", "Pred_biomass")],
    data_20[, c("Iteration", "Dataset", "Transect", "Patch", "Total_biomass", "Pred_biomass")]
  )
}


# Summary statistics of MAE and RSME
cv_mod$model_total[1] <- AIC(model_total)
cv_mod$model_total[2] <- summary(model_total)$r.sq
cv_mod$model_total[3] <- summary(model_total)$dev.expl
cv_mod$model_total[4] <- as.numeric(summary(model_total)$sp.criterion)
cv_mod$model_total[5] <- summary(model_total)$n
cv_mod$model_total[6] <- length(summary(model_total)$edf)
cv_mod$model_total[7] <- sum(influence(model_total))
cv_mod$model_total[8] <- df.residual(model_total)
cv_mod$model_total[9] <- mean(mae_80_list,na.rm=T)
cv_mod$model_total[10] <- mean(mae_20_list,na.rm=T)
cv_mod$model_total[11] <- mean(rsme_80_list,na.rm=T)
cv_mod$model_total[12] <- mean(rsme_20_list,na.rm=T)
cv_mod$model_total[13] <- mean(bias_80_list,na.rm=T)
cv_mod$model_total[14] <- mean(bias_20_list,na.rm=T)
cv_mod$model_total[15] <- mean(rsq_80_list, na.rm = T)
cv_mod$model_total[16] <- mean(rsq_20_list, na.rm = T)

print(cv_mod)

training <- subset(combined_results_total,Dataset == "Training")
testing <- subset(combined_results_total,Dataset == "Testing")

ggplot(combined_results_total, aes(x = Total_biomass, y = Pred_biomass, color = Dataset)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(data = training, aes(x = Total_biomass, y = Pred_biomass), method = "lm", color = "red", se = FALSE) +
  geom_smooth(data = testing, aes(x = Total_biomass, y = Pred_biomass), method = "lm", color = "blue", se = FALSE) +
  xlim(0,1000) +
  ylim(0,1000) +
  labs(title = "Predicted vs Actual Biomass",
       x = "Total Biomass",
       y = "Predicted Biomass") +
  scale_color_manual(
    values = c("Training" = "red", "Testing" = "blue"),  # Colors for the legend
    labels = c("Training Set", "Test Set")  # Custom legend text
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14), 
    axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    legend.position = c(0.12, 0.93),
    legend.title = element_blank(),  
    legend.text = element_text(size = 14)) +
  geom_text(
    data = data.frame(x = 840, y = 50, label = paste("Dev.Exp =", round(cv_mod$model_total[3],2),"\nRSME =", round(cv_mod$model_total[12],2),"\nBias =",round(cv_mod$model_total[14],2), "\nR²_20 =", round(cv_mod$model_total[16],2))),
    aes(x = x, y = y, label = label, hjust = 0),
    size = 5, color = "black", fontface = "italic")


##### Shrub Biomass #####

model_shrub <- gam(Shrub_biomass ~ s(Elevation_mean, k=5) + s(Slope_mean, k=5) +
                      s(Aspect_mean, bs = "cc", k = 5) +
                      s(NDVI_mean, by = LULC_RF, k = 5) + s(VHeight_min, by = LULC_RF, k = 5) + 
                      s(VHeight_max, by = LULC_RF, k = 5) + s(CC_min, k = 5) + s(vrd_0_mean, k = 5) + 
                      s(vrd_0.5_mean, k = 5) + s(vrd_1_mean, k = 5) + s(vrd_2_mean, k = 5) + s(vrd_2.5_mean, k = 5) +
                      s(Transect, bs = "re"),
                    data = biomass_data,
                    method = 'REML',
                    family = tw(),select=TRUE)

summary(model_shrub)
par(mfrow=c(2,2))
gam.check(model_shrub)

model_shrub_lulc <- gam(Shrub_biomass ~ LULC_RF + s(Elevation_mean, k=5) + s(Slope_mean, k=5) +
                           s(Aspect_mean, bs = "cc", k = 5) +
                           s(NDVI_mean, by = LULC_RF, k = 5) + s(VHeight_min, by = LULC_RF, k = 5) + 
                           s(VHeight_max, by = LULC_RF, k = 5) + s(CC_min, k = 5) + s(vrd_0_mean, k = 5) + 
                           s(vrd_0.5_mean, k = 5) + s(vrd_1_mean, k = 5) + s(vrd_2_mean, k = 5) + s(vrd_2.5_mean, k = 5) +
                           s(Transect, bs = "re"),
                    data = biomass_data,
                    method = 'REML',
                    family = tw(),select=TRUE)

summary(model_shrub_lulc)
par(mfrow=c(2,2))
gam.check(model_shrub_lulc)

AIC(model_shrub, model_shrub_lulc)
anova(model_shrub, model_shrub_lulc, test = "Chisq")

par(mfrow=c(3,2))
visreg(model_shrub_lulc, "NDVI_mean", by = "LULC_RF", type = "conditional",main="NDVI_mean")
visreg(model_shrub_lulc, "VHeight_min", by = "LULC_RF", type = "conditional",main="Vegetation Height (min)")
visreg(model_shrub_lulc, "VHeight_max", by = "LULC_RF", type = "conditional",main="Vegetation height (max)")
visreg(model_shrub_lulc, "LULC_RF", type = "conditional",main="Land-use")
visreg(model_shrub_lulc, "Elevation_mean", type = "conditional",main="Elevation")
visreg(model_shrub_lulc, "CC_min", type = "conditional",main="Canopy Cover (min)")
visreg(model_shrub_lulc, "vrd_1_mean", type = "conditional",main="VRD_1")
visreg(model_shrub_lulc, "vrd_2_mean", type = "conditional",main="VRD_2")
visreg(model_shrub_lulc, "vrd_2.5_mean", type = "conditional",main="VRD_2.5")

#### Cross-validation ####

combined_results_shrub <- data.frame(
  Iteration = integer(),
  Dataset = character(),
  Transect = integer(),
  Patch = factor(),
  Shrub_biomass = numeric(),
  Pred_biomass = numeric(),
  stringsAsFactors = FALSE
)

# Loop over the number of iterations
for (i in 1:n_iterations) {
  train_transects <- list()
  test_transects <- list()
  # Loop over each LULC_RF class to perform stratified sampling
  for (lulc in lulc_classes) {
    # Subset data for the current LULC_RF class
    class_data <- biomass_data[biomass_data$LULC_RF == lulc, ]
    # Get the unique transects for the current class
    unique_transects <- unique(class_data$Transect)
    # Perform stratified sampling: select 80% of the transects for training
    selected_train_transects <- sample(unique_transects, size = round(0.8 * length(unique_transects)), replace = FALSE)
    # Store the selected transects for training and testing
    train_transects[[lulc]] <- selected_train_transects
    test_transects[[lulc]] <- setdiff(unique_transects, selected_train_transects)
  }
  # Combine the selected transects across all LULC_RF classes for training and testing
  train_transects_combined <- unlist(train_transects)
  test_transects_combined <- unlist(test_transects)
  # Split the data into 80% training and 20% testing based on the selected transects
  data_80 <- biomass_data[biomass_data$Transect %in% train_transects_combined, ]
  data_20 <- biomass_data[biomass_data$Transect %in% test_transects_combined, ]
  # Ensure the levels of Patch are consistent in both training and testing data
  levels(data_80$Patch) <- levels(biomass_data$Patch)
  levels(data_20$Patch) <- levels(biomass_data$Patch)
  #   # Fit the model on the 80% training data
  model_rf_80 <- gam(model_shrub_lulc$formula,
                     data = data_80,
                     method = 'REML',
                     family = tw(),select=T)
  # Predict on the 80% training data
  data_80$Pred_biomass <- predict(model_rf_80, newdata = data_80, type = "response", exclude = "s(Transect)")
  # Calculate the prediction error (e.g., Mean Absolute Error)
  mae_80 <- mean(abs(data_80$Shrub_biomass - data_80$Pred_biomass), na.rm = TRUE)
  mae_80_list[i] <- mae_80
  rsme_80 <- sqrt(mean((data_80$Shrub_biomass - data_80$Pred_biomass)^2,na.rm = T))
  rsme_80_list[i] <- rsme_80
  # Calcualte bias
  bias_80 <- mean(data_80$Pred_biomass - data_80$Shrub_biomass,na.rm=T)
  bias_80_list[i] <- bias_80
  # Predict on the 20% training data
  data_20$Pred_biomass <- predict(model_rf_80, newdata = data_20, type = "response", exclude = "s(Transect)")
  # Calculate the prediction error (e.g., Mean Absolute Error)
  mae_20 <- mean(abs(data_20$Shrub_biomass - data_20$Pred_biomass), na.rm = TRUE)
  mae_20_list[i] <- mae_20
  rsme_20 <- sqrt(mean((data_20$Shrub_biomass - data_20$Pred_biomass)^2,na.rm = T))
  rsme_20_list[i] <- rsme_20
  # Calcualte bias
  bias_20 <- mean(data_20$Pred_biomass - data_20$Shrub_biomass,na.rm=T)
  bias_20_list[i] <- bias_20
  # # Plot fit
  # plot(data_20$Shrub_biomass~data_20$Pred_biomass,col="blue")
  # points(data_80$Shrub_biomass~data_80$Pred_biomass,col="red")
  # abline(lm(data_20$Shrub_biomass~data_20$Pred_biomass),col="blue")
  # abline(lm(data_80$Shrub_biomass~data_80$Pred_biomass),col="red")
  # abline(1,1,lty=2)
  rsq80 <- summary(lm(data_80$Shrub_biomass~data_80$Pred_biomass))$adj.r.squared
  rsq20 <- summary(lm(data_20$Shrub_biomass~data_20$Pred_biomass))$adj.r.squared
  # legend("bottomleft",legend = c(paste("R?(20):", round(rsq20, 2)),
  #                                paste("R?(80):", round(rsq80, 2))),
  #        col = c("blue", "red"),lty = 1)
  rsq_80_list[i] <- rsq80
  rsq_20_list[i] <- rsq20
  # Add a column for iteration and dataset type to data_80 and data_20
  data_80$Iteration <- i
  data_80$Dataset <- "Training"
  data_20$Iteration <- i
  data_20$Dataset <- "Testing"
  
  # Append the predicted and original values to combined_results
  combined_results_shrub <- rbind(
    combined_results_shrub,
    data_80[, c("Iteration", "Dataset", "Transect", "Patch", "Shrub_biomass", "Pred_biomass")],
    data_20[, c("Iteration", "Dataset", "Transect", "Patch", "Shrub_biomass", "Pred_biomass")]
  )
}

# Summary statistics of MAE and RSME
cv_mod$model_shrub_lulc[1] <- AIC(model_shrub_lulc)
cv_mod$model_shrub_lulc[2] <- summary(model_shrub_lulc)$r.sq
cv_mod$model_shrub_lulc[3] <- summary(model_shrub_lulc)$dev.expl
cv_mod$model_shrub_lulc[4] <- as.numeric(summary(model_shrub_lulc)$sp.criterion)
cv_mod$model_shrub_lulc[5] <- summary(model_shrub_lulc)$n
cv_mod$model_shrub_lulc[6] <- length(summary(model_shrub_lulc)$edf)
cv_mod$model_shrub_lulc[7] <- sum(influence(model_shrub_lulc))
cv_mod$model_shrub_lulc[8] <- df.residual(model_shrub_lulc)
cv_mod$model_shrub_lulc[9] <- mean(mae_80_list,na.rm=T)
cv_mod$model_shrub_lulc[10] <- mean(mae_20_list,na.rm=T)
cv_mod$model_shrub_lulc[11] <- mean(rsme_80_list,na.rm=T)
cv_mod$model_shrub_lulc[12] <- mean(rsme_20_list,na.rm=T)
cv_mod$model_shrub_lulc[13] <- mean(bias_80_list,na.rm=T)
cv_mod$model_shrub_lulc[14] <- mean(bias_20_list,na.rm=T)
cv_mod$model_shrub_lulc[15] <- mean(rsq_80_list, na.rm = T)
cv_mod$model_shrub_lulc[16] <- mean(rsq_20_list, na.rm = T)

print(cv_mod)

training <- subset(combined_results_shrub,Dataset == "Training")
testing <- subset(combined_results_shrub,Dataset == "Testing")

ggplot(combined_results_shrub, aes(x = Shrub_biomass, y = Pred_biomass, color = Dataset)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(data = training, aes(x = Shrub_biomass, y = Pred_biomass), method = "lm", color = "red", se = FALSE) +
  geom_smooth(data = testing, aes(x = Shrub_biomass, y = Pred_biomass), method = "lm", color = "blue", se = FALSE) +
  xlim(0,1000) +
  ylim(0,1000) +
  labs(title = "Predicted vs Actual Biomass (All Iterations)",
       x = "Shrub Biomass",
       y = "Predicted Biomass") +
  scale_color_manual(
    values = c("Training" = "red", "Testing" = "blue"),  # Colors for the legend
    labels = c("Training Set", "Test Set")  # Custom legend text
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14), 
    axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    legend.position = c(0.15, 0.93),
    legend.title = element_blank(),  
    legend.text = element_text(size = 14)) +
  geom_text(
    data = data.frame(x = 840, y = 50, label = paste("Dev.Exp =", round(cv_mod$model_shrub_lulc[3],2),"\nRSME =", round(cv_mod$model_shrub_lulc[14],2),"\nBias =",round(cv_mod$model_shrub_lulc[18],2), "\nR²_20 =", round(cv_mod$model_shrub_lulc[22],2))),
    aes(x = x, y = y, label = label, hjust = 0),
    size = 5, color = "black", fontface = "italic")


# Plot
ggplot(na.omit(combined_results_shrub), aes(x = Biomass_Interval, y = Pred_biomass, fill = Dataset)) +
  geom_boxplot(alpha = 0.6, position = position_dodge(width = 0.8)) +  # Boxplots for each interval
  geom_abline(slope = 100, intercept = 0, linetype = "dashed", color = "black") + # Reference line
  ylim(0,1100) +
  labs(
    title = "Predicted vs Actual Biomass (Grouped by Intervals)",
    x = "Shrub Biomass Interval",
    y = "Predicted Shrub Biomass"
  ) +
  scale_fill_manual(
    values = c("Training" = "red", "Testing" = "blue"),  # Colors for the legend
    labels = c("Training Set", "Test Set")  # Custom legend text
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    legend.position = c(0.11, 0.93),
    legend.title = element_blank(),
    legend.text = element_text(size = 13)
  ) +
  geom_text(
    data = data.frame(x = 0.85, y = 900, label = paste("Dev.Exp = 0.73", "\nRSME = 167.30 g/m²", "\nBias = -15.84 g/m²","\nR²_20 = 0.35")),
    aes(x = x, y = y, label = label, hjust = 0),
    size = 4.5, color = "black",
    inherit.aes = F
  )



##### Forb Biomass #####

model_forb <- gam(Forb_biomass ~ s(Slope_mean, k=5) +
                        s(Aspect_mean, bs = "cc", k = 5) +
                        s(NDVI_mean, by = LULC_RF, k = 5) + s(VHeight_min, by = LULC_RF, k = 5) + 
                        s(VHeight_max, by = LULC_RF, k = 5) + s(CC_min, k = 5) + s(vrd_0_mean, k = 5) + 
                        s(vrd_0.5_mean, k = 5) + s(vrd_1_mean, k = 5) + s(vrd_2_mean, k = 5) + s(vrd_2.5_mean, k = 5) +
                        s(Transect, bs = "re"),
                      data = biomass_data,
                      method = 'REML',
                      family = tw(),select=TRUE)

summary(model_forb)
par(mfrow=c(2,2))
gam.check(model_forb)

model_forb_lulc <- gam(Forb_biomass ~ LULC_RF + s(Slope_mean, k=5) +
                             s(Aspect_mean, bs = "cc", k = 5) +
                             s(NDVI_mean, by = LULC_RF, k = 5) + s(VHeight_min, by = LULC_RF, k = 5) + 
                             s(VHeight_max, by = LULC_RF, k = 5) + s(CC_min, k = 5) + s(vrd_0_mean, k = 5) + 
                             s(vrd_0.5_mean, k = 5) + s(vrd_1_mean, k = 5) + s(vrd_2_mean, k = 5) + s(vrd_2.5_mean, k = 5) +
                             s(Transect, bs = "re"),
                      data = biomass_data,
                      method = 'REML',
                      family = tw(),select=TRUE)

summary(model_forb_lulc)
par(mfrow=c(2,2))
gam.check(model_forb_lulc)

AIC(model_forb,model_forb_lulc)
anova(model_forb,model_forb_lulc,test="Chisq")

par(mfrow=c(3,2))
visreg(model_forb_lulc, "NDVI_mean", by = "LULC_RF", type = "conditional",main="NDVI_mean")
visreg(model_forb_lulc, "VHeight_min", by = "LULC_RF", type = "conditional",main="Vegetation Height (min)")
visreg(model_forb_lulc, "VHeight_max", by = "LULC_RF", type = "conditional",main="Vegetation height (max)")
visreg(model_forb_lulc, "LULC_RF", type = "conditional",main="Land-use")
visreg(model_forb_lulc, "Slope_mean", type = "conditional",main="Slope")
visreg(model_forb_lulc, "Aspect_mean", type = "conditional",main="Aspect")
visreg(model_forb_lulc, "vrd_1_mean", type = "conditional",main="VRD_1")
visreg(model_forb_lulc, "vrd_2_mean", type = "conditional",main="VRD_2")
visreg(model_forb_lulc, "vrd_2.5_mean", type = "conditional",main="VRD_2.5")


#### Cross-validation ####

combined_results_forb <- data.frame(
  Iteration = integer(),
  Dataset = character(),
  Transect = integer(),
  Patch = factor(),
  Forb_biomass = numeric(),
  Pred_biomass = numeric(),
  stringsAsFactors = FALSE
)


# Loop over the number of iterations
for (i in 1:n_iterations) {
  train_transects <- list()
  test_transects <- list()
  # Loop over each LULC_RF class to perform stratified sampling
  for (lulc in lulc_classes) {
    # Subset data for the current LULC_RF class
    class_data <- biomass_data[biomass_data$LULC_RF == lulc, ]
    # Get the unique transects for the current class
    unique_transects <- unique(class_data$Transect)
    # Perform stratified sampling: select 80% of the transects for training
    selected_train_transects <- sample(unique_transects, size = round(0.8 * length(unique_transects)), replace = FALSE)
    # Store the selected transects for training and testing
    train_transects[[lulc]] <- selected_train_transects
    test_transects[[lulc]] <- setdiff(unique_transects, selected_train_transects)
  }
  # Combine the selected transects across all LULC_RF classes for training and testing
  train_transects_combined <- unlist(train_transects)
  test_transects_combined <- unlist(test_transects)
  # Split the data into 80% training and 20% testing based on the selected transects
  data_80 <- biomass_data[biomass_data$Transect %in% train_transects_combined, ]
  data_20 <- biomass_data[biomass_data$Transect %in% test_transects_combined, ]
  # Ensure the levels of Patch are consistent in both training and testing data
  levels(data_80$Patch) <- levels(biomass_data$Patch)
  levels(data_20$Patch) <- levels(biomass_data$Patch)
  #   # Fit the model on the 80% training data
  model_rf_80 <- gam(model_forb_lulc$formula,
                     data = data_80,
                     method = 'REML',
                     family = tw(),select=T)
  # Predict on the 80% training data
  data_80$Pred_biomass <- predict(model_rf_80, newdata = data_80, type = "response", exclude = "s(Transect)")
  # Calculate the prediction error (e.g., Mean Absolute Error)
  mae_80 <- mean(abs(data_80$Forb_biomass - data_80$Pred_biomass), na.rm = TRUE)
  mae_80_list[i] <- mae_80
  rsme_80 <- sqrt(mean((data_80$Forb_biomass - data_80$Pred_biomass)^2,na.rm = T))
  rsme_80_list[i] <- rsme_80
  # Calcualte bias
  bias_80 <- mean(data_80$Pred_biomass - data_80$Forb_biomass,na.rm=T)
  bias_80_list[i] <- bias_80
  # Predict on the 20% training data
  data_20$Pred_biomass <- predict(model_rf_80, newdata = data_20, type = "response", exclude = "s(Transect)")
  # Calculate the prediction error (e.g., Mean Absolute Error)
  mae_20 <- mean(abs(data_20$Forb_biomass - data_20$Pred_biomass), na.rm = TRUE)
  mae_20_list[i] <- mae_20
  rsme_20 <- sqrt(mean((data_20$Forb_biomass - data_20$Pred_biomass)^2,na.rm = T))
  rsme_20_list[i] <- rsme_20
  # Calcualte bias
  bias_20 <- mean(data_20$Pred_biomass - data_20$Forb_biomass,na.rm=T)
  bias_20_list[i] <- bias_20
  # # Plot fit
  # plot(data_20$Forb_biomass~data_20$Pred_biomass,col="blue")
  # points(data_80$Forb_biomass~data_80$Pred_biomass,col="red")
  # abline(lm(data_20$Forb_biomass~data_20$Pred_biomass),col="blue")
  # abline(lm(data_80$Forb_biomass~data_80$Pred_biomass),col="red")
  # abline(1,1,lty=2)
  rsq80 <- summary(lm(data_80$Forb_biomass~data_80$Pred_biomass))$adj.r.squared
  rsq20 <- summary(lm(data_20$Forb_biomass~data_20$Pred_biomass))$adj.r.squared
  # legend("bottomleft",legend = c(paste("R?(20):", round(rsq20, 2)),
  #                                paste("R?(80):", round(rsq80, 2))),
  #        col = c("blue", "red"),lty = 1)
  rsq_80_list[i] <- rsq80
  rsq_20_list[i] <- rsq20
  # Add a column for iteration and dataset type to data_80 and data_20
  data_80$Iteration <- i
  data_80$Dataset <- "Training"
  data_20$Iteration <- i
  data_20$Dataset <- "Testing"
  
  # Append the predicted and original values to combined_results
  combined_results_forb <- rbind(
    combined_results_forb,
    data_80[, c("Iteration", "Dataset", "Transect", "Patch", "Forb_biomass", "Pred_biomass")],
    data_20[, c("Iteration", "Dataset", "Transect", "Patch", "Forb_biomass", "Pred_biomass")]
  )
}

training <- subset(combined_results_forb,Dataset == "Training")
testing <- subset(combined_results_forb,Dataset == "Testing")


# Summary statistics of MAE and RSME
cv_mod$model_forb_lulc[1] <- AIC(model_forb_lulc)
cv_mod$model_forb_lulc[2] <- summary(model_forb_lulc)$r.sq
cv_mod$model_forb_lulc[3] <- summary(model_forb_lulc)$dev.expl
cv_mod$model_forb_lulc[4] <- as.numeric(summary(model_forb_lulc)$sp.criterion)
cv_mod$model_forb_lulc[5] <- summary(model_forb_lulc)$n
cv_mod$model_forb_lulc[6] <- length(summary(model_forb_lulc)$edf)
cv_mod$model_forb_lulc[7] <- sum(influence(model_forb_lulc))
cv_mod$model_forb_lulc[8] <- df.residual(model_forb_lulc)
cv_mod$model_forb_lulc[9] <- mean(mae_80_list,na.rm=T)
cv_mod$model_forb_lulc[10] <- mean(mae_20_list,na.rm=T)
cv_mod$model_forb_lulc[11] <- mean(rsme_80_list,na.rm=T)
cv_mod$model_forb_lulc[12] <- mean(rsme_20_list,na.rm=T)
cv_mod$model_forb_lulc[13] <- mean(bias_80_list,na.rm=T)
cv_mod$model_forb_lulc[14] <- mean(bias_20_list,na.rm=T)
cv_mod$model_forb_lulc[15] <- mean(rsq_80_list, na.rm = T)
cv_mod$model_forb_lulc[16] <- mean(rsq_20_list, na.rm = T)

print(cv_mod)


ggplot(combined_results_forb, aes(x = Forb_biomass, y = Pred_biomass, color = Dataset)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(data = training, aes(x = Forb_biomass, y = Pred_biomass), method = "lm", color = "red", se = FALSE) +
  geom_smooth(data = testing, aes(x = Forb_biomass, y = Pred_biomass), method = "lm", color = "blue", se = FALSE) +
  facet_wrap(~ Transect) +
  xlim(0,100) +
  ylim(0,100) +
  labs(title = "Predicted vs Actual Biomass (All Iterations)",
       x = "Forb Biomass",
       y = "Predicted Biomass") +
  scale_color_manual(
    values = c("Training" = "red", "Testing" = "blue"),  # Colors for the legend
    labels = c("Training Set", "Test Set")  # Custom legend text
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14), 
    axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    legend.position = c(0.15, 0.93),
    legend.title = element_blank(),  
    legend.text = element_text(size = 14)) +
  geom_text(
    data = data.frame(x = 84, y = 5, label = paste("Dev.Exp =", round(summary(model_forb_lulc)$dev.expl, 2), "\nRSME =", round(mean(rsme_20_list,na.rm=T), 2), "\nBias =", round(mean(bias_20_list,na.rm=T), 2), "\nR²_20 =", round(mean(rsq_20_list, na.rm = T),2))),
    aes(x = x, y = y, label = label, hjust = 0),
    size = 5, color = "black", fontface = "italic"
  )



##### Grass Biomass #####

model_grass <- gam(Grass_biomass ~ s(Slope_mean, k=5) +
                         s(Aspect_mean, bs = "cc", k = 5) +
                         s(NDVI_mean, by = LULC_RF, k = 5) + s(VHeight_min, by = LULC_RF, k = 5) + 
                         s(VHeight_max, by = LULC_RF, k = 5) + s(CC_min, k = 5) + s(vrd_0_mean, k = 5) + 
                         s(vrd_0.5_mean, k = 5) + s(vrd_1_mean, k = 5) + s(vrd_2_mean, k = 5) + s(vrd_2.5_mean, k = 5) +
                         s(Transect, bs = "re"),
                       data = biomass_data,
                       method = 'REML',
                       family = tw(),select=TRUE)

summary(model_grass)
par(mfrow=c(2,2))
gam.check(model_grass)

model_grass_lulc <- gam(Grass_biomass ~ LULC_RF + s(Slope_mean, k=5) +
                              s(Aspect_mean, bs = "cc", k = 5) +
                              s(NDVI_mean, by = LULC_RF, k = 5) + s(VHeight_min, by = LULC_RF, k = 5) + 
                              s(VHeight_max, by = LULC_RF, k = 5) + s(CC_min, k = 5) + s(vrd_0_mean, k = 5) + 
                              s(vrd_0.5_mean, k = 5) + s(vrd_1_mean, k = 5) + s(vrd_2_mean, k = 5) + s(vrd_2.5_mean, k = 5) +
                              s(Transect, bs = "re"),
                       data = biomass_data,
                       method = 'REML',
                       family = tw(),select=TRUE)

summary(model_grass_lulc)
par(mfrow=c(2,2))
gam.check(model_grass_lulc)

AIC(model_grass, model_grass_lulc)
anova(model_grass, model_grass_lulc, test = "Chisq")

par(mfrow=c(2,2))
visreg(model_grass_lulc, "NDVI_mean", by = "LULC_RF", type = "conditional",main="NDVI_mean")
#visreg(model_grass_lulc, "VHeight_min", by = "LULC_RF", type = "conditional",main="Vegetation Height (min)")
visreg(model_grass_lulc, "VHeight_max", by = "LULC_RF", type = "conditional",main="Vegetation height (max)")
visreg(model_grass_lulc, "LULC_RF", type = "conditional",main="Land-use")
visreg(model_grass_lulc, "Slope_mean", type = "conditional",main="Slope")
visreg(model_grass_lulc, "Aspect_mean", type = "conditional",main="Aspect")
visreg(model_grass_lulc, "vrd_0_mean", type = "conditional",main="VRD_0")

#### Cross-validation ####

combined_results_grass <- data.frame(
  Iteration = integer(),
  Dataset = character(),
  Transect = integer(),
  Patch = factor(),
  Grass_biomass = numeric(),
  Pred_biomass = numeric(),
  stringsAsFactors = FALSE
)


# Loop over the number of iterations
for (i in 1:n_iterations) {
  train_transects <- list()
  test_transects <- list()
  # Loop over each LULC_RF class to perform stratified sampling
  for (lulc in lulc_classes) {
    # Subset data for the current LULC_RF class
    class_data <- biomass_data[biomass_data$LULC_RF == lulc, ]
    # Get the unique transects for the current class
    unique_transects <- unique(class_data$Transect)
    # Perform stratified sampling: select 80% of the transects for training
    selected_train_transects <- sample(unique_transects, size = round(0.8 * length(unique_transects)), replace = FALSE)
    # Store the selected transects for training and testing
    train_transects[[lulc]] <- selected_train_transects
    test_transects[[lulc]] <- setdiff(unique_transects, selected_train_transects)
  }
  # Combine the selected transects across all LULC_RF classes for training and testing
  train_transects_combined <- unlist(train_transects)
  test_transects_combined <- unlist(test_transects)
  # Split the data into 80% training and 20% testing based on the selected transects
  data_80 <- biomass_data[biomass_data$Transect %in% train_transects_combined, ]
  data_20 <- biomass_data[biomass_data$Transect %in% test_transects_combined, ]
  # Ensure the levels of Patch are consistent in both training and testing data
  levels(data_80$Patch) <- levels(biomass_data$Patch)
  levels(data_20$Patch) <- levels(biomass_data$Patch)
  #   # Fit the model on the 80% training data
  model_rf_80 <- gam(model_grass_lulc$formula,
                     data = data_80,
                     method = 'REML',
                     family = tw(),select=T)
  # Predict on the 80% training data
  data_80$Pred_biomass <- predict(model_rf_80, newdata = data_80, type = "response", exclude = "s(Transect)")
  # Calculate the prediction error (e.g., Mean Absolute Error)
  mae_80 <- mean(abs(data_80$Grass_biomass - data_80$Pred_biomass), na.rm = TRUE)
  mae_80_list[i] <- mae_80
  rsme_80 <- sqrt(mean((data_80$Grass_biomass - data_80$Pred_biomass)^2,na.rm = T))
  rsme_80_list[i] <- rsme_80
  # Calcualte bias
  bias_80 <- mean(data_80$Pred_biomass - data_80$Grass_biomass,na.rm=T)
  bias_80_list[i] <- bias_80
  # Predict on the 20% training data
  data_20$Pred_biomass <- predict(model_rf_80, newdata = data_20, type = "response", exclude = "s(Transect)")
  # Calculate the prediction error (e.g., Mean Absolute Error)
  mae_20 <- mean(abs(data_20$Grass_biomass - data_20$Pred_biomass), na.rm = TRUE)
  mae_20_list[i] <- mae_20
  rsme_20 <- sqrt(mean((data_20$Grass_biomass - data_20$Pred_biomass)^2,na.rm = T))
  rsme_20_list[i] <- rsme_20
  # Calcualte bias
  bias_20 <- mean(data_20$Pred_biomass - data_20$Grass_biomass,na.rm=T)
  bias_20_list[i] <- bias_20
  # # Plot fit
  # plot(data_20$Grass_biomass~data_20$Pred_biomass,col="blue")
  # points(data_80$Grass_biomass~data_80$Pred_biomass,col="red")
  # abline(lm(data_20$Grass_biomass~data_20$Pred_biomass),col="blue")
  # abline(lm(data_80$Grass_biomass~data_80$Pred_biomass),col="red")
  # abline(1,1,lty=2)
  rsq80 <- summary(lm(data_80$Grass_biomass~data_80$Pred_biomass))$adj.r.squared
  rsq20 <- summary(lm(data_20$Grass_biomass~data_20$Pred_biomass))$adj.r.squared
  # legend("bottomleft",legend = c(paste("R?(20):", round(rsq20, 2)),
  #                                paste("R?(80):", round(rsq80, 2))),
  #        col = c("blue", "red"),lty = 1)
  rsq_80_list[i] <- rsq80
  rsq_20_list[i] <- rsq20
  # Add a column for iteration and dataset type to data_80 and data_20
  data_80$Iteration <- i
  data_80$Dataset <- "Training"
  data_20$Iteration <- i
  data_20$Dataset <- "Testing"
  
  # Append the predicted and original values to combined_results
  combined_results_grass <- rbind(
    combined_results_grass,
    data_80[, c("Iteration", "Dataset", "Transect", "Patch", "Grass_biomass", "Pred_biomass")],
    data_20[, c("Iteration", "Dataset", "Transect", "Patch", "Grass_biomass", "Pred_biomass")]
  )
}

training <- subset(combined_results_grass,Dataset == "Training")
testing <- subset(combined_results_grass,Dataset == "Testing")

# Summary statistics of MAE and RSME
cv_mod$model_grass_lulc[1] <- AIC(model_grass_lulc)
cv_mod$model_grass_lulc[2] <- summary(model_grass_lulc)$r.sq
cv_mod$model_grass_lulc[3] <- summary(model_grass_lulc)$dev.expl
cv_mod$model_grass_lulc[4] <- as.numeric(summary(model_grass_lulc)$sp.criterion)
cv_mod$model_grass_lulc[5] <- summary(model_grass_lulc)$n
cv_mod$model_grass_lulc[6] <- length(summary(model_grass_lulc)$edf)
cv_mod$model_grass_lulc[7] <- sum(influence(model_grass_lulc))
cv_mod$model_grass_lulc[8] <- df.residual(model_grass_lulc)
cv_mod$model_grass_lulc[9] <- mean(mae_80_list,na.rm=T)
cv_mod$model_grass_lulc[10] <- mean(mae_20_list,na.rm=T)
cv_mod$model_grass_lulc[11] <- mean(rsme_80_list,na.rm=T)
cv_mod$model_grass_lulc[12] <- mean(rsme_20_list,na.rm=T)
cv_mod$model_grass_lulc[13] <- mean(bias_80_list,na.rm=T)
cv_mod$model_grass_lulc[14] <- mean(bias_20_list,na.rm=T)
cv_mod$model_grass_lulc[15] <- mean(rsq_80_list, na.rm = T)
cv_mod$model_grass_lulc[16] <- mean(rsq_20_list, na.rm = T)

print(cv_mod)

ggplot(combined_results_grass, aes(x = Grass_biomass, y = Pred_biomass, color = Dataset)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(data = training, aes(x = Grass_biomass, y = Pred_biomass), method = "lm", color = "red", se = FALSE) +
  geom_smooth(data = testing, aes(x = Grass_biomass, y = Pred_biomass), method = "lm", color = "blue", se = FALSE) +
  xlim(0,100) +
  ylim(0,100) +
  labs(title = "Predicted vs Actual Biomass (All Iterations)",
       x = "Grass Biomass",
       y = "Predicted Biomass") +
  scale_color_manual(
    values = c("Training" = "red", "Testing" = "blue"),  # Colors for the legend
    labels = c("Training Set", "Test Set")  # Custom legend text
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14), 
    axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    legend.position = c(0.15, 0.93),
    legend.title = element_blank(),  
    legend.text = element_text(size = 14)) +
  geom_text(
    data = data.frame(x = 84, y = 5, label = paste("Dev.Exp =", round(summary(model_grass_lulc)$dev.expl, 2), "\nRSME =", round(mean(rsme_20_list,na.rm=T), 2), "\nBias =", round(mean(bias_20_list,na.rm=T), 2), "\nR²_20 =", round(mean(rsq_20_list, na.rm = T),2))),
    aes(x = x, y = y, label = label, hjust = 0),
    size = 5, color = "black", fontface = "italic"
  )



##### Herbaceous Biomass ####

model_herbs <- gam(Herb_biomass ~ s(Slope_mean, k=5) +
                          s(Aspect_mean, bs = "cc", k = 5) +
                          s(NDVI_mean, by = LULC_RF, k = 5) + s(VHeight_min, by = LULC_RF, k = 5) + 
                          s(VHeight_max, by = LULC_RF, k = 5) + s(CC_min, k = 5) + s(vrd_0_mean, k = 5) + 
                          s(vrd_0.5_mean, k = 5) + s(vrd_1_mean, k = 5) + s(vrd_2_mean, k = 5) + s(vrd_2.5_mean, k = 5) +
                          s(Transect, bs = "re"),
                        data = biomass_data,
                        method = 'REML',
                        family = tw(),select=TRUE)

summary(model_herbs)
par(mfrow=c(2,2))
gam.check(model_herbs)

model_herbs_lulc <- gam(Herb_biomass ~ LULC_RF + s(Slope_mean, k=5) +
                               s(Aspect_mean, bs = "cc", k = 5) +
                               s(NDVI_mean, by = LULC_RF, k = 5) + s(VHeight_min, by = LULC_RF, k = 5) + 
                               s(VHeight_max, by = LULC_RF, k = 5) + s(CC_min, k = 5) + s(vrd_0_mean, k = 5) + 
                               s(vrd_0.5_mean, k = 5) + s(vrd_1_mean, k = 5) + s(vrd_2_mean, k = 5) + s(vrd_2.5_mean, k = 5) +
                               s(Transect, bs = "re"),
                        data = biomass_data,
                        method = 'REML',
                        family = tw(),select=TRUE)

summary(model_herbs_lulc)
par(mfrow=c(2,2))
gam.check(model_herbs_lulc)

AIC(model_herbs, model_herbs_lulc)
anova(model_herbs, model_herbs_lulc, test = "Chisq")

par(mfrow=c(2,2))
visreg(model_herbs_lulc, "NDVI_mean", by = "LULC_RF", type = "conditional",main="NDVI_mean")
#visreg(model_herbs_lulc, "VHeight_min", by = "LULC_RF", type = "conditional",main="Vegetation Height (min)")
visreg(model_herbs_lulc, "VHeight_max", by = "LULC_RF", type = "conditional",main="Vegetation height (max)")
visreg(model_herbs_lulc, "LULC_RF", type = "conditional",main="Land-use")
visreg(model_herbs_lulc, "Slope_mean", type = "conditional",main="Slope")
visreg(model_herbs_lulc, "Aspect_mean", type = "conditional",main="Aspect")
visreg(model_herbs_lulc, "CC_min", type = "conditional",main="Canopy Cover")
visreg(model_herbs_lulc, "vrd_0_mean", type = "conditional",main="VRD_0")
visreg(model_herbs_lulc, "vrd_0.5_mean", type = "conditional",main="VRD_0.5")
visreg(model_herbs_lulc, "vrd_2.5_mean", type = "conditional",main="VRD_2.5")

#### Cross-validation ####

combined_results_herb <- data.frame(
  Iteration = integer(),
  Dataset = character(),
  Transect = integer(),
  Patch = factor(),
  Herb_biomass = numeric(),
  Pred_biomass = numeric(),
  stringsAsFactors = FALSE
)


# Loop over the number of iterations
for (i in 1:n_iterations) {
  train_transects <- list()
  test_transects <- list()
  # Loop over each LULC_RF class to perform stratified sampling
  for (lulc in lulc_classes) {
    # Subset data for the current LULC_RF class
    class_data <- biomass_data[biomass_data$LULC_RF == lulc, ]
    # Get the unique transects for the current class
    unique_transects <- unique(class_data$Transect)
    # Perform stratified sampling: select 80% of the transects for training
    selected_train_transects <- sample(unique_transects, size = round(0.8 * length(unique_transects)), replace = FALSE)
    # Store the selected transects for training and testing
    train_transects[[lulc]] <- selected_train_transects
    test_transects[[lulc]] <- setdiff(unique_transects, selected_train_transects)
  }
  # Combine the selected transects across all LULC_RF classes for training and testing
  train_transects_combined <- unlist(train_transects)
  test_transects_combined <- unlist(test_transects)
  # Split the data into 80% training and 20% testing based on the selected transects
  data_80 <- biomass_data[biomass_data$Transect %in% train_transects_combined, ]
  data_20 <- biomass_data[biomass_data$Transect %in% test_transects_combined, ]
  # Ensure the levels of Patch are consistent in both training and testing data
  levels(data_80$Patch) <- levels(biomass_data$Patch)
  levels(data_20$Patch) <- levels(biomass_data$Patch)
  #   # Fit the model on the 80% training data
  model_rf_80 <- gam(model_herbs_lulc $formula,
                     data = data_80,
                     method = 'REML',
                     family = tw(),select=T)
  # Predict on the 80% training data
  data_80$Pred_biomass <- predict(model_rf_80, newdata = data_80, type = "response", exclude = "s(Transect)")
  # Calculate the prediction error (e.g., Mean Absolute Error)
  mae_80 <- mean(abs(data_80$Herb_biomass - data_80$Pred_biomass), na.rm = TRUE)
  mae_80_list[i] <- mae_80
  rsme_80 <- sqrt(mean((data_80$Herb_biomass - data_80$Pred_biomass)^2,na.rm = T))
  rsme_80_list[i] <- rsme_80
  # Calcualte bias
  bias_80 <- mean(data_80$Pred_biomass - data_80$Herb_biomass,na.rm=T)
  bias_80_list[i] <- bias_80
  # Predict on the 20% training data
  data_20$Pred_biomass <- predict(model_rf_80, newdata = data_20, type = "response", exclude = "s(Transect)")
  # Calculate the prediction error (e.g., Mean Absolute Error)
  mae_20 <- mean(abs(data_20$Herb_biomass - data_20$Pred_biomass), na.rm = TRUE)
  mae_20_list[i] <- mae_20
  rsme_20 <- sqrt(mean((data_20$Herb_biomass - data_20$Pred_biomass)^2,na.rm = T))
  rsme_20_list[i] <- rsme_20
  # Calcualte bias
  bias_20 <- mean(data_20$Pred_biomass - data_20$Herb_biomass,na.rm=T)
  bias_20_list[i] <- bias_20
  # # Plot fit
  # plot(data_20$Herb_biomass~data_20$Pred_biomass,col="blue")
  # points(data_80$Herb_biomass~data_80$Pred_biomass,col="red")
  # abline(lm(data_20$Herb_biomass~data_20$Pred_biomass),col="blue")
  # abline(lm(data_80$Herb_biomass~data_80$Pred_biomass),col="red")
  # abline(1,1,lty=2)
  rsq80 <- summary(lm(data_80$Herb_biomass~data_80$Pred_biomass))$adj.r.squared
  rsq20 <- summary(lm(data_20$Herb_biomass~data_20$Pred_biomass))$adj.r.squared
  # legend("bottomleft",legend = c(paste("R?(20):", round(rsq20, 2)),
  #                                paste("R?(80):", round(rsq80, 2))),
  #        col = c("blue", "red"),lty = 1)
  rsq_80_list[i] <- rsq80
  rsq_20_list[i] <- rsq20
  # Add a column for iteration and dataset type to data_80 and data_20
  data_80$Iteration <- i
  data_80$Dataset <- "Training"
  data_20$Iteration <- i
  data_20$Dataset <- "Testing"
  
  # Append the predicted and original values to combined_results
  combined_results_herb <- rbind(
    combined_results_herb,
    data_80[, c("Iteration", "Dataset", "Transect", "Patch", "Herb_biomass", "Pred_biomass")],
    data_20[, c("Iteration", "Dataset", "Transect", "Patch", "Herb_biomass", "Pred_biomass")]
  )
}

training <- subset(combined_results_herb,Dataset == "Training")
testing <- subset(combined_results_herb,Dataset == "Testing")

# Summary statistics of MAE and RSME
cv_mod$model_herbs_lulc [1] <- AIC(model_herbs_lulc )
cv_mod$model_herbs_lulc [2] <- summary(model_herbs_lulc )$r.sq
cv_mod$model_herbs_lulc [3] <- summary(model_herbs_lulc )$dev.expl
cv_mod$model_herbs_lulc [4] <- as.numeric(summary(model_herbs_lulc )$sp.criterion)
cv_mod$model_herbs_lulc [5] <- summary(model_herbs_lulc )$n
cv_mod$model_herbs_lulc [6] <- length(summary(model_herbs_lulc )$edf)
cv_mod$model_herbs_lulc [7] <- sum(influence(model_herbs_lulc ))
cv_mod$model_herbs_lulc [8] <- df.residual(model_herbs_lulc )
cv_mod$model_herbs_lulc [9] <- mean(mae_80_list,na.rm=T)
cv_mod$model_herbs_lulc [10] <- mean(mae_20_list,na.rm=T)
cv_mod$model_herbs_lulc [11] <- mean(rsme_80_list,na.rm=T)
cv_mod$model_herbs_lulc [12] <- mean(rsme_20_list,na.rm=T)
cv_mod$model_herbs_lulc [13] <- mean(bias_80_list,na.rm=T)
cv_mod$model_herbs_lulc [14] <- mean(bias_20_list,na.rm=T)
cv_mod$model_herbs_lulc [15] <- mean(rsq_80_list, na.rm = T)
cv_mod$model_herbs_lulc [16] <- mean(rsq_20_list, na.rm = T)

print(cv_mod)


ggplot(combined_results_herbs, aes(x = Herb_biomass, y = Pred_biomass, color = Dataset)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(data = training, aes(x = Herb_biomass, y = Pred_biomass), method = "lm", color = "red", se = FALSE) +
  geom_smooth(data = testing, aes(x = Herb_biomass, y = Pred_biomass), method = "lm", color = "blue", se = FALSE) +
  xlim(0,100) +
  ylim(0,100) +
  labs(title = "Predicted vs Actual Biomass (All Iterations)",
       x = "Herbal Biomass",
       y = "Predicted Biomass") +
  scale_color_manual(
    values = c("Training" = "red", "Testing" = "blue"),  # Colors for the legend
    labels = c("Training Set", "Test Set")  # Custom legend text
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14), 
    axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    legend.position = c(0.15, 0.93),
    legend.title = element_blank(),  
    legend.text = element_text(size = 14)) +
  geom_text(
    data = data.frame(x = 84, y = 5, label = paste("Dev.Exp =", round(summary(model_herbs_lulc )$dev.expl, 2), "\nRSME =", round(mean(rsme_20_list,na.rm=T), 2), "\nBias =", round(mean(bias_20_list,na.rm=T), 2), "\nR²_20 =", round(mean(rsq_20_list, na.rm = T),2))),
    aes(x = x, y = y, label = label, hjust = 0),
    size = 5, color = "black", fontface = "italic"
  )


cv_mod %>% 
  mutate_if(is.numeric, round,digits=2)

####### Predict Biomass across study area #######

LULC_color <- c("yellow","darkgoldenrod","darkgoldenrod1","chartreuse2","chartreuse4","gray62")

grid_data <- read.csv("grid_data.csv",header=T,sep=";")
grid_data$LULC_RF <- factor(grid_data$LULC_RF, levels = levels(biomass_data$LULC_RF))
grid_data$Transect <- "A1_2024"
grid_data$Transect <- factor(grid_data$Transect, levels(biomass_data$Transect))

grid_data <- na.omit(grid_data)

#### Total Biomass #####

# Predict biomass

grid_predict_total <- predict.gam(model_total,newdata=grid_data,se.fit=T,type="response", exclude = "s(Transect)")
grid_data$Predict_total <- grid_predict_total$fit
grid_data$SE_total <- grid_predict_total$se.fit
grid_data$CV_total <- grid_predict_total$se.fit/grid_predict_total$fit

# Get summary statistics 

mean(biomass_data$Total_biomass,na.rm = TRUE)
sd(biomass_data$Total_biomass, na.rm = TRUE) / sqrt(nrow(biomass_data))
biomass_data %>%
  group_by(LULC_RF) %>%
  summarise(mean_biomass = format(mean(Total_biomass, na.rm = TRUE), nsmall = 2),
            SE_biomass = format(sd(Total_biomass, na.rm = TRUE) / sqrt(n()), nsmall = 2),
            SD_biomass = format(sd(Total_biomass, na.rm = TRUE), nsmall = 2))

mean(grid_data$Predict_total)
sd(grid_data$Predict_total) / sqrt(nrow(grid_data))
grid_data %>%
  group_by(LULC_RF) %>%
  summarise(mean_biomass = format(mean(Predict_total), nsmall = 2),
            SE_biomass = format(sd(Predict_total) / sqrt(n()), nsmall = 2),
            SD_biomass = format(sd(Predict_total), nsmall = 2))

mean(grid_data$CV_total)
sd(grid_data$CV_total) / sqrt(nrow(grid_data))
grid_data %>%
  group_by(LULC_RF) %>%
  summarise(mean_CV = format(mean(CV_total), nsmall = 2),
            SE_CV = format(sd(CV_total) / sqrt(n()), nsmall = 2),
            SD_CV = format(sd(CV_total), nsmall = 2))

# Plot predictions 

p1 <- ggplot(biomass_data, aes(x = LULC_RF, y=Total_biomass,fill=LULC_RF)) +
  geom_boxplot() +
  scale_fill_manual(values=LULC_color) +
  ylim(0,1000) +
  labs(title="Observed Total Biomass ",x="Land use", y="Total biomass") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.position="none")

p2 <- ggplot(grid_data, aes(x = LULC_RF, y=Predict_total,fill=LULC_RF)) +
  geom_boxplot() +
  scale_fill_manual(values=LULC_color) +
  ylim(0,1000) +
  labs(title="Predicted Total Biomass",x="Land use", y="Total biomass") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.position="none")

print(grid.arrange(p1, p2, nrow=1,widths=c(1/2,1/2)))

p1 <- ggplot(grid_data, aes(x = LULC_RF, y=Predict_total,fill=LULC_RF)) +
  geom_boxplot() +
  scale_fill_manual(values=LULC_color) +
  ylim(0,1000) +
  labs(title="Predicted Total Biomass ",x="Land use", y="Predicted total biomass") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.position="none")

p2 <- ggplot(grid_data, aes(x = LULC_RF, y=CV_total,fill=LULC_RF)) +
  geom_boxplot() +
  scale_fill_manual(values=LULC_color) +
  ylim(0,1) +
  labs(title="Coefficient of Variation (CV) - Total Biomass",x="Land use", y="CV") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.position="none")

print(grid.arrange(p1, p2, nrow=1,widths=c(1/2,1/2)))

#### Shrub Biomass #####

grid_predict_shrub <- predict.gam(model_shrub_lulc,newdata=grid_data,se.fit=T,type="response", exclude = "s(Transect)")
grid_data$Predict_shrub <- grid_predict_shrub$fit
grid_data$SE_shrub <- grid_predict_shrub$se.fit
grid_data$CV_shrub <- grid_predict_shrub$se.fit/grid_predict_shrub$fit

# Get summary statistics 

mean(biomass_data$Shrub_biomass,na.rm = TRUE)
sd(biomass_data$Shrub_biomass, na.rm = TRUE) / sqrt(nrow(biomass_data))
biomass_data %>%
  group_by(LULC_RF) %>%
  summarise(mean_biomass = format(mean(Shrub_biomass, na.rm = TRUE), nsmall = 2),
            SE_biomass = format(sd(Shrub_biomass, na.rm = TRUE) / sqrt(n()), nsmall = 2),
            SD_biomass = format(sd(Shrub_biomass, na.rm = TRUE), nsmall = 2))

mean(grid_data$Predict_shrub)
sd(grid_data$Predict_shrub) / sqrt(nrow(grid_data))
grid_data %>%
  group_by(LULC_RF) %>%
  summarise(mean_biomass = format(mean(Predict_shrub), nsmall = 2),
            SE_biomass = format(sd(Predict_shrub) / sqrt(n()), nsmall = 2),
            SD_biomass = format(sd(Predict_shrub), nsmall = 2))

mean(grid_data$CV_shrub)
sd(grid_data$CV_shrub) / sqrt(nrow(grid_data))
grid_data %>%
  group_by(LULC_RF) %>%
  summarise(mean_CV = format(mean(CV_shrub), nsmall = 2),
            SE_CV = format(sd(CV_shrub) / sqrt(n()), nsmall = 2),
            SD_CV = format(sd(CV_shrub), nsmall = 2))

# Plot predictions 

p1 <- ggplot(biomass_data, aes(x = LULC_RF, y=Shrub_biomass,fill=LULC_RF)) +
  geom_boxplot() +
  scale_fill_manual(values=LULC_color) +
  ylim(0,1000) +
  labs(title="Observed Shrub Biomass ",x="Land use", y="Shrub biomass") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.position="none")

p2 <- ggplot(grid_data, aes(x = LULC_RF, y=Predict_shrub,fill=LULC_RF)) +
  geom_boxplot() +
  scale_fill_manual(values=LULC_color) +
  ylim(0,1000) +
  labs(title="Predicted Shrub Biomass",x="Land use", y="Shrub biomass") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.position="none")

print(grid.arrange(p1, p2, nrow=1,widths=c(1/2,1/2)))

p1 <- ggplot(grid_data, aes(x = LULC_RF, y=Predict_shrub,fill=LULC_RF)) +
  geom_boxplot() +
  scale_fill_manual(values=LULC_color) +
  ylim(0,1000) +
  labs(title="Predicted Shrub Biomass ",x="Land use", y="Shrub biomass") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.position="none")

p2 <- ggplot(grid_data, aes(x = LULC_RF, y=CV_shrub,fill=LULC_RF)) +
  geom_boxplot() +
  scale_fill_manual(values=LULC_color) +
  ylim(0,1) +
  labs(title="Coefficient of Variation (CV) - Shrub Biomass",x="Land use", y="CV") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.position="none")

print(grid.arrange(p1, p2, nrow=1,widths=c(1/2,1/2)))

#### Forb Biomass ####

grid_predict_forb <- predict.gam(model_forb_lulc,newdata=grid_data,se.fit=T,type="response", exclude = "s(Transect)")
grid_data$Predict_forb <- grid_predict_forb$fit
grid_data$SE_forb <- grid_predict_forb$se.fit
grid_data$CV_forb <- grid_predict_forb$se.fit/grid_predict_forb$fit

# Get summary statistics 

mean(biomass_data$Forb_biomass,na.rm = TRUE)
sd(biomass_data$Forb_biomass, na.rm = TRUE) / sqrt(nrow(biomass_data))
biomass_data %>%
  group_by(LULC_RF) %>%
  summarise(mean_biomass = format(mean(Forb_biomass, na.rm = TRUE), nsmall = 2),
            SE_biomass = format(sd(Forb_biomass, na.rm = TRUE) / sqrt(n()), nsmall = 2),
            SD_biomass = format(sd(Forb_biomass, na.rm = TRUE), nsmall = 2))

mean(grid_data$Predict_forb)
sd(grid_data$Predict_forb) / sqrt(nrow(grid_data))
grid_data %>%
  group_by(LULC_RF) %>%
  summarise(mean_biomass = format(mean(Predict_forb), nsmall = 2),
            SE_biomass = format(sd(Predict_forb) / sqrt(n()), nsmall = 2),
            SD_biomass = format(sd(Predict_forb), nsmall = 2))

mean(grid_data$CV_forb)
sd(grid_data$CV_forb) / sqrt(nrow(grid_data))
grid_data %>%
  group_by(LULC_RF) %>%
  summarise(mean_CV = format(mean(CV_forb), nsmall = 2),
            SE_CV = format(sd(CV_forb) / sqrt(n()), nsmall = 2),
            SD_CV = format(sd(CV_forb), nsmall = 2))

# Plot predictions 

p1 <- ggplot(biomass_data, aes(x = LULC_RF, y=Forb_biomass,fill=LULC_RF)) +
  geom_boxplot() +
  scale_fill_manual(values=LULC_color) +
  ylim(0,70) +
  labs(title="Observed Forb Biomass ",x="Land use", y="Forb biomass") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.position="none")

p2 <- ggplot(grid_data, aes(x = LULC_RF, y=Predict_forb,fill=LULC_RF)) +
  geom_boxplot() +
  scale_fill_manual(values=LULC_color) +
  ylim(0,70) +
  labs(title="Predicted Forb Biomass",x="Land use", y="Forb biomass") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.position="none")

print(grid.arrange(p1, p2, nrow=1,widths=c(1/2,1/2)))

p1 <- ggplot(grid_data, aes(x = LULC_RF, y=Predict_forb,fill=LULC_RF)) +
  geom_boxplot() +
  scale_fill_manual(values=LULC_color) +
  ylim(0,70) +
  labs(title="Predicted Forb Biomass",x="Land use", y="Forb biomass") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.position="none")

p2 <- ggplot(grid_data, aes(x = LULC_RF, y=CV_forb,fill=LULC_RF)) +
  geom_boxplot() +
  scale_fill_manual(values=LULC_color) +
  ylim(0,1) +
  labs(title="Coefficient of Variation (CV) - Forb Biomass",x="Land use", y="CV") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.position="none")

print(grid.arrange(p1, p2, nrow=1,widths=c(1/2,1/2)))


#### Grass Biomass ####

grid_predict_grass <- predict.gam(model_grass_lulc,newdata=grid_data,se.fit=T,type="response", exclude = "s(Transect)")
grid_data$Predict_grass <- grid_predict_grass$fit
grid_data$SE_grass <- grid_predict_grass$se.fit
grid_data$CV_grass <- grid_predict_grass$se.fit/grid_predict_grass$fit


# Get summary statistics 

mean(biomass_data$Grass_biomass, na.rm = TRUE)
sd(biomass_data$Grass_biomass, na.rm = TRUE) / sqrt(nrow(biomass_data))
biomass_data %>%
  group_by(LULC_RF) %>%
  summarise(mean_biomass = format(mean(Grass_biomass, na.rm = TRUE), nsmall = 2),
            SE_biomass = format(sd(Grass_biomass, na.rm = TRUE) / sqrt(n()), nsmall = 2),
            SD_biomass = format(sd(Grass_biomass, na.rm = TRUE), nsmall = 2))

mean(grid_data$Predict_grass)
sd(grid_data$Predict_grass) / sqrt(nrow(grid_data))
grid_data %>%
  group_by(LULC_RF) %>%
  summarise(mean_biomass = format(mean(Predict_grass), nsmall = 2),
            SE_biomass = format(sd(Predict_grass) / sqrt(n()), nsmall = 2),
            SD_biomass = format(sd(Predict_grass), nsmall = 2))

mean(grid_data$CV_grass)
sd(grid_data$CV_grass) / sqrt(nrow(grid_data))
grid_data %>%
  group_by(LULC_RF) %>%
  summarise(mean_CV = format(mean(CV_grass), nsmall = 2),
            SE_CV = format(sd(CV_grass) / sqrt(n()), nsmall = 2),
            SD_CV = format(sd(CV_grass), nsmall = 2))

# Plot predictions 

p1 <- ggplot(biomass_data, aes(x = LULC_RF, y=Grass_biomass,fill=LULC_RF)) +
  geom_boxplot() +
  scale_fill_manual(values=LULC_color) +
  ylim(0,70) +
  labs(title="Observed Grass Biomass ",x="Land use", y="Grass biomass") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.position="none")

p2 <- ggplot(grid_data, aes(x = LULC_RF, y=Predict_grass,fill=LULC_RF)) +
  geom_boxplot() +
  scale_fill_manual(values=LULC_color) +
  ylim(0,70) +
  labs(title="Predicted Grass Biomass",x="Land use", y="Grass biomass") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.position="none")

print(grid.arrange(p1, p2, nrow=1,widths=c(1/2,1/2)))

p1 <- ggplot(grid_data, aes(x = LULC_RF, y=Predict_grass,fill=LULC_RF)) +
  geom_boxplot() +
  scale_fill_manual(values=LULC_color) +
  ylim(0,70) +
  labs(title="Predicted Grass Biomass",x="Land use", y="Grass biomass") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.position="none")

p2 <- ggplot(grid_data, aes(x = LULC_RF, y=CV_grass,fill=LULC_RF)) +
  geom_boxplot() +
  scale_fill_manual(values=LULC_color) +
  ylim(0,1) +
  labs(title="Coefficient of Variation (CV) - Grass Biomass",x="Land use", y="CV") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.position="none")

print(grid.arrange(p1, p2, nrow=1,widths=c(1/2,1/2)))

#### Herbaceous Biomass #### 

grid_predict_herb <- predict.gam(model_herbs_lulc,newdata=grid_data,se.fit=T,type="response", exclude = "s(Transect)")
grid_data$Predict_herb <- grid_predict_herb$fit
grid_data$SE_herb <- grid_predict_herb$se.fit
grid_data$CV_herb <- grid_predict_herb$se.fit/grid_predict_herb$fit

# Get summary statistics 

mean(biomass_data$Herb_biomass, na.rm = TRUE)
sd(biomass_data$Herb_biomass, na.rm = TRUE) / sqrt(nrow(biomass_data))
biomass_data %>%
  group_by(LULC_RF) %>%
  summarise(mean_biomass = format(mean(Herb_biomass, na.rm = TRUE), nsmall = 2),
            SE_biomass = format(sd(Herb_biomass, na.rm = TRUE) / sqrt(n()), nsmall = 2),
            SD_biomass = format(sd(Herb_biomass, na.rm = TRUE), nsmall = 2))

mean(grid_data$Predict_herb)
sd(grid_data$Predict_herb) / sqrt(nrow(grid_data))
grid_data %>%
  group_by(LULC_RF) %>%
  summarise(mean_biomass = format(mean(Predict_herb), nsmall = 2),
            SE_biomass = format(sd(Predict_herb) / sqrt(n()), nsmall = 2),
            SD_biomass = format(sd(Predict_herb), nsmall = 2))

mean(grid_data$CV_herb,na.rm=T)
sd(grid_data$CV_herb) / sqrt(nrow(grid_data))
grid_data %>%
  group_by(LULC_RF) %>%
  summarise(mean_CV = format(mean(CV_herb), nsmall = 2),
            SE_CV = format(sd(CV_herb) / sqrt(n()), nsmall = 2),
            SD_CV = format(sd(CV_herb), nsmall = 2))

# Plot predictions 

p1 <- ggplot(biomass_data, aes(x = LULC_RF, y=Herb_biomass,fill=LULC_RF)) +
  geom_boxplot() +
  scale_fill_manual(values=LULC_color) +
  ylim(0,150) +
  labs(title="Observed Herbaceous Biomass ",x="Land use", y="Herbaceous biomass") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.position="none")

p2 <- ggplot(grid_data, aes(x = LULC_RF, y=Predict_herb,fill=LULC_RF)) +
  geom_boxplot() +
  scale_fill_manual(values=LULC_color) +
  ylim(0,150) +
  labs(title="Predicted Herbaceous Biomass",x="Land use", y="Herbaceous biomass") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.position="none")

print(grid.arrange(p1, p2, nrow=1,widths=c(1/2,1/2)))

p1 <- ggplot(grid_data, aes(x = LULC_RF, y=Predict_herb,fill=LULC_RF)) +
  geom_boxplot() +
  scale_fill_manual(values=LULC_color) +
  ylim(0,150) +
  labs(title="Predicted Herbaceous Biomass",x="Land use", y="Herbaceous biomass") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.position="none")

p2 <- ggplot(grid_data, aes(x = LULC_RF, y=CV_herb,fill=LULC_RF)) +
  geom_boxplot() +
  scale_fill_manual(values=LULC_color) +
  ylim(0,1) +
  labs(title="Coefficient of Variation (CV) - Herbaceous Biomass",x="Land use", y="CV") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.position="none")

print(grid.arrange(p1, p2, nrow=1,widths=c(1/2,1/2)))

###### Summary #####

# Save results 
write.table(grid_data,"Grid_predict.csv",row.names=F,sep=",")

list(
  total_biomass = sum(grid_data$Predict_total*100/ 1e6),
  shrub_biomass = sum(grid_data$Predict_shrub*100/ 1e6),
  forb_biomass = sum(grid_data$Predict_forb*100/ 1e6),
  grass_biomass = sum(grid_data$Predict_grass*100/ 1e6),
  herb_biomass = sum(grid_data$Predict_herb*100/ 1e6)
)

biomass_by_LULC <- grid_data %>%
  group_by(LULC_RF) %>%
  summarise(
    Area_ha = (n() * 100) / 10000,
    Total_biomass_tons = sum(Predict_total*100/ 1e6),
    Total_biomass_tons_ha = sum(Predict_total * 100/ 1e6) / Area_ha,
    Shrub_biomass_tons = sum(Predict_shrub * 100/ 1e6),
    Shrub_biomass_tons_ha = sum(Predict_shrub * 100/ 1e6) / Area_ha,
    Forb_biomass_tons = sum(Predict_forb * 100/ 1e6),
    Forb_biomass_tons_ha = sum(Predict_forb * 100 / 1e6) / Area_ha,
    Grass_biomass_tons = sum(Predict_grass * 100/ 1e6),
    Grass_biomass_tons_ha = sum(Predict_grass * 100/ 1e6) / Area_ha,
    Herb_biomass_tons = sum(Predict_herb * 100/ 1e6),
    Herb_biomass_tons_ha = sum(Predict_herb * 100/ 1e6) / Area_ha
  )

grid_data %>%
  group_by(LULC_RF) %>%
  summarise(
    Area = format(round((n() * 100) / 10000, 2),nsmall=2),
    Total_biomass = format(round(sum(Predict_total) * 10, 2), nsmall = 2),
    Shrub_biomass = format(round(sum(Predict_shrub) * 10, 2), nsmall = 2),
    Forb_biomass = format(round(sum(Predict_forb) * 10, 2), nsmall = 2),
    Grass_biomass = format(round(sum(Predict_grass) * 10, 2), nsmall = 2),
    Herb_biomass = format(round(sum(Predict_herb) * 10, 2), nsmall = 2)
  )


par(mfrow=c(1,3))
hist(grid_data$Predict_total,xlim=c(0,1000),breaks=50)
hist(grid_data$SE_total,xlim=c(0,200),breaks=200)
hist(grid_data$CV_total*100,xlim=c(0,100),breaks=20)

par(mfrow=c(1,3))
hist(grid_data$Predict_shrub,xlim=c(0,1000),breaks=100)
hist(grid_data$SE_shrub,xlim=c(0,200),breaks=1000)
hist(grid_data$CV_shrub*100,xlim=c(0,100),breaks=500)

par(mfrow=c(1,3))
hist(grid_data$Predict_forb,xlim=c(0,100),breaks=5)
hist(grid_data$SE_forb,xlim=c(0,50),breaks=10)
hist(grid_data$CV_forb*100,xlim=c(0,120),breaks=1000)

par(mfrow=c(1,3))
hist(grid_data$Predict_grass,xlim=c(0,100),breaks=10)
hist(grid_data$SE_grass,xlim=c(0,50),breaks=10)
hist(grid_data$CV_grass*100,xlim=c(0,120),breaks=5000)

par(mfrow=c(1,3))
hist(grid_data$Predict_herb,xlim=c(0,100),breaks=10)
hist(grid_data$SE_herb,xlim=c(0,50),breaks=10)
hist(grid_data$CV_herb*100,xlim=c(0,120),breaks=2000)



##### Alternative plots ##### 

combined_results_total$Biomass_Interval <- cut(
  combined_results_total$Total_biomass,
  breaks = seq(0, 1200, by = 100),  # Define the intervals (0-100, 100-200, ...)
  include.lowest = TRUE,            # Include the lowest interval
  labels = paste(seq(0, 1100, by = 100), seq(100, 1200, by = 100), sep = "-") # Label intervals
)

combined_results_shrub$Biomass_Interval <- cut(
  combined_results_shrub$Shrub_biomass,
  breaks = seq(0, 1200, by = 100),  # Define the intervals (0-100, 100-200, ...)
  include.lowest = TRUE,            # Include the lowest interval
  labels = paste(seq(0, 1100, by = 100), seq(100, 1200, by = 100), sep = "-") # Label intervals
)

combined_results_herb$Biomass_Interval <- cut(
  combined_results_herb$Herb_biomass,
  breaks = seq(0, 250, by = 25),  # Define the intervals (0-100, 100-200, ...)
  include.lowest = TRUE,            # Include the lowest interval
  labels = paste(seq(0, 225, by = 25), seq(25, 250, by = 25), sep = "-") # Label intervals
)


# Plot
total_plot <- ggplot(na.omit(combined_results_total), aes(x = Biomass_Interval, y = Pred_biomass, fill = Dataset)) +
  geom_boxplot(alpha = 0.6, position = position_dodge(width = 0.8)) +  # Boxplots for each interval
  geom_abline(slope = 100, intercept = 0, linetype = "dashed", color = "black") + # Reference line
  ylim(0,1100) +
  labs(
    title = "",
    x = "Total Biomass",
    y = "Predicted Total Biomass"
  ) +
  scale_fill_manual(
    values = c("Training" = "red", "Testing" = "blue"),  # Colors for the legend
    labels = c("Training Set", "Test Set")  # Custom legend text
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.title.x = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 18)
  ) +
  geom_text(
    data = data.frame(x = 0.8, y = 900, label = paste("Dev.Exp = 0.80", "\nRSME = 177.76 g/m²", "\nBias = -23.83 g/m²","\nR²_20 = 0.30")),
    aes(x = x, y = y, label = label, hjust = 0),
    size = 4.5, color = "red",
    inherit.aes = F
  )

# Plot
shrub_plot <- ggplot(na.omit(combined_results_shrub), aes(x = Biomass_Interval, y = Pred_biomass, fill = Dataset)) +
  geom_boxplot(alpha = 0.6, position = position_dodge(width = 0.8)) +  # Boxplots for each interval
  geom_abline(slope = 100, intercept = 0, linetype = "dashed", color = "black") + # Reference line
  ylim(0,1100) +
  labs(
    title = "",
    x = "Shrub Biomass",
    y = "Predicted Shrub Biomass"
  ) +
  scale_fill_manual(
    values = c("Training" = "red", "Testing" = "blue"),  # Colors for the legend
    labels = c("Training Set", "Test Set")  # Custom legend text
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.title.x = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    legend.position = "none") +
  geom_text(
    data = data.frame(x = 0.85, y = 900, label = paste("Dev.Exp = 0.73", "\nRSME = 167.30 g/m²", "\nBias = -15.84 g/m²","\nR²_20 = 0.35")),
    aes(x = x, y = y, label = label, hjust = 0),
    size = 4.5, color = "red",
    inherit.aes = F
  )


# Plot
herb_plot <- ggplot(na.omit(combined_results_herb), aes(x = Biomass_Interval, y = Pred_biomass, fill = Dataset)) +
  geom_boxplot(alpha = 0.6, position = position_dodge(width = 0.8)) +  # Boxplots for each interval
  geom_abline(slope = 25, intercept = 0, linetype = "dashed", color = "black") + # Reference line
  ylim(0,250) +
  labs(
    title = "",
    x = "Herbaceous Biomass",
    y = "Predicted Herbaceous Biomass"
  ) +
  scale_fill_manual(
    values = c("Training" = "red", "Testing" = "blue"),  # Colors for the legend
    labels = c("Training Set", "Test Set")  # Custom legend text
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.title.x = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    legend.position = "none") +
  geom_text(
    data = data.frame(x = 0.85, y = 210, label = paste("Dev.Exp = 0.73", "\nRSME = 167.30 g/m²", "\nBias = -15.84 g/m²","\nR²_20 = 0.35")),
    aes(x = x, y = y, label = label, hjust = 0),
    size = 4.5, color = "red",
    inherit.aes = F
  )


grid.arrange(total_plot, shrub_plot, herb_plot, layout_matrix = rbind(c(1,1), c(2,3)))



# Plot
total_plot <- ggplot(na.omit(combined_results_total), aes(x = Biomass_Interval, y = Pred_biomass, fill = Dataset)) +
  geom_violin(alpha = 0.6, position = position_dodge(width = 0.8)) +  # Boxplots for each interval
  stat_summary(
    fun = median, geom = "point", position = position_dodge(width = 0.8),
    size = 2, color = "black", shape = 16  # Black points for medians
  ) +
  geom_abline(slope = 200, intercept = 0, linetype = "dashed", color = "black") + # Reference line
  ylim(0,1100) +
  labs(
    title = "",
    x = "Total Biomass",
    y = "Predicted Total Biomass"
  ) +
  scale_fill_manual(
    values = c("Training" = "red", "Testing" = "blue"),  # Colors for the legend
    labels = c("Training Set", "Test Set")  # Custom legend text
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.title.x = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 18)
  ) +
  geom_text(
    data = data.frame(x = 0.5, y = 950, label = paste("Dev.Exp = 0.80", "\nRSME = 177.76 g/m²", "\nBias = -23.83 g/m²","\nR²_20 = 0.30")),
    aes(x = x, y = y, label = label, hjust = 0),
    size = 4.5, color = "red",
    inherit.aes = F
  )

# Plot
shrub_plot <- ggplot(na.omit(combined_results_shrub), aes(x = Biomass_Interval, y = Pred_biomass, fill = Dataset)) +
  geom_violin(alpha = 0.6, position = position_dodge(width = 0.8)) +  # Boxplots for each interval
  stat_summary(
    fun = median, geom = "point", position = position_dodge(width = 0.8),
    size = 2, color = "black", shape = 16  # Black points for medians
  ) +
  geom_abline(slope = 200, intercept = 0, linetype = "dashed", color = "black") + # Reference line
  ylim(0,1100) +
  labs(
    title = "",
    x = "Shrub Biomass",
    y = "Predicted Shrub Biomass"
  ) +
  scale_fill_manual(
    values = c("Training" = "red", "Testing" = "blue"),  # Colors for the legend
    labels = c("Training Set", "Test Set")  # Custom legend text
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.title.x = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    legend.position = "none") +
  geom_text(
    data = data.frame(x = 0.5, y = 950, label = paste("Dev.Exp = 0.73", "\nRSME = 167.30 g/m²", "\nBias = -15.84 g/m²","\nR²_20 = 0.35")),
    aes(x = x, y = y, label = label, hjust = 0),
    size = 4.5, color = "red",
    inherit.aes = F
  )


# Plot
herb_plot <- ggplot(na.omit(combined_results_herb), aes(x = Biomass_Interval, y = Pred_biomass, fill = Dataset)) +
  geom_violin(alpha = 0.6, position = position_dodge(width = 0.8)) +  # Boxplots for each interval
  stat_summary(
    fun = median, geom = "point", position = position_dodge(width = 0.8),
    size = 2, color = "black", shape = 16  # Black points for medians
  ) +
  geom_abline(slope = 50, intercept = 0, linetype = "dashed", color = "black") + # Reference line
  ylim(0,250) +
  labs(
    title = "",
    x = "Herbaceous Biomass",
    y = "Predicted Herbaceous Biomass"
  ) +
  scale_fill_manual(
    values = c("Training" = "red", "Testing" = "blue"),  # Colors for the legend
    labels = c("Training Set", "Test Set")  # Custom legend text
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.title.x = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    legend.position = "none") +
  geom_text(
    data = data.frame(x = 0.5, y = 220, label = paste("Dev.Exp = 0.73", "\nRSME = 167.30 g/m²", "\nBias = -15.84 g/m²","\nR²_20 = 0.35")),
    aes(x = x, y = y, label = label, hjust = 0),
    size = 4.5, color = "red",
    inherit.aes = F
  )


grid.arrange(total_plot, shrub_plot, herb_plot, layout_matrix = rbind(c(1,1), c(2,3)))
