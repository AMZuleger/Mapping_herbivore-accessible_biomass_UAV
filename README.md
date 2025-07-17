# Mapping herbivore-accessible biomass using UAV data
Data and Code used for analysis of the manuscript ***Mapping herbivore-accessible biomass across a heterogeneous mountain landscape using multisensor high-resolution UAV data***

**Abstract:**

Herbivore-accessible biomass (HAB), defined as aboveground biomass under 2 m, including leaves and soft branches, is a key metric for understanding ecosystem function, but remains poorly quantified. We estimated HAB across diverse habitats in the Peneda-Gerês National Park using high-resolution NDVI, LiDAR, topography and field data. Generalized Additive Mixed Models (GAMMs) revealed habitat-specific effects of NDVI and vegetation height, as well as terrain, and structural metrics across plant types. Models were evaluated using hold-out cross-validation on a 20% subset of the field data. The total HAB model performed well (Deviance Explained = 0.77, RMSE20 = 172.38 g/m²), while the shrub model performed slightly worse (Deviance Explained = 0.71, RMSE20 = 410.21 g/m²), and the herbaceous model exhibited a moderate fit and accuracy (Deviance Explained = 0.69 , RMSE20 = 34.25 g/m²). Average total HAB was 1.31 ± 0.83 tons/ha, dominated by shrubs (1.02 tons/ha) compared to herbaceous HAB (0.14 tons/ha). HAB density varied by habitat, highest in shrublands (up to 1.83 ton/ha) and lowest in oak forests (0.85 tons/ha), while agricultural areas supported the most herbaceous HAB (0.68 tons/ha). These values are substantially lower than shrub biomass estimates reported in other studies (e.g., up to 30 tons/ha), reflecting our focus on live biomass <2 m. Prediction uncertainty was low (CV: 22 - 34%), improving on other studies reporting up to 190%, and highlighting the strength of combining spectral and structural data for fine-scale forage estimation. This study provides the first spatially explicit HAB estimates for the area, supporting herbivore ecology and management.

**Data structure:**

- *Mapping_HAB.Rmd* - R Notebook document for the entire analysis.
  - Code runs all models and produces all results and figures (except for maps, which were created using QGIS).
- *Mapping_HAB.html* - HTML document created from R *Mapping_HAB.Rmd* 
  
- *Data.zip* - Contains datasets required for the analysis:

  - *biomass_data_2024.csv* - Field measurements of biomass collected in May, 2024 and remote sensing metrics obtained through Zonal Statistics in QGIS on the plot level
          
  - *grid_data.csv* - Remote sensing metrics across the entire study area on a 10 x 10 m grid for predictions
 
- *Figures.zip* - Contains all figures for the publication (including Supplemental):
  - *Figure_1_Study_Location_and_Sampling_Design.png* - Study Location and Sampling Design.
  - *Figure_2_Variable_Importance.png* - F-Statistics per variable obtained from GAMMs for total, shrub and herbaceous HAB
  - *Figure_3_Cross_validation_results.png* - Predicted vs. observed biomass across total (top), shrub (bottom left), and herbaceous HAB (bottom right) using hold-out cross-validation. 
  - *Figure_4_HAB_Predictions_CV_Maps.png* - HAB predictions and coefficient of variation of the predictions across the entire study area for total, shrub and herbaceous HAB.
  - *Figure_S2_Height_distribution.png* - Histogramms of mean and maximum height measurements from field data per vegetation type.
  - *Figure_S4_Remote_Sensing_Metrics.png* - Maps of all remote sensing metrics used in the final models.
  - *Figure_S5.1_Correlation_plot.png* - Correlation plot between all variables.
  - *Figure_S2.2_Dissimilarity_tree.png* - Dissimilarity tree of all variables.
  - *Figure_S6_Field_measurements.png* - Field measurements of herbivore-accessible AGB per biomass type and habitat class.
  - *Figure_S8_Variable_importance_forb_grass.png* - F-Statistics per variable obtained from GAMMs for forb and grass HAB.
  - *Figure_S9_Cross_validation_results_forb_grass.png* - Predicted vs. observed biomass across forb (left), and grass HAB (right) using hold-out cross-validation. 
  - *Figure_S11_HAB_Predictions_CV_Maps_herb.png* - Biomass prediction maps and Coefficient of variation (CV) across the entire study area for forbs and grasses
  - *Figure_S12.1_Observed_pred_biomass_distribution.png* - Observed and predicted biomass distribution per habitat and vegetation type. 
  - *Figure_S12.2_Pred_biomass_CV_distribution.png* - Predictions and Coefficient of Variation (CV) per habitat and vegetation type. 
  - *Figure_S12.3_Historgramms_HAB_distribution.png* - Histograms of HAB predictions, standard errors and coefficients of variation across vegetation types.
    
Remote Sensing (GIS) data will be made available in a separate repository in the near future and linked here.
          
