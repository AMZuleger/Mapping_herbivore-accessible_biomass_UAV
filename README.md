# Mapping herbivore-accessible biomass using UAV data
Data and Code used for analysis of the manuscript ***Mapping herbivore-accessible biomass across a heterogeneous mountain landscape using multisensor high-resolution UAV data***

**Abstract:**

Accurate estimation of herbivore-accessible biomass is essential for understanding ecosystem dynamics. The advent of drones has facilitated the acquisition of high-resolution imagery, such as multispectral and LiDAR data, opening new possibilities for estimating carbon stocks or fuel loads. However, the estimation of forage for herbivores has received less attention. This study develops a modeling approach to estimate herbivore-accessible biomass, integrating NDVI, LiDAR and topographic information with field measurements through GAMMs to predict herbaceous and shrub biomass across different habitats. We tested it in the Peneda-Gerês National Park across diverse habitats, including pastures, forests, and shrublands with varying biomass densities. Habitat-specific effects of NDVI, above-ground height, topography, and other LiDAR metrics significantly influenced biomass predictions across plant categories. The total biomass model showed a strong fit (R² = 0.74) but moderate precision, while the shrub model showed similar performance (R² = 0.73), and the herbaceous model exhibited a moderate fit with variable accuracy. Total biomass in the study area averaged 1.2 tons/ha, with shrubs contributing the largest share (0.92 tons/ha), and herbaceous biomass accounting for 0.14 tons/ha. Biomass density was highest in high shrub areas (2.10 tons/ha) and lowest in oak forests (0.65 tons/ha), with agricultural areas showing the highest herbaceous biomass (0.68 tons/ha). A key challenge was the high biomass variability within habitats, highlighting the need to consider spatial heterogeneity and include habitat-specific variables in biomass modeling. The integration of spectral and structural remote sensing data enabled us to provide the first plant biomass estimates for the Peneda-Gerês National Park.

**Data structure:**

- *Mapping_heterogeneous_biomass_UAV.R* - R-Code for the entire analysis.
  - Code automatically downloads repository from GitHub and produces all Results and Figures (except for maps, which were created using QGIS).
  
- *Data.zip* - Contains datasets required for the analysis:

  - *biomass_data_2024.csv* - Field measurements of biomass collected in May, 2024 and remote sensing metrics obtained through Zonal Statistics in QGIS on the plot level
          
  - *grid_data.csv* - Remote sensing metrics across the entire study area on a 10 x 10 m grid for predictions
 
- *Figures.zip* - Contains all figures for the publication (including Supplemental):
  - *Figure_1_Study_Location_and_Sampling_Design.png* - Study Location and Sampling Design.
  - *Figure_2_Cross_validation_results.png* - Predicted vs. observed biomass across total biomass (top), shrub biomass (bottom left), and herbaceous biomass (bottom right) using hold-out cross-validation. 
  - *Figure_3_Biomass_predictions.png* - Biomass predictions and coefficient of variation of the predictions across the entire study area. 
  - *Figure_S2.1_Variable_Correlation.png* - Correlation plot between all variables.
  - *Figure_S2.2_Variable_Dissimilarity.png* - Dissimilarity tree of all variables.
  - *Figure_S3_Field_measurements.png* - Field measurements of herbivore-accessible AGB per biomass type and habitat class.
  - *Figure_S6_Biomass_predictions_forb_grass.png* - Biomass prediction maps  and Coefficient of variation (CV) across the entire study area for forbs and grasses 
  - *Figure_S7_Prediction_CV_Biomass.png* - Predictions and Coefficient of Variation (CV) per habitat type for Total, Shrub and Herbaceous biomass.
  - *Figure_S8_Predictions_SE_CV_all.png* - Histograms of AGB predictions, standard errors and coefficients of variation across all biomass categories.
 
GIS data will be made available in a separate repository in the near future and linked here.
          
