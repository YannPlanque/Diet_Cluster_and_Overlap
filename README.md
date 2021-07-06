# Identifying diet clusters and interspecific overlap
### **CASE STUDY**: Harbour and grey seal diets from scat content analysis (in the *Baie de Somme*, France)
• **OBJECTIVES** : Analyse diet content data (here scat contents) to compare the diets of predator species (here harbour and grey seals). We aimed here at **quantifying the diet** (with the associated confidence interval) and **identifying different diet clusters** based on the composition of scats in functional prey (identified from the analysis of diagnostic hard parts). The distribution of samples in clusters is expected to inform on the typical diet contents found. We also aimed at using the **Pianka index** to measure the **diet overlap** between different species.

• **METHOD** : 
<br>(1) Quantify the **diet composition in % by mass of prey**, and measure the associated uncertainty (CI) using a bootstrap procedure (considering the number of samples). 
<br>(2) Identify **diet clusters** based on the composition of each sample (here scat) in % by mass of prey (here functionnal prey). Here we used an **agglomerative hierarchical cluster analysis** (using a Euclidian distance procedure and employing the Ward.D2 algorithm). Used function: 'clustCoDa' in *robCompositions* package. We also presented another way to cluster the data using this package, that was not chosen here. **Year-season variations** in the diet were also studied here, by looking at the % of samples in each cluster during year-season periods.
<br>(3) Identify the **diet overlap** by: (A) comparing the percentage of samples of different species in each diet clusters (several plots are proposed), and (B) calculating the Pianka index (Pianka 1974) between two diet data sets. 

• **Required data** : 
A data set summarising the (reconstructed and/or measured) **mass (e.g. g, kg)** of different **prey species/types (*in columns*)** found in **several samples (*in rows*)**. An example of data that could be used in the script: *Columns* = prey species/types (e.g. benthic flatfish, pelagic fish); *Rows* = samples used to analyse the diet content (e.g. scats, stomachs). Columns with the sample ID and the consumer/predator species are welcoming. Column with the sampling date can also be useful to study potential year-season variations.

• **LANGUAGE** : R (version 4.0.2).

• **CASE STUDY** : The harbour seals and grey seal diets were analysed based on the hard prey remains contained in scat samples collected in the *Baie de Somme* (France) from 2002 to 2019. Here the used data set included 193 harbour seal and 77 grey seal scats.

## Script developped as part of: 
Planque Y, Spitz J, Authier M, Guillou G, Vincent C, Caurant F. Trophic niche overlap between sympatric harbour seals (Phoca vitulina) and grey seals (Halichoerus grypus) at the southern limit of their European range (Eastern English Channel). Ecol Evol. 2021;00:1– 22. https://doi.org/10.1002/ece3.7739

## Data used in this script are freely available on SEANOE:
Planque Yann, Vincent Cécile, Caurant Florence, Spitz Jérôme (2020). Content of harbour seal (Phoca vitulina) and grey seal (Halichoerus grypus) scats in prey classified by functional groups (samples collected in the *Baie de Somme*, France, from 2002 to 2019). SEANOE. https://doi.org/10.17882/76780

Before running the script, please download the file repository in the ZIP file and place it on your desktop. Place the data previously dowloaded in the subfolder "Input".

First publication on GitHub : 2020-11-05

Last update : 2021-04-21 (Version 1.3)

## Author : Yann Planque<sup>(1)</sup>
 Affiliation :
    (1) Centre d'Etudes Biologiques de Chizé (CEBC, UMR 7372 CNRS - La Rochelle Université), La Rochelle, France

## Contact : yann.planque@univ-lr.fr ; yann.planque@hotmail.fr
