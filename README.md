# project-5.0
This repository allows the Project 5.0 Shiny application to be run locally. In
order to do this, the user must have both R and the Shiny library installed. To
run this application, call the following command from your R console
`shiny::runGitHub("project-5.0", "VictorKong94")`.

## Instructions
This application processes data exported from the SDS Real-Time PCR Software
developed by Applied Biosystems. Begin by exporting your qPCR results as a
"Tab-delimited Text (*.txt)" file. The application can now take this file as an
input. If you would like to rename your samples using a qPCR plate template
file, check `Submit qPCR Plate Template` (see guidelines below for template
requirements). If you would like to analyze only a subset of your plate for the
time being, check `Enter String to Subset Data` and enter the string used in
sample names to indicate they belong to a subset. Select your housekeeping
gene(s) and the quantification algorithm you want to use. If choosing `Relative
(ΔΔCt)`, select your control treatment(s). Change the selection under `Select
Data Output` to view your desired data output.

#### qPCR Plate Template
This application automatically assumes three replicates per sample on a 384-well
plate, and therefore requires any plate template to have 16 rows or less as well
as 8 columns or less. Include only sample names, without headers, row names,
etc. Save templates as "Comma-separated Values (*.csv)" files.

#### Special Case of Relative Quantification With a Single Housekeeping Gene
As requested by my lab group, I have modified the column for a lone housekeeping
gene to contain a pseudo-ΔΔCt and a pseudo-fold change for the purpose of
estimating error between samples. Note that these values are defined differently
from others in the table.
