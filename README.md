# project-4.0
This repository allows the Project 4.0 Shiny application to be run locally. In
order to do this, the user must have R, RStudio, as well as the Shiny library
installed. Run by calling the command
`shiny::runGitHub("project-4.0", "VictorKong94")`.

## Instructions
This application processes data exported from the SDS Real-Time PCR Software
developed by Applied Biosystems. To begin, select `File` > `Export`, and save
your data as a "Tab-delimited Text (*.txt)" file. The application can now take
this file as an input. If you would like to rename your samples using a qPCR
plate template file, check `Submit qPCR Plate Template` (see guidelines below
for template requirements). If you would like to sort your samples by replicates
of each treatment, check `Enter String to Sort by Replicates` and enter the
string used in sample names to indicate they are replicates (see guidelines
below for sort by replicates). Select your housekeeping gene and the
quantification algorithm you want to use. If choosing
`Relative ($\delta \delta$Ct)`, select your control treatment. Change the
selection under `Select Data Output`

#### qPCR Plate Template
This application automatically assumes three replicates per sample on a 384-well
plate, and therefore requires any plate template to have 16 rows or less as well
as 8 columns or less. Include only sample names, without headers, row names,
etc. Save templates as "Comma-separated Values(*.csv)" files.

#### Sort by Replicates

--------------------------------------------------------------------------------
