# BMI-project master's dissertation thesis Jef Moerman

For generating all analysis results only one file must be run. The file `calculateAll.R` generates the results of the statistical analyses and writes the figures in pdf format and tables in .tex format to the working directory. The purpose of this way of working is to obtain error reduction in the settings used for analysis.

Other R-files found in this repository contain functions used in `calculateAll.R`. File `convertToTexTable.R` contains a function that takes a matrix or data frame and converts it to a .tex file ready to include in the main document. 
File `dataExploration.R` contains functions generating figures and tables for the data exploration of this project. 
File `linearRegression.R` contains functions generating figures and tables for the regression analysis part of this project. 
File `metabolicAnalysis.R` contains functions generating figures and tables for the analysis on metabolic outliers part of this project. 

All the functions are called in file `calculateAll.R` for generating tracable output for reporting.
