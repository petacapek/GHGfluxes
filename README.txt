Important analysis-specific information:
1. GHG flux data files in csv are in project. Nevertheless, you will need your excel file "Top_Spectra_selected_and_variables130_JR.xlsx". Copy it to the folder please so the R script runs without issues. The .gitignore file is set to ignore any xlsx files.
2. "GHGflux_analysis.R" file contains the script to data analysis
3. I normalized CO2 flux only because there is no consistent temperature and moisture effect on CH4 and N2O fluxes that would bias FTIR predictions


Jessica scripts
1. The script 1.GHGflux_Petr_JR_V2.R is the V2 of your script, I made coments in all the included parts. This script generated the output output_R10_canopy. 
2.Script 2.Fraction_from_HR_plots.R is the fraction of hetero respiration from the plots were that was measured (CC=2,3; CO=31,42; CCF=52). Then it multiplies that fraction by each total respiration of the control plots normalized from output_R10_canopy/DATA_spectra_R10_canopy.csv
