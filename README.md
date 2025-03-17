analyses codes and computational models relevant to the bart implementation on jatos, summary statistics, and computational modeling



data processing: 

  filter_participants: deprecated; has some original code for exlusion criteria 
  
  jsonprocessor: turn the bart and questionnaire data into a csv 
  
  two_csv: compares common rows in two csv files 



summary scripts: 

  bart_summary_stats.rmd: plots for trial and block comparisons etc 
  
  questionnaire_summary: some plots for questionnaire data 
  
  questionnaire_pca: code to extract pc1 and pc2 and plot questionnaire subscales alongside 



stl folder: 

  stl_model_fit.rmd: take csv file and fit to stl model 
  
  stl_hbm_color.stan: stl model with changes made to incorporate different color balloons 
  
  stl_hbm_questionnaires: same as the color file but with regressors for spq, pdi, and caps 
  
  stl_plots: plotting functions from stl_model_fit.rmd 
  


ewmv folder: 
  coming soon - ewmv model
