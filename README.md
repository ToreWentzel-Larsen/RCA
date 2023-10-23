# RCA
Response congruence assessment

This is about response congruence assessment. The aim is to make an alternative to the field of response surface analysis (RSA, references will be included). Some main points of deviation

Congruence is not a dichotomy, in particular there is no p-value for congruence or not. Congruence should be assessed more than analysed.

Graphical tools for congruence assessment are emphasized.
The assessment of congruence should not be limited to specific model types, the assessment should be possible together with any regression model.

In particular, polynomial regression models have well-known limitations (references will be included). Quadratic models have limitations also discussed in the RSA literature. Cubic regression models also have limitations. The assessment procedures may, however, be used together with quadratic and cubic polynomial regression models if desired. 

Models with restricted cubic splines are recommended, and the code development will start with such models.

The development of RCA will start with models for data with independent observations, without cluster structures. Possible extension to clustered data may follow, but may not be easy or possible.

Instructions for downloading the (so far) main example data file mba_csv.csv, taken from supplementary material to the article Cubic Response Surface Analysis: Investigating Asymmetric and Level-Dependent Congruence Effects With Third-Order Polynomial Models, by Humberg et al. Also see the article The study of congruence in organizational behavior research: Critique and a proposed alternative) by 
Jeffrey R. Edwards, 1994. Download files for this article in http://public.kenan-flagler.unc.edu/faculty/edwardsj/downloads.htm , and save the Excel file mba as a csv file mba_csv.csv with semicolon as deliminator, in the same folder as the R script script rca.r made available in this repo. The script is at present preliminary and under work.
