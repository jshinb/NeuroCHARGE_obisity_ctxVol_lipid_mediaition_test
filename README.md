# To-fix:
Add - fasting duration in the base model for metabolites: 
# To-think
1. Make presentation files
2. Covariates for metabolites
3. Models: check!!
4. Which sets of analyses should we run? Potentially, it can be computationally expensive (e.g., for Rhineland) - at the moment analyses are performed sequentially

2026-11-05
* The previous versions - I used functions to run each set of analyses; however, it was hard to debug the code when errors occurred, so I have decided to create a separate script for each analysis and source it.

Models:
Lipid ~ Obesity + base covariates + (additional covariates)
Brain ~ Obesity + Lipid + base covariates + (additional covariates)

Metabolite-specific covariates
* Age at blood sampling 
* Fasting duration: hours between the last meal and blood draw for metabolomic analysis (if available)
* Cohort specific covariates: e.g.,batch, etc. 

Additional covariates:
* Statin treatment (yes/no)
* Other lipid lowering medication (yes/no)
* Estimated glomerular filtration rate, eGFR (mL/min per 1.73 m2)
https://www.niddk.nih.gov/health-information/communication-programs/nkdep/laboratory-evaluation/glomerular-filtration-rate-calculators 

Base covariates: 
* All considered in stage 1 + metabolomic-specific covariates

Additional covariates: 
* All considered in stage 1 + metabolomic-specific covariates

<img width="842" height="223" alt="image" src="https://github.com/user-attachments/assets/dc3972da-5d5c-4fc0-bac9-a7756fe7dc8e" />
<img width="846" height="133" alt="image" src="https://github.com/user-attachments/assets/f090393c-1a6e-4119-8923-26ec96e83932" />

