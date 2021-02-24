# Covid Zero charts for B117 variant

See https://www.endcoronavirus.org/ for more information on #covidzero

Creates charts of alternative scenarios for development of covid given the new variant B117.

Current charts are for Feb. 14 2021 
Current states supported are CA, AZ, TX, GA, NC, LA, MI, IN, PA, MA, FL
Other states have insufficient B117 data from Helix as of 2/23/2021

# Usage
python new_variant_estimate.py

# R Values:
  Initial R0 is given by data from Epiforecasts
  Old variants vs. new variant R assumes that new variant has 50% higher R than old variants
  #covidzero R values are always 0.6 and 0.9, calculated from status where new variant is on a downward trajectory.

# Data Sources:
- Historical rates: The covid tracking project (https://covidtracking.com/data) (using JSON API for daily rates)
- R0 projectio based on historical rates: Epiforecasts: https://epiforecasts.io/covid/posts/national/united-states/ 
  - Paper for epiforecasts method: "Estimating the time-varying reproduction number of SARS-CoV-2 using national and subnational case counts", Abbot et. al, https://wellcomeopenresearch.org/articles/5-112/v1  
- B117 variant percent: Helix: https://public.tableau.com/profile/helix6052#!/vizhome/SGTFDashboard/SGTFDashboard 
- Growth rate for new variant (~50% higher transmission): Martina Reichmuth et al 2021, "Transmission of SARS-CoV-2 variants in Switzerland", https://ispmbern.github.io/covid-19/variants/ 

Do not use data CSVs directly for other projects. Please download fresh data from these sources.
