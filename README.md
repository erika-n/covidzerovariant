B117 variant charts
Covid Zero charts for B117 variant

See https://www.endcoronavirus.org/ for more information on #covidzero

Creates charts of alternative scenarios for development of covid given the new variant B117.

Usage: python new_variant_estimate.py

Current charts are for Feb. 14 2021 

Current states supported are CA, AZ, TX, GA, NC, LA, MI, IN, PA, MA, FL
Other states have insufficient B117 data from Helix as of 2/23/2021

R Values:
  Initial R0 is given by data from Epiforecasts
  Old variants vs. new variant R assumes that new variant has 50% higher R than old variants
  #covidzero R values are always 0.6 and 0.9, calculated from status where new variant is on a downward trajectory.

Data Sources:
  The covid tracking project (https://covidtracking.com/data) for historical rates (using JSON API for daily rates)
	Epiforecasts: https://epiforecasts.io/covid/posts/national/united-states/ for R for projection for old variant
  Helix: https://public.tableau.com/profile/helix6052#!/vizhome/SGTFDashboard/SGTFDashboard for new variant percent
  Growth rate for new variant: Martina Reichmuth et al 2021, "Transmission of SARS-CoV-2 variants in Switzerland", https://ispmbern.github.io/covid-19/variants/ (~50% higher transmission)


