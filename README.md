# covidzerovariant
Covid Zero charts for B117 variant

Creates charts of alternative scenarios for development of covid given the new variant B117.

NEEDS REVIEW BEFORE PUBLISHING CHARTS

Usage: python new_variant_estimate.py

Change state and days history/projection in the code.

Current data is for Feb. 14 2021

Current states supported are MA, FL, and CA

Data Sources:
  The covid tracking project (https://covidtracking.com/data) for historical rates (using JSON API for daily rates)
	Epiforecasts: https://epiforecasts.io/covid/posts/national/united-states/ for R for projection for old variant
  Helix: https://public.tableau.com/profile/helix6052#!/vizhome/SGTFDashboard/SGTFDashboard for new variant percent
  Growth rate for new variant: Martina Reichmuth et al 2021, "Transmission of SARS-CoV-2 variants in Switzerland", https://ispmbern.github.io/covid-19/variants/ (~50% higher transmission)


