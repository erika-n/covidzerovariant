import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as ticker
import matplotlib.dates as mdates
from matplotlib import rc,rcParams
from matplotlib.dates import DateFormatter
import matplotlib.ticker as mticker
import seaborn as sns
import matplotlib.font_manager
import datetime
import urllib.request, json

import os

from matplotlib import rcParams
from labellines import labelLines
from scipy import ndimage


# Estimate projections for old, B117, and 
# combined variants per U.S. state

# Copyright 2021 Erika Nesse and NECSI.org with MIT license 
# (see LICENSE)


state_abbrevs = {
  'AK': 'Alaska',
  'AL': 'Alabama',
  'AR': 'Arkansas',
  'AS': 'American Samoa',
  'AZ': 'Arizona',
  'CA': 'California',
  'CO': 'Colorado',
  'CT': 'Connecticut',
  'DC': 'District of Columbia',
  'DE': 'Delaware',
  'FL': 'Florida',
  'GA': 'Georgia',
  'GU': 'Guam',
  'HI': 'Hawaii',
  'IA': 'Iowa',
  'ID': 'Idaho',
  'IL': 'Illinois',
  'IN': 'Indiana',
  'KS': 'Kansas',
  'KY': 'Kentucky',
  'LA': 'Louisiana',
  'MA': 'Massachusetts',
  'MD': 'Maryland',
  'ME': 'Maine',
  'MI': 'Michigan',
  'MN': 'Minnesota',
  'MO': 'Missouri',
  'MP': 'Northern Mariana Islands',
  'MS': 'Mississippi',
  'MT': 'Montana',
  'NA': 'National',
  'NC': 'North Carolina',
  'ND': 'North Dakota',
  'NE': 'Nebraska',
  'NH': 'New Hampshire',
  'NJ': 'New Jersey',
  'NM': 'New Mexico',
  'NV': 'Nevada',
  'NY': 'New York',
  'OH': 'Ohio',
  'OK': 'Oklahoma',
  'OR': 'Oregon',
  'PA': 'Pennsylvania',
  'PR': 'Puerto Rico',
  'RI': 'Rhode Island',
  'SC': 'South Carolina',
  'SD': 'South Dakota',
  'TN': 'Tennessee',
  'TX': 'Texas',
  'UT': 'Utah',
  'VA': 'Virginia',
  'VI': 'Virgin Islands',
  'VT': 'Vermont',
  'WA': 'Washington',
  'WI': 'Wisconsin',
  'WV': 'West Virginia',
  'WY': 'Wyoming'
}



# Make a projection for cases over n_days
# R: reproduction rate
# init_cases: initial number of cases
# n_days: number of days to project
# pad: prepend blanks
# generation: number of generation days that R is based on
def projectCases(R, init_cases, n_days=30, pad=0, generation=4):

  generation_time = 3.6 # generation time from https://epiforecasts.io/covid/methods.html
  days = np.arange(n_days)
  projection = init_cases*R**(days/generation_time) 
  if pad > 0:
    pad = np.full(pad, np.nan)
    projection = np.concatenate(( pad, projection))

  return projection



# Get historical data for this state
def getDataForState(data, state):

  state_data = data[data['state']==state].drop('state', axis=1)
  state_data.index = np.arange(len(state_data))
  window = 7

  state_data['average'] = state_data['positiveIncrease'].iloc[::-1].rolling(window=window, min_periods=1, center=False).mean().iloc[::-1]

  state_data.fillna(0, inplace=True)
  return state_data

# Pad out a column to have the correct number of x axis values
def padColumn(data, column, n_pad):
  pad = np.full(n_pad, np.nan) # ignore the pad values in the chart
  average = np.concatenate((pad, data[column]))[::-1]

  return average

# Set up each matplotlib axis
def setUpAxis(axis, y_cutoff, labels=False):
  # DATE FORMATTER #
  date_form = DateFormatter("%#m/%#d")


  # SETTING LABEL TO TOP & FORMATTING DATE #
  
  axis.xaxis.set_label_position('bottom')
  axis.xaxis.set_major_formatter(date_form)    


  # SETTING TICK PARAMATERS #
  axis.tick_params(axis='y',labelcolor= 'grey',labelsize = 18,width=0, length=8)
  axis.tick_params(axis='x',labelcolor= 'grey',labelsize = 18,labelrotation=0, width=1.5, length=8)
  ticks_loc = axis.get_yticks().tolist()
  axis.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
  axis.set_yticklabels(['{:,}'.format(int(x)) for x in ticks_loc])

  # SET Y LIMIT #
  axis.set_ylim((0, int(y_cutoff)))


  # SETTING GRID #
  axis.grid(axis='y')

  # LABEL AXES #
  if labels:
    axis.set_xlabel('Date', fontsize=18, color='grey', labelpad=8)
    axis.set_ylabel('Cases', fontsize=18, color='grey', labelpad=8)



        
# Matplotlib setup
def setUpPlots(n_plots):
  rc('font',weight='light')
  rc('axes', lw=0.01)
  rcParams.update({'figure.autolayout': True})

  with plt.style.context("seaborn-white"):
    fig, axis= plt.subplots(1, n_plots,figsize=(17,7) , squeeze=False)
  axes = axis.flatten()
  # SETTING WIDTH BETWEEN THE PLOTS #
  plt.subplots_adjust(wspace=0.65)

  plt.tight_layout()
  rc('font',weight='light')
  rc('axes', lw=0.01)


    # SETTING NOTE BELOW FIGURE #
  plt.figtext(0.5, -0.2, """Data sources:
The Covid Tracking Project (historical data)
Epiforecasts (projected R)
Helix (B117 percentage)
M. Reichmuth et al 2021 (~50% increase for new variant)"""
      , ha="center", fontsize=18, bbox={ "facecolor":"white", "pad":5})

  return axes


# Make a sub-chart of historical and projection data
def makeSubChart(df_historical, df_r_estimates, state, R_covid, R_variant, 
    axis, header_text, n_days_data, n_days_projection, legend=False, current=True, pct_variant=1.0, overlap=False):

    historical_start_date = df_historical['date'][0]

    if(current):
      projection_start_date = datetime.datetime(2021, 2, 15)
    else:
      projection_start_date = datetime.datetime(2021, 2, 23)

    projection_overlap = int((historical_start_date - projection_start_date).days)



    if current and overlap:
      data=df_historical[:n_days_data]
    else:
      data=df_historical[projection_overlap:n_days_data + projection_overlap]
      projection_overlap = 0
    

    data = data.reset_index(drop=True)



    # get dates for x axis
    dates = pd.date_range(start=data['date'][len(data) -1 ],periods=(len(data) + n_days_projection - projection_overlap), freq='D')


    # PROJECTIONS #

    # 3.6 days generation time (from Epiforecasts model https://epiforecasts.io/covid/methods)
    generation = 3.6 

    # initial cases per day for original strain

    init_covid_cases = data['average'].iat[projection_overlap]

    # initial cases per day for B117 variant

    init_variant_cases = pct_variant*init_covid_cases
    init_covid_cases -= init_variant_cases



    # get projections
    projection_covid = projectCases(R_covid, init_covid_cases, n_days_projection  + 1,generation=generation, pad=len(data) - projection_overlap - 1)
    projection_variant = projectCases(R_variant, init_variant_cases, n_days_projection  + 1, generation=generation, pad=len(data) - projection_overlap - 1)
    projection_total = projection_covid + projection_variant

    # cutoff where greater than y_cutoff
    y_cutoff = data['average'].max()
    x_cutoff = np.where(projection_total > y_cutoff)

    if len(x_cutoff) > 0 and x_cutoff[0].size > 0:
      x_cutoff = x_cutoff[0][0] 

      if x_cutoff > 1:
        projection_total[x_cutoff + 1:] = np.nan

        projection_variant[x_cutoff + 1:] = np.nan
        #projection_covid[x_cutoff:] = np.nan
    else:
      x_cutoff = -1

    # HISTORICAL DATA
    average = padColumn(data, 'average', n_days_projection - projection_overlap)


    # HISTORICAL: LINE PLOT #
    axis.plot(dates, average, color = 'black', lw=3, label='7 day avg')


    # PROJECTION PLOTS #

    old_var_str = r"Old: R=%.2f"% (R_covid)
    new_var_str = r"B117: R=%.2f"%(R_variant)
    axis.plot(dates, projection_covid, '--', color='blue', lw=1, label=old_var_str)
    axis.plot(dates, projection_variant, '--', color='red', lw=2, label=new_var_str)
    axis.plot(dates, projection_total, color='blue', lw=1, label="Total")




    # SUBPLOT TITLE#
    axis.set_title(header_text,fontsize=18, y=1.03, alpha=0.9, fontweight='bold')

    axis_cutoff = 1.55*y_cutoff
    
    axis_cutoff = int(1500*axis_cutoff)/1500 -1
    setUpAxis(axis, axis_cutoff, current)
    
    if legend:
      axis.legend(fontsize=16, loc='best')        



# Make a two paneled chart for the given state and values
def makeChart(state, n_days_projection, n_days_data, data_folder, update_data=False, legend=False, fig_folder='fig', overlap=False):

  state_full_name = str(state_abbrevs[state])
  


  if not os.path.exists(data_folder):
    os.makedirs(data_folder)
  if not os.path.exists(fig_folder):
    os.makedirs(fig_folder)

  
  # DOWNLOAD HISTORICAL DATA #
  # automatic update of covidtracking.com data
  filepath_historical = data_folder + '/covid_daily.csv'  
  if update_data or not os.path.exists(filepath_historical):
    print('downloading daily data')
    with urllib.request.urlopen("https://covidtracking.com/api/v1/states/daily.json") as url:
        historical = json.loads(url.read().decode())
    df_historical = pd.DataFrame(historical)
    df_historical.to_csv(filepath_historical)
  else:
    df_historical = pd.read_csv(filepath_historical, index_col=0)

  # SET UP DATA
  # historical
  df_historical['date'] = pd.to_datetime(df_historical.date, format='%Y%m%d')
  df_historical = getDataForState(df_historical, state)
  # last date of historical data
  date = df_historical['date'][0]
  date_str = date.strftime("%#m-%#d-%Y")

  # r estimates
  df_r_estimates = pd.read_csv(data_folder + '/Covid-19 National and Subnational estimates for the United States of America.csv', index_col=0)
  df_r_estimates.rename(columns=lambda x: x.strip(), inplace=True)
  df_r_estimates['R'] = df_r_estimates['Effective reproduction no.'].str.split(" ").str.get(0)

  # data on B117 rates
  df_b117 = pd.read_csv(data_folder + '/Helix_B117.csv', index_col=0)

  # Estimated R0
  R0 = float(df_r_estimates['R'][state_full_name])

  # calculate pct_variant = pct cases SGTF*pct of SGTF cases B117 (5 day moving avg)
  pct_SGTF = df_b117['pct_SGTF'][state]/100.0 
  pct_SGTF_B117 = df_b117['pct_SGTF_B117'][state]/100.0
  pct_variant = pct_SGTF*pct_SGTF_B117

  # calculate R from estimated R0 for all cases
  R_covid = R0/(1.5*pct_variant + (1-pct_variant))
  R_variant = 1.5*R_covid # 50% increas in reproduction rate

  # lockdown R calculated from requirement to reduce total cases (variant R = 0.9)
  R_covid_lockdown = 0.6 
  R_variant_lockdown = 1.5*R_covid_lockdown


  # CHART SETUP
  title_current = r'%s projections for new variants'% (state_full_name.upper())
  title_covidzero = 'Projections with #CovidZero policies in place'

  axes = setUpPlots(2)
  makeSubChart(df_historical, df_r_estimates, state, R_covid, R_variant, axes[0], title_current, n_days_data, n_days_projection, legend=legend, current=True, pct_variant=pct_variant, overlap=overlap)
  makeSubChart(df_historical, df_r_estimates, state, R_covid_lockdown, R_variant_lockdown, axes[1], title_covidzero, n_days_data, n_days_projection, legend=legend, current=False, pct_variant=pct_variant, overlap=overlap)

  filename = r'%s\projection_%s_%s.png'% (fig_folder, state, date_str)
  plt.savefig(filename, dpi=150, bbox_inches='tight', pad_inches=1)
  
  print('Saved ', filename)
  #plt.show()


  
if __name__ == '__main__':
  # SETTINGS #
  n_days_projection = 90
  n_days_data = 60
  data_folder = 'data'
  overlap = True
  update_data = False
  legend = True
  state = 'ALL'

  if overlap:
    fig_folder='figs_overlap'
  else:
    fig_folder = 'figs'

  df_b117 = pd.read_csv(data_folder + '/Helix_B117.csv')
  states = df_b117['state']
  if state == 'ALL':
    for state in states:
      makeChart(state, n_days_projection, n_days_data, data_folder, update_data, legend=legend, fig_folder=fig_folder, overlap=overlap)
  else:
    makeChart(state, n_days_projection, n_days_data, data_folder, update_data, legend=legend, fig_folder=fig_folder, overlap=overlap)

# Data Sources:
#   Historical rates: The covid tracking project (https://covidtracking.com/data) (using JSON API for daily rates)
# 	R0 projectio based on historical rates: Epiforecasts: https://epiforecasts.io/covid/posts/national/united-states/ 
#     Paper for epiforecasts method: "Estimating the time-varying reproduction number of SARS-CoV-2 using national and subnational case counts", Abbot et. al, https://wellcomeopenresearch.org/articles/5-112/v1  
#   B117 variant percent: Helix: https://public.tableau.com/profile/helix6052#!/vizhome/SGTFDashboard/SGTFDashboard 
#   Growth rate for new variant (~50% higher transmission): Martina Reichmuth et al 2021, "Transmission of SARS-CoV-2 variants in Switzerland", https://ispmbern.github.io/covid-19/variants/

