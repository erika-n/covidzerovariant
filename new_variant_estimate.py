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

# Estimate projections for old, new, and combined variants per U.S. state


states = {
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

  state_data['average'] = state_data['positiveIncrease'].rolling(window=window, min_periods=1, center=False).mean()

  state_data.fillna(0, inplace=True)
  return state_data

# Pad out a column to have the correct number of x axis values
def padColumn(data, column, n_pad):
  pad = np.full(n_pad, np.nan) # ignore the pad values in the chart
  average = np.concatenate((pad, data[column]))[::-1]

  return average

# Set up each matplotlib axis
def setUpAxis(axis, y_cutoff):
  # DATE FORMATTER #
  date_form = DateFormatter("%b-%y")


  # SETTING LABEL TO TOP & FORMATTING DATE #
  
  axis.xaxis.set_label_position('top')
  axis.xaxis.set_major_formatter(date_form)    


  # SETTING TICK PARAMATERS #
  axis.tick_params(axis='y',labelcolor= 'grey',labelsize = 12,width=0, length=8)
  axis.tick_params(axis='x',labelcolor= 'grey',labelsize = 12,labelrotation=0, width=1.5, length=8)

  # SETTING TICK LIMITS #
  #axis.set_ylim(bottom=0)

  
  # FOR PUTTING COMMAS IN NUMBERS #
  #axis.set_yticklabels(['{:,}'.format(int(x)) for x in axis.get_yticks().tolist()])      

  ticks_loc = axis.get_yticks().tolist()
  axis.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
  axis.set_yticklabels(['{:,}'.format(int(x)) for x in ticks_loc])

  axis.set_ylim((0, int(y_cutoff)))


  axis.legend()        


  # SETTING GRID #
  axis.grid(axis='y')



        
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
  plt.figtext(0.5, -0.15, """Data sources:
The Covid Tracking Project (historical data)
Epiforecasts (projected R for old variant)
Helix (B117 percentage)
M. Reichmuth et al 2021 (~50% infection rate increase for new variant)"""
      , ha="center", fontsize=12, bbox={ "facecolor":"white", "pad":5})

  return axes


# Make a sub-chart of historical and projection data
def makeSubChart(df_historical, df_r_estimates, df_emerging_variants, state, R_covid, R_variant, 
    axis, header_text, n_days_data, n_days_projection):

    data=df_historical[:n_days_data]


    dates = pd.date_range(start=data['date'][len(data) -1 ],periods=(len(data) + n_days_projection), freq='D')
    date = data['date'][0]
    date = date.strftime("%#m-%#d-%Y")

    # PROJECTIONS #

    generation = 3.6 # 3.6 days generation time (from Epiforecasts model https://epiforecasts.io/covid/methods)

    # initial cases per day for original strain
    init_covid_cases = data['average'][0]

    y_cutoff = data['average'].max()

    # from https://public.tableau.com/profile/helix6052#!/vizhome/SGTFDashboard/SGTFDashboard
    # (@helix)
    # Collcted on 2/20/2021
    # 5 day moving average for last day counted (February 13 or 14)
    percent_cases_variant = {
      'FL': 0.1604,
      'MA': 0.0455,
      'CA': 0.06706
    }


    # initial cases per day for B117 variant
    if state in percent_cases_variant:
      init_variant_cases = percent_cases_variant[state]*init_covid_cases
    else:
      print('This state does not have variant data')
      exit()
    # overlap one day between projection and original data
    projection_covid = projectCases(R_covid, init_covid_cases, n_days_projection + 1,generation=generation, pad=len(data) - 1)
    projection_variant = projectCases(R_variant, init_variant_cases, n_days_projection + 1, generation=generation, pad=len(data) - 1 )
    projection_total = projection_covid + projection_variant
    projection_diff = np.abs(projection_covid - projection_variant)[-n_days_projection:]

    # crossover point where original and variant strains are equal (if exists)
    use_crossover_point = True
    crossover_point = np.argmin(projection_diff)
    crossover_point += len(data)
    crossover_date = dates[crossover_point]
    crossover_cases = projection_total[crossover_point]


  

    if crossover_point >= len(projection_total) - 1 or projection_total[crossover_point + 1] < projection_total[crossover_point]:
      use_crossover_point = False



    # cutoff greater than y_cutoff
    x_cutoff = np.where(projection_total > y_cutoff)

    if len(x_cutoff) > 0 and x_cutoff[0].size > 0:
      x_cutoff = x_cutoff[0][0] - 1

      if x_cutoff > 1:
        projection_total[x_cutoff:] = np.nan

        projection_variant[x_cutoff:] = np.nan
        #projection_covid[x_cutoff:] = np.nan
    else:
      x_cutoff = -1


    # total cases at end of projection window

    if x_cutoff > 1:
      end_cases = projection_total[x_cutoff - 1]
    else:
      end_cases = projection_total[-1]



    # HISTORICAL DATA
    average = padColumn(data, 'average', n_days_projection)
    positiveIncrease = padColumn(data, 'positiveIncrease', n_days_projection)



    # PROJECTION PLOTS #
    axis.plot(dates, projection_covid, '--', color='blue', lw=1, label=r"Projected new cases: Old variants, R=%.2f"% (R_covid))
    axis.plot(dates, projection_variant, '--', color='red', lw=2, label=r"Projected new cases: Variant B.1.1.7, R=%.2f (estimated 50%% > baseline)"%(R_variant))
    axis.plot(dates, projection_total, color='blue', lw=1, label="Projected new cases: Sum")
 


    # HISTORICAL: LINE PLOT #
    axis.plot(dates, average, color = 'black', lw=3, label='New cases (historical): 7 day average')


    # HISTORICAL: SCATTER PLOT #
    axis.scatter(x=dates, y=positiveIncrease, color ='grey', label="New cases (historical)")

    # PLOT LINE TO KEEP X AXIS VALUES #
    axis.plot(dates, np.zeros(dates.size), color='black')

    # VERTICAL LINES
    if use_crossover_point:
      crossover_line_height = 0.5*positiveIncrease[:-n_days_projection - 1].max()

      crossover_line_str = r"New variant dominant on %s"%(crossover_date.strftime("%m-%d-%Y"))
      axis.vlines(crossover_date, 0, crossover_line_height, color='pink', lw=3, label=crossover_line_str)

    # if x_cutoff > 0:
    #   cutoff_line_height = 0.7*positiveIncrease[:-n_days_projection - 1].max()
    #   cutoff_date = dates[x_cutoff - 1]
    #   crossover_line_str = cutoff_date.strftime("%m-%d-%Y")
    #   axis.vlines(cutoff_date, 0, cutoff_line_height, color='grey', lw=1, label=crossover_line_str)



    # ANNOTATIONS

    # round numbers
    init_covid_cases_str = str(int(10*round(init_covid_cases/10))) 
    init_variant_cases_str = str(int(10*round(init_variant_cases/10)))
    crossover_cases_str = str(int(10*round(crossover_cases/10)))
    end_cases_str = str(int(10*round(end_cases/10)))
    end_cases_str += ': ' + dates[x_cutoff - 1].strftime("%m-%d-%Y")
    # annotate numbers for initial conditions
    axis.annotate(init_covid_cases_str, 
      xy=(dates[n_days_data - 1], init_covid_cases + 100),
      xytext=(dates[n_days_data + 2], init_covid_cases + 1500),
      arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color='blue'),
      color='blue')
    axis.annotate(init_variant_cases_str, 
      xy=(dates[n_days_data -1 ], init_variant_cases + 100),
      xytext=(dates[n_days_data + 2], init_variant_cases + 1500),
      arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color='red'),
      color='red')
    # annotate crossover point
    if use_crossover_point:
      axis.annotate(crossover_cases_str, 
        xy=(crossover_date, crossover_cases + 100),
        xytext=(crossover_date, crossover_cases + 1500),
        arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color='blue'),
        color='blue')    
    # annotate end value and date
    axis.annotate(end_cases_str, 
      xy=(dates[x_cutoff - 1], end_cases + 0),
      xytext=(dates[x_cutoff - 6], end_cases + 1500),
      arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color='blue'),
      color='blue')  



    


    # SUBPLOT TITLE#
    axis.set_title(header_text,fontsize=18, y=1.03, alpha=0.5)

    setUpAxis(axis, y_cutoff + 0.7*y_cutoff)
 


# Make a two paneled chart for the given state and values
def makeChart(state, n_days_projection, n_days_data, data_folder, update_data=False):

  state_full_name = str(states[state])
  


  if not os.path.exists(data_folder):
    os.makedirs(data_folder)
  if not os.path.exists('figs'):
    os.makedirs('figs')

  
  # DOWNLOAD DATA #
  # automatic update of covidtracking.com data
  filepath_historical = data_folder + '/covid_daily.csv'  
  if update_data or not os.path.exists(filepath_historical):
    print('downloading daily data')
    with urllib.request.urlopen("https://covidtracking.com/api/v1/states/daily.json") as url:
        historical = json.loads(url.read().decode())
    df_historical = pd.DataFrame(historical)
    df_historical.to_csv(filepath_historical)
  else:
    df_historical = pd.read_csv(filepath_historical)

  # SET UP DATA
  # histoical
  df_historical['date'] = pd.to_datetime(df_historical.date, format='%Y%m%d')
  df_historical = getDataForState(df_historical, state)


  # r estimates
  df_r_estimates = pd.read_csv(data_folder + '/Covid-19 National and Subnational estimates for the United States of America.csv', index_col=0)
  df_r_estimates.rename(columns=lambda x: x.strip(), inplace=True)
  df_r_estimates['R'] = df_r_estimates['Effective reproduction no.'].str.split(" ").str.get(0)

  # emerging variants
  df_emerging_variants = pd.read_csv(data_folder + '/Emerging Variant Cases in the United States.csv', index_col=0)
  df_emerging_variants = df_emerging_variants[df_emerging_variants['filter'] == 'Variant B.1.1.7']
  df_emerging_variants.rename(columns=lambda x: x.strip(), inplace=True)


  axes = setUpPlots(2)

  
  R_covid = float(df_r_estimates['R'][state_full_name])
  R_variant = 1.5*R_covid # 50% increas in reproduction rate

  R_covid_lockdown = 0.6 # example value
  R_variant_lockdown = 1.5*R_covid_lockdown


  date = df_historical['date'][0]
  date = date.strftime("%#m-%#d-%Y")

  title_current = r'SAMPLE Current %d day projections for %s on %s'% (n_days_projection, state_full_name, date)
  title_covidzero = 'SAMPLE Solution: #COVIDZero'
  makeSubChart(df_historical, df_r_estimates, df_emerging_variants, state, R_covid, R_variant, axes[0], title_current, n_days_data, n_days_projection)
  makeSubChart(df_historical, df_r_estimates, df_emerging_variants, state, R_covid_lockdown, R_variant_lockdown, axes[1], title_covidzero, n_days_data, n_days_projection)

  filename = r'figs\projection_%s_%s.png'% (state, date)

  plt.savefig(filename, dpi=150, bbox_inches='tight', pad_inches=1)
  print('Saved ', filename)
  plt.show()


  
if __name__ == '__main__':
  # SETTINGS #
  state = 'MA'
  n_days_projection = 60
  n_days_data = 60
  data_folder = 'data_02_15'
  update_data = False

  makeChart(state, n_days_projection, n_days_data, data_folder, update_data)


# Data Sources:
#   The covid tracking project (https://covidtracking.com/data) for historical rates (using JSON API for daily rates, see below)
# 	Epiforecasts: https://epiforecasts.io/covid/posts/national/united-states/ for R for projection for old variant
#   Helix: https://public.tableau.com/profile/helix6052#!/vizhome/SGTFDashboard/SGTFDashboard for new variant percent
# Growth rate for new variant: Martina Reichmuth et al 2021, "Transmission of SARS-CoV-2 variants in Switzerland", https://ispmbern.github.io/covid-19/variants/ (~50% higher transmission)

