import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as ticker
import matplotlib.dates as mdates
from matplotlib import rc, rcParams
from matplotlib.dates import DateFormatter
import matplotlib.ticker as mticker
import seaborn as sns
import matplotlib.font_manager
import matplotlib.dates as mdates
import datetime
import urllib.request
import json

import os

from matplotlib import rcParams
from labellines import labelLines
from scipy import ndimage
import argparse


# Estimate projections for old, B117, and
# combined variants per U.S. state

# Copyright 2021 Erika Nesse and NECSI.org with MIT license
# (see LICENSE)


# Data Sources:
#   Historical rates: New York Times (https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv)
# 	R0 projection based on historical rates: Epiforecasts: https://epiforecasts.io/covid/posts/national/united-states/
#     Paper for epiforecasts method: "Estimating the time-varying reproduction number of SARS-CoV-2 using national and subnational case counts", Abbot et. al, https://wellcomeopenresearch.org/articles/5-112/v1
#   B117 variant percent: Helix: https://public.tableau.com/profile/helix6052#!/vizhome/SGTFDashboard/SGTFDashboard
#   Growth rate for new variant (~50% higher transmission): Martina Reichmuth et al 2021, "Transmission of SARS-CoV-2 variants in Switzerland", https://ispmbern.github.io/covid-19/variants/


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

R_ratio = {}


# Make a projection for cases over n_days
# R: reproduction rate
# init_cases: initial number of cases
# n_days: number of days to project
# pad: prepend blanks
# generation: number of generation days that R is based on
def projectCases(R, init_cases, n_days=30, pad=0, generation_time=3.6):

  days = np.arange(n_days)
  projection = init_cases*R**(days/generation_time)
  if pad > 0:
    pad = np.full(pad, np.nan)
    projection = np.concatenate((pad, projection))

  return projection


# Get historical data for this state
def getDataForState(data, state):
    state = str(state_abbrevs[state])
    state_data = data[data['state'] == state].drop('state', axis=1)
    state_data.index = np.arange(len(state_data))
    window = 7
    state_data['case_diff'] = state_data['cases'].iloc[::-1].diff().iloc[::-1]
    state_data['average'] = state_data['case_diff'].iloc[::-
        1].rolling(window=window, min_periods=1, center=False).mean().iloc[::-1]

    state_data.fillna(0, inplace=True)
    return state_data

# Pad out a column to have the correct number of x axis values
def padColumn(data, column, n_pad):
  pad = np.full(n_pad, np.nan)  # ignore the pad values in the chart
  padded = np.concatenate((pad, data[column]))[::-1]

  return padded

# Set up each matplotlib axis


def setUpAxis(axis, y_cutoff, labels=False):
  # DATE FORMATTER #
  date_form = DateFormatter("%#m/%#d")

  # SETTING LABEL TO TOP & FORMATTING DATE #

  axis.xaxis.set_label_position('bottom')
  axis.xaxis.set_major_formatter(date_form)

  # SETTING Y TICK PARAMETERS #
  axis.tick_params(axis='y', labelcolor='grey',
                   labelsize=18, width=0, length=8)
  axis.tick_params(axis='x', labelcolor='grey', labelsize=18,
                   labelrotation=0, width=1.5, length=8)
  ticks_loc = axis.get_yticks().tolist()
  axis.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
  axis.set_yticklabels(['{:,}'.format(int(x)) for x in ticks_loc])

  # SETTNG DATE TICKS #
  months = mdates.MonthLocator()  # every month
  axis.xaxis.set_major_locator(months)
  months_fmt = mdates.DateFormatter('%#m/%#d')
  axis.xaxis.set_major_formatter(months_fmt)

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
  rc('font', weight='light')
  rc('axes', lw=0.01)
  rcParams.update({'figure.autolayout': True})

  with plt.style.context("seaborn-white"):
    fig, axis = plt.subplots(1, n_plots, figsize=(17, 7), squeeze=False)
  axes = axis.flatten()
  # SETTING WIDTH BETWEEN THE PLOTS #
  plt.subplots_adjust(wspace=0.65)

  plt.tight_layout()
  rc('font', weight='light')
  rc('axes', lw=0.01)

    # SETTING NOTE BELOW FIGURE #
  plt.figtext(0.5, -0.2, """Data sources:
The Covid Tracking Project (historical data)
Epiforecasts (projected R)
Helix (B117 percentage)
M. Reichmuth et al 2021 (~50% increase for new variant)""", ha="center", fontsize=18, bbox={"facecolor": "white", "pad": 5})

  return axes




# Make a sub-chart of historical and projection data
# def makeSubChart(df_historical, df_r_estimates, state, R_covid, R_variant,
#     axis, header_text, n_days_data, n_days_projection, current=True, pct_variant=1.0, pct_variant_update=0.0, overlap=False, generation_time=4):

def makeSubChart(state, df_historical, df_r_estimates, df_b117, axis, title='Sample Title', n_days_data=60, n_days_projection=60,
        current=True):

    historical_start_date = df_historical['date'][0]
    R0 = float(df_r_estimates['R'][state_abbrevs[state]])

    if(current):
        projection_start_date = df_b117.iloc[0]['collection_date']
    else:
        projection_start_date = historical_start_date
 
      

    projection_start_date_str = projection_start_date.strftime("%#m/%#d")

    if current:
        projection_overlap = int(
            (historical_start_date - projection_start_date).days)
    else:
        projection_overlap = 0
    

    data = df_historical[:n_days_data]
    data = data.reset_index(drop=True)

    # get dates for x axis
    dates = pd.date_range(start=data['date'][len(
        data) - 1], periods=(len(data) + n_days_projection - projection_overlap), freq='D')

    # PROJECTIONS #

    # initial cases

    init_covid_cases = data['average'].iat[projection_overlap]

    init_variant_cases = df_b117.loc[0]['ratio_B117']*init_covid_cases
    init_covid_cases = init_covid_cases - init_variant_cases


    # HISTORICAL DATA
    average = padColumn(
        data, 'average', n_days_projection - projection_overlap)


    # Generation time for R
    generation_time = 5.2

    if current:

        # VARIANT DATA POINTS
        
        variant_data_points = np.full(len(dates), np.nan)
        x_vals = []
        for index, row in df_b117.iterrows():

            x = int((row['collection_date'] - dates[0]).days)

        
            x_vals.append(x)
            y = average[x]*row['ratio_B117']
            variant_data_points[x] = y

        # VARIANT R
        ratio_B117 = df_b117['ratio_B117'][0]
        R_covid = R0/(1.5*ratio_B117 + (1-ratio_B117))

        v0 = variant_data_points[x_vals[0]]
        v1 = variant_data_points[x_vals[-1]]

        variant_days = int(
            (df_b117.iloc[-1]['collection_date'] - df_b117.iloc[0]['collection_date']).days)

        R_variant = (v1/v0)**(generation_time/variant_days)

        R_ratio[state] = R_variant/R_covid


    else:
        # R target for lockdown
        R_covid = 0.6
        R_variant = 1.5*R_covid



    # PROJECTIONS

    projection_covid = projectCases(R_covid, init_covid_cases, n_days_projection + 1,
                                    generation_time=generation_time, pad=len(data) - projection_overlap - 1)
    projection_variant = projectCases(R_variant, init_variant_cases, n_days_projection + 1,
                                      generation_time=generation_time, pad=len(data) - projection_overlap - 1)

    projection_total = projection_covid + projection_variant

    # cutoff where greater than y_cutoff
    y_cutoff = data['average'].max()
    x_cutoff = np.where(projection_total > y_cutoff)

    if len(x_cutoff) > 0 and x_cutoff[0].size > 0:
      x_cutoff = x_cutoff[0][0]

      if x_cutoff > 1:
        projection_total[x_cutoff + 1:] = np.nan

        projection_variant[x_cutoff + 1:] = np.nan
        # projection_covid[x_cutoff:] = np.nan
    else:
      x_cutoff = -1

    # PLOT DATA #
    axis.plot(dates, average, color='black', lw=3, label='7 day avg')

    old_var_str = r"Old: R=%.2f" % (R_covid)
    new_var_str = r"B117: R=%.2f" % (R_variant)
    axis.plot(dates, projection_covid, '--',
              color='blue', lw=1, label=old_var_str)
    axis.plot(dates, projection_variant, '--',
              color='red', lw=2, label=new_var_str)
    axis.plot(dates, projection_total, color='blue', lw=1, label="Total")

    if current:
      axis.scatter(dates, variant_data_points, color='green',
                   label='B117 data points')

    # VERTICAL LINE FOR PROJECTION START #

    # if current:
    axis.vlines(dates[n_days_data - projection_overlap],
                0, 0.8*y_cutoff, color='grey')
    axis.text(dates[n_days_data - projection_overlap + 3], 0.6*y_cutoff,
              "Projection\nstart (" + projection_start_date_str + ")", fontsize=16, color='grey')

    # SUBPLOT TITLE#
    axis.set_title(title, fontsize=18, y=1.03,
                   alpha=0.9, fontweight='bold')

    axis_cutoff = 1.65*y_cutoff

    axis_cutoff = int(1500*axis_cutoff)/1500 - 1
    setUpAxis(axis, axis_cutoff, current)

 
    axis.legend(fontsize=16, loc='best')


def makeChartsForState(state, df_historical, df_r_estimates, df_b117, fig_folder='new_figs', n_days_data=90, n_days_projection=90):
    state_full_name = str(state_abbrevs[state])
    date = df_historical['date'][0]
    date_str = date.strftime("%#m-%#d-%Y")

    title_variant = r'%s projections for new variants' % (
        state_full_name.upper())
    title_covidzero = 'Projections with #CovidZero policies in place'
    
    axes = setUpPlots(2)

    makeSubChart(state, df_historical, df_r_estimates, df_b117, axes[0], title=title_variant, n_days_data=n_days_data, n_days_projection=n_days_projection,
        current=True)
    makeSubChart(state, df_historical, df_r_estimates, df_b117, axes[1], title=title_covidzero, n_days_data=n_days_data, n_days_projection=n_days_projection,
        current=False)
    

    filename = r'%s\projection_%s_%s.png' % (fig_folder, state, date_str)
    plt.savefig(filename, dpi=150, bbox_inches='tight', pad_inches=1)

    print('Saved ', filename)


# Make a two paneled chart for the given state and values
def makeAllCharts(states, data_folder, fig_folder):

    if not os.path.exists(data_folder):
        os.makedirs(data_folder)
    if not os.path.exists(fig_folder):
        os.makedirs(fig_folder)

    filepath_historical = data_folder + '/covid_daily.csv'
    df_historical = pd.read_csv(filepath_historical, index_col=0)

    # SET UP DATA
    # historical
    df_historical['date'] = pd.to_datetime(df_historical.date, format='%Y-%m-%d')

    # r estimates
    df_r_estimates = pd.read_csv(data_folder + '/R_estimates.csv', index_col=0)
    # strip whitespace from column headers
    df_r_estimates.rename(columns=lambda x: x.strip(), inplace=True)
    df_r_estimates['R'] = df_r_estimates['Effective reproduction no.'].str.split(
        " ").str.get(0)

    # data on B117 rates
    df_b117 = pd.read_csv(data_folder + '/Helix_B117.csv')
    df_b117.rename(columns=lambda x: x.strip(), inplace=True)
    df_b117['collection_date'] = pd.to_datetime(
        df_b117.collection_date, format='%m/%d/%Y')

    # Get total estimated ratio of B117 cases
    df_b117['ratio_SGTF'] = df_b117['pct_SGTF']/100.0
    df_b117['ratio_SGTF_B117'] = df_b117['pct_SGTF_B117']/100.0
    df_b117['ratio_B117'] = df_b117['ratio_SGTF']*df_b117['ratio_SGTF_B117']

    for state in states:
        df_historical_state = getDataForState(df_historical, state)
        df_b117_state = df_b117[df_b117['state'] == state].copy()
        df_b117_state = df_b117_state.set_index('date_code', drop=True)
        makeChartsForState(state, df_historical_state, df_r_estimates, df_b117_state, fig_folder=args.fig_folder, n_days_data=args.n_days_data, n_days_projection=args.n_days_projection)
    

if __name__ == '__main__':
  # SETTINGS #

    parser = argparse.ArgumentParser()
    parser.add_argument("--n_days_projection", type=int,
                        default="60", help="number of days to project")
    parser.add_argument("--n_days_data", type=int, default=90,
                        help="number of days of historical data")
    parser.add_argument("--update_data", default=False,
                        type=bool, help="update historical data")
    parser.add_argument("--state", default='ALL',
                        type=str, help="specify a state")
    parser.add_argument("--fig_folder", default='new_figs',
                        type=str, help="specify a state")
    parser.add_argument("--data_folder", default='data',
                        type=str, help="specify a state")

    args = parser.parse_args()


  # DOWNLOAD HISTORICAL DATA #
  # automatic update of covidtracking.com data
    filepath_historical = args.data_folder + '/covid_daily.csv'
    if args.update_data or not os.path.exists(filepath_historical):
        print('downloading daily data')
 
        df = pd.read_csv('https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv')
        df = df[::-1].reset_index(drop=True)
        df.to_csv(filepath_historical)  


  # get states
    df_b117 = pd.read_csv(args.data_folder + '/Helix_B117.csv')
    states = df_b117['state'].tolist()
    states = list(dict.fromkeys(states))
    if args.state != 'ALL':
        if not args.state in states:
            print('state not in set of available states (see Helix_B117 data)')
            exit()
        states = [args.state]

    # make charts
    makeAllCharts(states, args.data_folder, args.fig_folder)

    print('Ratio of R_variant to R_covid as measured')
    print(R_ratio)


