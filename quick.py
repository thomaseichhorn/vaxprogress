# Script to plot the vaccination data available at:
# https://impfdashboard.de/static/data/germany_vaccinations_timeseries_v2.tsv
#
# Invoke with: python3 quick.py germany_vaccinations_timeseries_v2.tsv

import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.dates as dates
from scipy import optimize
from scipy.optimize import curve_fit
from statistics import mean
import csv
import numpy as np
import sys
import scipy.fftpack

# fit function: pol5 + pol3 * sin + const
def func ( x, a0, a1, a2, a3, a4, a5, b2, b3, b4, c1, c2, c3, d1 ) :
    return ( a5 * np.array ( x ) ** 5 + a4 * np.array ( x ) ** 4 + a3 * np.array ( x ) ** 3 + a2 * np.array ( x ) ** 2 + a1 * np.array ( x ) + a0 ) #+ ( b1 * np.array ( x ) ** 3 + b2 * np.array ( x ) ** 2 + b3 * np.array ( x ) + b4 ) * ( c1 * np.sin ( c2 * np.array ( x ) ) + c3 ) + d1

# max line count to read in a file
myrange = 1000000

# population size
full_pop = 83200000

# days for rolling average
rollingdays = 7

#0 date
#1 dosen_kumulativ
#2 dosen_differenz_zum_vortag
#3 dosen_erst_differenz_zum_vortag
#4 dosen_zweit_differenz_zum_vortag
#5 dosen_biontech_kumulativ
#6 dosen_moderna_kumulativ
#7 dosen_astrazeneca_kumulativ
#8 personen_erst_kumulativ
#9 personen_voll_kumulativ
#10 impf_quote_erst
#11 impf_quote_voll
#12 indikation_alter_dosen
#13 indikation_beruf_dosen
#14 indikation_medizinisch_dosen
#15 indikation_pflegeheim_dosen
#16 indikation_alter_erst
#17 indikation_beruf_erst
#18 indikation_medizinisch_erst
#19 indikation_pflegeheim_erst
#20 indikation_alter_voll
#21 indikation_beruf_voll
#22 indikation_medizinisch_voll
#23 indikation_pflegeheim_voll
#24 dosen_dim_kumulativ
#25 dosen_kbv_kumulativ
#26 dosen_johnson_kumulativ
#27 dosen_biontech_erst_kumulativ
#28 dosen_biontech_zweit_kumulativ
#29 dosen_moderna_erst_kumulativ
#30 dosen_moderna_zweit_kumulativ
#31 dosen_astrazeneca_erst_kumulativ
#32 dosen_astrazeneca_zweit_kumulativ 

# the arrays we read into
datearray = []
datearray_ext = []
personen_erst_kumulativ = []
personen_voll_kumulativ = []
dosen_kumulativ = []
dosen_dim_kumulativ = []
dosen_kbv_kumulativ = []

# temp var to get daily from cummulative
temp_count = []
temp_count = [0 for i in range ( 33 )]

# arrays we fill - daily
dosen_biontech_daily = []
dosen_moderna_daily = []
dosen_astrazeneca_daily = []
dosen_johnson_daily = []
dosen_biontech_erst_daily = []
dosen_biontech_zweit_daily = []
dosen_moderna_erst_daily = []
dosen_moderna_zweit_daily = []
dosen_astrazeneca_erst_daily = []
dosen_astrazeneca_zweit_daily = []

dosen_biontech_rolling = []
dosen_moderna_rolling = []
dosen_astrazeneca_rolling = []
dosen_johnson_rolling = []
dosen_biontech_erst_rolling = []
dosen_biontech_zweit_rolling = []
dosen_moderna_erst_rolling = []
dosen_moderna_zweit_rolling = []
dosen_astrazeneca_erst_rolling = []
dosen_astrazeneca_zweit_rolling = []

# counter for read lines
i = 0

# show plot window?
gui = 1

# input file to plot
if ( len ( sys.argv ) == 1 ) :
	inputfile = input ( "Please input file name " )
elif ( len ( sys.argv ) == 2 ) :
	inputfile = sys.argv[1]
else :
	print ( "Too many arguments! Usage: python", sys.argv[0], " <optional file to use>" )
	sys.exit ( )

# extend the dates read from file to extrapolate into the future
finaldate = -1

# reading
with open ( inputfile, 'r' ) as csvfile :
	print ( )
	print ( "Reading", inputfile, ":" )
	plots = csv.reader ( csvfile, delimiter = '\t' )
	for row in plots :
		# only plot every n-th value
		n = 1
		if ( i % n == 0 ) :
			# skip header line
			if ( i > 0 ) :
				try :
					# date
					datearray.append ( dates.datestr2num ( row[0] ) )
					datearray_ext.append ( dates.datestr2num ( row[0] ) )
					# save the last date
					finaldate = dates.datestr2num ( row[0] )

					# people
					personen_erst_kumulativ.append ( int ( row[8] ) )
					personen_voll_kumulativ.append ( int ( row[9] ) )

					# doses total
					dosen_kumulativ.append ( int ( row[1] ) )

					# doses by location
					dosen_dim_kumulativ.append ( int ( row[24] ) )
					dosen_kbv_kumulativ.append ( int ( row[25] ) )

					# calculate daily doses by type, subtract previous
					dosen_biontech_daily.append ( int ( row[5] ) - temp_count[5] )
					dosen_moderna_daily.append ( int ( row[6] ) - temp_count[6] )
					dosen_astrazeneca_daily.append ( int ( row[7] ) - temp_count[7] )
					dosen_johnson_daily.append ( int ( row[26] ) - temp_count[26] )

					# calculate daily first and second doses by type, again subtract previous
					dosen_biontech_erst_daily.append ( int ( row[27] ) - temp_count[27] )
					dosen_biontech_zweit_daily.append ( int ( row[28] ) - temp_count[28] )
					dosen_moderna_erst_daily.append ( int ( row[29] ) - temp_count[29] )
					dosen_moderna_zweit_daily.append ( int ( row[30] ) - temp_count[30] )
					dosen_astrazeneca_erst_daily.append ( int ( row[31]) - temp_count[31] )
					dosen_astrazeneca_zweit_daily.append ( int ( row[32]) - temp_count[32] )

					# save values in temp array
					for k in [5, 6, 7, 26, 27, 28, 29, 30, 31, 32] :
						temp_count[k] = int ( row[k] )
					
				except :
					break

				# rolling average
				if ( i >= rollingdays ) :
					temp_astrazeneca = 0
					temp_biontech = 0
					temp_johnson = 0
					temp_moderna = 0
					# sum up and divide by N... range starts at 0 so add 1
					for k in range ( rollingdays ) :
						temp_astrazeneca += dosen_astrazeneca_daily[ -1 * ( k + 1 ) ]
						temp_biontech += dosen_biontech_daily[ -1 * ( k + 1 ) ]
						temp_johnson += dosen_johnson_daily[ -1 * ( k + 1 ) ]
						temp_moderna += dosen_moderna_daily[ -1 * ( k + 1 ) ]
					dosen_astrazeneca_rolling.append ( temp_astrazeneca / ( rollingdays * 1.0 ) )
					dosen_biontech_rolling.append ( temp_biontech / ( rollingdays * 1.0 ) )
					dosen_johnson_rolling.append ( temp_johnson / ( rollingdays * 1.0 ) )
					dosen_moderna_rolling.append ( temp_moderna / ( rollingdays * 1.0 ) )
				# don't calculate
				else :
					dosen_astrazeneca_rolling.append ( 0 )
					dosen_biontech_rolling.append ( 0 )
					dosen_johnson_rolling.append ( 0 )
					dosen_moderna_rolling.append ( 0 )

		# quit early
		if i == myrange :
			print ( "Reached max reading limit!" )
			break
		i += 1

	# extend the date into the future
	for j in range ( 200 ) :
		datearray_ext.append ( finaldate + j )
	
	# func fit to population graphs
	popt1, pcov1 = curve_fit ( func, datearray, personen_erst_kumulativ )
	popt2, pcov2 = curve_fit ( func, datearray, personen_voll_kumulativ )
	intercept11 = np.interp ( full_pop * 0.75, func ( datearray_ext, *popt1 ), datearray_ext )
	intercept12 = np.interp ( full_pop, func ( datearray_ext, *popt1 ), datearray_ext )
	intercept21 = np.interp ( full_pop * 0.75, func ( datearray_ext, *popt2 ), datearray_ext )
	intercept22 = np.interp ( full_pop, func ( datearray_ext, *popt2 ), datearray_ext )
	print ( )
	print ( "Func Fit:" )
	print ( "First: 70%: ", dates.num2date  ( intercept11 ) .strftime ( '%Y-%m-%d' ), " Full: ", dates.num2date  ( intercept12 ) .strftime ( '%Y-%m-%d' ) )
	print ( "Second: 70%: ", dates.num2date  ( intercept21 ) .strftime ( '%Y-%m-%d' ), " Full: ", dates.num2date  ( intercept22 ) .strftime ( '%Y-%m-%d' ) )

	# func fit to the daily dose
	popt3, pcov3 = curve_fit ( func, datearray, dosen_kumulativ )

	# only plot if user selection
	if gui == 1 :
		# date format for visualisation
		plt.gca ( ) .xaxis.set_major_formatter ( dates.DateFormatter ( '%Y-%m-%d' ) )
		# don't plot every day
		plt.gca ( ) .xaxis.set_major_locator ( dates.DayLocator ( interval = 14 ) )
		# x axis range for dates
		plt.gca ( ) .set_xlim ( [dt.date ( 2020, 12, 1 ), dt.date ( 2021, 12, 31 ) ] )
		# y axis limit
		plt.gca ( ) .set_ylim ( 0, 100000000 )

		# horizontal line at 100%
		plt.axhline ( y = full_pop, color = 'r', linestyle = '-' )
		plt.axhline ( y = full_pop * 0.75, color = 'r', linestyle = '--' )

		plt.plot_date ( datearray_ext, func ( datearray_ext, *popt1 ), xdate = True, fmt = ':g', label = 'Fit 1 Func' )
		plt.plot_date ( datearray, personen_erst_kumulativ, tz = None, xdate = True, fmt = '+--b', label = 'Data 1st' )

		plt.plot_date ( datearray_ext, func ( datearray_ext, *popt2 ), xdate = True, fmt = ':m', label = 'Fit 2 Func' )
		plt.plot_date ( datearray, personen_voll_kumulativ, tz = None, xdate = True, fmt = '+--r', label = 'Data 2nd' )

		plt.gcf ( ) .autofmt_xdate ( )
		plt.xlabel ( 'Date' )
		plt.ylabel ( 'N' )
		plt.title ( 'Population Vaccinated' )
		plt.legend ( loc = 'upper left' )
		plt.ticklabel_format ( useOffset = False, style = 'plain', axis = 'y' )
		plt.grid ( )
		plt.figure ( )

		# date format for visualisation
		plt.gca ( ) .xaxis.set_major_formatter ( dates.DateFormatter ( '%Y-%m-%d' ) )
		# don't plot every day
		plt.gca ( ) .xaxis.set_major_locator ( dates.DayLocator ( interval = 14 ) )
		# x axis range for dates
		plt.gca ( ) .set_xlim ( [dt.date ( 2020, 12, 1 ), dt.date ( 2021, 12, 1 ) ] )
		# y axis limit
		plt.gca ( ) .set_ylim ( 0, 300000000 )

		plt.plot_date ( datearray_ext, func ( datearray_ext, *popt3 ), xdate = True, fmt = ':m', label = 'Fit 3 Func' )
		plt.plot_date ( datearray, dosen_kumulativ, tz = None, xdate = True, fmt = '+-r', label = 'Total Dose' )
		plt.gcf ( ) .autofmt_xdate ( )
		plt.xlabel ( 'Date' )
		plt.ylabel ( 'N' )
		plt.title ( 'Total Doses' )
		plt.legend ( loc = 'upper left' )
		plt.ticklabel_format ( useOffset = False, style = 'plain', axis = 'y' )
		plt.grid ( )

		plt.figure ( )
		# date format for visualisation
		plt.gca ( ) .xaxis.set_major_formatter ( dates.DateFormatter ( '%Y-%m-%d' ) )
		# don't plot every day
		plt.gca ( ) .xaxis.set_major_locator ( dates.DayLocator ( interval = 14 ) )
		# x axis range for dates
		plt.gca ( ) .set_xlim ( [dt.date ( 2020, 12, 1 ), dt.date ( 2021, 12, 1 ) ] )
		# y axis limit
		#plt.gca ( ) .set_ylim ( 0, 3000000 )

		plt.plot_date ( datearray, dosen_astrazeneca_daily, tz = None, xdate = True, fmt = '+-r', label = 'Astra' )
		plt.plot_date ( datearray, dosen_biontech_daily, tz = None, xdate = True, fmt = '+-c', label = 'BioNTech' )
		plt.plot_date ( datearray, dosen_johnson_daily, tz = None, xdate = True, fmt = '+-b', label = 'Johnson' )
		plt.plot_date ( datearray, dosen_moderna_daily, tz = None, xdate = True, fmt = '+-g', label = 'Moderna' )

		plt.gcf ( ) .autofmt_xdate ( )
		plt.xlabel ( 'Date' )
		plt.ylabel ( 'N' )
		plt.title ( 'Daily Doses by Type' )
		plt.legend ( loc = 'upper left' )
		plt.ticklabel_format ( useOffset = False, style = 'plain', axis = 'y' )
		plt.grid ( )

		plt.figure ( )
		# date format for visualisation
		plt.gca ( ) .xaxis.set_major_formatter ( dates.DateFormatter ( '%Y-%m-%d' ) )
		# don't plot every day
		plt.gca ( ) .xaxis.set_major_locator ( dates.DayLocator ( interval = 14 ) )
		# x axis range for dates
		plt.gca ( ) .set_xlim ( [dt.date ( 2020, 12, 1 ), dt.date ( 2021, 12, 1 ) ] )
		# y axis limit
		#plt.gca ( ) .set_ylim ( 0, 3000000 )

		plt.plot_date ( datearray, dosen_astrazeneca_erst_daily, tz = None, xdate = True, fmt = '+-r', label = 'Astra Erst' )
		plt.plot_date ( datearray, dosen_astrazeneca_zweit_daily, tz = None, xdate = True, fmt = '+:r', label = 'Astra Zweit' )
		plt.plot_date ( datearray, dosen_biontech_erst_daily, tz = None, xdate = True, fmt = '+-c', label = 'BioNTech Erst' )
		plt.plot_date ( datearray, dosen_biontech_zweit_daily, tz = None, xdate = True, fmt = '+:c', label = 'BioNTech Zweit' )
		plt.plot_date ( datearray, dosen_johnson_daily, tz = None, xdate = True, fmt = '+-b', label = 'Johnson' )
		plt.plot_date ( datearray, dosen_moderna_erst_daily, tz = None, xdate = True, fmt = '+-g', label = 'Moderna Erst' )
		plt.plot_date ( datearray, dosen_moderna_zweit_daily, tz = None, xdate = True, fmt = '+:g', label = 'Moderna Zweit' )

		plt.gcf ( ) .autofmt_xdate ( )
		plt.xlabel ( 'Date' )
		plt.ylabel ( 'N' )
		plt.title ( 'Daily Doses by Type' )
		plt.legend ( loc = 'upper left' )
		plt.ticklabel_format ( useOffset = False, style = 'plain', axis = 'y' )
		plt.grid ( )

		plt.figure ( )
		# date format for visualisation
		plt.gca ( ) .xaxis.set_major_formatter ( dates.DateFormatter ( '%Y-%m-%d' ) )
		# don't plot every day
		plt.gca ( ) .xaxis.set_major_locator ( dates.DayLocator ( interval = 14 ) )
		# x axis range for dates
		plt.gca ( ) .set_xlim ( [dt.date ( 2020, 12, 1 ), dt.date ( 2021, 12, 1 ) ] )
		# y axis limit
		#plt.gca ( ) .set_ylim ( 0, 3000000 )

		plt.plot_date ( datearray, dosen_dim_kumulativ, tz = None, xdate = True, fmt = '+-r', label = 'Impfzentrum' )
		plt.plot_date ( datearray, dosen_kbv_kumulativ, tz = None, xdate = True, fmt = '+-c', label = 'Hausärzte' )

		plt.gcf ( ) .autofmt_xdate ( )
		plt.xlabel ( 'Date' )
		plt.ylabel ( 'N' )
		plt.title ( 'Total Doses by Location' )
		plt.legend ( loc = 'upper left' )
		plt.ticklabel_format ( useOffset = False, style = 'plain', axis = 'y' )
		plt.grid ( )

		plt.figure ( )
		# date format for visualisation
		plt.gca ( ) .xaxis.set_major_formatter ( dates.DateFormatter ( '%Y-%m-%d' ) )
		# don't plot every day
		plt.gca ( ) .xaxis.set_major_locator ( dates.DayLocator ( interval = 14 ) )
		# x axis range for dates
		plt.gca ( ) .set_xlim ( [dt.date ( 2020, 12, 1 ), dt.date ( 2021, 12, 1 ) ] )
		# y axis limit
		#plt.gca ( ) .set_ylim ( 0, 3000000 )

		plt.plot_date ( datearray, dosen_astrazeneca_rolling, tz = None, xdate = True, fmt = '+-r', label = 'Astra' )
		plt.plot_date ( datearray, dosen_biontech_rolling, tz = None, xdate = True, fmt = '+-c', label = 'BioNTech' )
		plt.plot_date ( datearray, dosen_johnson_rolling, tz = None, xdate = True, fmt = '+-b', label = 'Johnson' )
		plt.plot_date ( datearray, dosen_moderna_rolling, tz = None, xdate = True, fmt = '+-g', label = 'Moderna' )

		plt.gcf ( ) .autofmt_xdate ( )
		plt.xlabel ( 'Date' )
		plt.ylabel ( 'N' )
		plt.title ( 'Rolling Average over {} Days'.format ( rollingdays ) )
		plt.legend ( loc = 'upper left' )
		plt.ticklabel_format ( useOffset = False, style = 'plain', axis = 'y' )
		plt.grid ( )

		plt.show ( )


#