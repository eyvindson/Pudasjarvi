# SpaFHyWeather
Get weather data from Luke 'weather' database. 
Get weather forecast from Norway database.

Though being independent module SpaFHyWeather is intended to be used with
SpaFHy (the spatial hydrology model) and its applications like DjangoSpaFHy (wetnessdemo www-interface).

SpaFHyWeather can fetch data from two
'weather' database tables: 'grid10_day' and 'grid_day'.

SpaFHyWeather requires python3 with 'scipy', 'psycopg2' (PostgreSQL adapter) 
and 'pandas' python packages. 

SpaFHyWeather has three python scripts: fmidata.py 
(access Luke 'weather' database)i, checkfmidata.py and experimental
yrforecast.py (access Norway forecast database).

A: fmidata.py: The  '-h' command line option gives  help for usage.

Sample command line for SpapFHy wetness demo Vihti (grid_day):
python fmidata.py -e 3356593 -n 6703104 -l Vihti -f VihtiFMIData -u user -p Luke_network_password

The -e (eastern) and -n (northern) are the coordinates for Vihti. In short, fmidata.py first finds 
the closest coordinate point available in the 'weather' database in the 
'grid_day' table (Euclidian distance) and then using the closest point retrieves weather 
data needed for SpaFHy. The results will appear in two files:
VihtiFMIData.csv (ascii, ';' separeted) and VihtiFMIData.xlsx (excel). 
The database queries may take a while.

PLEASE NOTE: The coordinate system is KKJ(YKJ) and please, please double check that you do not 
accidentally mix eastern and northern coordinates.

1. To tailor database queries for your needs with this script here are some functions to build your own scripts:

            1.1 fetch_fmidata(host_server,database,user,password,sql_query).
            Returns the weather data in a list.
            The function 'fetch_fmidata' shows the basic techniques 
            how to connect and pass SQL queries to 'weather' database.
            This is probably you what you will use at some point

            1.2 find_closest_point(x,y,ls): returns 3-tuple (x_closest,y_closest,date)

            1.3 find_closest_point2 2-tuple(x,y,ls): return 2-tuple (x_closest,y_closest).
      
            The functions 1.2 and 1.3 find the closest point in the the list of (x,y) coordinate points 
            queried from the weather database. Note: the closest distance is 2D Euclidian. 
            If you need more advanced you probably know what to do. The yrforcast.py 
            contains algorithm for Haversine distance that takes the curvature of the 
            earth into account. 
            
            See functions 2.1 and 2.2 how they first query at least one day of data covering whole country from their 
            respective tables before quering requested data closest to the point needed.

2. Utility functions for SpaFHy:
    
       2.1 location_to_fmidata(x,y,user,password): Fetches all weather data (from 1961) from 'grid_day' for wetnessdemo. 
            Returns a list of weather data.
    
       2.2 location_to_fmidata_grid10_day(x,y,user,password): Fetches long term weather data from 01.01.1981 from grid10_day for Sompa. Returns weather data in a pandas dataframe. Uses computational time to determie the closest point to (x,y) in the table.
    
       2.3 location_to_fmidata_grid10_day_right_x_y(x,y,user,password). As 2.2 but assumes that x and y exist in grid10_day


B: checkfmidata.py

The 'weather' database is not perfect. It may have missing dates
(grid10_day certainly has).  

1. check_fmi_data(df) and check_and_correct_fmi_data(df)

Both functions have weather data as a dataframe argument. If you have used
e.g. location_to_fmidata that return a list, simply change it
to dataframe: df = pd.DataFrame(ls).
   
      1.1 check_fmi_data: returns dictionary (lookup table) for missing dates.
   
      1.2 check_and_correct_fmi_data(df): returns 2-tuple (dictionary_missing_dates,corrected dataframe)

The function 'check_fmi' returns dictionary for missing dates. 
The keys are row indices (df.index) to the beginning of missing dates 
and the values date differences between beginning and end of missing dates 
as returned by datetime package. For example difference between 
2019.01.01 and 2019.01.02 is one and no missing dates (YYYY.MM.DD notation). 
The difference between 2019.1.1 and 2019.1.3 is two and there is one missing day.

The function 'check_and_correct_fmi_data' returns both the dictionary
for missing dates and the corrected dataframe. The weather data for
missing dates is simply linearly interpolated.

Unfortunately grid10_day and grid_day use different names for same
things in the two tables. But I guess there is a reason for that.

fmidata.py script will run daily for SpaFHy and DjangoSpaFHy 
(unix crontab) to get the latest weather data as input
for SpaFHy in wetness.luke.fi/wetnessdemo.  

The host server for Luke 'weather' database is  lukedb1.ns.luke.fi.
Your Luke network user name and password are fine but please note
you may need explicitely granted read permissions for the databse.
Contact Arto Aalto at Luke.

You can browse 'weather' database with PgAdmin program.
'weather_kannan_taulukuvaukset.docx' will become outdated 
as soon as the structure of the database changes.
Contact Anne Tapola (Luke) for further information
regarding Luke's 'weather' database. 


