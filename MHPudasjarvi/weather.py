import os
import argparse
import datetime
import pickle
import pathlib
import glob
import numpy as np
import scipy.interpolate as interpolate
import pandas as pd 
import SpaFHyWeather.fmidata as fmi
import SpaFHyWeather.checkfmidata as cfmi



daydaydate1="'"+'1961-01-01'+"'"+"::date"
daydaydate2="'"+'1961-01-01'+"'"+"::date"
query = "SELECT longitude,latitude FROM grid10_day WHERE pvm BETWEEN "+daydaydate1+" AND "+daydaydate2
N=7253044.583
E=3500934.506
def fetch_data(user,passwd):
    ls1 = fmi.fetch_fmidata(host="lukedb1.ns.luke.fi",database='weather',user=user,
                            password=passwd,query=query)
    print('Number of coordinates as a list',len(ls1))
    ls1 = list(set(ls1))
    print('Number of coordinates as a set',len(ls1),ls1[0])
    t=fmi.find_closest_point2(E,N,ls1)
    ls = fmi.location_to_fmidata_grid10_day_right_x_y_start_date(x,y,user,password,start_date
    return ls
