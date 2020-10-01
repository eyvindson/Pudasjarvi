#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 15:03:11 2019

@author: eyvindso
"""
import PIL.Image as im
import PIL.ImageOps as imops
import numpy as np
import pandas as pd
import xarray
import xarray as xr


import os
os.chdir('/proj/project_2000611/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/results/')

def leap_year(year, calendar='standard'):
    """Determine if year is a leap year"""
    leap = False
    if ((calendar in ['standard', 'gregorian',
        'proleptic_gregorian', 'julian']) and
        (year % 4 == 0)):
        leap = True
        if ((calendar == 'proleptic_gregorian') and
            (year % 100 == 0) and
            (year % 400 != 0)):
            leap = False
        elif ((calendar in ['standard', 'gregorian']) and
                 (year % 100 == 0) and (year % 400 != 0) and
                 (year < 1583)):
            leap = False
    return leap

def get_dpm(time, calendar='standard'):
    """
    return a array of days per month corresponding to the months provided in `months`
    """
    month_length = np.zeros(len(time), dtype=np.int)

    cal_days = dpm[calendar]

    for i, (month, year) in enumerate(zip(time.month, time.year)):
        month_length[i] = cal_days[month]
        if leap_year(year, calendar=calendar) and month == 2:
            month_length[i] += 1
    return month_length

dpm = {'noleap': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '365_day': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'standard': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'gregorian': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'proleptic_gregorian': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'all_leap': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '366_day': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '360_day': [0, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]}




#Alue = 0
#lengths = [100,100,100,100]
#for area in ['amaa','pkarjala','lansisuomi','itasuomi','lappi']:
#for area in ['lappi']:
def Mapping(TT):    
    print(TT)
    area = TT[0]
    lengths = TT[1]
    ds={}
    for x in range(1,TT[1]):
        ds= xr.open_mfdataset(TT[0]+str(x)+'.nc')
        
        # Make a DataArray with the number of days in each month, size = len(time)            
        month_length = xarray.DataArray(get_dpm(ds['soil_ground_water_level'].date.to_index(), calendar='noleap'),
                                    coords=[ds['soil_ground_water_level'].date], name='month_length')
        
        # Calculate the weights by grouping by 'time.season'.
        # Conversion to float type ('astype(float)') only necessary for Python 2.x
        weights = month_length.groupby('date.season') / month_length.astype(float).groupby('date.season').sum()
        
        # Test that the sum of the weights for each season is 1.0
        np.testing.assert_allclose(weights.groupby('date.season').sum().values, np.ones(4))
    
        ds= ds.assign(mean_summer_sgwl = (ds['soil_ground_water_level'] * weights).groupby('date.season').sum('date').sel(season = "JJA") )


        print("TEST")
        t_lon, t_lat = [round(elm, 4) for elm in list(ds['parameters_lon'][0].values)] ,[round(elm, 4) for elm in list(ds['parameters_lat'][0].values)]
        parameter = list((ds['soil_ground_water_level'] * weights).groupby('date.season').sum('date').sel(season = "JJA").values[0] )
    
        set_lon = set(t_lon)
        set_lat = set(t_lat)
        df =pd.DataFrame(index = set_lat, columns = set_lon) 
        for k in range(0,len(t_lon)):
            df[t_lon[k]][t_lat[k]]=parameter[k]
        
        df = df.fillna(-99999)
        df= df.sort_index()
        df= df.sort_index(axis = 1)
        
        df2 = pd.np.array(df).astype(np.uint8)
        #df2.to_csv('/proj/project_2000611/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/results/TEST'+str(TT[0])+"_"+str(x)+'.csv')
        imarray= im.fromarray(df2)
    #imarray.save('/proj/project_2000611/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/results/TEST.tif')
    
    #imops.invert(
    #imarray.rotate(180))
    #    
        imops.invert(imarray).transpose(im.FLIP_TOP_BOTTOM).save('/proj/project_2000611/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/results/TEST'+str(TT[0])+"_"+str(x)+'.tif')

Mapping(["amaa",5])
    
    print(sum(df[df > 0].count()))
    print(len(t_lon))


if __name__ == '__main__':
    import multiprocessing as mp
    
    Variables = [['amaa',10],['pkarjala',100],['lansisuomi',100],['itasuomi',100],['lappi',100]]
    
    pool = mp.Pool(5)
    pool.map(Mapping,[i for i in Variables])
