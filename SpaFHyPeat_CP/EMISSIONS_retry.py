#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 12:17:14 2019

@author: eyvindso
"""

import geopandas as gpd
import pickle
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import pandas as pd
import xarray
import xarray as xr
import os
import math

os.chdir('/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/')
import koordinaattimuunnos      
os.chdir('/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/results/')

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

def Mapping_old(TT):    
    print(TT)
    area = TT[0]
    fpath = "/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/"+str(TT[0])+"/parameters/"
    lengths = TT[1]
    ds={}
    
    x = TT[1]

    ds= xr.open_dataset(TT[0]+str(x)+'.nc')            
    
    fpath = "/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/"+str(TT[0])+"/parameters/"
    
    # Make a DataArray with the number of days in each month, size = len(time)            
    month_length = xarray.DataArray(get_dpm(ds['soil_ground_water_level'].date.to_index(), calendar='noleap'),
                                coords=[ds['soil_ground_water_level'].date], name='month_length')
    
    # Calculate the weights by grouping by 'time.season'.
    # Conversion to float type ('astype(float)') only necessary for Python 2.x
    weights = month_length.groupby('date.season') / month_length.astype(float).groupby('date.season').sum()
    
    # Test that the sum of the weights for each season is 1.0
    np.testing.assert_allclose(weights.groupby('date.season').sum().values, np.ones(4))

    ds= ds.assign(mean_summer_sgwl = (ds['soil_ground_water_level'] * weights).groupby('date.season').sum('date').sel(season = "JJA") )
    
    t_lon, t_lat = [round(elm, 5) for elm in list(ds['parameters_lon'][0].values)] ,[round(elm, 5) for elm in list(ds['parameters_lat'][0].values)]
    parameter = list((ds['soil_ground_water_level'] * weights).groupby('date.season').sum('date').sel(season = "JJA").values[0] )
    
    
    data = {'x': t_lon, 'y' : t_lat, 'z':parameter} 
    
    df =pd.DataFrame(data = data) 
                    
    return df

def Mapping(TT):    
    print(TT)
    def is_year(year,yr):
        return (year == yr)
    area = TT[0]
    fpath = "/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/"+str(TT[0])+"/parameters/"
    lengths = TT[1]
    ds={}
    x = TT[1]
    year = TT[2]
    ds= xr.open_dataset("/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/results/"+ TT[0]+str(year)+'.nc')            
    fpath = "/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/"+str(TT[0])+"/parameters/"
    # Make a DataArray with the number of days in each monthöls, size = len(time)            
    month_length = xarray.DataArray(get_dpm(ds['soil_ground_water_level'].date.to_index(), calendar='noleap'),coords=[ds['soil_ground_water_level'].date], name='month_length')
    
    num_years = len(list(ds['forcing_air_temperature'].groupby('date.year').groups.keys()))
    incl_months = [5,6,7,8,9,10]
    # Calculate the weights by grouping by 'time.season'.
    # Conversion to float type ('astype(float)') only necessary for Python 2.x
    weights = month_length.groupby('date.month') / month_length.astype(float).groupby('date.month').sum()
    weights = weights.month.isin(incl_months)*weights
    t_lon, t_lat = [round(elm, 5) for elm in list(ds['parameters_lon'][0].values)] ,[round(elm, 5) for elm in list(ds['parameters_lat'][0].values)]
    Mean_LAI_dec = list((ds['parameters_lai_decid_max'] * weights).sum('date').values[0]/(len(incl_months) ))
    MEAN_LAI_conf = list((ds['parameters_lai_conif'] * weights).sum('date').values[0]/(len(incl_months) ))

    Mean_AirTemp = []
    Mean_precip = []
    for i in list(ds['forcing_air_temperature'].groupby('date.year').groups.keys()):
        weights1 = weights * is_year(ds['forcing_air_temperature']['date.year'],i)*(num_years)
        if len(Mean_AirTemp) == 0:
            Mean_AirTemp = np.array((ds['forcing_air_temperature'] * weights1).sum('date').values[0]/(len(incl_months) ))
            Mean_precip = np.array((ds['forcing_precipitation'] * weights1).sum('date').values[0]/(len(incl_months) ))
        else:
            Mean_AirTemp = Mean_AirTemp + np.array((ds['forcing_air_temperature'] * weights1).sum('date').values[0]/(len(incl_months) ))
            Mean_precip = Mean_precip + np.array((ds['forcing_precipitation'] * weights1).sum('date').values[0]/(len(incl_months) ))
    Mean_AirTemp = Mean_AirTemp/len(list(ds['forcing_air_temperature'].groupby('date.year').groups.keys()))
    Mean_precip= Mean_precip/len(list(ds['forcing_precipitation'].groupby('date.year').groups.keys()))

    Mean_AirTemp_10 = []
    Mean_precip_10 = []
    for i in list(ds['forcing_air_temperature'].groupby('date.year').groups.keys())[-10:]:
        weights1 = weights * is_year(ds['forcing_air_temperature']['date.year'],i)*(num_years)
        if len(Mean_AirTemp_10) == 0:
            Mean_AirTemp_10 = np.array((ds['forcing_air_temperature'] * weights1).sum('date').values[0]/(len(incl_months) ))
            Mean_precip_10 = np.array((ds['forcing_precipitation'] * weights1).sum('date').values[0]/(len(incl_months) ))
        else:
            Mean_AirTemp_10 = Mean_AirTemp_10 + np.array((ds['forcing_air_temperature'] * weights1).sum('date').values[0]/(len(incl_months) ))
            Mean_precip_10 = Mean_precip_10 + np.array((ds['forcing_precipitation'] * weights1).sum('date').values[0]/(len(incl_months) ))
    Mean_AirTemp_10 = Mean_AirTemp_10/len(list(ds['forcing_air_temperature'].groupby('date.year').groups.keys())[-10:])  
    Mean_precip_10= Mean_precip_10/len(list(ds['forcing_air_temperature'].groupby('date.year').groups.keys())[-10:])
    data = {'x': t_lon, 'y' : t_lat, 'DEC_LAI': Mean_LAI_dec, 'CON_LAI': MEAN_LAI_conf,'Air_Temp': Mean_AirTemp,'Precip': Mean_precip,'Air_Temp_10': Mean_AirTemp_10,'Precip_10': Mean_precip_10}
    k = 0
    Mean_GWT = []
    for i in list(ds['soil_ground_water_level'].groupby('date.year').groups.keys()):
        weights1 = weights * is_year(ds['soil_ground_water_level']['date.year'],i)*(num_years)
        data1 = {"Mean_GWT" +str(i):list((ds['soil_ground_water_level'] * weights1).sum('date').values[0]/(len(incl_months) ))}
        data.update(data1)
        k = k+1
    df =pd.DataFrame(data = data)
    df.to_csv("/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/test.csv")
    return df

def N2O_CH4_C2O_EMISSIONS(GTW,peat_loc,region,file_start,DD,UPFERT,MPU,year):
    path = "/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/"
    with open(path + "results/"+region+".pkl", 'rb') as pickle_file:
        GTW = pickle.load(pickle_file)
    folder = "/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/"+region+"/Data/"
    name_data_MPU = "V5/PUDASJARVI_P3005"+MPU+"_SKE4_4_V5.MPU"
    data_MPU = pd.read_csv(folder+ name_data_MPU,delim_whitespace=True,header = None)
    kasvu = data_MPU[3]
    paa = data_MPU[2]    
    #kasvu = pd.read_csv("/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/"+region+"/parameters/UPDATED_Fertility.dat",skiprows=6,header = None)
    #pixels = np.loadtxt(path + region + "/segmaj1.pseudotype.hits")
    #np.savetxt(path + region + "/"+region+"_segmaj1.pseudotype.hits",pixels)
    #pixels = pd.read_csv(path + region + "/"+region+"_segmaj1.pseudotype.hits",header = None, sep = " ")
    #paa = np.loadtxt(path + region + "/segmaj1.paatyyppi_vmi1x_1216.hits")
    #np.savetxt(path + region + "/"+region+"_segmaj1.paatyyppi_vmi1x_1216.hits",paa)
    #paa = pd.read_csv(path + region + "/"+region+"_segmaj1.paatyyppi_vmi1x_1216.hits",header = None,sep=" ")
    peat_loc = pd.read_csv('/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/'+region+'/'+file_start+'_peat_loc.csv',header = None)   
    Drainage = pd.read_csv("/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/"+region+"/parameters/ditch_spacing.dat",skiprows=6,header = None)
    Depth =  pd.read_csv("/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/"+region+"/parameters/ditch_depth.dat",skiprows=6,header = None)
    #Canopy_height =  pd.read_csv("/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/"+region+"/parameters/hc.dat",skiprows=6,header = None)
    #Canopy_fraction =  pd.read_csv("/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/"+region+"/parameters/cf.dat",skiprows=6,header = None)
    East =  pd.read_csv("/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/"+region+"/parameters/Ecoord.dat",skiprows=6,header = None)
    North =  pd.read_csv("/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/"+region+"/parameters/Ncoord.dat",skiprows=6,header = None)
    
    df1=pd.DataFrame(data=peat_loc)
    df2=pd.DataFrame(data=kasvu) 
    #df3=pd.DataFrame(data=pixels)
    df4=pd.DataFrame(data=paa)    
    df5=pd.DataFrame(data=Drainage)
    df6=pd.DataFrame(data=Depth)
    #df7=pd.DataFrame(data=Canopy_height)
    #df8=pd.DataFrame(data=Canopy_fraction)
    df9=pd.DataFrame(data=East)
    df10=pd.DataFrame(data=North)
    segm = data_MPU[0]    
    df11=pd.DataFrame(data=segm)
    osuus = data_MPU[193]    
    df12=pd.DataFrame(data=osuus)
        
    def merge_dfs(df1,df_loc,dup,x = True):
        DF1 = df1[0].isin(df_loc[0])
        if x == True:
            m_DF1 = DF1* df1[2]
        else:
            m_DF1 = DF1* df1[1]
        m_DF1 = list(m_DF1)
        m_DF1 = list(filter(lambda a: a != 0, m_DF1))
        t_DF1 = pd.DataFrame(m_DF1, columns=["Segment"])
        t_DF1 =t_DF1.drop(dup.values)
        return t_DF1
    #if UPFERT == False:
    #    t_KASVU =merge_dfs(df2,df1,duplicate_part)
    #else:
    #    t_KASVU = df2
    #    t_KASVU.columns = ["Segment"]
    #t_PIXELS =merge_dfs(df3,df1,duplicate_part,x = False)
    t_PAA =df4
    t_PAA.columns= ["Segment"]
    t_KASVU = df2
    t_KASVU.columns = ["Segment"]
    t_seg = df11
    t_seg.columns = ["Segment"]
    t_osuus = df12
    t_osuus.columns = ["Segment"]
    t_Drain = df5
    t_Drain.columns = ["Segment"]
    t_Depth = df6
    t_Depth.columns = ['Segment']
    #t_Canopy_height = df7
    #t_Canopy_height.columns = ['Segment']
    #t_Canopy_fraction = df8
    #t_Canopy_fraction.columns = ['Segment']
    t_East = df9
    t_East.columns = ['Segment']
    t_North = df10
    t_North.columns = ['Segment']
    
    t_KASVU_N2O = t_KASVU.replace(1,20).replace(2,20).replace(3,25).replace(4,29).replace(5,39).replace(6,48).replace(7,48).replace(8,48).replace(9,48).replace(10,48)
    Nut_rich = t_KASVU.replace(1,1).replace(2,1).replace(3,1).replace(4,0).replace(5,0).replace(6,0).replace(7,0).replace(8,0).replace(9,0).replace(10,0)
    Nut_poor = Nut_rich*-1+1
    #Calculate CH4 ## calculated based on Ojanen et al. 2010. Soil–atmosphere CO2, CH4 and N2O fluxesinborealforestry-drainedpeatlands. Forest Ecology and Management.
    CH4_y0 =  -0.378
    CH4_a = 12.3
    CH4_b = 0.121
    #CH4_flux = CH4_y0 +CH4_a*np.exp(-CH4_b*GTW['z'])
    CH4_flux = []
    key = list(GTW.keys())
    remove_keys = ['x','y','DEC_LAI','CON_LAI','Air_Temp', 'Precip', 'Air_Temp_10','Precip_10']
    for rem in remove_keys:
        if rem in key:
            key.remove(rem)
    
    key = list(key)
    for i in key:
        CH4_flux = CH4_flux + [CH4_y0 +CH4_a*np.exp(-CH4_b*(-1*GTW[i]*100))]
    #Calculate C02:
    Nut_rich = t_KASVU.replace(1,1).replace(2,1).replace(3,1).replace(4,0).replace(5,0).replace(6,0).replace(7,0).replace(8,0).replace(9,0).replace(10,0)
    Nut_poor = Nut_rich*-1+1
    NetSoil_CO2 = []
    for i in key:
        NetSoil_CO2_NR_t = -115+12*(-1*GTW[i]*100)
        NetSoil_CO2_NP_t = -259+6*(-1*GTW[i]*100)
        NetSoil_CO2_t = Nut_rich['Segment']*NetSoil_CO2_NR_t+Nut_poor['Segment']*NetSoil_CO2_NP_t
        NetSoil_CO2= NetSoil_CO2+[NetSoil_CO2_t]        
    #Calculate N20   -- BASED ON equation provided by Ojanen "Lähde/viite tällä erää: Minkkinen, K., Ojanen, P., Koskinen, M. & Penttilä, T. Nitrous oxide emissions of undrained, forestry-drained and rewetted boreal peatlands. Submitted."
    N2O_a = 0.699
    N2O_b = 0.00990
    N2O_c = 0.211
    N2O_d = 0.0794
    N2O = []
    for i in key:
        N2O = N2O+[(N2O_a+N2O_b*(-100*GTW[i])**2)*np.exp(-N2O_c*t_KASVU_N2O["Segment"])+N2O_d]
    #Create dataframe to return
    EM_data = {"Segment":df1[0],"ORIG_seg":t_seg["Segment"],"Osuus":t_osuus["Segment"],"LAI_DEC":GTW['DEC_LAI'], "LAI_CON":GTW['CON_LAI'],"Depth":t_Depth['Segment'],"Kasvu":t_KASVU['Segment'],"Ditch":t_Drain['Segment'],"East":t_East['Segment'],"North":t_North['Segment'],"Air_Temp_Mean_ALL":GTW['Air_Temp'],"Precip_Mean_ALL":GTW['Precip'],"Air_Temp_Mean_Last_10":GTW['Air_Temp_10'],"Precip_Mean_Last_10":GTW['Precip_10'],"Soil_type":t_KASVU['Segment']}
   
    k = 0
    for i in key:
        EM_data_1 = {"N2O"+str(i[8:]):N2O[k],"CH4"+str(i[8:]): CH4_flux[k],"CO2"+str(i[8:]): NetSoil_CO2[k],"GWT"+str(i[8:]): GTW[i]}
        EM_data.update(EM_data_1)
        k = k+1
    EM_data_pd = pd.DataFrame(EM_data)
    EM_data_pd.to_csv(path+region+"/"+region+"EMISSIONS_output"+year+".csv")
    return key


if __name__ == '__main__':

    import argparse
    import multiprocessing as mp
    import glob
    import pickle
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', help='parameter folder', type=str)
    parser.add_argument('--MPU', help ='which MPU data', type=str)
    parser.add_argument('--year', help ='which year data', type=str)
    run_map = False
    
    args = parser.parse_args()
            
    pool = mp.Pool(mp.cpu_count())
    end = 100
    if args.folder in["amaa", "Test_case"]:
        end = 10
    DATA =  pool.map(Mapping,[[args.folder,0,args.year]])
    
    DATA1 = DATA[0]
    #for dat in range(1,len(DATA)):
    #    DATA1 = pd.concat([DATA1,DATA[dat]],ignore_index = True)
    with open(args.folder + '.pkl', 'wb') as f:
        pickle.dump(DATA1, f)
    import time
    time.sleep(6)
    file_start = "amaa"
    rr = "amaa"
    if args.folder != "amaa":
        if args.folder == "lappi":
            file_start = "la"
            rr = "la"
        elif args.folder == "itasuomi":
            file_start = "is"
            rr = "is2"
        elif args.folder == "lansisuomi":
            file_start = "ls"
            rr = "ls"
        elif args.folder == "pkarjala":
            file_start = "pk"
            rr = "pk"
        elif args.folder == "Puhdasjarvi_data_51":
            file_start = "Pd_51"
            rr = "Pd_51"
        elif args.folder == "Puhdasjarvi_data_52":
            file_start = "Pd_52"
            rr = "Pd_53"
        elif args.folder == "Puhdasjarvi_data_53":
            file_start = "Pd_52"
            rr = "Pd_53"
     
    #Mapping WTD
    #all_loc_data = pd.read_csv("/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/"+args.folder+"/"+file_start+"_segmean_XY_ERTS-YKJ-latlongEPSG2458.txt",sep=" ",header = None)
    Peat_loc_data_E = pd.read_csv("/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/"+args.folder+"/parameters/Ecoord.dat",skiprows=6,header = None)
    Peat_loc_data_E = Peat_loc_data_E.rename(columns={0: "E"})
    Peat_loc_data_N = pd.read_csv("/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/"+args.folder+"/parameters/Ncoord.dat",skiprows=6,header = None)
    Peat_loc_data_N = Peat_loc_data_N.rename(columns={0: "N"})
    Peat_loc_data = pd.concat([Peat_loc_data_N, Peat_loc_data_E], axis=1, sort=False)
    tt = []

    #all_combined = all_loc_data[4]+all_loc_data[5]*100000000
    all_combined = Peat_loc_data["E"]+Peat_loc_data["N"]*100000000
    
    #duplicate_part = part_combined[part_combined.duplicated()].index
    #duplicate_all = all_combined[all_combined.duplicated()].index

    df1=pd.DataFrame(data=all_combined)
    #df2=pd.DataFrame(data=part_combined)
    #TEST = df1[0].isin(df2[0])
    #all_loc = all_loc_data[0]
    #all_loc.drop(duplicate_all.values)
    #tt = TEST* all_loc#.values
    #tt = list (tt)
    #tt = list(filter(lambda a: a != 0, tt))
    peat_loc = df1#pd.DataFrame(tt, columns=["Segment"])
    #peat_loc =peat_loc.drop(duplicate_part.values)
    if os.path.isfile('/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/'+args.folder+'/'+file_start+'_peat_loc.csv') == False:
        peat_loc.to_csv('/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/'+args.folder+'/'+file_start+'_peat_loc.csv', header = None, index=False)
    Drainage = pd.read_csv("/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/"+args.folder+"/parameters/ditch_spacing.dat",skiprows=6,header = None)
    Depth =  pd.read_csv("/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/"+args.folder+"/parameters/ditch_depth.dat",skiprows=6,header = None)
    D2 = Drainage
    D2['Depth'] = Depth
    D2.columns = ["Ditch","Depth"]
    
    key = N2O_CH4_C2O_EMISSIONS(DATA1,peat_loc,args.folder,file_start,D2,True,args.MPU,args.year)
    