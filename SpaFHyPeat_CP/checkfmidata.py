import argparse
import matplotlib
matplotlib.use('Agg')
import os
import datetime
import numpy as np
import scipy.interpolate as interpolate
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import fmidata as fmi

aika_header =['aika']
date_header =['vuosi','kk','paiva']
weather_header = ['t_mean','t_max','t_min','rainfall','radiation','hpa','lamposumma_v','rainfall_v']

def plot_weather(df,years:int,file_name:str):
    """
    Plot rainfall as barplot and temperature as line 
    starting from 'from_date' as pd.Timestamp(YYYY,MM,DD)
    """
    dpi_value=600
    with PdfPages(file_name) as pdf:
        plt.rc('text', usetex=True)
        (fig, ax1) = plt.subplots()
        plt.rc('font', family='serif')
        plt.rc('mathtext',fontset='stix')
        plt.title('Daily Weather',fontsize=20)
        df['aika'] = df['aika'].astype(str)
        df['aika'] = df['aika'].apply(pd.to_datetime)
        days = years*365
        rainfall_ls = list(df['rainfall'])
        slice_days = max(len(rainfall_ls)-days,0)
        rainfall_ls = rainfall_ls[slice_days:]
        tmean_ls = list(df['t_mean'])
        tmean_ls = tmean_ls[slice_days:]
        date_ls = list(df['aika'])
        date_ls = date_ls[slice_days:]
        datemin = date_ls[0]
        #datemin = np.datetime64(datemin,'D')
        #print(datemin)
        datemax = date_ls[-1]
        #print(datemax)
        #datemax = np.datetime64(datemax,'D')+np.timedelta64(1,'D')
        print(datemin,datemax)
        ax1.set_xlim(datemin,datemax)
        ax1.bar(date_ls,rainfall_ls,align='center',color='b',alpha=0.8)
        ax1.set_ylabel('Daily Precipitation (mm)',color='b',fontsize=16)
        ax2 = ax1.twinx()
        cold_limit= -5.0
        ax2.axhline(y=cold_limit,lw=1)
        tmean_array = np.asarray(tmean_ls)
        t_cold_array = np.ma.masked_where(tmean_array >= cold_limit,tmean_array)
        t_warm_array = np.ma.masked_where(tmean_array <= cold_limit,tmean_array)
        #ax2.plot(date_ls, t_cold_array,date_ls,t_warm_array,lw=0.5)
        #ax2.plot(date_ls,tmean_array,lw=0.5,ls=':',color='r')
        ax2.plot(date_ls,tmean_array,'ro',markersize=1)
        ax2.set_ylabel('Daily Mean Temperature ($^\circ$C)', color='r',fontsize=16)
        #Rotate x-axis labels (dates) automatically if needed 
        fig.autofmt_xdate()
        fig.tight_layout()
        pdf.savefig()
        png_file=os.path.splitext(file_name)[0]
        #png suites better for www-pages, it seems e.g. that
        #Chrome does not show pdf files by default 
        png_file=png_file+'.png'
        plt.savefig(png_file,dpi=dpi_value)
        plt.close()
            
def insert_row(df,row_number,row):
    start_upper=0
    end_upper=row_number
    start_lower = row_number+1
    end_lower = df.shape[0]+1
    upper_index_ls = list(range(start_upper,end_upper))
    lower_index_ls = list(range(start_lower,end_lower))
    new_index_ls = upper_index_ls+lower_index_ls
    df.index = new_index_ls
    df.loc[row_number]=row
    df = df.sort_index()
    return df

def insert_row_by_transpose(df,row_number,row_ls):
    """
    Insert row by first transposing data frame,
    then inserting the row as a column and finally
    transpose the data frame back to original form.
    It is assumed that row index is with numbers
    """
    df=df.transpose()
    #insert happens "in place"
    df.insert(row_number,row_number,row_ls,allow_duplicates=True)
    df.columns = range(0,len(df.columns))
    df=df.transpose()
    return df
      
def check_fmidata(df):
    #shape = np.shape(df)
    #rows = shape[0]
    rows = len(df.index)
    print('Rows',rows)
    delta_dict=dict()
    print('Checking duplicate rows') 
    series=df[df.duplicated()]
    print('Done')
    if len(series) > 0:
        print('Duplicate rows',len(series))
    for row in range(rows-1):
        r1 = df.loc[row,['vuosi','kk','paiva']]
        r2 = df.loc[row+1,['vuosi','kk','paiva']]
        d1 = datetime.date(r1.vuosi,r1.kk,r1.paiva)
        d2 = datetime.date(r2.vuosi,r2.kk,r2.paiva)
        delta = d2-d1
        days = delta.days
        if days > 1:
            print('Missing days error',row,d1,d2,days)
            delta_dict[row]=([r1.paiva,r1.kk,r1.vuosi,r2.paiva,r2.kk,r2.vuosi,days])
        if days == 0:
            print('Two same dates error',row,row+1,d1,d2,days)
        if days < 0:
            print('Earlier date error',row,row+1,d2,days)
    return delta_dict

def check_and_correct_weather_data(df,drop_duplicates=True):
    """
    Correct weather data: interpolate missing values
    df: weather data dataframe, return by for example read_weather_data_grid10_day_right_x_y in fmidata.
    return corrected weather data
    """  
    delta_dict = check_fmidata(df)
    keys = sorted(delta_dict.keys())
    increment=0
    if drop_duplicates==True:
        df.drop_duplicates(inplace=True)
    print("Missing days",delta_dict.keys())
    for key in keys:
        #Beginning of the missing dates
        row1 = df.loc[key+increment,]
        #End of missing dates
        row2 = df.loc[key+increment+1,]
        #print('Day1',key+increment,row1)
        #print('Day2',key+increment+1,row2)
        dates = delta_dict[key]
        d1=datetime.date(dates[2],dates[1],dates[0])
        #Weather data from the two dates
        weather_ls1=row1[weather_header]
        weather_ls2=row2[weather_header]
        #Create bounding y-values 
        t_mean_ls=[weather_ls1[0],weather_ls2[0]]
        t_max_ls=[weather_ls1[1],weather_ls2[1]]
        t_min_ls=[weather_ls1[2],weather_ls2[2]]
        rainfall_ls=[weather_ls1[3],weather_ls2[3]]
        radiation_ls=[weather_ls1[4],weather_ls2[4]]
        hpa_ls=[weather_ls1[5],weather_ls2[5]]
        lamposumma_v_ls=[weather_ls1[6],weather_ls2[6]]
        rainfall_v_ls=[weather_ls1[7],weather_ls2[7]]
        d2=datetime.date(dates[5],dates[4],dates[3])
        time_delta = d2-d1
        #Number of days (bounding x-values)
        x_ls = [0,time_delta.days]
        #Interpolation functions
        f0 = interpolate.interp1d(x_ls,t_mean_ls)
        f1 = interpolate.interp1d(x_ls,t_max_ls)
        f2 = interpolate.interp1d(x_ls,t_min_ls)
        f3 = interpolate.interp1d(x_ls,rainfall_ls)
        f4 = interpolate.interp1d(x_ls,radiation_ls)
        f5 = interpolate.interp1d(x_ls,hpa_ls)
        f6 = interpolate.interp1d(x_ls,lamposumma_v_ls)
        f7 = interpolate.interp1d(x_ls,rainfall_v_ls)
        for i in range(1,time_delta.days):
            print(key,increment,time_delta.days,i)
            weatherinterp1d_ls = [f0(i),f1(i),f2(i),f3(i),f4(i),f5(i),f6(i),f7(i)]
            d2=d1+datetime.timedelta(days=i)
            row1[date_header]=[d2.year,d2.month,d2.day]
            row1[weather_header] = weatherinterp1d_ls
            aika = fmi.date_to_aika(d2)
            print('Corrected aika',aika)
            row1[aika_header] = aika
            df=insert_row_by_transpose(df,key+increment+(i),row1)
        #For example, difference 1 day has no days in between,
        #difference 2 days has 1 day in between, 3 days 2 etc.
        increment=increment+(time_delta.days-1)
    return (delta_dict,df)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f',type=str,dest='f',help="FMI weather data, csv file ';' separated")
    parser.add_argument('-p',dest='p',action='store_true',help="Plot weather data, pdf file")
    parser.add_argument('-y',type=int,dest='y',default=0,help="Display number of years")
    args=parser.parse_args()
    if not args.f:
        print("No FMI weather data")
        quit()
    df = pd.read_csv(args.f,sep=';')
    #delta_dict = check_fmidata(df)
    #if delta_dict:
    #    print(delta_dict)
    #else:
    #    print("No missing days")
    if args.p:
        fname=args.f.replace('.csv','.pdf')
        plot_weather(df,args.y,fname)
