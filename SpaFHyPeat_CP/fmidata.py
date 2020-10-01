#!python
import psycopg2
from scipy.spatial import distance
import pandas as pd
import datetime
import argparse
import checkfmidata as cfmi

#index for data frame 
ildatacols=['OmaTunniste','OmaIta','OmaPohjoinen','Kunta','siteid','aika','vuosi','kk','paiva','longitude','latitude','t_mean',\
            't_max','t_min','rainfall','radiation','hpa','lamposumma_v','rainfall_v']
            
def compare_distance(x):
    return x[0]

#Find closest point, euclidian distance
def find_closest_point(x,y,ls):
    """
    The list ls consists of 3-tuple of coordinates and (x,y,date)
    return 2-tuple (distance,(x,y)) 
    """
    ls1 = [(distance.euclidean([x,y],[x1,y1]),(x1,y1)) for (x1,y1,z) in ls]
    ls1 = sorted(ls1,key=compare_distance)
    return ls1[0]

def find_closest_point2(x,y,ls):
    """
    The list ls consists only of coordinate tuples (x,y).
    Return 2-tuple (distance,(x,y))
    """
    ls1 = [(distance.euclidean([x,y],[x1,y1]),(x1,y1)) for (x1,y1) in ls]
    ls1 = sorted(ls1,key=compare_distance)
    return ls1[0]

#column 'aika' = str(year)+str(month)+str(day)
def date_to_aika(date):
    year=str(date.year)
    month=str(date.month)
    day=str(date.day)
    if date.month < 10:
        month=str(0)+month
    if date.day < 10:
        day=str(0)+day
    date=year+month+day
    return date

def add_aika_vuosi_kk_paiva(ls):
    ls1 = [(x,y,date_to_aika(date),date.year,date.month,date.day,x,y,temp_avg,temp_max,temp_min,prec,global_rad,vapour_press) for\
               (x,y,date,x,y,temp_avg,temp_max,temp_min,prec,global_rad,vapour_press)  in ls]
    return ls1

def add_aika_vuosi_kk_paiva_grid10_day(ls):
    ls1 = [(x,y,siteid,date_to_aika(date),date.year,date.month,date.day,x,y,temp_avg,temp_max,temp_min,prec,global_rad,vapour_press) for\
               (x,y,siteid,date,x,y,temp_avg,temp_max,temp_min,prec,global_rad,vapour_press)  in ls]
    return ls1
def add_fill_missing_column(df,col_num,name,initial_value):
    df.insert(col_num,name,initial_value)
    return df

def calculate_annual_rainfall(df):
    annual_rainfall = 0.0
    for i in range(len(df.index)):
        day = df.loc[i,'paiva']
        month = df.loc[i,'kk']
        if day == 1 and month == 1:
            annual_rainfall=0.0
        rainfall = df.loc[i,'rainfall']
        annual_rainfall=annual_rainfall+rainfall
        df.loc[i,'rainfall_v'] = annual_rainfall
    return df
        
def calculate_day_degree(df):
    day_degree=0.0
    dd_limit=5.0
    for i in range(len(df.index)):
        day = df.loc[i,'paiva']
        month = df.loc[i,'kk']
        if day == 1 and month == 1:
            day_degree=0.0
        t_mean=df.loc[i,'t_mean']
        dd_tmp=t_mean-dd_limit
        day_degree = day_degree + max(0,dd_tmp)
        df.loc[i,'lamposumma_v']=day_degree
    return df
    
def fetch_fmidata(host,database,user,password,query): 
    conn = psycopg2.connect(host=host,database=database,user=user,password=password) 
    cur = conn.cursor()
    #cur.execute('SELECT version()')                                                                      
    #version=cur.fetchone()
    #print(version)
    cur.execute(query)
    ls = cur.fetchall()
    conn.close()
    return ls

#This is for Wetness demo with SpaFHy, grid_day table
def location_to_fmidata(x,y,user,password):
    """"
    x: longitude
    y: latitude
    """
    ls1 = []
    time_delta = 10
    while ls1 ==[]:
        dayday=(datetime.date.today()-datetime.timedelta(time_delta)).strftime('%Y-%m-%d')
        #SQL syntax requires the use of single quotes (i.e. the "'"-character)
        daydaydate="'"+dayday+"'"+"::date"
        #This weather database/table should update daily. But life is. Fetch first a couple of days of data
        ls1 = fetch_fmidata(host="lukedb1.ns.luke.fi",database='weather',user=user,password=password,
                           query="SELECT x,y,date FROM grid_day WHERE date >="+daydaydate)
        #then try again if necessary
        time_delta = time_delta+10
    t = find_closest_point(x,y,ls1)
    x1=t[1][0]
    y1=t[1][1]
    #With the closest point select all data
    #print('Closest point','N',y1,'E',x1)    
    ls = fetch_fmidata(host="lukedb1.ns.luke.fi",database='weather',user=user,password=password,
                      query="SELECT x,y,date,x,y,temp_avg,temp_max,temp_min,prec,global_rad,vapour_press FROM grid_day WHERE x="+str(x1)+" AND y="+str(y1)+" AND date >= '1961-01-01'::date")
    return ls

def location_to_fmidata_grid10_day(x,y,user,password):
    """"
    x: longitude
    y: latitude
    """
    #n=x+1000
    #s=x-1000
    ls1 = []
    #SQL syntax requires the use of single quotes (i.e. the "'"-character)
    #grid10_day has data up to 2016 and BETWEEN query/command seems to be faster than e.g.
    #arithmetic operations
    daydaydate="'"+'2016-01-01'+"'"+"::date"
    daydaydate2="'"+'2016-01-01'+"'"+"::date"
    #Fetch first one day of data, grid10_day table seems to have data for 2016-01-01
    ls1 = fetch_fmidata(host="lukedb1.ns.luke.fi",database='weather',user=user,password=password,
                        query="SELECT longitude,latitude,pvm FROM grid10_day WHERE pvm BETWEEN "+daydaydate+" AND "+daydaydate2)
    t = find_closest_point(x,y,ls1)
    x1=t[1][0]
    y1=t[1][1]
    #With the closest point select all data
    #print('Closest point','N',y1,'E',x1)    
    ls = fetch_fmidata(host="lukedb1.ns.luke.fi",database='weather',user=user,password=password,
                       query="SELECT longitude,latitude,siteid,pvm,longitude,latitude,t_mean,t_max,t_min,rainfall,radiation,hpa FROM grid10_day WHERE longitude="+str(x1)+" AND latitude="+str(y1)+" AND pvm >= '1981-01-01'::date")
    return ls

def location_to_fmidata_grid10_day_right_x_y(x,y,user,password):
    #With the closest point select all data
    #print('Closest point','N',y1,'E',x1)    
    ls = fetch_fmidata(host="lukedb1.ns.luke.fi",database='weather',user=user,password=password,
                       query="SELECT longitude,latitude,siteid,pvm,longitude,latitude,t_mean,t_max,t_min,rainfall,radiation,hpa FROM grid10_day WHERE longitude="+str(x)+" AND latitude="+str(y)+" AND pvm >= '1981-01-01'::date")
    return ls
    
def write_ls_to_excel(file_name,ls,column_ls):
    writer = pd.ExcelWriter(file_name,engine='xlsxwriter')
    df = pd.DataFrame(ls)
    df.columns=column_ls
    df.to_excel(writer,sheet_name='ILData')
    writer.save()

def write_df_to_excel(file_name,sheet_name,df,write_index=True):
    writer = pd.ExcelWriter(file_name,engine='xlsxwriter')
    df.to_excel(writer,index=write_index,sheet_name=sheet_name)
    writer.save()

def write_df_to_csv(file_path,df,sep=';'):
    df.to_csv(file_path,sep)

#This is for Wetness demo, SpaFHY using grid_day table.
def read_and_write_weather_data(x,y,loc,file_path,user,passwd):
    #database query returns the list of found database entries 
    ls = location_to_fmidata(x,y,user,passwd)
    #add date related columns SpaFHy needs
    ls1 = add_aika_vuosi_kk_paiva(ls)
    #create data frame 
    df = pd.DataFrame(ls1)
    #add rest of the columns SpaFHy needs with default values
    df = add_fill_missing_column(df,0,'OmaTunniste',1)
    df = add_fill_missing_column(df,3,'Kunta',loc)
    df = add_fill_missing_column(df,4,'siteid',9999)
    df = add_fill_missing_column(df,len(df.columns),'lamposumma_v',0)
    df = add_fill_missing_column(df,len(df.columns),'rainfall_v',0)
    df.columns=ildatacols
    df['OmaIta']=x
    df['OmaPohjoinen']=y
    #Sort Values!!!
    df.sort_values(by='aika',kind='mergesort',inplace=True)
    df.index = range(0,len(df.index))
    #calculate real values for day degree and annual rainfall
    df=calculate_day_degree(df)
    df=calculate_annual_rainfall(df)
    #write csv file
    fcsv=file_path
    if fcsv.endswith('.xlsx'):
        fcsv = fcsv.replace('.xlsx','.csv')
    if not fcsv.endswith('.csv'):
        fcsv = fcsv+'.csv'
    df.to_csv(fcsv,sep=';',index=False)
    #write excel file
    fxlsx = fcsv.replace('.csv','.xlsx')
    write_df_to_excel(fxlsx,'FMIData'+loc,df)
    return df

def read_and_write_weather_data_grid10_day(x,y,loc,file_path,user,passwd):
    #database query returns the list of found database entries 
    ls = location_to_fmidata_grid10_day(x,y,user,passwd)
    #add date related columns SpaFHy needs
    ls1 = add_aika_vuosi_kk_paiva_grid10_day(ls)
    #create data frame 
    df = pd.DataFrame(ls1)
    #add rest of the columns SpaFHy needs with default values
    df = add_fill_missing_column(df,0,'OmaTunniste',1)
    df = add_fill_missing_column(df,3,'Kunta',loc)
    df = add_fill_missing_column(df,len(df.columns),'lamposumma_v',0)
    df = add_fill_missing_column(df,len(df.columns),'rainfall_v',0)
    df.columns=ildatacols
    df['OmaIta']=x
    df['OmaPohjoinen']=y
    df.sort_values(by='aika',kind='mergesort',inplace=True)
    df.index = range(0,len(df.index))
    #calculate real values for day degree and annual rainfall
    df=calculate_day_degree(df)
    df=calculate_annual_rainfall(df)
    #write csv file
    fcsv=file_path
    if fcsv.endswith('.xlsx'):
        fcsv = fcsv.replace('.xlsx','.csv')
    if not fcsv.endswith('.csv'):
        fcsv = fcsv+'.csv'
    df.to_csv(fcsv,sep=';',index=False)
    #write excel file
    fxlsx = fcsv.replace('.csv','.xlsx')
    write_df_to_excel(fxlsx,'FMIData'+loc,df)
    return df

def read_weather_data_grid10_day(x,y,loc,user,passwd):
    #database query returns the list of found database entries 
    ls = location_to_fmidata_grid10_day(x,y,user,passwd)
    #add date related columns SpaFHy needs
    ls1 = add_aika_vuosi_kk_paiva_grid10_day(ls)
    #create data frame 
    df = pd.DataFrame(ls1)
    #add rest of the columns SpaFHy needs with default values
    df = add_fill_missing_column(df,0,'OmaTunniste',1)
    df = add_fill_missing_column(df,3,'Kunta',loc)
    df = add_fill_missing_column(df,len(df.columns),'lamposumma_v',0)
    df = add_fill_missing_column(df,len(df.columns),'rainfall_v',0)
    df.columns=ildatacols
    df['OmaIta']=x
    df['OmaPohjoinen']=y
    df.sort_values(by='aika',kind='mergesort',inplace=True)
    df.index = range(0,len(df.index))
    #calculate real values for day degree and annual rainfall
    df=calculate_day_degree(df)
    df=calculate_annual_rainfall(df)
    return df

def read_weather_data_grid10_day_right_x_y(x,y,own_x,own_y,loc,user,passwd):
    #database query returns the list of found database entries 
    ls = location_to_fmidata_grid10_day_right_x_y(x,y,user,passwd)
    #add date related columns SpaFHy needs
    ls1 = add_aika_vuosi_kk_paiva_grid10_day(ls)
    #print('Weather',ls1[-1])
    #create data frame 
    df = pd.DataFrame(ls1)
    #add rest of the columns SpaFHy needs with default values
    df = add_fill_missing_column(df,0,'OmaTunniste',1)
    df = add_fill_missing_column(df,3,'Kunta',loc)
    df = add_fill_missing_column(df,len(df.columns),'lamposumma_v',0)
    df = add_fill_missing_column(df,len(df.columns),'rainfall_v',0)
    df.columns=ildatacols
    df['OmaIta']=own_x
    df['OmaPohjoinen']=own_y
    #calculate real values for day degree and annual rainfall
    df.sort_values(by='aika',kind='mergesort',inplace=True)
    df.index = range(0,len(df.index))
    #print("After sorting checking data")
    #cfmi.check_fmidata(df)
    #print("Done")
    df=calculate_day_degree(df)
    df=calculate_annual_rainfall(df)
    return df
    
#Vihti research site: (x=3356593 AND y=6703104)
def test(x=3356593,y=6703104):
    ls = location_to_ildata(x,y)
    ls1 = add_aika_vuosi_kk_paiva(ls)
    df = pd.DataFrame(ls1)
    df = add_fill_missing_column(df,0,'OmaTunniste',1)
    df = add_fill_missing_column(df,3,'Kunta','Vihti')
    df = add_fill_missing_column(df,4,'siteid',9999)
    df = add_fill_missing_column(df,len(df.columns),'lamposumma_v',0)
    df = add_fill_missing_column(df,len(df.columns),'rainfall_v',0)
    df.columns=ildatacols
    df['OmaIta']=x
    df['OmaPohjoinen']=y
    df=calculate_day_degree(df)
    df=calculate_annual_rainfall(df)
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-e',type=int,dest='e',help='East coordinate, KKJ(YKJ) (Finnish National Coordinate System)')
    parser.add_argument('-n',type=int,dest='n',help='North coordinate, KKJ(YKJ) (Finnish National Coordinate System)')
    parser.add_argument('-l',type=str,dest='l',default='xxxx',help='Location name')
    parser.add_argument('-f',type=str,dest='f',help='Output file name (type CSV)')
    parser.add_argument('-u','--user',type=str,dest='u',help='User name (Luke network) for weather database')
    parser.add_argument('-p','--passwd',type=str,dest='p',help='Password for weather database')
    parser.add_argument('-t','--table',type=str,dest='t',default='grid_day',help="Weather database table")
    args=parser.parse_args()
    #print(args.x,args.y,args.f,args.u,args.p)
    if args.t=='grid_day':
        read_and_write_weather_data(args.e,args.n,args.l,args.f,args.u,args.p)
    else:
        read_and_write_weather_data_grid10_day(args.e,args.n,args.l,args.f,args.u,args.p)
