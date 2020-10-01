import os
import argparse
import datetime
import pickle
import pathlib
import glob
import numpy as np
import scipy.interpolate as interpolate
import pandas as pd 
import fmidata as fmi
import checkfmidata as cfmi

date_header =['vuosi','kk','paiva']
weather_header = ['t_mean','t_max','t_min','rainfall','radiation','hpa','lamposumma_v','rainfall_v']
segment_file_header = ['Segment','Pixels','ETRSx','ETRSy','YKJx','YKJy',
                       'Latitude','Longitude']
pseudotype_file_header = ["Segment","Pixels","Pseudotype"]
gis_file_header1 = "ncols 1\n"
gis_file_header2 = "xllcorner 0\nyllcorner 0\ncellsize 16\nNODATA_value -9999\n"
segment_file_missing_file_columns=['Weather_id','Line','Day1','Month1','Year1','Day2','Month2','Year2','DeltaDays']
weather_file_prefix = "Weather_id_"
regions_ls = ['amaa','itasuomi','lansisuomi','pkarjala','lappi']
class WeatherId:
    id=0
    def __init__(self):
        self.id = -1
    def next_id(self):
        self.id = self.id+1
        return self.id
    def reset(self):
        self.id = -1

class CoordinateDict:
    def __init__(self):
        self.current_line=0
        self.coord_dict={}
        self.weather_id = WeatherId()
    def insert(self,longitude:int,latitude:int):
        if (longitude,latitude) in self.coord_dict:
            #No insertion, return False
            return False
        else:
            #New value, return True
            value= self.weather_id.next_id()
            self.coord_dict[(longitude,latitude)]=value
            return True
    #Returns exception if value not found
    def lookup(self,longitude:int,latitude:int):
        return self.coord_dict[(longitude,latitude)]
    #Like lookup but return None if value not found
    def get(self,longitude:int,latitude:int):
        return self.coord_dict.get(longitude,latitude)
    def reset_weather_id(self):
        self.weather_id.reset()
    def reset(self):
        self.weather_id.reset()
        self.coord_dict={}
        self.current_line=0
    def update_current_line(self):
        self.current_line=self.current_line+1
    
class ForcingId:
    id=-1
    def next_id():
        ForcingId.id = ForcingId.id+1
        return ForcingId.id
    def reset():
        ForcingId.id = -1

def write_gis_header(fname:str,header:str):
    f=open(fname,'w')
    f.write(header)
    return f

def create_forcing_id():
    id=ForcingId.next_id()
    return id

def read_segment_file(file_name:str):
    df = pd.read_csv(file_name,names=segment_file_header,delim_whitespace=True,header=None)
    return df

def select_YKJ(df,row:int):
    return df.loc[row,['YKJx','YKJy']].values.astype(int)


def foliage_mass_to_LAI(df_pseudotype,fm_file:str,canopy_closure_file:str,species:int,out_file:str,unit_change1:float,unit_change2:float):
    #Mg to kg
    pine_limit = 10.0*unit_change1
    spruce_limit = 20.0*unit_change1
    deciduous_limit = 5.0*unit_change1
    cf_limit=0
    fm_limit=0
    cf_specific1=0
    cf_specific2=0
    if species == 1:
        fm_limit=pine_limit
        cf_limit=0.5
        cf_specific1=12.5
        cf_specific2=15.0
    if species == 2:
        fm_limit =spruce_limit
        cf_limit=0.7
        cf_specific1=8.0
        cf_specific2=11.0
    if species > 2:
        fm_limit=deciduous_limit
        cf_limit=0.5
        cf_specific1=28.0
        cf_specific2=0
    #Give column names, easier to follow
    df_fm_all = pd.read_csv(fm_file,names=['Segment','Pixels','FM'],delim_whitespace=True,header=None,dtype={'FM':float})
    #print('ALL',df_fm_all.head())
    df_fm = merge_by_segment_number(df_pseudotype,df_fm_all)
    #print('FM',df_fm.head())
    #It is now assumed that canopy closure file has been written out
    df_cf = pd.read_csv(canopy_closure_file,names=['CF'],delim_whitespace=True,header=None,skiprows=6)
    #print('CF',df_cf.head())
    #10kg to kg
    #print("FM kg")
    df_fm['FM']=df_fm['FM'].multiply(unit_change2)
    #Limit max foliage mass
    #print("Limit fm")
    df_fm['FM']=df_fm['FM'].apply(lambda x:min(fm_limit,x))
    #concat df_cf and df_fm
    #df_cf is one column and df_fm is three column wide
    df_cf_fm_lai = pd.concat([df_cf,df_fm],axis=1)
    df_cf_fm_lai['LAI']=0
    #Hectars to m2 with 10000
    #print("LAI")
    #Pine or Spruce 
    if species == 1 or species ==2:
        df_cf_fm_lai['LAI'] = df_cf_fm_lai.apply(lambda row: 0.5*row['FM']*cf_specific1/10000.0 if row['CF'] < cf_limit else row['LAI'],axis=1)
        df_cf_fm_lai['LAI'] = df_cf_fm_lai.apply(lambda row: 0.5*row['FM']*(0.5*cf_specific1+0.5*cf_specific2)/10000.0 if row['CF'] >= cf_limit else row['LAI'],axis=1)
     #Deciduous
    else:
        df_cf_fm_lai['LAI'] = df_cf_fm_lai.apply(lambda row: (0.5*row['FM']*cf_specific1)/10000.0,axis=1)
    #print("To file")
    #print(df_cf_fm_lai.head())
    df_cf_fm_lai.to_csv(out_file,mode='a',columns=['LAI'],header=False,index=False)
    p=pathlib.Path(out_file)
    stem=str(p.stem)
    parent = str(p.parent)
    all_data_file = parent+'/'+stem+'LAIData.csv'
    df_cf_fm_lai.to_csv(all_data_file,sep = ' ',index=False)
    excel_file= parent+'/'+stem+'LAIData.xlsx'
    if len(df_cf_fm_lai.index) <= 1048576:
        fmi.write_df_to_excel(excel_file,'LAIData',df_cf_fm_lai)
    else:
        print('Too many rows for excel',len(df_cf_fm_lai.index),'(Max 1048576)')
    
def soiltype_id(df_pseudotype,input_file:str,output_file:str,gis_file_header:str):
    print("Soiltype")
    f=write_gis_header(output_file,gis_file_header)
    f.close()
    df_all = pd.read_csv(input_file,names=['Segment','Pixels','ID'],delim_whitespace=True,header=None)
    df = merge_by_segment_number(df_pseudotype,df_all)
    print(df.head())
    #print(len(df.index))
    #print(output_file)
    df.to_csv(output_file,mode='a',sep=' ',columns=['ID'],header=False,index=False)

def generate_soil_id(regions_ls:list):
    for region in regions_ls:
        df_pseudotype = read_and_select_pseudotype(region+'/'+'segmaj1.pseudotype.hits',400.0)
        soil_id_file = region+'/'+'segmaj1.paatyyppi_vmi1x_1216.hits'
        pine_file = glob.glob(region+'/'+'out_'+'*'+'_segmean_manty_vmi1x_1216.txt')
        spruce_file = glob.glob(region+'/'+'out_'+'*'+'_segmean_kuusi_vmi1x_1216.txt')
        birch_file = glob.glob(region+'/'+'out_'+'*'+'_segmean_koivu_vmi1x_1216.txt')
        decid_file = glob.glob(region+'/'+'out_'+'*'+'_segmean_muulp_vmi1x_1216.txt')
        veg_class_file = region+'/'+'segmaj1.kasvupaikka_vmi1x_1216.hits'
        output_file=region+'_soil_id.dat'
        print(soil_id_file, pine_file[0],spruce_file[0],birch_file[0],decid_file[0],veg_class_file)
        generate_soiltype_id_stem_volume_veg_class(df_pseudotype,soil_id_file,pine_file[0],spruce_file[0],birch_file[0],decid_file[0],
                                                   veg_class_file,output_file,"")
        print("Done")
def generate_soiltype_id_stem_volume_veg_class(df_pseudotype,soil_id_file:str,pine_file:str,spruce_file:str,birch_file:str,decid_file:str,
                                               veg_class_file:str,output_file:str,gis_file_header:str):
    print("Soiltype")
    nrows = len(df_pseudotype.index)
    gis_file_header=gis_file_header1+"nrows "+str(nrows)+"\n"+gis_file_header2
    f=write_gis_header(output_file,gis_file_header)
    f.close()
    df_soil_id_all = pd.read_csv(soil_id_file,names=['Segment','Pixels','ID'],delim_whitespace=True,header=None)
    df = merge_by_segment_number(df_pseudotype,df_soil_id_all)
    df_pine = pd.read_csv(pine_file,names=['Segment','Pixels','StemVolPine'],delim_whitespace=True,header=None)
    df_spruce = pd.read_csv(spruce_file,names=['Segment','Pixels','StemVolSpruce'],delim_whitespace=True,header=None)
    df_birch = pd.read_csv(birch_file,names=['Segment','Pixels','StemVolBirch'],delim_whitespace=True,header=None)
    df_decid = pd.read_csv(decid_file,names=['Segment','Pixels','StemVolDecid'],delim_whitespace=True,header=None)
    df_vegclass = pd.read_csv(veg_class_file,names=['Segment','Pixels','VegClass'],delim_whitespace=True,header=None)
    df1 = merge_by_segment_number(df,df_pine)
    df2 = merge_by_segment_number(df1,df_spruce)
    df3 = merge_by_segment_number(df2,df_birch)
    df4 = merge_by_segment_number(df3,df_decid)
    df_soil_id = merge_by_segment_number(df4,df_vegclass)
    df_soil_id['DecidVolTot'] = df_soil_id['StemVolBirch']+df_soil_id['StemVolDecid']
    df_soil_id['Soil_ID']=0
    #This is straightforward for Soil_ID
    #VegClass > 4, i.e. >=5
    df_soil_id['Soil_ID']=np.where(df_soil_id['VegClass']>4,1,df_soil_id['Soil_ID'])
    #VegClass <5, i.e. <=4
    df_soil_id['Soil_ID']=np.where(df_soil_id['VegClass']<5,2,df_soil_id['Soil_ID'])
    df_soil_id['MainSpecies']=0
    #Final part depends on main tree species and vegetation class (kasvillisuusluokka)
    #Assign main tree species one at a time, easier to read. Note the increasingly '>' comparison.
    #Pine wins
    df_p = df_soil_id.assign(MainSpecies=np.where((df_soil_id['StemVolPine'] >= df_soil_id['StemVolSpruce']) & (df_soil_id['StemVolPine'] >= df_soil_id['DecidVolTot']),1,df_soil_id['MainSpecies']))
    #Spruce wins
    df_ps = df_p.assign(MainSpecies=np.where((df_p['StemVolSpruce'] > df_p['StemVolPine']) & (df_p['StemVolSpruce'] >= df_p['DecidVolTot']),2,df_p['MainSpecies']))
    #Deciduous wins 
    df_psd = df_ps.assign(MainSpecies=np.where((df_ps['DecidVolTot'] > df_ps['StemVolPine']) & (df_ps['DecidVolTot'] > df_ps['StemVolSpruce']),3,df_ps['MainSpecies']))
    #Finally the Soil_ID: (VegClass == 3 or VegClass == 4) and MainSpecies<3 (i.e., MainSpecies <=2) 
    df_soil_id = df_psd.assign(Soil_ID=np.where(((df_psd['VegClass'] == 3) | (df_psd['VegClass'] == 4)) & (df_psd['MainSpecies'] < 3),1,df_psd['Soil_ID']))
    df_soil_id.to_csv(output_file,mode='a',sep=' ',columns=['Soil_ID'],header=False,index=False)
    print(len(df_soil_id.index))

def dominant_height_hc(df_pseudotype,input_file:str,output_file:str,gis_file_header:str,unit_change:float):
    f=write_gis_header(output_file,gis_file_header)
    f.close()
    df_all = pd.read_csv(input_file,names=['Segment','Pixels','Hc'],delim_whitespace=True,header=None)
    #print('Shape',np.shape(df))
    #print(df.head())
    #change dm to m
    df = merge_by_segment_number(df_pseudotype,df_all)
    df['Hc']=df['Hc'].multiply(unit_change)
    #print('Shape',np.shape(df))
    #print(df.head())
    #print(output_file)
    df.to_csv(output_file,mode='a',sep=' ',columns=['Hc'],header=False,index=False)

def canopy_cover_cf(df_pseudotype,input_file:str,output_file:str,gis_file_header:str,unit_change:float):
    f = write_gis_header(output_file,gis_file_header)
    f.close()
    df_all = pd.read_csv(input_file,names=['Segment','Pixels','Cf'],delim_whitespace=True,header=None)
    df = merge_by_segment_number(df_pseudotype,df_all)
    #change percent [0:100] to decimal [0:1]representation
    df['Cf']=df['Cf'].multiply(unit_change)
    df.to_csv(output_file,mode='a',columns=['Cf'],header=False,index=False)

def ditch_data_as_constant(output_file:str,gis_file_header:str,c:float,nrows:int):
    f = write_gis_header(output_file,gis_file_header)
    for i in range(0,nrows):
        f.write(str(c)+'\n')
    f.close()
    
def read_and_select_pseudotype(input_file,pseudo_type:int):
    df_pseudotype = pd.read_csv(input_file,names=pseudotype_file_header,delim_whitespace=True,header=None)
    df=df_pseudotype.loc[df_pseudotype["Pseudotype"]==pseudo_type]
    return df

def merge_by_segment_number(df_pseudotype,df_segments):
    #pseudotype has the line numbers we need, less that number of segments
    df = pd.merge(df_pseudotype,df_segments,how='left',on=['Segment'])
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f',type=str,dest='f',help='Segment coordinate input file')
    parser.add_argument('-p',type=str,dest='p',help='Output directory for parameters')
    parser.add_argument('-w',type=str,dest='w',help='Output directory for weather')
    parser.add_argument('-q',type=float,dest='q',help='Segment pseudotype needed')
    parser.add_argument('--pickle',type=str,dest='pickle',help='Closest points in pickle file')
    parser.add_argument('-r',type=str,dest='r',help='Region symbol')
    parser.add_argument('--user',type=str,dest='user',help='User name Luke weather database')
    parser.add_argument('--passwd',type=str,dest='passwd',help='Password Luke weather database')
    parser.add_argument('--cd',type=str,dest='cd',help='Coordinate dictionary pickle file\nto continue after database hangup')
    args = parser.parse_args()
    closest_points = None
    if args.pickle:
        print('Closest points from pickle file',args.pickle)
        closest_points = pd.read_pickle(args.pickle)
    if args.cd == None:
        print("No Cooradinate dictionary file, quitting program")
        quit()
    df_all_segments = read_segment_file(args.f)
    #print(gis_file_header)
    parameter_dir = args.p+'/'
    weather_dir = args.w+'/'
    p=pathlib.Path(args.p)
    region = str(p.parent)+'/'
    df_pseudotype = read_and_select_pseudotype(region+'segmaj1.pseudotype.hits',args.q)
    df_segments =  merge_by_segment_number(df_pseudotype,df_all_segments)
    index_length = len(df_segments.index)
    gis_file_header = gis_file_header1+"nrows "+str(index_length)+"\n"+gis_file_header2
    #Parameter files
    print("Writing out canopy cover",parameter_dir+'cf.dat')
    canopy_cover_cf(df_pseudotype,region+'out_'+args.r+'_segmean_latvuspeitto_vmi1x_1216.txt',parameter_dir+'cf.dat',gis_file_header,1/100.0)
    print("Writing out dominant height",parameter_dir+'hc.dat')
    dominant_height_hc(df_pseudotype,region+'out_'+args.r+'_segmean_keskipituus_vmi1x_1216.txt',parameter_dir+'hc.dat',gis_file_header,1.0/10.0)
    print("Writing out soiltype",parameter_dir+'soil_id.dat')
    soiltype_id(df_pseudotype,region+'segmaj1.paatyyppi_vmi1x_1216.hits',parameter_dir+'soil_id.dat',gis_file_header)
    print("Writing ditch depth",parameter_dir+'ditch_depth.dat',0.287)
    ditch_data_as_constant(parameter_dir+'ditch_depth.dat',gis_file_header,0.287,index_length)
    print("Writing ditch spacing",parameter_dir+'ditch_spacing.dat',36)
    ditch_data_as_constant(parameter_dir+'ditch_spacing.dat',gis_file_header,36.0,index_length)
    #File close mandatory: data frames for LAI will append data later
    lai_decid = write_gis_header(parameter_dir+'LAI_decid.dat',gis_file_header)
    lai_decid.close()
    lai_pine = write_gis_header(parameter_dir+'LAI_pine.dat',gis_file_header)
    lai_pine.close()
    lai_spruce = write_gis_header(parameter_dir+'LAI_spruce.dat',gis_file_header)
    lai_spruce.close()
    print("Calculating LAI")
    print("Manty")
    foliage_mass_to_LAI(df_pseudotype,region+'out_'+args.r+'_segmean_bm_manty_neulaset_vmi1x_1216.txt',parameter_dir+'cf.dat',
                        1,parameter_dir+'LAI_pine.dat',1000.0,10.0)
    print("Kuusi")
    foliage_mass_to_LAI(df_pseudotype,region+'out_'+args.r+'_segmean_bm_kuusi_neulaset_vmi1x_1216.txt',parameter_dir+'cf.dat',
                        2,parameter_dir+'LAI_spruce.dat',1000.0,10.0)
    print("Lehtipuu")
    foliage_mass_to_LAI(df_pseudotype,region+'out_'+args.r+'_segmean_bm_lehtip_neulaset_vmi1x_1216.txt',parameter_dir+'cf.dat',
                        3,parameter_dir+'LAI_decid.dat',1000.0,10.0)
    #Weather data
    print("Getting and correcting weather data")
    coordinate_dict = CoordinateDict()
    #Collect all missing dates here
    delta_date_ls=[]
    print('Max index',index_length)
    daydaydate1="'"+'2014-01-01'+"'"+"::date"
    daydaydate2="'"+'2014-01-01'+"'"+"::date"
    ls1 = []
    if not closest_points:
        print("Fetching all coordinate points for weather data for a single day")
        ls1 = fmi.fetch_fmidata(host="lukedb1.ns.luke.fi",database='weather',user=args.user,password=args.passwd,
                                query="SELECT longitude,latitude FROM grid10_day WHERE pvm BETWEEN "+daydaydate1+" AND "+daydaydate2)
        print('Number of coordinates as a list',len(ls1))
        ls1 = list(set(ls1))
        print('Number of coordinates as a set',len(ls1),ls1[0])
    #Coordinates and forcing id's will be written when found
    print('Starting writing out weather')
    e_coord=None
    n_coord=None
    forcing_id=None
    if os.path.isfile(args.cd):
        #After database hangup continue with help saved coordinate dict
        print("Using existing Coordinate dictionary file for weather data")
        with open(args.cd,'rb') as f_pickle:
            coordinate_dict = pickle.load(f_pickle)
        e_coord = open(parameter_dir+'Ecoord.dat','a')
        n_coord = open(parameter_dir+'Ncoord.dat','a')
        forcing_id = open(parameter_dir+'forcing_id.dat','a')
    else:
        #Start from the beginning
        print("Starting weather data from the beginning") 
        e_coord=write_gis_header(parameter_dir+'Ecoord.dat',gis_file_header)
        n_coord=write_gis_header(parameter_dir+'Ncoord.dat',gis_file_header)
        forcing_id=write_gis_header(parameter_dir+'forcing_id.dat',gis_file_header)
    print('Coordinate dictionary current line',coordinate_dict.current_line)
    for i in range(coordinate_dict.current_line,index_length):
        coord = select_YKJ(df_segments,i)
        #returns euclidian distance 'd' and the point '(x1,y1)' as (d,(x1,y1))
        p1=0
        p2=0
        if not closest_points:
            dp = fmi.find_closest_point2(coord[0],coord[1],ls1)
            p1 = dp[1][0]
            p2 = dp[1][1]
        else:
            (p1,p2) = closest_points[(coord[0],coord[1])]
        #Update coordinate dictionary
        if coordinate_dict.insert(p1,p2):
            print('Reading weather data',i,p1,p2,coord[0],coord[1])
            df_w = fmi.read_weather_data_grid10_day_right_x_y(p1,p2,coord[0],coord[1],args.r,args.user,args.passwd)
            #'longitude' and 'latitude' denote also the closest point
            coord_array = df_w.loc[0,['longitude','latitude']]
            print(coord[0],coord[1],p1,p2,coord_array[0],coord_array[1])
            #If a new coordinate point, create a new weather (i.e. forcing) file
            print('Checking and correcting weather data')
            (delta_dict,df_wc) = cfmi.check_and_correct_weather_data(df_w)
            weather_id = coordinate_dict.lookup(coord_array[0],coord_array[1])
            weather_file = weather_dir+weather_file_prefix+str(weather_id)+'.csv'
            df_wc.to_csv(weather_file,sep=';',index=False)
            #Write original file to excel
            fmi.write_df_to_excel(weather_dir+weather_file_prefix+str(weather_id)+'_orig.xlsx','FMIData',df_w)
            keys_ls = sorted(delta_dict.keys())
            #print('Printing missing dates')
            #print(keys_ls)
            for key in keys_ls:
                date_ls = delta_dict[key]
                date_ls = [weather_id,key]+date_ls
                delta_date_ls.append(date_ls)
                #print(key,date_ls)
            print('Done')
        #Update parameter files
        e_coord.write(str(coord[0])+'\n')
        n_coord.write(str(coord[1])+'\n')
        weather_id = coordinate_dict.lookup(p1,p2)
        forcing_id.write(str(weather_id)+'\n')
        #Keep track of the current line in case of failure with weather database
        #Take the saved coordinate_dict in use and continue from the current line
        coordinate_dict.update_current_line()
        #fpickle = open(region+'closest_point_weather_id_1.pickle','wb')
        with  open(str(args.cd),'wb') as fpickle:
            print("Writing Coordinate dictionary pickle file",str(args.cd))
            pickle.dump(coordinate_dict,fpickle, pickle.HIGHEST_PROTOCOL)
            print("Done")
    #Close parameter files, note thet after database hangup operating system will close files
    #Improvement: in case of database hangup catch psycopg2 DatabaseError and quit the program in more controller manner 
    print("Weather data done, closing files and print missing weather data to excel")
    e_coord.close()
    n_coord.close()
    forcing_id.close()
    lai_decid.close()
    lai_pine.close()
    lai_spruce.close()    
    #Write all missing dates to excel
    df_delta_dates = pd.DataFrame(delta_date_ls)
    df_delta_dates.columns=segment_file_missing_file_columns
    stem=pathlib.Path(args.f).stem
    parent=pathlib.Path(args.f).parent
    fmi.write_df_to_excel(str(parent)+'/'+weather_dir+str(stem)+'Missing.xlsx','Missing',df_delta_dates)
    
