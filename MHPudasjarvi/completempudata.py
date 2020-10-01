import os
import argparse
import glob
import pathlib
import numpy as np
from scipy import interpolate
import pandas as pd

waste1="Hukkapuu massat (1000 t DM)\n"
waste2="Hukkapuu massat (1000 t DM);;;;;;\n"
def complete_areas(df_area,col_name,value_ls):
    df_area.insert(1,col_name,value_ls)
    return df_area

def complete_volumes(df_volume,year,years):
    y0 = df_volume.loc[0:0,1:]
    y1 = df_volume.loc[1:1,1:]
    y2 = df_volume.loc[2:2,1:]
    y3 = df_volume.loc[3:3,1:]
    x = years
    f0 = interpolate.interp1d(x, y0, fill_value = "extrapolate")
    f1 = interpolate.interp1d(x, y1, fill_value = "extrapolate")
    f2 = interpolate.interp1d(x, y2, fill_value = "extrapolate")
    f3 = interpolate.interp1d(x, y3, fill_value = "extrapolate")
    [v0] = f0(year)  
    [v1] = f1(year)
    [v2] = f2(year)
    [v3] = f3(year)
    
    df_volume.insert(1,year,0)
    #Do not let volume go below zero
    df_volume.loc[0:0,year:year] = np.max(v0,0)
    df_volume.loc[1:1,year:year] = np.max(v1,0)
    df_volume.loc[2:2,year:year] = np.max(v2,0)
    df_volume.loc[3:3,year:year] = np.max(v3,0)
    return df_volume

def complete_zeros(df_dm,year):
    """So far we use zero for waste wood 
       and the second part (North Finland) is zero
    """
    df_dm.insert(1,year,0)
    return df_dm

def complete_one_file(df1,df2,df3,df4,df5,df6,df7,df8):
    df1 = complete_areas(df1,2015,[2015,1.0,0.0,0.0,0.0])
    df1 = complete_areas(df1,2010,[2010,1.0,0.0,0.0,0.0])
    df2 = complete_volumes(df2,2015,list(range(2020,2070,10)))
    df2 = complete_volumes(df2,2010,list(range(2020,2070,10)))
    df3 = complete_volumes(df3,2015,list(range(2020,2060,10)))
    df3 = complete_volumes(df3,2010,list(range(2020,2060,10)))
    df4 = complete_zeros(df4,2015)
    df4 = complete_zeros(df4,2010)
    df5 = complete_zeros(df5,2015)

    df5 = complete_zeros(df5,2010)
    df6 = complete_zeros(df6,2015)
    df6 = complete_zeros(df6,2010)
    df7 = complete_zeros(df7,2015)
    df7 = complete_zeros(df7,2010)
    df8 = complete_zeros(df8,2015)
    df8 = complete_zeros(df8,2010)
    return (df1,df2,df3,df4,df5,df6,df7,df8)
                             
def read_yasso_mpu_data(file_name):
    #Areas, 5 columns
    df1 = pd.read_csv(file_name,sep=';',decimal=',',header=None,nrows=5)
    #Volume, 5 columns
    df2 = pd.read_csv(file_name,sep=';',decimal=',',header=None,skiprows=6,nrows=10)
    #Natural thinning, 4 columns
    df3 = pd.read_csv(file_name,sep=';',decimal=',',header=None,skiprows=17,nrows=10)
    #Wastewood, 4 columns, 1(title)+45(data) rows,no data after title
    #Skip title and add explicitely when writing 
    df4 =  pd.read_csv(file_name,sep=';',decimal=',',header=None,skiprows=30,nrows=45)
    #Now the same for the second part
    #Areas
    df5 = pd.read_csv(file_name,sep=';',decimal=',',header=None,skiprows=75,nrows=4)
    #Volume
    df6 = pd.read_csv(file_name,sep=';',decimal=',',header=None,skiprows=81,nrows=10)
    #Natural thinning
    df7 = pd.read_csv(file_name,sep=';',decimal=',',header=None,skiprows=92,nrows=10)
    #Wastewood
    #Skip the title and add explicitely with ";;;;;;" as with Risto
    df8 = pd.read_csv(file_name,sep=';',decimal=',',header=None,skiprows=104,nrows=45)
    return (df1,df2,df3,df4,df5,df6,df7,df8)

def complete_yasso_mpu_file(file_name):
    (df1,df2,df3,df4,df5,df6,df7,df8) = read_yasso_mpu_data(file_name)
    (df1,df2,df3,df4,df5,df6,df7,df8) = complete_one_file(df1,df2,df3,df4,df5,df6,df7,df8)
    return (df1,df2,df3,df4,df5,df6,df7,df8)

def write_new_yasso_mpu_file(df1,df2,df3,df4,df5,df6,df7,df8,f):
    df1 = df1.astype(np.str)
    df_0 = df1.loc[0:0,:]
    df_0 = df_0.apply(lambda x: x.str.replace('\.0',''))
    df_0 = df_0.apply(lambda x: x.str.replace('nan',''))
    df_str = df1.astype(np.str)
    df_str = df_str.apply(lambda x:x.str.replace('.',','))
    df_str = df_str.apply(lambda x:x.str.replace('nan',''))
    df_concat = pd.concat([df_0,df_str.loc[1:,:]])
    df_concat.to_csv(f,sep=';',index=False,header=False,mode='a')
    f.write(os.linesep)
    df2_str = df2.astype(np.str)
    df2_str = df2_str.apply(lambda x:x.str.replace('.',','))
    df2_str.to_csv(f,sep=';',index=False,header=False,mode='a')
    f.write(os.linesep)
    df3_str = df3.astype(np.str)
    df3_str = df3_str.apply(lambda x:x.str.replace('.',','))
    df3_str.to_csv(f,sep=';',index=False,header=False,mode='a')
    f.write(os.linesep)
    f.write(waste1)
    df4_str = df4.astype(np.str)
    df4_str = df4_str.apply(lambda x:x.str.replace('.',','))
    df4_str.to_csv(f,sep=';',index=False,header=False,mode='a')
    f.write(os.linesep)
    df5_str = df5.astype(np.str)
    df5_str = df5_str.apply(lambda x:x.str.replace('.',','))
    df5_str.to_csv(f,sep=';',index=False,header=False,mode='a')
    f.write(os.linesep)
    df6_str = df6.astype(np.str)
    df6_str = df6_str.apply(lambda x:x.str.replace('.',','))
    df6_str.to_csv(f,sep=';',index=False,header=False,mode='a')
    f.write(os.linesep)
    df7_str = df7.astype(np.str)
    df7_str = df7_str.apply(lambda x:x.str.replace('.',','))
    df7_str.to_csv(f,sep=';',index=False,header=False,mode='a')
    f.write(os.linesep)
    f.write(waste2)
    df8_str = df8.astype(np.str)
    df8_str = df8_str.apply(lambda x:x.str.replace('.',','))
    df8_str.to_csv(f,sep=';',index=False,header=False,mode='a')
    
def update_yasso_mpu_file(input_file_name:str,output_dir:str):
     (df1,df2,df3,df4,df5,df6,df7,df8) = complete_yasso_mpu_file(input_file_name)
     p = pathlib.Path(input_file_name)
     name = p.name
     dir_out = pathlib.Path(output_dir)
     output_file = dir_out.joinpath(name)
     if not dir_out.exists():
         print(dir_out)
         dir_out.mkdir()
     f = open(str(output_file),'w')
     write_new_yasso_mpu_file(df1,df2,df3,df4,df5,df6,df7,df8,f)
     f.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f',type=str,dest='f',help='MELA Yasso MPU data files (wild card search)')
    parser.add_argument('-o',type=str,dest='o',help='Output directory for updated Yasso initial files')
    args = parser.parse_args()
    if args.f == None:
        print("No input directory")
        quit()
    if args.o == None:
        print("No output directory")
        quit()
    
    ls = glob.glob(args.f)
    print("No FILES", len(ls))
    f = ls[0]
    for f in ls:
        update_yasso_mpu_file(f,args.o)
