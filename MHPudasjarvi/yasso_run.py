import os
import argparse
import glob
import pathlib
import numpy as np
from scipy import interpolate
import pandas as pd

def read_co2_file(file_name):
    df = pd.read_csv(file_name,sep=' ')
    return df

def extract_grid_id(file_name):
    p = pathlib.Path(file_name)
    name = p.name
    name_str = str(name)
    name_ls = name_str.split('_')
    #First part is the ID
    grid_id = name_ls[0]
    return grid_id

def add_id(file_name,df):
    grid_id = extract_grid_id(file_name)
    df.insert(0,'ID',grid_id)
    return df

def concat_df(df_new,df):
    df = pd.concat([df,df_new])
    return df

def write_df(df,file_name,mode):
    df.to_csv(file_name,sep=' ',mode=mode)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f',type=str,dest='f',help='CO2_result file after single yasso run')
    parser.add_argument('-o',type=str,dest='o',help='Combined output file after all yasso runs')
    args = parser.parse_args()
    if args.f == None:
        print("No CO2result input file name")
        quit()
    if args.o == None:
        print("No final CO2result output file name")
    args = parser.parse_args()
    ls = ['794812_PUDASJARVI_P30061_SKE4_4.txt','832917_PUDASJARVI_P30061_SKE4_4.txt','807276_PUDASJARVI_P30063_SKE4_4.txt']
    fname = ls[0]
    ls.pop(0)
    #run yasso here once
    df_final_co2 = read_co2_file(args.f)
    df_final_co2 = add_id(fname,df_final_co2)
    for fname in ls:
        #run yasso here
        df_co2 = read_co2_file(args.f)
        df_co2 = add_id(fname,df_co2)
        df_final_co2 = concat_df(df_co2,df_final_co2)
    write_df(df_final_co2,args.o,'w')
