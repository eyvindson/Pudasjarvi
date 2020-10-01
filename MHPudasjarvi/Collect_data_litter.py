import os
import argparse
import glob
import pathlib
import numpy as np
from scipy import interpolate
import pandas as pd
import subprocess
from subprocess import PIPE, run
from os import listdir
from os.path import isfile, join
import multiprocessing as mp
import pandas
import psutil
import time
import os
from pathlib import Path
from multiprocessing import Pool


def read_co2_file(file_name):
    df = pd.read_csv(file_name,sep=',')
    return df

def extract_grid_id(file_name):
    p = pathlib.Path(file_name)
    name = p.name
    name_str = str(name)
    name_ls = name_str.split('_')
    #First part is the ID
    grid_id = name_ls[0]
    if name_ls[5] == "litter.csv":
        osuus = 1
    else:
        osuus = name_ls[5]
    return grid_id, osuus

def add_id(file_name,df):
    grid_id,osuus = extract_grid_id(file_name)
    df.insert(0,'ID',grid_id)
    df.insert(1,'Osuus',osuus)
    return df

def concat_df(df_new,df):
    df = pd.concat([df,df_new])
    return df

def write_df(df,file_name,mode):
    if mode == "a":
        df.to_csv(file_name,sep=';',mode=mode,index = False,header = False)
    else:
        df.to_csv(file_name,sep=';',mode=mode,index = False)

  
mypath = "/scratch/project_2003225/GIT/Yasso/Litter_output_files_yasso/"#PUDASJARVI_P30061_2010_2060/"
ls_all = [f for f in listdir(mypath) if isfile(join(mypath, f)) if f.endswith(".csv")]

fname = ls_all[0]
ls_all.pop(0)
df_final_co2 = read_co2_file(mypath+fname)
df_final_co2 = df_final_co2.drop(["Unnamed: 0"],axis = 1)
df_final_co2 = add_id(fname,df_final_co2)
write_df(df_final_co2,'/scratch/project_2003225/temp_litter.txt','w')
ll = 0
for fname in ls_all:#[0:100]:
    df_co2 = read_co2_file(mypath+fname)
    df_co2 = df_co2.drop(["Unnamed: 0"],axis = 1)
    df_co2 = add_id(fname,df_co2)
    write_df(df_co2,'/scratch/project_2003225/temp_litter.txt','a')