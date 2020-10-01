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
    parser.add_argument('--f',type=str,dest='f',help='CO2_result file after single yasso run')
    parser.add_argument('--o',type=str,dest='o',help='Combined output file after all yasso runs')
    parser.add_argument('--s',type=int,dest='s',help='Where to start analyses -- from 0 to number of plots in folder')
    #parser.add_argument('--e',type=int,dest='e',help='Where to end analyses -- from 0 to number of plots in folder')
    args = parser.parse_args()
    if args.f == None:
        print("Input folder name")
        quit()
    if args.o == None:
        print("Output folder name")
    if args.s == None:
        print("No start number")
        quit()
    #if args.e == None:
    #    print("No end number")
    #    quit()
        
        
    
    args = parser.parse_args()

    os.chdir('/scratch/project_2003225/GIT/MELASoil/R/')

    outfolder = args.o#'/scratch/project_2003225/GIT/Yasso/Output_files_Yasso/'

    
    mypath = args.f#"/scratch/project_2003225/GIT/Yasso/Input_files_Yasso/PUDASJARVI_P30061_2010_2060/"
    ls_all = [f for f in listdir(mypath) if isfile(join(mypath, f)) if f.endswith(".txt")]

    #ls = ls[0:5]#['794812_PUDASJARVI_P30061_SKE4_4.txt','832917_PUDASJARVI_P30061_SKE4_4.txt','807276_PUDASJARVI_P30063_SKE4_4.txt']
    folder = args.f#'/scratch/project_2003225/GIT/Yasso/Input_files_Yasso/PUDASJARVI_P30061_2010_2060/'
    ls = ls_all[args.s:args.s+200]#interval of 200.
    fname = ls[0]
    ls.pop(0)
    #run yasso here once
    #subprocess.call(['Rscript', '/scratch/project_2003225/GIT/MELASoil/R/MELAToC_main.R', fname, folder])
    command = ['Rscript', '/scratch/project_2003225/GIT/MELASoil/R/MELAToC_main.R', fname, folder,outfolder]
    run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True)

    #Edit to move outside of R loop
    #df_final_co2 = read_co2_file(outfolder+fname[:-4]+"_out.csv")
    #df_final_co2 = add_id(fname,df_final_co2)


    #cmds_list = [['./bunzip2_file.py', file_name] for file_name in f_list]
    cmds_list =[['Rscript', '/scratch/project_2003225/GIT/MELASoil/R/MELAToC_main.R', fname, folder,outfolder] for fname in ls]
    procs_list = [run(cmd, stdout=PIPE, stderr=PIPE) for cmd in cmds_list]
    i = 0
    for proc in procs_list:
        proc
        print(i)
        i=i+1

    #for fname in ls:

    #    command = ['Rscript', '/scratch/project_2003225/GIT/MELASoil/R/MELAToC_main.R', fname, folder,outfolder]
    #    run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True)

    #    df_co2 = read_co2_file(outfolder+fname[:-4]+"_out.csv")
    #    df_co2 = add_id(fname,df_co2)
    #    df_final_co2 = concat_df(df_co2,df_final_co2)
    #write_df(df_final_co2,'/scratch/project_2003225/temp2.txt','w')
