import argparse
import datetime
import pathlib
import glob
import numpy as np
import pandas as pd

def check_weather_file(file_name,sep=';'):
    df = pd.read_csv(file_name,sep=sep)
    rows = len(df.index)
    delta_dict=dict()
    same_day_dict=dict()
    previous_day_dict=dict()
    print('Checking duplicate rows')
    series = df[df.duplicated()]
    print('Number of duplicate rows',len(series))
    if len(series) > 0:
        print(series)
    for row in range(rows-1):
        r1 = df.loc[row,['vuosi','kk','paiva']]
        r2 = df.loc[row+1,['vuosi','kk','paiva']]
        d1 = datetime.date(r1.vuosi,r1.kk,r1.paiva)
        d2 = datetime.date(r2.vuosi,r2.kk,r2.paiva)
        delta = d2-d1
        days = delta.days
        if days > 1:
            print('Missing days',d1,d2)
            delta_dict[row]=([r1.paiva,r1.kk,r1.vuosi,r2.paiva,r2.kk,r2.vuosi,days])
        if days == 0:
            print('Two same dates error',d1,d2)
            same_day_dict[row]=(d1,d2)
        if days < 0:
            print('Earlier date error',d1,d2)
            previous_day_dict[row]=(d1,d2)
    s = df['aika']-df['aika'].shift(1)
    for i in range(len(s)):
        if s[i] != 1:
            print(i,s[i])
    return (delta_dict,same_day_dict,previous_day_dict)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f',type=str,dest='f',help='Segment coordinate input files (with wild card search use double quotes)')
    args = parser.parse_args()
    file_ls = glob.glob(args.f)
    for file_name in sorted(file_ls):
        print('Checking File',file_name)
        (delta_dict,same_day_dict,previous_day_dict)=check_weather_file(file_name)
        df1 = pd.DataFrame.from_dict(delta_dict)
        df2 = pd.DataFrame.from_dict(same_day_dict)
        df3 = pd.DataFrame.from_dict(previous_day_dict)
        p = pathlib.Path(file_name)
        stem=pathlib.Path(file_name).stem
        parent=pathlib.Path(file_name).parent
        file_excel = str(parent)+'/'+str(stem)+'Checked.xlsx'
        excel_writer = pd.ExcelWriter(file_excel,engine='xlsxwriter')
        print('Printing File',file_excel)
        df1.to_excel(excel_writer,sheet_name='Missing days')
        df2.to_excel(excel_writer,sheet_name='Double days')
        df3.to_excel(excel_writer,sheet_name='False order days')
        excel_writer.save()
        
