import argparse
import pickle
import pathlib
import pandas as pd
import SpaFHyWeather.fmidata as fmi
import createsimulationdata as csd
coord_dict = {}
coord_ls = []

#Tried this with dataframe.apply() but no faster
def closest_point(x,ls):
    e = x['YKJx']
    n = x['YKJy']
    s = x['Segment']
    dp = fmi.find_closest_point2(e,n,ls)
    p1 = dp[1][0]
    p2 = dp[1][1]
    if s % 100 == 0:
        print(len(ls),s,e,n,p1,p2)
    coord_dict[(e,n)]=(p1,p2)
    #coord_ls.append([e,n,p1,p2])
    return x

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f',type=str,dest='f',help='Segment coordinate input file')
    parser.add_argument('-d',type=str,dest='d',help='Directory for segment pseudotype file')
    parser.add_argument('-o',type=str,dest='o',help='Weather data closest point pickle file')
    parser.add_argument('-u',type=str,dest='u',help='User name')
    parser.add_argument('-p',type=str,dest='p',help='Password')
    parser.add_argument('-q',type=float,dest='q',help='Segment pseudotype')
    args = parser.parse_args()
    region = args.d+'/'
    fout = region+args.o
    df_all_segments = csd.read_segment_file(args.f)
    df_pseudotype = csd.read_and_select_pseudotype(region+'segmaj1.pseudotype.hits',args.q)
    df_segments =  csd.merge_by_segment_number(df_pseudotype,df_all_segments)
    print(df_segments.head())
    daydaydate="'"+'2014-01-01'+"'"+"::date"
    daydaydate2="'"+'2014-01-01'+"'"+"::date"
    #Fetch first a couple of days of data, grid10_day table seems to have
    print("Fetching all coordinate points for a single day")
    ls1 = fmi.fetch_fmidata(host="lukedb1.ns.luke.fi",database='weather',user=args.u,password=args.p,
                            query="SELECT longitude,latitude FROM grid10_day WHERE pvm BETWEEN "+daydaydate+" AND "+daydaydate2)
    print('Number of coordinates',len(ls1))
    #df.apply(closest_point,axis=1,args=(ls1,))
    weather_id_dict = csd.CoordinateDict()
    coord_dict = {}
    coord_ls = []
    index_length=len(df_segments.index)
    print('Number of pseudotype',args.q,'lines',index_length)
    for i in range(0,index_length):
        coord = csd.select_YKJ(df_segments,i)
        dp = fmi.find_closest_point2(coord[0],coord[1],ls1)
        p1 = dp[1][0]
        p2 = dp[1][1]
        coord_dict[(coord[0],coord[1])]=(p1,p2)
        coord_ls.append([coord[0],coord[1],p1,p2])
        weather_id_dict.insert(p1,p2)
        if i % 1000 == 0:
            print('Line',i)
    with open(fout,'wb') as f:
        pickle.dump(coord_dict,f, pickle.HIGHEST_PROTOCOL)
    p = pathlib.Path(region+args.o)
    stem = p.stem
    new_name = stem+'weather_id.pickle'
    p_new = p.with_name(new_name)
    with open(str(p_new),'wb') as f:
        pickle.dump(weather_id_dict,f, pickle.HIGHEST_PROTOCOL)
    if len(coord_ls) <= 1000000:
        file_excel = p.with_suffix('.xlsx')
        fmi.write_ls_to_excel(str(file_excel),coord_ls,['East','North','ClosestEast','ClosestNorth'])
    else:
        print('Too many lines for excel',len(coord_ls))
