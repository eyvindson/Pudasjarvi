import pandas as pd
import numpy as np
import csv
import numpy as np


#FUNCTIONS:
def foliage_mass_to_LAI(df_fm,df_cf,species):
    unit_change1= 1000
    unit_change2 = 1
    df_cf = df_cf.multiply(100)

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
        sp_name = "pine"
    if species == 2:
        fm_limit =spruce_limit
        cf_limit=0.7
        cf_specific1=8.0
        cf_specific2=11.0
        sp_name = "spruce"
    if species > 2:
        fm_limit=deciduous_limit
        cf_limit=0.5
        cf_specific1=28.0
        cf_specific2=0
        sp_name = "dec"

    unit_change2=10
    bm_name = "BM_fol_"+sp_name
    df_fm=df_fm.multiply(100)
    #Limit max foliage mass
    df_fm=df_fm.apply(lambda x:min(fm_limit,x))
    #df_cf is one column and df_fm is three column wide

    df_cf_fm_lai = pd.concat([df_cf,df_fm],axis=1)
    df_cf_fm_lai['LAI']=0
    df_cf_fm_lai.columns = ["CF","FM","LAI"]
    
    #Hectars to m2 with 10000
    #Pine or Spruce 
    df_cf_fm_lai = df_cf_fm_lai.rename(columns={"cf":"CF",bm_name:"FM"})
    
    if species == 1 or species ==2:
        df_cf_fm_lai['LAI'] = df_cf_fm_lai.apply(lambda row: 0.5*row['FM']*cf_specific1/10000.0 if row['CF'] < cf_limit else row['LAI'],axis=1)
        df_cf_fm_lai['LAI'] = df_cf_fm_lai.apply(lambda row: 0.5*row['FM']*(0.5*cf_specific1+0.5*cf_specific2)/10000.0 if row['CF'] >= cf_limit else row['LAI'],axis=1)
    else:
        df_cf_fm_lai['LAI'] = df_cf_fm_lai.apply(lambda row: (0.5*row['FM']*cf_specific1)/10000.0,axis=1)
    df = df_cf_fm_lai['LAI']
    df[df <0 ] = 0
    return(df)

def write_gis_header(fname:str,header:str):
    f=open(fname,'w')
    f.write(header)
    return f

def create_all_input(data_loc):
    base = "/scratch/project_2003225/GIT/KE/SpaFHyPeat_CP/"
    folder = base + "MELA_Data/"
    folder_lisa = base + "MELA_Data/Lisatiedot/"


    name_data_51 = "V5/PUDASJARVI_P30051_SKE4_4_V5.MPU"
    name_data_52 = "V5/PUDASJARVI_P30052_SKE4_4_V5.MPU"
    name_data_53 = "V5/PUDASJARVI_P30053_SKE4_4_V5.MPU"

    start_data_51 = "30051_MuodostaRSU_HILAPISTE_PR11.txt"
    start_data_52 = "30052_MuodostaRSU_HILAPISTE_PR12.txt"
    start_data_53 = "30053_MuodostaRSU_HILAPISTE_PR13.txt"

    header_file = "V3/header.csv"

    with open(folder+header_file, newline="\r\n") as f:
        reader = csv.reader(f)
        header = list(reader)
    header = header[0]

    if data_loc == "data_51":
        data_site = pd.read_csv(folder+name_data_51,delim_whitespace=True,header = None)
        s_data = pd.read_csv(folder_lisa+start_data_51,delim_whitespace=True,header = None,names = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32])
    elif data_loc == "data_52":
        data_site = pd.read_csv(folder+name_data_52,delim_whitespace=True,header = None)
        s_data = pd.read_csv(folder_lisa+start_data_52,delim_whitespace=True,header = None,names = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32])    
    elif data_loc == "data_53":
        data_site = pd.read_csv(folder+name_data_53,delim_whitespace=True,header = None)
        s_data = pd.read_csv(folder_lisa+start_data_53,delim_whitespace=True,header = None,names = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32])

    s_data.index= s_data[1]

    i = 0
    a1 = data_site[15+i]/data_site[[15+i,20+i,25+i]].sum(axis=1)
    a2 = data_site[20+i]/data_site[[15+i,20+i,25+i]].sum(axis=1)
    a3 = data_site[25+i]/data_site[[15+i,20+i,25+i]].sum(axis=1)

    a1 = a1.fillna(0)
    a2 = a2.fillna(0)
    a3 = a3.fillna(0)

    a1v = data_site[15+i]
    a2v = data_site[20+i]
    a3v = data_site[25+i]


    a4 = data_site[[15+i,20+i,25+i]].sum(axis=1)
    pine_ppa = data_site[183]*a1
    spruce_ppa = data_site[183]*a2
    dec_ppa = data_site[183]*a3

    dat = {"pine_ppa": pine_ppa, "spruce_ppa":spruce_ppa,"dec_ppa":dec_ppa,"vol":a4 ,"pine_vol":a1v,"spruce_vol":a2v,"dec_vol":a3v}
    dat = pd.DataFrame(dat)
    cf = (dat["pine_ppa"]*3.2234+dat["spruce_ppa"]*2.1661+dat["dec_ppa"]*5.0358)
    cf_vol =0.549*dat["pine_vol"]+.1387*dat["spruce_vol"]+.7537*dat["dec_vol"]

    x_y = {"x":cf,"x1": cf_vol,"y" :(dat["pine_ppa"]+dat["spruce_ppa"]+dat["dec_ppa"]),"y1": dat["vol"]}
    x_y = pd.DataFrame(x_y)

    i = 0
    df_fm_pine_2020 = data_site[[158+i,163+i,168+i]].sum(axis = 1)*(data_site[15+i]/data_site[[15+i,20+i,25+i]].sum(axis=1))
    i = 1
    df_fm_pine_2030 = data_site[[158+i,163+i,168+i]].sum(axis = 1)*(data_site[15+i]/data_site[[15+i,20+i,25+i]].sum(axis=1))
    i = 2
    df_fm_pine_2040 = data_site[[158+i,163+i,168+i]].sum(axis = 1)*(data_site[15+i]/data_site[[15+i,20+i,25+i]].sum(axis=1))
    i = 3
    df_fm_pine_2050 = data_site[[158+i,163+i,168+i]].sum(axis = 1)*(data_site[15+i]/data_site[[15+i,20+i,25+i]].sum(axis=1))
    i = 4
    df_fm_pine_2060 = data_site[[158+i,163+i,168+i]].sum(axis = 1)*(data_site[15+i]/data_site[[15+i,20+i,25+i]].sum(axis=1))
    i = 0
    df_fm_spruce_2020 = data_site[[158+i,163+i,168+i]].sum(axis = 1)*(data_site[20+i]/data_site[[15+i,20+i,25+i]].sum(axis=1))
    i =1
    df_fm_spruce_2030 = data_site[[158+i,163+i,168+i]].sum(axis = 1)*(data_site[20+i]/data_site[[15+i,20+i,25+i]].sum(axis=1))
    i = 2
    df_fm_spruce_2040 = data_site[[158+i,163+i,168+i]].sum(axis = 1)*(data_site[20+i]/data_site[[15+i,20+i,25+i]].sum(axis=1))
    i = 3
    df_fm_spruce_2050 = data_site[[158+i,163+i,168+i]].sum(axis = 1)*(data_site[20+i]/data_site[[15+i,20+i,25+i]].sum(axis=1))
    i = 4
    df_fm_spruce_2060 = data_site[[158+i,163+i,168+i]].sum(axis = 1)*(data_site[20+i]/data_site[[15+i,20+i,25+i]].sum(axis=1))
    i = 0
    df_fm_dec_2020 = data_site[[158+i,163+i,168+i]].sum(axis = 1)*(data_site[25+i]/data_site[[15+i,20+i,25+i]].sum(axis=1))
    i = 1
    df_fm_dec_2030 = data_site[[158+i,163+i,168+i]].sum(axis = 1)*(data_site[25+i]/data_site[[15+i,20+i,25+i]].sum(axis=1))
    i = 2
    df_fm_dec_2040 = data_site[[158+i,163+i,168+i]].sum(axis = 1)*(data_site[25+i]/data_site[[15+i,20+i,25+i]].sum(axis=1))
    i = 3
    df_fm_dec_2050 = data_site[[158+i,163+i,168+i]].sum(axis = 1)*(data_site[25+i]/data_site[[15+i,20+i,25+i]].sum(axis=1))
    i = 4
    df_fm_dec_2060 = data_site[[158+i,163+i,168+i]].sum(axis = 1)*(data_site[25+i]/data_site[[15+i,20+i,25+i]].sum(axis=1))

    i=0

    cf_2020 = (data_site[15+i]*0.5490+data_site[20+i]*0.1387+data_site[25+i]*0.7537)/100
    cf_2020 = pd.Series(np.where(cf_2020 > .88,.88,cf_2020))
    i=1
    cf_2030 = (data_site[15+i]*0.5490+data_site[20+i]*0.1387+data_site[25+i]*0.7537)/100
    cf_2030 = pd.Series(np.where(cf_2030 > .88,.88,cf_2030))
    i=2
    cf_2040 = (data_site[15+i]*0.5490+data_site[20+i]*0.1387+data_site[25+i]*0.7537)/100
    cf_2040 = pd.Series(np.where(cf_2040 > .88,.88,cf_2040))
    i=3
    cf_2050 = (data_site[15+i]*0.5490+data_site[20+i]*0.1387+data_site[25+i]*0.7537)/100
    cf_2050 = pd.Series(np.where(cf_2050 > .88,.88,cf_2050))
    i=4
    cf_2060 = (data_site[15+i]*0.5490+data_site[20+i]*0.1387+data_site[25+i]*0.7537)/100
    cf_2060 = pd.Series(np.where(cf_2060 > .88,.88,cf_2060))

        
    ALL_results = {"Soil":data_site[3],"CF_2020":cf_2020,"CF_2030":cf_2030,"CF_2040":cf_2040,"CF_2050":cf_2050,"CF_2060":cf_2060}
    ALL_results = pd.DataFrame(data = ALL_results)

    ALL_results["LAI_pine_2020"] = foliage_mass_to_LAI(df_fm_pine_2020,cf_2020,1)
    ALL_results["LAI_pine_2030"] = foliage_mass_to_LAI(df_fm_pine_2030,cf_2030,1)
    ALL_results["LAI_pine_2040"] = foliage_mass_to_LAI(df_fm_pine_2040,cf_2040,1)
    ALL_results["LAI_pine_2050"] = foliage_mass_to_LAI(df_fm_pine_2050,cf_2050,1)
    ALL_results["LAI_pine_2060"] = foliage_mass_to_LAI(df_fm_pine_2060,cf_2060,1)

    ALL_results["LAI_spruce_2020"] = foliage_mass_to_LAI(df_fm_spruce_2020,cf_2020,2)
    ALL_results["LAI_spruce_2030"] = foliage_mass_to_LAI(df_fm_spruce_2030,cf_2030,2)
    ALL_results["LAI_spruce_2040"] = foliage_mass_to_LAI(df_fm_spruce_2040,cf_2040,2)
    ALL_results["LAI_spruce_2050"] = foliage_mass_to_LAI(df_fm_spruce_2050,cf_2050,2)
    ALL_results["LAI_spruce_2060"] = foliage_mass_to_LAI(df_fm_spruce_2060,cf_2060,2)

    ALL_results["LAI_dec_2020"] = foliage_mass_to_LAI(df_fm_dec_2020,cf_2020,3)
    ALL_results["LAI_dec_2030"] = foliage_mass_to_LAI(df_fm_dec_2030,cf_2030,3)
    ALL_results["LAI_dec_2040"] = foliage_mass_to_LAI(df_fm_dec_2040,cf_2040,3)
    ALL_results["LAI_dec_2050"] = foliage_mass_to_LAI(df_fm_dec_2050,cf_2050,3)
    ALL_results["LAI_dec_2060"] = foliage_mass_to_LAI(df_fm_dec_2060,cf_2060,3)

    ALL_results["height_2020"] = data_site[[173,178]].max(axis=1)
    ALL_results["height_2030"] = data_site[[173+1,178+1]].max(axis=1)
    ALL_results["height_2040"] = data_site[[173+2,178+2]].max(axis=1)
    ALL_results["height_2050"] = data_site[[173+3,178+3]].max(axis=1)
    ALL_results["height_2060"] = data_site[[173+4,178+4]].max(axis=1)

    ALL_results.index = data_site[0]

    ALL_results["N_Coord"] = s_data[s_data[2] == 1][10]
    ALL_results["E_Coord"] = s_data[s_data[2] == 1][11]+3000000

    site = "pkarjala" #Taking data from pkarjala
    E = pd.read_csv(base + site+"/Ecoord.dat",header = None,skiprows=6) 
    N = pd.read_csv(base + site+"/Ncoord.dat",header = None,skiprows=6)  
    forcing = pd.read_csv(base +site+"/forcing_id.dat",header = None,skiprows=6) 
    depth = pd.read_csv(base +site+"/ditch_depth.dat",header = None,skiprows=6) 
    spacing = pd.read_csv(base +site+"/ditch_spacing.dat",header = None,skiprows=6) 
    soil = pd.read_csv(base +site+"/soil_id.dat",header = None,skiprows=6) 

    forcing_id = []
    depth_val = []
    spacing_val = []
    soil_val = []
    t={"Ecoord":list(E[0]),"Ncoord":list(N[0]),"Forcing":list(forcing[0]),"depth":list(depth[0]),"spacing":list(spacing[0]),"soil":list(soil[0])}
    E_N_data = pd.DataFrame(data = t)

    for i in ALL_results["E_Coord"].keys():
        E_N_data["TEST"] =  abs(E_N_data["Ecoord"]-(ALL_results["E_Coord"][i]))+abs(E_N_data["Ncoord"]-ALL_results["N_Coord"][i])
        forcing_id = forcing_id + [int(E_N_data.iloc[E_N_data["TEST"].idxmin()]["Forcing"])]
        depth_val = depth_val + [E_N_data.iloc[E_N_data["TEST"].idxmin()]["depth"]]
        spacing_val = spacing_val + [E_N_data.iloc[E_N_data["TEST"].idxmin()]["spacing"]]
        soil_val = soil_val + [E_N_data.iloc[E_N_data["TEST"].idxmin()]["soil"]]
    ALL_results["Forcing_id"] = forcing_id
    ALL_results["Ditch_Spacing"] = spacing_val
    ALL_results["Ditch_depth"] = depth_val
    ALL_results["soil_check"] = soil_val

    depth_val = [i/10 for i in depth_val]
    ALL_results["Ditch_depth"] = depth_val
    ALL_results["soil_check"] = soil_val

    ALL_results["SC"] = list(data_site[3])
    ALL_results['Soil_ID']=0

    ALL_results['Soil_ID'] = np.where(ALL_results['SC']>4,1,ALL_results['Soil_ID'])
    ALL_results['Soil_ID'] = np.where(ALL_results['SC']<5,2,ALL_results['Soil_ID'])
    ALL_results['V_Sp_2020'] = list(data_site[15])
    ALL_results['V_Ns_2020'] = list(data_site[20])
    ALL_results['V_Dec_2020'] = list(data_site[25])

    ALL_results['MainSpecies']=0

    df_p = ALL_results.assign(MainSpecies=np.where((ALL_results['V_Sp_2020'] >= ALL_results['V_Ns_2020']) & (ALL_results['V_Sp_2020'] >= ALL_results['V_Dec_2020']),1,ALL_results['MainSpecies']))

    df_ps = df_p.assign(MainSpecies=np.where((df_p['V_Ns_2020'] > df_p['V_Sp_2020']) & (df_p['V_Ns_2020'] >= df_p['V_Dec_2020']),2,df_p['MainSpecies']))
    df_psd = df_ps.assign(MainSpecies=np.where((df_ps['V_Dec_2020'] > df_ps['V_Sp_2020']) & (df_ps['V_Dec_2020'] > df_ps['V_Ns_2020']),3,df_ps['MainSpecies']))
    ALL_results = df_psd.assign(Soil_ID=np.where(((df_psd['SC'] == 3) | (df_psd['SC'] == 4)) & (df_psd['MainSpecies'] == 2),1,df_psd['Soil_ID']))

    for j in range(2,6):
        for i in range(1,10):
            ALL_results["CF_20"+str(j)+str(i)] = (ALL_results["CF_20"+str(j)+"0"]*(10-i)/10+ALL_results["CF_20"+str(j+1)+"0"]*i/10)
            ALL_results["height_20"+str(j)+str(i)] = (ALL_results["height_20"+str(j)+"0"]*(10-i)/10+ALL_results["height_20"+str(j+1)+"0"]*i/10)
            ALL_results["LAI_pine_20"+str(j)+str(i)] = (ALL_results["LAI_pine_20"+str(j)+"0"]*(10-i)/10+ALL_results["LAI_pine_20"+str(j+1)+"0"]*i/10)
            ALL_results["LAI_spruce_20"+str(j)+str(i)] = (ALL_results["LAI_spruce_20"+str(j)+"0"]*(10-i)/10+ALL_results["LAI_spruce_20"+str(j+1)+"0"]*i/10)
            ALL_results["LAI_dec_20"+str(j)+str(i)] = (ALL_results["LAI_dec_20"+str(j)+"0"]*(10-i)/10+ALL_results["LAI_dec_20"+str(j+1)+"0"]*i/10)

    #SINGLE VALUES
    gis_file_header1 = "ncols 1\n"
    gis_file_header2 = "xllcorner 0\nyllcorner 0\ncellsize 16\nNODATA_value -9999\n"
    gis_file_header = gis_file_header1+"nrows "+str(len(ALL_results))+"\n"+gis_file_header2

    D_spacing = ALL_results["Ditch_Spacing"]#cf_dat.index.get_level_values(level=2) - get from earlier data
    D_Depth = ALL_results["Ditch_depth"]#vcf_dat.index.get_level_values(level=3) - get from earlier data
    X1 = ALL_results["E_Coord"]#Ditch_cf_dat.index.get_level_values(level=4) - get from initial data
    Y1 = ALL_results["N_Coord"]#cf_dat.index.get_level_values(level=5) - get from initial data
    Forcing_id = ALL_results["Forcing_id"]#cf_dat.index.get_level_values(level=6) - link to Xand Y
    Soil_ID = ALL_results["Soil_ID"]#data_site[4]cf_dat.index.get_level_values(level=7) # column 4

    dat = {"D_spacing":list(D_spacing),"D_Depth":list(D_Depth),"X1":list(X1),"Y1":list(Y1),"Forcing_id":list(Forcing_id),"Soil_ID":list(Soil_ID)}
    dat = pd.DataFrame(dat)

    extract = ['X1', 'Y1', 'Forcing_id', 'D_spacing','D_Depth', "Soil_ID"]
    names = ['Ecoord.dat','Ncoord.dat','forcing_id.dat','ditch_spacing.dat','ditch_depth.dat','soil_id.dat']

    i = 0
    for ex in extract:
        output_file = base + "Puhdasjarvi_"+data_loc +"/parameters/"+names[i]
        f=write_gis_header(output_file,gis_file_header)
        f.close()
        dat.to_csv(output_file,mode='a',columns=[ex],header=False,index=False,sep =" ")
        i = i+1

    #Multiple VALUES
    gis_file_header1 = "ncols 1\n"
    gis_file_header2 = "xllcorner 0\nyllcorner 0\ncellsize 16\nNODATA_value -9999\n"
    gis_file_header = gis_file_header1+"nrows "+str(len(ALL_results))+"\n"+gis_file_header2

    CF = ALL_results[["CF_20"+str(i) for i in range(20,61)]]
    LAI_pine = ALL_results[["LAI_pine_20"+str(i) for i in range(20,61)]]
    LAI_spruce = ALL_results[["LAI_spruce_20"+str(i) for i in range(20,61)]]
    LAI_dec = ALL_results[["LAI_dec_20"+str(i) for i in range(20,61)]]
    hc = ALL_results[["height_20"+str(i) for i in range(20,61)]]


    dat = {"cf":list(CF),"LAI_ns":list(LAI_spruce),"LAI_sp":list(LAI_pine),"LAI_dec":list(LAI_dec),"hc":list(hc)}
    dat = pd.DataFrame(dat)

    extract = ['cf', 'LAI_ns', 'LAI_sp', 'LAI_dec','H_gM']
    names = ['cf.dat','LAI_spruce.dat','LAI_pine.dat','LAI_decid.dat','hc.dat']
    i = 0
    years = [i for i in range(2020,2061)]
    for ex in extract:
        for y in years:
            output_file = base + "Puhdasjarvi_"+data_loc +"/parameters/"+names[i][:-4] + "_"+str(y)+".dat"
            f=write_gis_header(output_file,gis_file_header)
            f.close()
            if ex == "cf":
                CF["CF_"+str(y)].to_csv(output_file,mode='a',header=False,index=False,sep =" ")
            elif ex == "LAI_ns":
                LAI_spruce["LAI_spruce_"+str(y)].to_csv(output_file,mode='a',header=False,index=False,sep =" ")
            elif ex == "LAI_sp":
                LAI_pine["LAI_pine_"+str(y)].to_csv(output_file,mode='a',header=False,index=False,sep =" ")
            elif ex == "LAI_dec":
                LAI_dec["LAI_dec_"+str(y)].to_csv(output_file,mode='a',header=False,index=False,sep =" ")
            elif ex == "H_gM":
                hc["height_"+str(y)].to_csv(output_file,mode='a',header=False,index=False,sep =" ")
        i = i+1
        
        
#Check start of function create_all_input, as there are folders set at the beginning of the function.
create_all_input("data_53")