import os
import argparse
import pathlib
import pandas as pd
import numpy as np
column_titles = ["No","Maaluokka","Alaryhma","Kasvupaikka","Ojitus2020","Ojitus2030","Ojitus2040","Ojitus2050","Ojitus2060",
                 "Kasittely2020","Kasittely2030","Kasittely2040","Kasittely2050","Kasittely2060","ha", #Ha = col 15
                 "manty_m32020","manty_m32030","manty_m32040","manty_m32050","manty_m32060",           #manty_m32060 = col 20
                 "kuusi_m32020","kuusi_m32030","kuusi_m32040","kuusi_m32050","kuusi_m32060",           #kuusi_m32060 = col 25
                 "lepu_m32020","lepu_m32030","lepu_m32040","lepu_m32050","lepu_m32060",#lepu_m32060 = col 30
                 "manty_poistuma_m3/v2020_2029","manty_poistuma_m3/v2030_2039","manty_poistuma_m3/v2040_2049","manty_poistuma_m3/v2050_2059", #col 34
                 "kuusi_poistuma_m3/v2020_2029","kuusi_poistuma_m3/v2030_2039","kuusi_poistuma_m3/v2040_2049","kuusi_poistuma_m3/v2050_2059", #col 38
                 "lepu_poistuma_m3/v2020_2029","lepu_poistuma_m3/v2030_2039","lepu_poistuma_m3/v2040_2049","lepu_poistuma_m3/v2050_2059",     #col 42
                 "manty_hukka_runko_t/v2020_2029","manty_hukka_runko_t/v2030_2039","manty_hukka_runko_t/v2040_2049","manty_hukka_runko_t/v2050_2059",
                 "manty_hukka_elavat_oksat_t/v2020_2029","manty_hukka_elavat_oksat_t/v2030_2039","manty_hukka_elavat_oksat_t/v2040_2049","manty_hukka_elavat_oksat_t/v2050_2059",
                 "manty_hukka_kuolleet_oksat_t/v2020_2029","manty_hukka_kuolleet_oksat_t/v2030_2039","manty_hukka_kuolleet_oksat_t/v2040_2049","manty_hukka_kuolleet_oksat_t/v2050_2059",
                 "manty_lehdet/neulaset_t/v2020_2029","manty_lehdet/neulaset_t/v2030_2039","manty_lehdet/neulaset_t/v2040_2049","manty_lehdet/neulaset_t/v2050_2059",
                 "manty_kanto_t/v2020_2029","manty_kanto_t/v2030_2039","manty_kanto_t/v2040_2049","manty_kanto_t/v2050_2059",
                 "manty_juuret_t/v2020_2029","manty_juuret_t/v2030_2039","manty_juuret_t/v2040_2049","manty_juuret_t/v2050_2059",
                 "kuusi_hukka_runko_t/v2020_2029","kuusi_hukka_runko_t/v2030_2039","kuusi_hukka_runko_t/v2040_2049","kuusi_hukka_runko_t/v2050_2059",
                 "kuusi_hukka_elavat_oksat_t/v2020_2029","kuusi_hukka_elavat_oksat_t/v2030_2039","kuusi_hukka_elavat_oksat_t/v2040_2049","kuusi_hukka_elavat_oksat_t/v2050_2059",
                 "kuusi_hukka_kuolleet_oksat_t/v2020_2029","kuusi_hukka_kuolleet_oksat_t/v2030_2039","kuusi_hukka_kuolleet_oksat_t/v2040_2049","kuusi_hukka_kuolleet_oksat_t/v2050_2059",
                 "kuusi_lehdet/neulaset_t/v2020_2029","kuusi_lehdet/neulaset_t/v2030_2039","kuusi_lehdet/neulaset_t/v2040_2049","kuusi_lehdet/neulaset_t/v2050_2059",
                 "kuusi_kanto_t/v2020_2029","kuusi_kanto_t/v2030_2039","kuusi_kanto_t/v2040_2049","kuusi_kanto_t/v2050_2059",
                 "kuusi_juuret_t/v2020_2029","kuusi_juuret_t/v2030_2039","kuusi_juuret_t/v2040_2049","kuusi_juuret_t/v2050_2059",
                 "lepu_hukka_runko_t/v2020_2029","lepu_hukka_runko_t/v2030_2039","lepu_hukka_runko_t/v2040_2049","lepu_hukka_runko_t/v2050_2059",
                 "lepu_hukka_elavat_oksat_t/v2020_2029","lepu_hukka_elavat_oksat_t/v2030_2039","lepu_hukka_elavat_oksat_t/v2040_2049","lepu_hukka_elavat_oksat_t/v2050_2059",
                 "lepu_hukka_kuolleet_oksat_t/v2020_2029","lepu_hukka_kuolleet_oksat_t/v2030_2039","lepu_hukka_kuolleet_oksat_t/v2040_2049","lepu_hukka_kuolleet_oksat_t/v2050_2059",
                 "lepu_lehdet/neulaset_t/v2020_2029","lepu_lehdet/neulaset_t/v2030_2039","lepu_lehdet/neulaset_t/v2040_2049","lepu_lehdet/neulaset_t/v2050_2059",
                 "lepu_kanto_t/v2020_2029","lepu_kanto_t/v2030_2039","lepu_kanto_t/v2040_2049","lepu_kanto_t/v2050_2059",
                 "lepu_juuret_t/v2020_2029","lepu_juuret_t/v2030_2039","lepu_juuret_t/v2040_2049","lepu_juuret_t/v2050_2059",
                 "Osuus"]

mela_years = ";2020;2030;2040;2050;2060\n\n"
turvemaat = "Turvemaat;0;0;0;0;0\n"
ojittamaton = " ojittamaton;0;0;0;0;0\n"
ojitettu = " ojitettu;0;0;0;0;0\n\n"
kankaat_zero = "Kankaat;0;0;0;0;0\n"

def check_share(df):
    share = df.loc[df.index[0],'Osuus']
    return share

def check_treatment(df):
    treatment =  df.loc[df.index[0],'Kasittely2020']
    return treatment

def check_no(df):
    No = df.loc[df.index[0],'No']
    return No

def create_empty_df(rows,cols):
    a=np.zeros(shape=(rows,cols))
    df = pd.DataFrame(a)
    return df

def read_mpu_file(file_name):
    """Read MELA MPU file (e.g.'PUDASJARVI_P30051_SKE4_4.MPU')"""
    df = pd.read_csv(file_name,header=None,delim_whitespace=True)
    df.columns = column_titles
    return df

def select_mineral_soils(df):
     df_mineral = df[df['Alaryhma'] == 1]
     return df_mineral
 
#Mineral: alaryhma==1; Maaluokka: 1==forest land, 2 == poorly productive, 3 == unproductive
def select_mineral_soil_areas(df,alaryhma):
    """Return area of the unit for MELA years, remember to multiply with 'Osuus' to not to over estimate.
       Area does not change in MELA years
    """
    df_mineral_1_2020 = df[((df['Maaluokka']==1) & df['Alaryhma']==alaryhma)]
    df_mineral_2_2020 = df[((df['Maaluokka']==2) & df['Alaryhma']==alaryhma)]
    df_mineral_3_2020 = df[((df['Maaluokka']==3) & df['Alaryhma']==alaryhma)]
    df_mineral_1_2020 = pd.concat([df_mineral_1_2020,df_mineral_2_2020,df_mineral_3_2020])
    df_area_1_2020 = df_mineral_1_2020['ha']
    df_area_1_2020 = df_area_1_2020*df_mineral_1_2020['Osuus']
    sum_area_1_2020 = df_area_1_2020.sum()
    df = pd.DataFrame([sum_area_1_2020,sum_area_1_2020,sum_area_1_2020,sum_area_1_2020,sum_area_1_2020])
    df = df.astype(np.str)
    df = df.apply(lambda x: x.str.replace('.',','))
    return df.transpose()


#Mineral: alaryhma==1; Maaluokka: 1==forest land, 2 == poorly productive, 3 == unproductive
def select_mineral_soil_m3(df,alaryhma,maaluokka):
    """Return stock m3 for mineral soils for each MELA year"""
    #Kasittely == 7 means outside timber production
    df_mineral_1_2020_2060 = df[(df['Alaryhma']==alaryhma) & (df['Maaluokka']==maaluokka)]
    #print("INDEX",len(df_mineral_1_2020_2060.index),df_mineral_1_2020_2060.index,df_mineral_1_2020_2060['Kasittely2020'] )
    if len(df_mineral_1_2020_2060) == 0:
        df_zero1=pd.DataFrame([[0,0,0,0,0]])
        df_zero1.columns = ['manty_m32020','manty_m32030','manty_m32040','manty_m32050','manty_m32060']
        df_zero2=pd.DataFrame([[0,0,0,0,0]])
        df_zero2.columns = ['kuusi_m32020','kuusi_m32030','kuusi_m32040','kuusi_m32050','kuusi_m32060']
        df_zero3=pd.DataFrame([[0,0,0,0,0]])
        df_zero3.columns = ['lepu_m32020','lepu_m32030','lepu_m32040','lepu_m32050','lepu_m32060']
        tuple3 = (df_zero1,df_zero2,df_zero3)
        return tuple3 
    df_ma_2020_2060 = df_mineral_1_2020_2060[['manty_m32020','manty_m32030','manty_m32040','manty_m32050','manty_m32060']]
    df_ku_2020_2060 = df_mineral_1_2020_2060[['kuusi_m32020','kuusi_m32030','kuusi_m32040','kuusi_m32050','kuusi_m32060']]
    df_lepu_2020_2060 = df_mineral_1_2020_2060[['lepu_m32020','lepu_m32030','lepu_m32040','lepu_m32050','lepu_m32060']]
    #Remember to multiply with 'Osuus' to not to over estimate.
    share =  df_mineral_1_2020_2060['Osuus']
    #Indexes are the same as in the original
    #This is important: we get a number
    share = share[share.index[0]]
    df_ma_2020_2060 = df_ma_2020_2060*share
    df_ku_2020_2060 = df_ku_2020_2060*share
    df_lepu_2020_2060 = df_lepu_2020_2060*share
    return (df_ma_2020_2060,df_ku_2020_2060,df_lepu_2020_2060)

#Mineral: alaryhma==1; Maaluokka: 1==forest land, 2 == poorly productive, 3 == unproductive
def select_mineral_soil_natural_loss_m3(df,alaryhma,maaluokka):
    """Return stock m3 natural loss for mineral soils for each MELA year"""
    #Kasittely == 7 means outside timber production
    df_mineral_1_2020_2060 = df[(df['Alaryhma']==alaryhma) & (df['Maaluokka']==maaluokka)]
    if len(df_mineral_1_2020_2060) == 0:
        df_zero1=pd.DataFrame([[0,0,0,0]])
        df_zero1.columns = ['manty_poistuma_m3/v2020_2029','manty_poistuma_m3/v2030_2039','manty_poistuma_m3/v20400_2049',
                            'manty_poistuma_m3/v2050_2059']
        df_zero2=pd.DataFrame([[0,0,0,0]])
        df_zero2.columns = ['kuusi_poistuma_m3/v2020_2029','kuusi_poistuma_m3/v2030_2039','kuusi_poistuma_m3/v20400_2049',
                            'kuusi_poistuma_m3/v2050_2059']
        df_zero3=pd.DataFrame([[0,0,0,0]])
        df_zero3.columns = ['lepu_poistuma_m3/v2020_2029','lepu_poistuma_m3/v2030_2039','lepu_poistuma_m3/v20400_2049',
                            'lepu_poistuma_m3/v2050_2059']
        tuple3 = (df_zero1,df_zero2,df_zero3)
        return tuple3
    df_ma_2020_2060 = df_mineral_1_2020_2060[['manty_poistuma_m3/v2020_2029','manty_poistuma_m3/v2030_2039',
                                              'manty_poistuma_m3/v2040_2049','manty_poistuma_m3/v2050_2059']]
    df_ku_2020_2060 = df_mineral_1_2020_2060[['kuusi_poistuma_m3/v2020_2029','kuusi_poistuma_m3/v2030_2039',
                                              'kuusi_poistuma_m3/v2040_2049','kuusi_poistuma_m3/v2050_2059']]
    df_lepu_2020_2060 = df_mineral_1_2020_2060[['lepu_poistuma_m3/v2020_2029','lepu_poistuma_m3/v2030_2039',
                                                'lepu_poistuma_m3/v2040_2049','lepu_poistuma_m3/v2050_2059']]
    share = df_mineral_1_2020_2060['Osuus']
    #This is important: we get the number
    share = share[share.index[0]]
    df_ma_2020_2060 = df_ma_2020_2060*share
    df_ku_2020_2060 = df_ku_2020_2060*share
    df_lepu_2020_2060 = df_lepu_2020_2060*share
    return (df_ma_2020_2060,df_ku_2020_2060, df_lepu_2020_2060)
#Mineral: alaryhma==1; Maaluokka: 1==forest land, 2 == poorly productive, 3 == unproductive,
#Waste: 1. stem, 2. living branches + dead branches, 3.foliage, 4. stump, 5. roots
#Note: we need to get living and  dead branches separately. Thus: 2.1 living branches and 2.2 dead branches
def select_mineral_waste_dm(df,alaryhma,maaluokka,waste):
    """Return wastewood (1000t DryMass)  for mineral soils for each MELA year"""
    df_mineral_1_2020_2060 = df[(df['Alaryhma']==alaryhma) & (df['Maaluokka']==maaluokka)]
    if waste == 1:
        waste_ma_2020_2060 = ["manty_hukka_runko_t/v2020_2029","manty_hukka_runko_t/v2030_2039","manty_hukka_runko_t/v2040_2049",
                              "manty_hukka_runko_t/v2050_2059"]
        waste_ku_2020_2060 = ["kuusi_hukka_runko_t/v2020_2029","kuusi_hukka_runko_t/v2030_2039","kuusi_hukka_runko_t/v2040_2049",
                              "kuusi_hukka_runko_t/v2050_2059"]
        waste_lepu_2020_2060 = ["lepu_hukka_runko_t/v2020_2029","lepu_hukka_runko_t/v2030_2039","lepu_hukka_runko_t/v2040_2049",
                                "lepu_hukka_runko_t/v2050_2059"]
    elif waste == 2.1:
        waste_ma_2020_2060 = ["manty_hukka_elavat_oksat_t/v2020_2029","manty_hukka_elavat_oksat_t/v2030_2039",
                               "manty_hukka_elavat_oksat_t/v2040_2049","manty_hukka_elavat_oksat_t/v2050_2059"]
        waste_ku_2020_2060 = ["kuusi_hukka_elavat_oksat_t/v2020_2029","kuusi_hukka_elavat_oksat_t/v2030_2039","kuusi_hukka_elavat_oksat_t/v2040_2049",
                              "kuusi_hukka_elavat_oksat_t/v2050_2059"]
        waste_lepu_2020_2060 = ["lepu_hukka_elavat_oksat_t/v2020_2029","lepu_hukka_elavat_oksat_t/v2030_2039","lepu_hukka_elavat_oksat_t/v2040_2049",
                                "lepu_hukka_elavat_oksat_t/v2050_2059"]
    elif waste == 2.2:
        waste_ma_2020_2060 = ["manty_hukka_kuolleet_oksat_t/v2020_2029","manty_hukka_kuolleet_oksat_t/v2030_2039",
                               "manty_hukka_kuolleet_oksat_t/v2040_2049","manty_hukka_kuolleet_oksat_t/v2050_2059"]
        waste_ku_2020_2060 = ["kuusi_hukka_kuolleet_oksat_t/v2020_2029","kuusi_hukka_kuolleet_oksat_t/v2030_2039","kuusi_hukka_kuolleet_oksat_t/v2040_2049",
                              "kuusi_hukka_kuolleet_oksat_t/v2050_2059"]
        waste_lepu_2020_2060 = ["lepu_hukka_kuolleet_oksat_t/v2020_2029","lepu_hukka_kuolleet_oksat_t/v2030_2039","lepu_hukka_kuolleet_oksat_t/v2040_2049",
                                "lepu_hukka_kuolleet_oksat_t/v2050_2059"]
    elif waste == 3:
        waste_ma_2020_2060 = ["manty_lehdet/neulaset_t/v2020_2029","manty_lehdet/neulaset_t/v2030_2039","manty_lehdet/neulaset_t/v2040_2049","manty_lehdet/neulaset_t/v2050_2059"]
        waste_ku_2020_2060 = ["kuusi_lehdet/neulaset_t/v2020_2029","kuusi_lehdet/neulaset_t/v2030_2039","kuusi_lehdet/neulaset_t/v2040_2049","kuusi_lehdet/neulaset_t/v2050_2059"]
        waste_lepu_2020_2060 = ["lepu_lehdet/neulaset_t/v2020_2029","lepu_lehdet/neulaset_t/v2030_2039","lepu_lehdet/neulaset_t/v2040_2049",
                                "lepu_lehdet/neulaset_t/v2050_2059"]   
    elif waste == 4:
        waste_ma_2020_2060 = ["manty_kanto_t/v2020_2029","manty_kanto_t/v2030_2039","manty_kanto_t/v2040_2049","manty_kanto_t/v2050_2059"]
        waste_ku_2020_2060 = ["kuusi_kanto_t/v2020_2029","kuusi_kanto_t/v2030_2039","kuusi_kanto_t/v2040_2049","kuusi_kanto_t/v2050_2059",]
        waste_lepu_2020_2060 = ["lepu_kanto_t/v2020_2029","lepu_kanto_t/v2030_2039","lepu_kanto_t/v2040_2049","lepu_kanto_t/v2050_2059"]
        
    elif waste == 5:
        waste_ma_2020_2060 = ["manty_juuret_t/v2020_2029","manty_juuret_t/v2030_2039","manty_juuret_t/v2040_2049","manty_juuret_t/v2050_2059"]
        waste_ku_2020_2060 = ["kuusi_juuret_t/v2020_2029","kuusi_juuret_t/v2030_2039","kuusi_juuret_t/v2040_2049","kuusi_juuret_t/v2050_2059"]
        waste_lepu_2020_2060 = ["lepu_juuret_t/v2020_2029","lepu_juuret_t/v2030_2039","lepu_juuret_t/v2040_2049","lepu_juuret_t/v2050_2059"]
    else:
        print("No matching waste type, quitting")
        quit()
    if len(df_mineral_1_2020_2060) == 0:
        df_zero1=pd.DataFrame([[0,0,0,0]])
        df_zero1.columns = waste_ma_2020_2060
        df_zero2=pd.DataFrame([[0,0,0,0]])
        df_zero2.columns = waste_ku_2020_2060
        df_zero3=pd.DataFrame([[0,0,0,0]])
        df_zero3.columns = waste_lepu_2020_2060
        tuple4 = (df_zero1,df_zero2,df_zero3)
        return tuple4
    df_ma_2020_2060 = df_mineral_1_2020_2060[waste_ma_2020_2060]
    df_ku_2020_2060 = df_mineral_1_2020_2060[waste_ku_2020_2060]
    df_lepu_2020_2060 = df_mineral_1_2020_2060[waste_lepu_2020_2060]
    #Remember to multiply with 'Osuus' to not to over estimate.
    share = df_mineral_1_2020_2060['Osuus']
    #This is important: we get the number
    share = share[share.index[0]]
    df_ma_2020_2060 = df_ma_2020_2060*share
    df_ku_2020_2060 = df_ku_2020_2060*share
    df_lepu_2020_2060 = df_lepu_2020_2060*share
    return (df_ma_2020_2060,df_ku_2020_2060, df_lepu_2020_2060)
 
            
def merge_living_dead_branches(living,dead):
    """Combine living and dead branches here"""
    df_fol_ma_1 = living[0]
    df_fol_ku_1 = living[1]
    df_fol_lepu_1 = living[2]
    df_fol_ma_2 = dead[0]
    df_fol_ku_2 = dead[1]
    df_fol_lepu_2 = dead[2]
    df_ma = pd.DataFrame(np.concatenate([df_fol_ma_1,df_fol_ma_2],axis=0))
    df_ku = pd.DataFrame(np.concatenate([df_fol_ku_1,df_fol_ku_2],axis=0))
    df_lepu = pd.DataFrame(np.concatenate([df_fol_lepu_1,df_fol_lepu_2],axis=0))
    df_ma_sum=df_ma.sum()
    df_ma = pd.DataFrame(df_ma_sum).transpose()
    df_ku_sum=df_ku.sum()
    df_ku = pd.DataFrame(df_ku_sum).transpose()
    df_lepu_sum=df_lepu.sum()
    df_lepu = pd.DataFrame(df_lepu_sum).transpose()
    return (df_ma,df_ku,df_lepu)

def select_waste_with_land_class(df,maaluokka):
    tuple3_dm_1_1 = select_mineral_waste_dm(df,1,maaluokka,1)
    tuple3_dm_1_21 = select_mineral_waste_dm(df,1,maaluokka,2.1)
    tuple3_dm_1_22 = select_mineral_waste_dm(df,1,maaluokka,2.2)
    tuple3_dm_1_3 = select_mineral_waste_dm(df,1,maaluokka,3)
    tuple3_dm_1_4 = select_mineral_waste_dm(df,1,maaluokka,4)
    tuple3_dm_1_5 = select_mineral_waste_dm(df,1,maaluokka,5)
    #The numbering corresponds Ristos input file
    df_manty_113 = tuple3_dm_1_1[0]
    df_kuusi_123 = tuple3_dm_1_1[1]
    df_lepu_133 = tuple3_dm_1_1[2]
    df_manty_111 = tuple3_dm_1_3[0]
    df_kuusi_121 = tuple3_dm_1_3[1]
    df_lepu_131 = tuple3_dm_1_3[2]
    df_manty_114 = tuple3_dm_1_4[0]
    df_kuusi_124 = tuple3_dm_1_4[1]
    df_lepu_134 = tuple3_dm_1_4[2]
    df_manty_115= tuple3_dm_1_5[0]
    df_kuusi_125 = tuple3_dm_1_5[1]
    df_lepu_135 = tuple3_dm_1_5[2]
    #Risto merges living and dead branches
    (df_manty_112,df_kuusi_122,df_lepu_132) = merge_living_dead_branches(tuple3_dm_1_21,tuple3_dm_1_22)
    df_all = pd.DataFrame(np.concatenate([df_manty_113.values,df_manty_112.values,df_manty_111.values],axis=0))
    df_all = pd.DataFrame(np.concatenate([df_all,df_manty_114.values,df_manty_115.values],axis=0))
    df_all = pd.DataFrame(np.concatenate([df_all,df_kuusi_123.values,df_kuusi_122.values],axis=0))
    df_all = pd.DataFrame(np.concatenate([df_all,df_kuusi_121.values,df_kuusi_124.values,df_kuusi_125.values],axis=0))
    df_all = pd.DataFrame(np.concatenate([df_all,df_lepu_133.values,df_lepu_132.values],axis=0))
    df_all = pd.DataFrame(np.concatenate([df_all,df_lepu_131.values,df_lepu_134.values,df_lepu_135.values],axis=0))
    #df_all_sum =  df_all.sum()
    #df_all_sum = pd.DataFrame(df_all_sum).transpose()
    #df_all = pd.concat([df_all_sum,df_all])
    return df_all

def sum_mineral_m3(df1,df2,df3):
    """Sum three volumes, usually three "maaluokka" classes for one species"""
    #df_all = pd.concat([df1,df2,df3],axis=0)
    df_all = pd.DataFrame(np.concatenate([df1.values,df2.values,df3.values],axis=0))
    #Series as result
    s_all_sum = df_all.sum()
    #Back to data frame, right shape
    df_all_sum = pd.DataFrame(s_all_sum).transpose()
    return df_all_sum

def species_stock_m3(tuple_stock_fl,tuple_stock_poor,tuple_stock_waste,index):
    """Calculate esentially "1 1", "1 2" and "1 3
       index: 0 pine, 1 spruce and 2 deciduous"""
    #print("Before")
    df1 = tuple_stock_fl[index]
    df2 = tuple_stock_poor[index]
    df3 = tuple_stock_waste[index]
    #print(df1)
    #print(df1.shape,df2.shape,df3.shape)
    #print("After")
    df = sum_mineral_m3(df1,df2,df3)
    #print(df.shape)
    return df

def mela_to_yasso(mpu_file,outdir,logfile):
    """Construct from MELA file the file that Risto uses for Yasso
       Construct a data frame that is then printed out.
    """
    p = pathlib.Path(mpu_file,)
    stem = p.stem
    df=read_mpu_file(mpu_file)
    df=select_mineral_soils(df)
    for i in range(0,len(df)):
        #print("Start ",i)
        df_current = pd.DataFrame(df.loc[i,])
        df_current = df_current.transpose()
        out_file = str(int(df_current.loc[i,'No']))
        treatment = check_treatment(df_current)
        No = check_no(df_current)
        if treatment == 7:
            print("TREATMENT",No,treatment)
        share = check_share(df_current)
        if share < 1.0:
            print("SHARE",No,share)
            out_file = outdir+"/"+out_file+"_"+stem+"_"+str(share)+".txt"
        else:
            out_file = outdir+"/"+out_file+"_"+stem+".txt"
        f = open(out_file,'w')
        logfile.write(str(i)+"Starting: "+ out_file+'\n')
        f.write(mela_years)
        df_1_fl = select_mineral_soil_areas(df_current,1)
        df_1_fl.insert(0,'Kankaat','Kankaat')
        df_1_fl.to_csv(f,mode='a',sep=';',index=False,header=False)
        f.write(turvemaat)
        f.write(ojittamaton)
        f.write(ojitettu)
        #Mineral soil m3
        tuple3_m3_1 = select_mineral_soil_m3(df_current,1,1)
        tuple3_m3_2 = select_mineral_soil_m3(df_current,1,2)
        tuple3_m3_3 = select_mineral_soil_m3(df_current,1,3)
        
        #Pine, "1 1"
        df_11 = species_stock_m3(tuple3_m3_1,tuple3_m3_2,tuple3_m3_3,0)
        #Spruce, "1 2" 
        df_12  = species_stock_m3(tuple3_m3_1,tuple3_m3_2,tuple3_m3_3,1)
        #Deciduous, "1 3"
        df_13 = species_stock_m3(tuple3_m3_1,tuple3_m3_2,tuple3_m3_3,2)
        #Peat, "3 1", we should get zeroes with df_current mineral only
        tuple_31 = select_mineral_soil_m3(df_current,2,1)
        #Peat, "3 2"
        tuple_32 = select_mineral_soil_m3(df_current,2,1)
        #Peat, "3 3"
        tuple_33 = select_mineral_soil_m3(df_current,2,1)
        #Peat, "2 1"
        tuple_21 = select_mineral_soil_m3(df_current,2,1)
        #Peat, "2 2"
        tuple_22 = select_mineral_soil_m3(df_current,2,1)
        #Peat, "2 3"
        tuple_23 = select_mineral_soil_m3(df_current,2,1) 
        df_all = pd.DataFrame(np.concatenate([df_11.values,df_12.values,df_13.values,tuple_31[0].values,
                                              tuple_32[0].values,tuple_33[0].values,tuple_21[0].values,
                                              tuple_22[0].values,tuple_23[0].values],axis=0))
        df_all_sum =  df_all.sum()
        df_all_sum = pd.DataFrame(df_all_sum).transpose()
        df_all = pd.concat([df_all_sum,df_all])
        df_min_m3_zero = df_all*0
        df_all.insert(0,"index",["Tilavuus (1000 m3)","1 1","1 2", "1 3","3 1","3 2","3 3","2 1", "2 2","2 3"])
        df_min_m3_zero.insert(0,"index",["Tilavuus (1000 m3)","1 1","1 2", "1 3","3 1","3 2","3 3","2 1", "2 2","2 3"])
        df_all = df_all.astype(np.str)
        df_all = df_all.apply(lambda x: x.str.replace('.',','))
        df_all.to_csv(f,mode='a',sep=';',index=False,header=False)
        f.write(os.linesep)
        #Natural loss m3"
        tuple3_m3_1 = select_mineral_soil_natural_loss_m3(df_current,1,1)
        tuple3_m3_2 = select_mineral_soil_natural_loss_m3(df_current,1,2)
        tuple3_m3_3 = select_mineral_soil_natural_loss_m3(df_current,1,3)
        #Pine, "1 1"
        df_11 = species_stock_m3(tuple3_m3_1,tuple3_m3_2,tuple3_m3_3,0)
        #Spruce, "1 2" 
        df_12  = species_stock_m3(tuple3_m3_1,tuple3_m3_2,tuple3_m3_3,1)
        #Deciduous, "1 3"
        df_13 = species_stock_m3(tuple3_m3_1,tuple3_m3_2,tuple3_m3_3,2)
        #Peat, "3 1", we should get zeroes with df_current mineral only
        #tuple_31 = select_mineral_soil_m3(df_current,2,1)
        #Peat, "3 2"
        #tuple_32 = select_mineral_soil_m3(df_current,2,1)
        #Peat, "3 3"
        #tuple_33 = select_mineral_soil_m3(df_current,2,1)
        #Peat, "2 1"
        #tuple_21 = select_mineral_soil_m3(df_current,2,1)
        #Peat, "2 2"
        #tuple_22 = select_mineral_soil_m3(df_current,2,1)
        #Peat, "2 3"
        peat_array = np.zeros(shape=df_11.shape)
        tuple_23 = select_mineral_soil_m3(df_current,2,1)
        tuple_23 = select_mineral_soil_m3(df_current,2,1)
        df_tmp=pd.DataFrame(np.concatenate([df_11.values,df_12.values,df_13.values],axis=0))
        #df_all = pd.DataFrame(np.concatenate([df_11.values,df_12.values,df_13.values,tuple_31[0].values,
        #                                      tuple_32[0].values,tuple_33[0].values,tuple_21[0].values,
        #                                      tuple_22[0].values,tuple_23[0].values],axis=0))
        df_all = pd.DataFrame(np.concatenate([df_11.values,df_12.values,df_13.values,peat_array,
                                              peat_array,peat_array,peat_array,
                                              peat_array,peat_array],axis=0))
        df_all_sum =  df_all.sum()
        df_all_sum = pd.DataFrame(df_all_sum).transpose()
        df_all = pd.concat([df_all_sum,df_all])
        df_natural_loss_zero = df_all*0
        df_all.insert(0,"index",["Luonnonpoistuma (1000 m3)","1 1","1 2", "1 3","3 1","3 2","3 3","2 1", "2 2","2 3"])
        df_natural_loss_zero.insert(0,"index",["Luonnonpoistuma (1000 m3)","1 1","1 2", "1 3","3 1","3 2","3 3","2 1", "2 2","2 3"])
        df_all = df_all.astype(np.str)
        df_all = df_all.apply(lambda x: x.str.replace('.',','))
        df_all.to_csv(f,mode='a',sep=';',index=False,header=False)
        f.write(os.linesep)
        #Waste, forest land, 3 land classes
        dfw1 = select_waste_with_land_class(df_current,1)
        dfw2 = select_waste_with_land_class(df_current,2)
        dfw3 = select_waste_with_land_class(df_current,3)
        df_all = dfw1+dfw2+dfw3
        df_all_sum =  df_all.sum()
        df_all_sum = pd.DataFrame(df_all_sum).transpose()
        df_all = pd.concat([df_all_sum,df_all])
        df_zeros = create_empty_df(30,4)
        df_all = pd.DataFrame(np.concatenate([df_all.values,df_zeros.values],axis=0))
        df_waste_zero = df_all*0
        col_title_ls = ["Hukkapuu massat (1000 t DM)",
                        "1 1 3","1 1 2","1 1 1","1 1 4","1 1 5",
                        "1 2 3","1 2 2","1 2 1","1 2 4","1 2 5",
                        "1 3 3","1 3 2","1 3 1","1 3 4","1 3 5",
                        "3 1 3","3 1 2","3 1 1","3 1 4","3 1 5",
                        "3 2 3","3 2 2","3 2 1","3 2 4","3 2 5",
                        "3 3 3","3 3 2","3 3 1","3 3 4","3 3 5",
                        "2 1 3","2 1 2","2 1 1","2 1 4","2 1 5",
                        "2 2 3","2 2 2","2 2 1","2 2 4","2 2 5",
                        "2 3 3","2 3 2","2 3 1","2 3 4","2 3 5"]
                        
        df_all.insert(0,"index",col_title_ls)
        df_waste_zero.insert(0,"index",col_title_ls)
        df_all = df_all.astype(np.str)
        f.write("Hukkapuu massat (1000 t DM)\n")
        #df_all.loc[0:] = ["","","","",""]
        df_all = df_all.apply(lambda x: x.str.replace('.',','))
        df_all.loc[1:,:].to_csv(f,mode='a',sep=';',index=False,header=False)
        f.write(os.linesep)
        f.write(kankaat_zero)
        f.write(turvemaat)
        f.write(ojittamaton)
        f.write(ojitettu)
        df_min_m3_zero = df_min_m3_zero.astype(np.str)
        df_min_m3_zero = df_min_m3_zero.apply(lambda x: x.str.replace('.',','))
        df_min_m3_zero.to_csv(f,mode='a',sep=';',index=False,header=False)
        f.write(os.linesep)
        df_natural_loss_zero = df_natural_loss_zero.astype(np.str)
        df_natural_loss_zero = df_natural_loss_zero.apply(lambda x: x.str.replace('.',','))
        df_natural_loss_zero.to_csv(f,mode='a',sep=';',index=False,header=False)
        f.write(os.linesep)
        df_waste_zero = df_waste_zero.astype(np.str)
        df_waste_zero = df_waste_zero.apply(lambda x: x.str.replace('.',','))
        f.write("Hukkapuu massat (1000 t DM);;;;;;\n")
        df_waste_zero.loc[1:,:].to_csv(f,mode='a',sep=';',index=False,header=False)
        f.close()
        logfile.write("Finished: " + out_file+'\n')
        #print("End ",i)
        #return df
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f',type=str,dest='f',help='MELA MPU data input file')
    parser.add_argument('-d',type=str,dest='d',help='Output directory Yasso initial files')
    parser.add_argument('-l',type=str,dest='l',help='Log file')
    args = parser.parse_args()
    if args.f == None:
        print("No MELA MPU file")
        quit()
    if args.d == None:
        print("No output directory")
        quit()
    if args.l == None:
        print("No logfile")
        quit()
    logfile = open(args.l,"w")
    mela_to_yasso(args.f,args.d,logfile)
    logfile.close()
