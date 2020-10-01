import pandas as pd
import geopandas as gpd
import numpy as np
#


def driver(name, nm):
    peat_loc = pd.read_csv('/scratch/project_2000611/SOMPASOILS/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/'+name+'/'+nm+'_peat_loc.csv', header = None)
    big_file = "/scratch/project_2000611/SOMPASOILS/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/image_files.gdb"
    if nm == "is":
        rr = "is2"
    else:
        rr = nm
    data = gpd.read_file(big_file, driver = "FileGDB",layer = rr)
    DATA1= peat_loc
    DATA1['TEST']= 1
    peat_loc = pd.read_csv('/scratch/project_2000611/SOMPASOILS/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/'+name+'/'+nm+'_peat_loc.csv', header = None)
    DATA1['gridcode'] = peat_loc
    j2 = data.merge(DATA1, on='gridcode')
    kunta_shp = gpd.read_file("/scratch/project_2000611/SOMPASOILS/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/Kuntajako/SuomenKuntajako_2020_10k.shp")
    shape_by_KUNTA = gpd.sjoin(j2, kunta_shp[['NAMEFIN', 'geometry','NATCODE']], how='left', op='intersects')
    tt = shape_by_KUNTA[['NAMEFIN','gridcode','NATCODE']].copy()
    kuntat = {"Aland": ["035","043","060","062","065","076","170","295","318","417","438","478","736","766","771","941"],"Varsinais-Suomi": ["019","202","284","304","322","400","423","430","445","480","481","503","529","538","561","577","631","636","680","704","734","738","761","833","853","895","918"],"Uusimaa": ["018","049","078","091","092","106","149","186","224","235","245","257","407","434","444","504","505","543","611","616","638","710","753","755","858","927"],"Satakunta": ["050","051","079","099","102","181","214","230","271","442","484","531","608","609","684","747","783","886"],"Pohjois-Savo": ["140","174","204","239","263","297","402","420","595","686","687","749","762","778","844","857","915","921","925"],"Pohjois-Pohjanmaa": ["009","069","071","072","139","208","244","305","317","425","436","483","494","535","563","564","615","625","626","630","678","691","746","748","785","791","832","859","889","977"],"Pohjois-Karjala": ["146","167","176","260","276","309","422","426","541","607","707","848","911"],"Pohjanmaa": ["152","231","280","287","288","399","440","475","499","545","598","599","893","905","946"],"Pirkanmaa": ["020","108","143","177","211","250","418","508","536","562","581","604","619","635","702","790","837","887","908","922","936","980"],"Keski-Pohjanmaa": ["074","217","236","272","421","584","849","924"],"Paijat-Hame": ["016","081","098","111","316","398","560","576","781"],"Lappi": ["047","148","240","241","261","273","320","498","583","614","683","698","732","742","751","758","845","851","854","890","976"],"Kymenlaakso": ["075","142","285","286","489","624","935"],"Kanta-Hame": ["061","082","086","103","109","165","169","433","694","834","981"],"Kainuu": ["105","205","290","578","620","697","765","777"],"Etela-Savo": ["046","090","097","171","178","213","491","507","588","593","623","681","740","768"],"Etela-Karjala": ["153","405","416","441","580","689","700","739","831"],"Etela-Pohjanmaa": ["005","010","052","145","151","218","232","233","300","301","403","408","743","759","846","934","989"], "Keski-Suomi": ["077","172","179","182","216","226","249","256","265","275","291","312","410","435","495","500","592","601","729","850","892","931","992"]} 
    for k in kuntat.keys():
        tt.loc[tt['NATCODE'].isin(kuntat[k]),"NATCODE"] = k
    tt = tt.drop_duplicates()
    tt.drop_duplicates(subset="gridcode",keep="first",inplace=True)
    drainage =  pd.read_csv('/scratch/project_2000611/SOMPASOILS/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/nuts3_ditch_depth.csv', sep = ",")  
    drainage_mean =  pd.read_csv('/scratch/project_2000611/SOMPASOILS/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/nuts3_ditch_depth_means.csv', sep = ",")  
    peat_loc = pd.read_csv('/scratch/project_2000611/SOMPASOILS/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/'+name+'/'+nm+'_peat_loc.csv',header = None)  
    kasvu = np.loadtxt('/scratch/project_2000611/SOMPASOILS/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/'+name+'/segmaj1.kasvupaikka_vmi1x_1216.hits')
    np.savetxt('/scratch/project_2000611/SOMPASOILS/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/'+name+'/'+name+"_segmaj1.kasvupaikka_vmi1x_1216.hits",kasvu)
    kasvu_paikka = pd.read_csv('/scratch/project_2000611/SOMPASOILS/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/'+name+'/'+name+'_segmaj1.kasvupaikka_vmi1x_1216.hits',header = None,sep = " ")  
    paatyppi = np.loadtxt('/scratch/project_2000611/SOMPASOILS/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/'+name+'/segmaj1.paatyyppi_vmi1x_1216.hits')
    np.savetxt('/scratch/project_2000611/SOMPASOILS/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/'+name+'/'+name+"_segmaj1.paatyyppi_vmi1x_1216.hits",paatyppi)
    paatyppi = pd.read_csv('/scratch/project_2000611/SOMPASOILS/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/'+name+'/'+name+'_segmaj1.paatyyppi_vmi1x_1216.hits',header = None,sep = " ") 
    '''kasvu_dict = {}
    paatyppi_dict = {}
    for i in range(0,len(kasvu_paikka)):
        kasvu_dict[kasvu_paikka[0][i]]=kasvu_paikka[2][i]
        paatyppi_dict[paatyppi[0][i]]=paatyppi[2][i]
    PP = []
    KUNT = []
    KUNT_MEAN = []
    typpi = []
    
    for i in range(0,len(peat_loc)):
        PP = PP+[kasvu_dict[peat_loc[0][i]]]
        SC = kasvu_dict[peat_loc[0][i]]
        if SC <= 5:
            KUNT = KUNT +[float(drainage.loc[(drainage['kptyy'] == SC) & (drainage['Nuts3Nimi']== peat_loc["NATCODE"][i])]['keski'])]
            KUNT_MEAN = KUNT_MEAN + [float(drainage_mean.loc[(drainage_mean['Nuts3Nimi'] == peat_loc["NATCODE"][i])]['keski'])]
        else:
            KUNT = KUNT +[float(drainage_mean.loc[(drainage_mean['Nuts3Nimi'] == peat_loc["NATCODE"][i])]['keski'])]
            KUNT_MEAN = KUNT_MEAN + [float(drainage_mean.loc[(drainage_mean['Nuts3Nimi'] == peat_loc["NATCODE"][i])]['keski'])]
        typpi = typpi + [paatyppi_dict[peat_loc[0][i]]]'''
    peat_loc =  tt.merge(peat_loc,left_on="gridcode",right_on=0)
    peat_loc["SC"] = list(paatyppi[kasvu_paikka[0].isin(peat_loc[0])][2])
    peat_loc["Paatyppi"] = list(paatyppi[paatyppi[0].isin(peat_loc[0])][2])
    peat_loc["depth"] = 0
    peat_loc["depth_mean"] = 0
    for SC in range(1,6):
        for NUTS in kuntat.keys():
            if SC <= 5:
                peat_loc.loc[(peat_loc["NATCODE"] == NUTS),"depth_mean"] = float(drainage_mean.loc[(drainage_mean['Nuts3Nimi'] == NUTS)]['keski'])
                if  len(drainage.loc[(drainage['kptyy'] == SC) & (drainage['Nuts3Nimi']== NUTS)])>0:
                    peat_loc.loc[(peat_loc["SC"]==SC) & (peat_loc["NATCODE"] == NUTS),"depth"] =float(drainage.loc[(drainage['kptyy'] == SC) & (drainage['Nuts3Nimi']== NUTS)]['keski'])
                else:
                    peat_loc.loc[(peat_loc["SC"]==SC) & (peat_loc["NATCODE"] == NUTS),"depth"] =float(drainage_mean.loc[(drainage_mean['Nuts3Nimi'] == NUTS)]['keski'])
            else:
                peat_loc.loc[(peat_loc["NATCODE"] == NUTS),"depth_mean"] = float(drainage_mean.loc[(drainage_mean['Nuts3Nimi'] == NUTS)]['keski'])
                peat_loc.loc[(peat_loc["SC"]==SC) & (peat_loc["NATCODE"] == NUTS),"depth"] =float(drainage_mean.loc[(drainage_mean['Nuts3Nimi'] == NUTS)]['keski'])
    header = pd.read_csv('/scratch/project_2000611/SOMPASOILS/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/'+name+'/parameters/cf.dat',header = None)
    with open('/scratch/project_2000611/SOMPASOILS/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/'+name+'/parameters/ditch_depth.dat','w') as csv_file:
        for i in range(0,6):
            csv_file.write(str(header[0][i])+"\r\n")
        for dep in peat_loc['depth']:
            csv_file.write(str(dep)+"\r\n")
            
    with open('/scratch/project_2000611/SOMPASOILS/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/'+name+'/parameters/ditch_depth_mean.dat','w') as csv_file:
        for i in range(0,6):
            csv_file.write(str(header[0][i])+"\r\n")
        for dep in peat_loc['depth_mean']:
            csv_file.write(str(dep)+"\r\n")
    #UPDATING DITCH SPACING
    ds = pd.read_csv('/scratch/project_2000611/SOMPASOILS/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/'+name+'/'+nm+'_seg_ds_data.csv')
    #DITCH = ds.loc[ds['gridcode'].isin(peat_loc['gridcode'])]['gridcode','DS']
    peat_loc2 = pd.read_csv('/scratch/project_2000611/SOMPASOILS/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/'+name+'/'+nm+'_peat_loc.csv', header = None)
    DITCH = ds.loc[ds['gridcode'].isin(peat_loc2[0])]
    DIT2 = peat_loc2
    DIT2.columns = ["gridcode"]
    DIT = DIT2.merge(DITCH,on='gridcode',how='left')
    DIT['DSlim2'].fillna("44",inplace = True)
    DIT['DS'].fillna("44",inplace = True)
    with open('/scratch/project_2000611/SOMPASOILS/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/'+name+'/parameters/ditch_spacing.dat','w') as csv_file:
        for i in range(0,6):
            csv_file.write(str(header[0][i])+"\r\n")
        for dit in list(DIT['DSlim2']):
            csv_file.write(str(dit)+"\r\n")
    #UPDATING soil class 
    fert = pd.read_csv('/scratch/project_2000611/SOMPASOILS/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/Updated_site_type_data/'+nm+'_mvmi_metsaan_fertility.csv',sep=";")
    #DITCH = ds.loc[ds['gridcode'].isin(peat_loc['gridcode'])]['gridcode','DS']
    peat_loc3 = pd.read_csv('/scratch/project_2000611/SOMPASOILS/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/'+name+'/'+nm+'_peat_loc.csv', header = None)
    FERTILITY = fert.loc[fert['gridcode'].isin(peat_loc3[0])]
    FERTIL2 = peat_loc3
    FERTIL2.columns = ["gridcode"]
    FERTIL = FERTIL2.merge(FERTILITY,on='gridcode',how='left')
    with open('/scratch/project_2000611/SOMPASOILS/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/'+name+'/parameters/UPDATED_Fertility.dat','w') as csv_file:
        for i in range(0,6):
            csv_file.write(str(header[0][i])+"\r\n")
        for fer in list(FERTIL['fclass.med']):
            csv_file.write(str(fer)+"\r\n")
    peat_loc.to_csv('/scratch/project_2000611/SOMPASOILS/SompaMaps/soil/SpaFHyPeat_Kyle/SpaFHyPeat/'+name+'/parameters/ditch_values.dat')


if __name__ == '__main__':

    import argparse
    import multiprocessing as mp
    import glob
    import time
    start_time = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('--name', help='parameter folder', type=str)
    parser.add_argument('--nm', help='name of subfiles start', type=str)
    args = parser.parse_args()
    print(args.name,args.nm)
    driver(name=args.name, nm=args.nm)
    print("--- %s seconds ---" % (time.time() - start_time))