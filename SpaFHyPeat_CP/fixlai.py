import pandas as pd

def df_to_LAI_file(df,file_name,nrows):
    #Header for LAI files
    header1 = "ncols 1\n"
    header2 = "nrows "+str(nrows)+"\n"
    header3 =  "xllcorner 0\nyllcorner 0\ncellsize 16\nNODATA_value -9999\n"
    f=open(file_name,'w')
    f.write(header1+header2+header3)
    f.close()
    df.to_csv(file_name,mode='a',columns=['LAI'],header=False,index=False)
    
def fix_lai(region:str):
    """NFI data may produce too big LAI values. Fix LAI to max 9.0"""
    f_pine = 'LAI_pine.dat'
    f_spruce = 'LAI_spruce.dat'
    f_decid = 'LAI_decid.dat'
    #Collect the fixed LAI values by region to LAI_fixed/lappi etc.
    outdir = 'LAI_fixed/'+region+'/'
    #Read files. Name columns, easier to follow indexing
    df_pine = pd.read_csv(region+'/parameters/'+f_pine,delim_whitespace=True,header=None,skiprows=6)
    df_pine.columns=['LAI']
    df_spruce =  pd.read_csv(region+'/parameters/'+f_spruce,delim_whitespace=True,header=None,skiprows=6)
    df_spruce.columns=['LAI']
    df_decid =  pd.read_csv(region+'/parameters/'+f_decid,delim_whitespace=True,header=None,skiprows=6)
    df_decid.columns=['LAI']
    #Remember the number of rows for LAI file header
    nrows = len(df_pine.index)
    #Calculate the sum and find LAI values >= 9.0
    df_lai = df_pine+df_spruce+df_decid
    df_big_lai = df_lai[df_lai['LAI']>=9.0]
    #Now, find the LAI propotion of the tree species
    #Note: df_big_lai.index has the row indexes that contain the big LAI values in each tree species.
    df_pine_share = df_pine.loc[df_big_lai.index,:]/df_big_lai
    df_spruce_share = df_spruce.loc[df_big_lai.index,:]/df_big_lai
    df_decid_share = df_decid.loc[df_big_lai.index,:]/df_big_lai
    #Just to check shares add to 1
    #df_tot_share = df_pine_share+df_spruce_share+df_decid_share
    #Now calculate the new LAI value assumiing LAI max == 9.0
    df_pine.loc[df_pine_share.index,:]=df_pine_share*9.0
    df_spruce.loc[df_spruce_share.index,:]=df_spruce_share*9.0
    df_decid.loc[df_decid_share.index,:]=df_decid_share*9.0
    #Write the files
    df_to_LAI_file(df_pine,outdir+f_pine,nrows)
    df_to_LAI_file(df_spruce,outdir+f_spruce,nrows)
    df_to_LAI_file(df_decid,outdir+f_decid,nrows)
    return (df_pine,df_spruce,df_decid)
