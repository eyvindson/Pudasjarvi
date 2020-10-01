# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 10:34:37 2016

@author: slauniai & khaahti

"""
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from soilprofile import gwl_Wsto, gwl_Ksat, nan_function
from koordinaattimuunnos import koordTG

eps = np.finfo(float).eps  # machine epsilon
workdir = os.getcwd()

def read_soil_gisdata(fpath, partitions = "", plotgrids=False):
    """
    reads gis-data grids and returns numpy 2d-arrays
    Args:
        fpath - relative path to data folder (str)
        plotgrids - True plots
    Returns:
        gis - dict of gis-data rasters
            cmask
            soil_id
            ditch_depth
            ditch_spacing
    """
    fpath = os.path.join(workdir, fpath)


    # soil classification
    soilclass, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'soil_id.dat'))

    # ditch depth and spacing
    ditch_depth, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'ditch_depth.dat'))
    ditch_spacing, _, _, _, _ = read_AsciiGrid(os.path.join(fpath,'ditch_spacing.dat'))

    # catchment mask cmask[i,j] == 1, np.NaN outside
    if os.path.isfile(os.path.join(fpath, 'cmask.dat')):
        cmask, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'cmask.dat'))
    else:
        cmask = np.ones(np.shape(soilclass))

    # dict of all rasters
    gis = {'cmask': cmask,
           'soilclass': soilclass,
           'ditch_depth': ditch_depth,
           'ditch_spacing': ditch_spacing
           }
    
    for key in gis.keys():
        gis[key] *= cmask

    if plotgrids is True:
        plt.figure()
        plt.subplot(311); plt.imshow(soilclass); plt.colorbar(); plt.title('soiltype')
        plt.subplot(312); plt.imshow(ditch_depth); plt.colorbar();
        plt.title('ditch depth (m)')
        plt.subplot(313); plt.imshow(ditch_spacing); plt.colorbar();
        plt.title('ditch spacing (m)')

    return gis

def read_cpy_gisdata(fpath, plotgrids=False,partitions=""):
    """
    reads gis-data grids and returns numpy 2d-arrays
    Args:
        fpath - relative path to data folder (str)
        plotgrids - True plots
    Returns:
        gis - dict of gis-data rasters
            cmask
            LAI_pine, LAI_spruce - pine and spruce LAI (m2m-2)
            LAI_conif - conifer total annual max LAI (m2m-2)
            LAI_dedid - deciduous annual max LAI (m2m-2)
            cf - canopy closure (-)
            hc - mean stand height (m)

    """
    fpath = os.path.join(workdir, fpath)

    # tree height [m]
    hc, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'hc_'+str(partitions)+'.dat'))

    # canopy closure [-]
    cf, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'cf_'+str(partitions)+'.dat'))

    # leaf area indices
    try:
        print(os.path.join(fpath, 'LAI_pine.dat'))
        LAI_pine, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'LAI_pine_'+str(partitions)+'.dat'))
        LAI_spruce, _, _, _, _ = read_AsciiGrid(os.path.join(fpath,'LAI_spruce_'+str(partitions)+'.dat'))
        print(np.shape(LAI_pine),np.shape(LAI_spruce),np.shape(cf))
        print('HEAD pine',LAI_pine[0,-10:])
        print('HEAD spruce',LAI_spruce[0,-10:])
        LAI_conif = LAI_pine + LAI_spruce
    except:
        LAI_conif, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'LAI_conif_'+str(partitions)+'.dat'))
    LAI_decid, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'LAI_decid_'+str(partitions)+'.dat'))

    # catchment mask cmask[i,j] == 1, np.NaN outside
    if os.path.isfile(os.path.join(fpath, 'cmask.dat')):
        cmask, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'cmask.dat'))
    else:
        cmask = np.ones(np.shape(hc))

    # dict of all rasters
    gis = {'cmask': cmask,
           'LAI_conif': LAI_conif,
           'LAI_decid': LAI_decid, 'hc': hc, 'cf': cf}

    for key in gis.keys():
        gis[key] *= cmask

    if plotgrids is True:

        plt.figure()
        plt.subplot(221); plt.imshow(LAI_pine+LAI_spruce); plt.colorbar();
        plt.title('LAI conif (m2/m2)')
        plt.subplot(222); plt.imshow(LAI_decid); plt.colorbar();
        plt.title('LAI decid (m2/m2)')
        plt.subplot(223); plt.imshow(hc); plt.colorbar(); plt.title('hc (m)')
        plt.subplot(224); plt.imshow(cf); plt.colorbar(); plt.title('cf (-)')

    return gis

def read_forcing_gisdata(fpath, partitions=""):
    """
    reads gis-data grids and returns numpy 2d-arrays
    Args:
        fpath - relative path to data folder (str)
    Returns:
        gis - dict of gis-data rasters
            cmask
            lat
            lon
            forcing_id

    """
    fpath = os.path.join(workdir, fpath)

    # latitude and longitude
    Ncoord, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'Ncoord.dat'))
    Ecoord, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'Ecoord.dat'))
    lat, lon = koordTG(Ncoord, Ecoord)

    forcing_id, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'forcing_id.dat'))

    # catchment mask cmask[i,j] == 1, np.NaN outside
    if os.path.isfile(os.path.join(fpath, 'cmask.dat')):
        cmask, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'cmask.dat'))
    else:
        cmask = np.ones(np.shape(Ncoord))

    # dict of all rasters
    gis = {'cmask': cmask, 'lat': lat, 'lon': lon, 'forcing_id': forcing_id}

    for key in gis.keys():
        gis[key] *= cmask

    return gis

def preprocess_soildata(psp, peatp, gisdata, spatial=True):
    """
    creates input dictionary for initializing SoilGrid
    Args:
        soil parameters
        soiltype parameters
        gisdata
            cmask
            soilclass
        spatial
    """
    # create dict for initializing soil profile.
    # copy pbu into sdata and make each value np.array(np.shape(cmask))
    data = psp.copy()
    data.update((x, y * gisdata['cmask']) for x, y in data.items())

    data.update({'soiltype': np.empty(np.shape(gisdata['cmask']),dtype=object),
                 'depth_id': np.empty(np.shape(gisdata['cmask']),dtype=int)})

    if spatial == False:
        data['soilclass'] = psp['soil_id'] * gisdata['cmask']
    else:
        data['soilclass'] = gisdata['soilclass']
        data['ditch_depth'] = gisdata['ditch_depth']
        data['ditch_spacing'] = gisdata['ditch_spacing']

    data['gwl_to_Ksat'] = np.full(
            len(peatp)*len(np.unique(np.round(data['ditch_depth'],2))),
            nan_function, dtype=object)

    i = 0
    for key, value in peatp.items():
        c = value['soil_id']
        ix = np.where(data['soilclass'] == c)
        data['soiltype'][ix] = key
        # interpolation function between wsto and gwl
        value.update(gwl_Wsto(value['z'], value['pF']))
        # interpolation function between root_wsto and gwl
        value.update(gwl_Wsto(value['z'][:2], {key: value['pF'][key][:2] for key in value['pF'].keys()}, root=True))
        for depth in np.unique(np.round(data['ditch_depth'][ix],2)):
            data['gwl_to_Ksat'][i] = gwl_Ksat(value['z'],
                    value['saturated_conductivity'], depth)
            ixx = np.where((np.round(data['ditch_depth'],2) == depth) &
                           (data['soiltype'] == key))
            data['depth_id'][ixx] = i
            i=i+1

    data['gwl_to_Ksat'] = data['gwl_to_Ksat'][:i]

    data['wtso_to_gwl'] = {soiltype: peatp[soiltype]['to_gwl'] for soiltype in peatp.keys()}
    data['gwl_to_wsto'] = {soiltype: peatp[soiltype]['to_wsto'] for soiltype in peatp.keys()}
    data['gwl_to_rootmoist'] = {soiltype: peatp[soiltype]['to_rootmoist'] for soiltype in peatp.keys()}

    return data

def preprocess_cpydata(pcpy, gisdata, spatial=True):
    """
    creates input dictionary for initializing CanopyGrid
    Args:
        canopy parameters
        gisdata
            cmask
            LAI_pine, LAI_spruce - pine and spruce LAI (m2m-2)
            LAI_conif - conifer total annual max LAI (m2m-2)
            LAI_dedid - deciduous annual max LAI (m2m-2)
            cf - canopy closure (-)
            hc - mean stand height (m)
            (lat, lon)
        spatial
    """
    # inputs for CanopyGrid initialization: update pcpy using spatial data
    cstate = pcpy['state'].copy()

    if spatial:
        cstate['lai_conif'] = gisdata['LAI_conif']
        cstate['lai_decid_max'] = gisdata['LAI_decid']
        cstate['cf'] = gisdata['cf']
        cstate['hc'] = gisdata['hc']
        for key in ['w', 'swe']:
            cstate[key] *= gisdata['cmask']
        if {'lat','lon'}.issubset(gisdata.keys()):
            pcpy['loc']['lat'] = gisdata['lat']
            pcpy['loc']['lon'] = gisdata['lon']
    else:
        for key in cstate.keys():
            cstate[key] *= gisdata['cmask']

    pcpy['state'] = cstate

    return pcpy

def read_FMI_weather(start_date, end_date, sourcefile, CO2=380.0, U=2.0, ID=0):
    """
    reads FMI interpolated daily weather data from file
    IN:
        ID - sve catchment ID. set ID=0 if all data wanted
        start_date - 'yyyy-mm-dd'
        end_date - 'yyyy-mm-dd'
        sourcefile - optional
        CO2 - atm. CO2 concentration (float), optional
        U - wind speed, optional
    OUT:
        fmi - pd.dataframe with datetimeindex
    """

    # OmaTunniste;OmaItä;OmaPohjoinen;Kunta;siteid;vuosi;kk;paiva;longitude;latitude;t_mean;t_max;t_min;
    # rainfall;radiation;hpa;lamposumma_v;rainfall_v;lamposumma;lamposumma_cum
    # -site number
    # -date (yyyy mm dd)
    # -latitude (in KKJ coordinates, metres)
    # -longitude (in KKJ coordinates, metres)
    # -T_mean (degrees celcius)
    # -T_max (degrees celcius)
    # -T_min (degrees celcius)
    # -rainfall (mm)
    # -global radiation (per day in kJ/m2)
    # -H2O partial pressure (hPa)

    sourcefile = os.path.join(sourcefile)

    ID = int(ID)

    # import forcing data
    try:
        fmi = pd.read_csv(sourcefile, sep=';', header='infer',
                          usecols=['OmaTunniste', 'Kunta', 'aika', 'longitude',
                          'latitude', 't_mean', 't_max', 't_min', 'rainfall',
                          'radiation', 'hpa', 'lamposumma_v', 'rainfall_v'],
                          parse_dates=['aika'],encoding="ISO-8859-1")

        fmi = fmi.rename(columns={'aika': 'date',
                                  'OmaTunniste': 'ID',
                                  't_mean': 'air_temperature',
                                  'rainfall': 'precipitation',
                                  'radiation': 'global_radiation',
                                  'hpa': 'h2o'})

        time = pd.to_datetime(fmi['date'], format='%Y%m%d')
    except:
        try:
            fmi = pd.read_csv(sourcefile, sep=';', header='infer',
                              usecols=['x','y','date','temp_avg','prec',
                              'wind_speed_avg','global_rad','vapour_press'],
                              parse_dates=['date'],encoding="ISO-8859-1")

            fmi = fmi.rename(columns={'temp_avg': 'air_temperature',
                                      'prec': 'precipitation',
                                      'global_rad': 'global_radiation',
                                      'vapour_press': 'h2o',
                                      'wind_speed_avg':'wind_speed'})
            time = pd.to_datetime(fmi['date'], format='%Y-%m-%d')
        except:
            raise ValueError('Problem reading forcing data')

    fmi.index = time
    # get desired period and catchment
    fmi = fmi[(fmi.index >= start_date) & (fmi.index <= end_date)]

    if ID > 0:
        fmi = fmi[fmi['ID'] == ID]

    fmi['h2o'] = 1e-1*fmi['h2o']  # hPa-->kPa
    fmi['global_radiation'] = 1e3 / 86400.0*fmi['global_radiation']  # kJ/m2/d-1 to Wm-2
    fmi['par'] = 0.45*fmi['global_radiation']

    # saturated vapor pressure
    esa = 0.6112*np.exp(
            (17.67*fmi['air_temperature']) / (fmi['air_temperature'] + 273.16 - 29.66))  # kPa
    vpd = esa - fmi['h2o']  # kPa
    vpd[vpd < 0] = 0.0
    rh = 100.0*fmi['h2o'] / esa
    rh[rh < 0] = 0.0
    rh[rh > 100] = 100.0

    fmi['RH'] = rh
    fmi['esa'] = esa
    fmi['vapor_pressure_deficit'] = vpd

    fmi['doy'] = fmi.index.dayofyear
    # replace nan's in prec with 0.0
    fmi['precipitation'] = fmi['precipitation'].fillna(0.0)

    # add CO2 and wind speed concentration to dataframe
    if 'CO2' not in fmi:
        fmi['CO2'] = float(CO2)
    if 'wind_speed' not in fmi:
        fmi['wind_speed'] = float(U)

    fmi['wind_speed'] = fmi['wind_speed'].fillna(U)

#    print("NaN values in forcing data:")
#    print(fmi.isnull().any())

    dates = pd.date_range(start_date, end_date).tolist()
    fmi = fmi.drop_duplicates(keep='first')
    if len(dates) != len(fmi):
        print(str(len(dates) - len(fmi)) + ' days missing from forcing file, interpolated')
    forcing = pd.DataFrame(index=dates, columns=[])
    forcing = forcing.merge(fmi, how='outer', left_index=True, right_index=True)
    forcing = forcing.fillna(method='ffill')

    return forcing

def initialize_netcdf(pgen, cmask, filepath, filename, description):
    """
    netCDF4 format output file initialization

    Args:
        variables (list): list of variables to be saved in netCDF4
        cmask
        filepath: path for saving results
        filename: filename
        description: description
    """
    from netCDF4 import Dataset, date2num
    from datetime import datetime

    # dimensions
    date_dimension = None
    i_dimension, j_dimension = np.shape(cmask)

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    ff = os.path.join(filepath, filename)

    # create dataset and dimensions
    ncf = Dataset(ff, 'w')
    ncf.description = 'SpaFHy results : ' + description
    ncf.history = 'created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    ncf.source = 'modified SpaFHy v.1.0'

    ncf.createDimension('date', date_dimension)
    ncf.createDimension('i', i_dimension)
    ncf.createDimension('j', j_dimension)

    date = ncf.createVariable('date', 'f8', ('date',))
    date.units = 'days since 0001-01-01 00:00:00.0'
    date.calendar = 'standard'
    tvec = pd.date_range(pgen['spinup_end'], pgen['end_date']).tolist()[1:]
    date[:] = date2num(tvec, units=date.units, calendar=date.calendar)

    for var in pgen['variables']:

        var_name = var[0]
        var_unit = var[1]

        if (var_name.split('_')[0] == 'forcing' and
            pgen['spatial_forcing'] == False):
            var_dim = ('date')
        elif var_name.split('_')[0] == 'parameters':
            var_dim = ('i', 'j')
        else:
            var_dim = ('date','i', 'j')

        variable = ncf.createVariable(
                var_name, 'f4', var_dim)

        variable.units = var_unit

    return ncf, ff

def write_ncf(results, ncf, steps=None):
    """
    Writes model simultaion results in netCDF4-file

    Args:
        index (int): model loop index
        results (dict): calculation results from group
        ncf (object): netCDF4-file handle
    """

    keys = results.keys()
    variables = ncf.variables.keys()

    for key in keys:

        if key in variables and key != 'date':
            if len(ncf[key].shape) > 2:
                if steps==None:
                    ncf[key][:,:,:] = results[key]
                else:
                    ncf[key][steps[0]:steps[1],:,:] = results[key][0:steps[1]-steps[0],:,:]
            elif len(ncf[key].shape) > 1:
                ncf[key][:,:] = results[key]
            else:
                if steps==None:
                    ncf[key][:] = results[key]
                else:
                    ncf[key][steps[0]:steps[1]] = results[key][0:steps[1]-steps[0]]

def read_AsciiGrid(fname, setnans=True):
    """
    reads AsciiGrid format in fixed format as below:
        ncols         750
        nrows         375
        xllcorner     350000
        yllcorner     6696000
        cellsize      16
        NODATA_value  -9999
        -9999 -9999 -9999 -9999 -9999
        -9999 4.694741 5.537514 4.551162
        -9999 4.759177 5.588773 4.767114
    IN:
        fname - filename (incl. path)
    OUT:
        data - 2D numpy array
        info - 6 first lines as list of strings
        (xloc,yloc) - lower left corner coordinates (tuple)
        cellsize - cellsize (in meters?)
        nodata - value of nodata in 'data'
    Samuli Launiainen Luke 7.9.2016
    """
    import numpy as np

    fid = open(fname, 'r')
    info = fid.readlines()[0:6]
    fid.close()

    # print info
    # conversion to float is needed for non-integers read from file...
    xloc = float(info[2].split(' ')[-1])
    yloc = float(info[3].split(' ')[-1])
    cellsize = float(info[4].split(' ')[-1])
    nodata = float(info[5].split(' ')[-1])

    # read rest to 2D numpy array
    data = np.loadtxt(fname, skiprows=6)

    if setnans is True:
        data[data == nodata] = np.NaN
        nodata = np.NaN

    data = np.array(data, ndmin=2)

    return data, info, (xloc, yloc), cellsize, nodata

def write_AsciiGrid(fname, data, info, fmt='%.18e'):
    """ writes AsciiGrid format txt file
    IN:
        fname - filename
        data - data (numpy array)
        info - info-rows (list, 6rows)
        fmt - output formulation coding

    Samuli Launiainen Luke 7.9.2016
    """
    import numpy as np

    # replace nans with nodatavalue according to info
    nodata = int(info[-1].split(' ')[-1])
    data[np.isnan(data)] = nodata
    # write info
    fid = open(fname, 'w')
    fid.writelines(info)
    fid.close()

    # write data
    fid = open(fname, 'a')
    np.savetxt(fid, data, fmt=fmt, delimiter=' ')
    fid.close()

def read_results(outputfile):
    """
    Opens simulation results netcdf4 dataset in xarray
    Args:
        outputfile (str): outputfilename
    Returns:
        results (xarray): simulation results from given outputfile
    """

    import xarray as xr

    result = xr.open_dataset(outputfile)
    result.coords['i'] = -np.arange(0,result.dims['i'])
    result.coords['j'] = np.arange(0,result.dims['j'])

    return result

def create_input_GIS(fpath, plotgrids=False):
    """

    """
    fpath = os.path.join(workdir, fpath)

    # specific leaf area (m2/kg) for converting leaf mass to leaf area
    SLA = {'pine': 6.8, 'spruce': 4.7, 'decid': 14.0}  # Härkönen et al. 2015 BER 20, 181-195

    # mask, cmask == 1, np.NaN outside
    cmask, _, pos, cellsize, _ = read_AsciiGrid(os.path.join(fpath, 'cmask.asc'))

    # latitude, longitude arrays
    nrows, ncols = np.shape(cmask)
    pos = [round(poss) for poss in pos]
    lon0 = np.arange(pos[0], pos[0] + cellsize*ncols, cellsize)
    lat0 = np.arange(pos[1], pos[1] + cellsize*nrows, cellsize)
    lat0 = np.flipud(lat0)  # indexes count from top left, pos is bottom left

    cmask[np.isfinite(cmask)] = 1.0
    ix = np.where(cmask == 1.0)
    rows = [min(ix[0]), max(ix[0])+2]
    cols = [min(ix[1]), max(ix[1])+2]

    cmask = cmask[rows[0]:rows[-1],cols[0]:cols[-1]]

    lat = lat0[rows[0]:rows[-1]]
    lon = lon0[cols[0]:cols[-1]]

    # peat (only peat areas simulated)
    peat, _, pos, cellsize, _ = read_AsciiGrid(os.path.join(fpath, 'suo_16m.asc'))
    r, c = np.shape(peat)
    pos = [round(poss) for poss in pos]
    ix = r - (lat0[rows] - pos[1]) / cellsize
    iy = (lon0[cols] - pos[0]) / cellsize
    peat = peat[int(ix[0]):int(ix[1]),int(iy[0]):int(iy[-1])]
    peat[np.isfinite(peat)] = 1.0

    cmask = cmask * peat

    # needle/leaf biomass to LAI
    bmleaf_pine, _, pos, cellsize, _ = read_AsciiGrid(os.path.join(fpath, 'bm_manty_neulaset.asc'))
    r, c = np.shape(bmleaf_pine)
    pos = [round(poss) for poss in pos]
    ix = r - (lat0[rows] - pos[1]) / cellsize
    iy = (lon0[cols] - pos[0]) / cellsize
    bmleaf_pine = bmleaf_pine[int(ix[0]):int(ix[1]),int(iy[0]):int(iy[-1])]

    bmleaf_spruce, _, pos, cellsize, _ = read_AsciiGrid(os.path.join(fpath, 'bm_kuusi_neulaset.asc'))
    r, c = np.shape(bmleaf_spruce)
    pos = [round(poss) for poss in pos]
    ix = r - (lat0[rows] - pos[1]) / cellsize
    iy = (lon0[cols] - pos[0]) / cellsize
    bmleaf_spruce = bmleaf_spruce[int(ix[0]):int(ix[1]),int(iy[0]):int(iy[-1])]

    bmleaf_decid, _, pos, cellsize, _ = read_AsciiGrid(os.path.join(fpath, 'bm_lehtip_neulaset.asc'))
    r, c = np.shape(bmleaf_decid)
    pos = [round(poss) for poss in pos]
    ix = r - (lat0[rows] - pos[1]) / cellsize
    iy = (lon0[cols] - pos[0]) / cellsize
    bmleaf_decid = bmleaf_decid[int(ix[0]):int(ix[1]),int(iy[0]):int(iy[-1])]

    LAI_pine = 1e-3*bmleaf_pine*SLA['pine']  # 1e-3 converts 10kg/ha to kg/m2
    LAI_spruce = 1e-3*bmleaf_spruce*SLA['spruce']
    LAI_decid = 1e-3*bmleaf_decid*SLA['decid']

    # tree height
    hc, _, pos, cellsize, _ = read_AsciiGrid(os.path.join(fpath, 'keskipituus.asc'))
    r, c = np.shape(hc)
    pos = [round(poss) for poss in pos]
    ix = r - (lat0[rows] - pos[1]) / cellsize
    iy = (lon0[cols] - pos[0]) / cellsize
    hc = hc[int(ix[0]):int(ix[1]),int(iy[0]):int(iy[-1])]
    hc = 0.1*hc # m

    # canopy closure
    cf, _, pos, cellsize, _ = read_AsciiGrid(os.path.join(fpath, 'latvuspeitto.asc'))
    r, c = np.shape(cf)
    pos = [round(poss) for poss in pos]
    ix = r - (lat0[rows] - pos[1]) / cellsize
    iy = (lon0[cols] - pos[0]) / cellsize
    cf = cf[int(ix[0]):int(ix[1]),int(iy[0]):int(iy[-1])]
    cf = 1e-2*cf

    # ditch depth
    ditch_depth = cmask * 0.8

    # ditch spacing
    ditch_spacing, _, pos, cellsize, _ = read_AsciiGrid(os.path.join(fpath, 'ditch_spacing.asc'))
    r, c = np.shape(ditch_spacing)
    pos = [round(poss) for poss in pos]
    ix = r - (lat0[rows] - pos[1]) / cellsize
    iy = (lon0[cols] - pos[0]) / cellsize
    ditch_spacing = ditch_spacing[int(ix[0]):int(ix[1]),int(iy[0]):int(iy[-1])]
    ditch_spacing = np.minimum(200., ditch_spacing)
    ditch_spacing = np.maximum(20., ditch_spacing)

    # soil_id
    soil_id = cmask * 2.0

    nrows, ncols = np.shape(cmask)
    info = ['ncols         ' + str(nrows) + '\n',
            'nrows         ' + str(ncols) + '\n',
            'xllcorner     ' + str(lon[0]) + '\n',
            'yllcorner     ' + str(lat[0]) + '\n',
            'cellsize      ' + str(cellsize) + '\n',
            'NODATA_value  -9999\n']

    # dict of all rasters
    GisData = {'cmask': cmask, 'ditch_spacing': ditch_spacing * cmask, 'ditch_depth': ditch_depth * cmask,
               'LAI_pine': LAI_pine * cmask, 'LAI_spruce': LAI_spruce * cmask, 'LAI_decid': LAI_decid * cmask,
               'hc': hc * cmask, 'cf': cf * cmask, 'soil_id': soil_id}

    if plotgrids is True:
        xy = np.meshgrid(lon, lat)

        plt.figure(99,figsize=(12, 12))
        i=1
        for key, gdata in GisData.items():
            if key != 'cmask':
                plt.subplot(3,3,i)
                plt.pcolor(xy[0], xy[1], gdata)
                plt.colorbar()
                plt.title(key)
                i+=1

    fpath = os.path.join(fpath, 'inputs')
    if os.path.isdir(fpath) == False:
        os.mkdir(fpath)

    for key, gdata in GisData.items():
        write_AsciiGrid(os.path.join(fpath, key + '.dat'), gdata, info, fmt='%.6e')


def rw_FMI_files(sourcefiles, out_path, plot=False):
    """
    reads and writes FMI interpolated daily weather data
    """
    frames = []
    for sourcefile in sourcefiles:
        sourcefile = os.path.join(sourcefile)

        # import forcing data
        try:
            fmi = pd.read_csv(sourcefile, sep=',', header='infer',index_col=False,
                              usecols=['pvm','latitude','longitude','t_mean','t_max','t_min',
                                       'rainfall','radiation','hpa','site'],
                              parse_dates=['pvm'],encoding="ISO-8859-1")

            fmi = fmi.rename(columns={'pvm': 'date',
                                      't_mean': 'temp_avg',
                                      't_max': 'temp_max',
                                      't_min': 'temp_min',
                                      'rainfall': 'prec',
                                      'radiation': 'global_rad',
                                      'hpa': 'vapour_press',
                                      'longitude':'x',
                                      'latitude':'y'})
            fmi = fmi[fmi['date']<'2016-07-03']
        except:
            try:
                fmi = pd.read_csv(sourcefile, sep=',', header='infer',index_col=False,
                                  usecols=['x','y','date','temp_avg','temp_min','temp_max',
                                           'prec', 'wind_speed_avg','global_rad','vapour_press',
                                           'snow_depth','pot_evap','site'],
                                  parse_dates=['date'],encoding="ISO-8859-1")

                fmi = fmi.rename(columns={})
            except:
                raise ValueError('Problem reading forcing data')

        time = pd.to_datetime(fmi['date'], format='%Y-%m-%d')
        fmi.index=time

        frames.append(fmi.copy())

    fmi = pd.concat(frames, sort=False)

    sites = list(set(fmi['site']))
    sites.sort()
    index = 0
    readme = 'Indices of weather files'
    for site in sites:
        fmi[fmi['site']==site].to_csv(path_or_buf=out_path + 'weather_id_' + str(index) + '.csv', sep=';', na_rep='NaN', index=False)
        readme += '\n'+ str(index) + ':' + site
        index+=1
        if plot:
            fmi[fmi['site']==site].plot(subplots=True)
    outF = open(out_path + "weather_readme.txt", "w")
    print(readme, file=outF)
    outF.close()
    return fmi

