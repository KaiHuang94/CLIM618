import xarray as xr 
import matplotlib.pyplot as plt 
import numpy as np 
from numpy import linalg as LA
from scipy.stats import linregress as lr

import cartopy.crs as ccrs 
import cartopy.mpl.ticker as cticker 
from cartopy.util import add_cyclic_point

def label_latlon(ax,lons,lats):
    
    # Longitude labels
    ax.set_xticks(lons, crs=ccrs.PlateCarree())
    lon_formatter = cticker.LongitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)

    # Latitude labels
    ax.set_yticks(lats, crs=ccrs.PlateCarree())
    lat_formatter = cticker.LatitudeFormatter()
    ax.yaxis.set_major_formatter(lat_formatter)
    
def MonthlyClimatology(ds):
    ds_clm    = ds.groupby('time.month').mean()
    return ds_clm

def DailyClimatology(ds):
    ds_clm    = ds.groupby('time.dayofyear').mean()
    return ds_clm

def DailyAnomaly(ds):
    ds_clm    = ds.groupby('time.dayofyear').mean()
    ds_anom   = ds.groupby('time.dayofyear')-ds_clm
    return ds_anom

def LinearRegression_Map(v,p,ndimlat,ndimlon):
    p_array = np.zeros((ndimlat,ndimlon))
    r_array = np.zeros((ndimlat,ndimlon))
    m_array = np.zeros((ndimlat,ndimlon))

    for i in range(ndimlon):
        for j in range(ndimlat):
            x=p
            y=v[:,j,i]
            
            m,b,r,p,e=lr(x,y)
        
            m_array[j,i]=m
            r_array[j,i]=r
            p_array[j,i]=p 

    return [m_array,r_array,p_array]

def rmm(ds_u850,ds_u200,ds_olr):
    latS = -15 
    latN =  15
    
    #compute the anomalies   
    ds_u850a = DailyAnomaly(ds_u850)
    ds_u850a = ds_u850a.sel(lat=slice(latS,latN))
    ds_u200a = DailyAnomaly(ds_u200)
    ds_u200a = ds_u200a.sel(lat=slice(latS,latN))
    ds_olra  = DailyAnomaly(ds_olr)
    ds_olra  = ds_olra.sel(lat=slice(latS,latN))
    
    #extract dimension and variable values
    ntime = len(ds_u850a['time'])
    nlat  = len(ds_u850a['lat'])
    nlon  = len(ds_u850a['lon'])
    
    u850a = ds_u850a['u'].values
    u200a = ds_u200a['u'].values 
    olra  = ds_olra['olr'].values
    
    # take out the running mean of previous 120 days
    u850sm = np.zeros((ntime,nlat,nlon))
    u200sm = np.zeros((ntime,nlat,nlon))
    olrsm  = np.zeros((ntime,nlat,nlon))
    
    u850sm[:120,:,:] = u850a[:120,:,:]
    u200sm[:120,:,:] = u200a[:120,:,:]
    olrsm[:120,:,:]  = olra[:120,:,:]
    
    for i in range(120,ntime):
        u850sm[i,:,:] = u850a[i,:,:]-np.mean(u850a[i-120:i-1,:,:],axis=0)
        u200sm[i,:,:] = u200a[i,:,:]-np.mean(u200a[i-120:i-1,:,:],axis=0)
        olrsm[i,:,:]  = olra[i,:,:]-np.mean(olra[i-120:i-1,:,:],axis=0)
    
    # meridional average
    u850 = np.mean(u850sm,axis=1)
    u200 = np.mean(u200sm,axis=1)
    olr  = np.mean(olrsm,axis=1)
    
    # compute the temporal variance at each longitude
    var_u850 = np.var(u850,axis=0)
    var_u200 = np.var(u200,axis=0)
    var_olr  = np.var(olr,axis=0)
    
    # compute the zonal mean of the temporal variance 
    zavg_var_u850 = np.mean(var_u850)
    zavg_var_u200 = np.mean(var_u200)
    zavg_var_olr  = np.mean(var_olr)
    
    #normalize by sqrt(avg_var)
    olr  = olr/np.sqrt(zavg_var_olr)
    u850 = u850/np.sqrt(zavg_var_u850)
    u200 = u200/np.sqrt(zavg_var_u200)
    
    #combine the normalized data into one variable 
    cdata = np.zeros((3*nlon,ntime))
    for i in range(nlon):
        cdata[i,:] = olr[:,i]
        cdata[i+nlon,:]=u850[:,i]
        cdata[i+2*nlon,:]=u200[:,i]
    
    #calculate EOFs 
    C=np.cov(cdata)
    eigenvalues, eigenvectors = LA.eig(C)
    
    #sort the eigenvalues and eigenvectors from high to low
    idx=eigenvalues.argsort()[::-1] 
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:,idx]
    
    #get the EOF patterns
    EOF=eigenvectors
    
    nvar = 3 
    neof = 2 
    ceof_bk = np.zeros((nvar,neof,nlon))

    for i in range(neof):
        ceof_bk[0,i,:] = EOF[0:nlon,i]
        ceof_bk[1,i,:] = EOF[nlon:2*nlon,i]
        ceof_bk[2,i,:] = EOF[2*nlon:,i]
    
    ceof = np.zeros((nvar,neof,nlon))
    ceof[:,0,:] = -ceof_bk[:,1,:]
    ceof[:,1,:] = ceof_bk[:,0,:]
    
    # PC Time Series
    PC=np.dot(cdata.T,EOF)
    PC1 = -PC[:,1]
    PC1 = PC1/np.std(PC1)
    PC2 = PC[:,0]
    PC2 = PC2/np.std(PC2)
    
    # MJO index
    MJO_index_bk = PC1**2+PC1**2
    nt = len(MJO_index_bk)
    MJO_index = np.zeros(nt)

    for i in range(45):
        MJO_index[i] = np.mean(MJO_index_bk[i:i+91])

    for i in range(nt-45,nt):
        MJO_index[i] = np.mean(MJO_index_bk[i-90:i+1])
    
    for i in range(45,nt-45):
        MJO_index[i] = np.mean(MJO_index_bk[i-45:i+46])
        
    # Explained Variance
    vexp=(eigenvalues/np.sum(eigenvalues))*100
    var1=vexp[0]
    var2=vexp[1]
    
    ds_ceof=xr.DataArray(ceof,
                         coords={'var':np.arange(1,nvar+1,1),
                                 'eof':np.arange(1,neof+1,1),
                                 'lon': ds_olra['lon']},
                         dims=['var','eof','lon'])        
    ds_ceof=ds_ceof.to_dataset(name='ceof')
    
    ds_PC1 =xr.DataArray(PC1,
                         coords={'time':ds_olra['time']},
                         dims=['time'])
    ds_PC1=ds_PC1.to_dataset(name='PC1')
    
    ds_PC2 =xr.DataArray(PC2,
                         coords={'time':ds_olra['time']},
                         dims=['time'])
    ds_PC2=ds_PC2.to_dataset(name='PC2')
    
    ds_MJO_index =xr.DataArray(MJO_index,
                         coords={'time':ds_olra['time']},
                         dims=['time'])
    ds_MJO_index=ds_MJO_index.to_dataset(name='MJO_INDEX')
    
    ds_rmm=xr.merge([ds_ceof,ds_PC1,ds_PC2,ds_MJO_index])
    
    return [ds_rmm,var1,var2]