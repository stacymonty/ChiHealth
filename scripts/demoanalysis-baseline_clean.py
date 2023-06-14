#demo_analysis-baseline_clean.py
# calculate exposure

#---------------------------------------------------------------------
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import netCDF4
import math
from scipy.interpolate import griddata
import scipy.stats as st
import cartopy.feature as cfeature 
from cartopy import crs as ccrs;
from shapely.ops import unary_union, cascaded_union
from geopandas.tools import sjoin
from shapely.geometry import Point, shape, Polygon
from cartopy import crs as ccrs;
import geopandas as gpd
from scipy import ndimage, misc
import seaborn as sns
from scipy import ndimage, misc
import matplotlib.path as mpath;
from cartopy.io.shapereader import Reader
import string
import geopandas as gpd
#---------------------------------------------------------------------

def mask_given_shapefile(lon,lat,shapefile):
   '''
   Make a mask given a shapefile
   lon - array of grid lons
   lat - array of grid lats
   shapefile - geopandas geodataframe of a shapefile, needs geometry column
   '''
   union=gpd.GeoSeries(unary_union(shapefile.geometry))
   mask=np.ones(lon.shape,dtype=bool)
   mask[:] = False
   for i in range(len(lon)):
      for j in range(len(lon[0])):
         pt = Point(lon[i][j],lat[i][j])
         mask[i][j] =  pt.within(union[0])
   #
   return mask
#---------------------------------------------------------------------

#------------------------------------------------------------------------
# Start of new analysis
#------------------------------------------------------------------------

# configure the data
#doe = gpd.read_file('/projects/b1045/montgomery/Paper2/DOE_J40_DAC_Shapefiles_v2022c/DOE_J40_DAC_Shapefiles_v2022c.shp')
#doe = doe.drop('geometry',axis=1)
#ac2 = 
#ac = gpd.read_file('/projects/b1045/montgomery/Paper2/baseline_demos_update_v4.shp')
#ac = gpd.GeoDataFrame(pd.merge(doe,ac,on = 'GEOID'),crs = ac.crs)
#$chi = chi.to_file('chi_with_mortality_asthma.geojson',driver='GeoJSON')

chi = gpd.read_file('chi_with_mortality_asthma.geojson')

chi['Over65'] = [int(chi.Over65[i]) for i in range(len(chi))]
chi['Under19'] = [int(chi.Under19[i]) for i in range(len(chi))]

# make some more features
chi['pOver65'] = chi['Over65']/chi['total_pop']
chi['pUnder19'] = chi['Under19']/chi['total_pop']
chi['pvacancy'] = (chi.total_hous-chi.total_occu)/chi.total_hous

popgroup = ['pop_30-35', 'pop_35-44', 'pop_45-54', 'pop_55-64', 'pop_65-74', 'pop_75-84', 'pop_85+']
mortgroup=[ 'mort_25-34', 'mort_35-44', 'mort_45-54', 'mort_55-64', 'mort_65-74', 'mort_75-84', 'mort_85-99']

def make_mortall_weighted(chi,popgroup,mortgroup,pre):
#for i in range(1):
   A = chi[popgroup]; B = chi[mortgroup]
   A.columns = B.columns
   chi[pre+'popweight-mortall'] = np.sum(A*B,axis=1)/np.sum(A,axis=1)
   return chi

chi = make_mortall_weighted(chi,popgroup,mortgroup,'')
#msa = make_mortall_weighted(msa,popgroup,mortgroup,'')

#make race-weighted mortality rates
chi = make_mortall_weighted(chi,Bpopgroup,mortgroup,'B')


#cook_shapefile = 

############################################################################################################################################
# Make categories
# what limits
ul = 0.9
ll = 0.1

def make_cats(shp,ul,ll):
   values = []
   # make cateogries to check highest concs
   c = ['total_pop', 'pWhite', 'pBlack', 'pLatino','pAsian','income_per','med_income','jobaccess','trnsptbrdn','no2','pm25_y','mdao3']
   c = ['pOver65', 'pUnder19', 'total_pop', 'pWhite', 'pBlack', 'pLatino','med_income','incpc_all','incpc_w','incpc_bl','incpc_lat','urban_flood_suscep', 'heatisl','trees_n','PubAs-SNAP','svi_pecentile','ndvi','pvacancy','med_age','no2','pm25','mdao3']
   c = c+['asthRT_no2','asthRT_pm25','mortRT_pm25','mortRT_no2','mortRT_o3','asthma_rt','no2','pm25_y','mdao3']
   #c = ['med_age']
   for i in range(len(c)):
      # for each category c
      hi = shp[c[i]] > np.nanquantile(shp[c[i]],ul)
      lo = shp[c[i]] < np.nanquantile(shp[c[i]],ll)
      print(np.nanquantile(shp[c[i]],ul))
      print(np.nanquantile(shp[c[i]],ll))
      hilo = []
      for t in range(len(shp)):
         if hi[t] == True: hilo.append('high')
         if lo[t] == True: hilo.append('low')
         if (hi[t] == False) & (lo[t] == False): hilo.append(np.nan)
      #
      shp['cat_'+c[i]] = hilo
      values.append([np.nanquantile(shp[c[i]],.9),np.nanquantile(shp[c[i]],.1),shp[shp['cat_'+c[i]]=='high'][c[i]].mean(),shp[shp['cat_'+c[i]]=='low'][c[i]].mean(),shp[c[i]].mean()])
      #
   values = pd.DataFrame(values)
   values.columns = ['90%ile','10%ile','high_avg','low_avg','avg']
   values['cat'] = c
   return shp

chi = make_cats(chi,0.9,0.1)

####################################################################################
# make multiple categories
lims = [0.2,0.4,0.6,0.8]

def make_cats_multi(shp,lims):
   vals = []
   # make cateogries to check highest concs
   #c = ['total_pop', 'pWhite', 'pBlack', 'pLatino','income_per','med_income','jobaccess','trnsptbrdn','no2','pm25_y','mdao3']
   #c= ['incpc_w','incpc_bl','incpc_lat','incpc_as','PubAs-SNAP','svi_pecentile','ndvi','pvacancy',]
   #c = ['pOver65', 'pUnder19', 'total_pop', 'pWhite', 'pBlack', 'pLatino','med_income','incpc_all','med_age','PubAs-SNAP','svi_pecentile','ndvi','no2','pm25','mdao3']
   c = ['no2','pm25','mdao3','mortRT_pm25','mortRT_no2','mortRT_o3','asthma_rt','med_income','incpc_all','med_age']+['Under19ER_no2','Under19ER_pm25','Under19ER_o3','Under19ER_no2_rt','Under19ER_pm25_rt','Under19ER_o3_rt']+['asth_no2','asth_pm25','asth_o3','mort_pm25','mort_no2','mort_o3','asthRT_o3','asthRT_pm25','asthRT_no2','pWhite', 'pBlack', 'pLatino','pAsian']
   #c = ['med_age']
   #c = ['no2','pm25','mdao3','mortRT_pm25','mortRT_no2','mortRT_o3']
   for i in range(len(c)):
      # for each category c
      cats = []
      cats.append((shp[c[i]] < np.nanquantile(shp[c[i]],lims[0])) & ~np.isnan(shp[c[i]]))
      for j in range(1,len(lims)):
         cats.append((shp[c[i]] < np.nanquantile(shp[c[i]],lims[j])) & (shp[c[i]] > np.nanquantile(shp[c[i]],lims[j-1])))
      #
      cats.append((shp[c[i]] > np.nanquantile(shp[c[i]],lims[len(lims)-1])))
      hilo = [cats[t]*t for t in range(len(cats))]
      #
      shp['mcat_'+c[i]] = np.sum(hilo,axis = 0)
      vals.append([c[i],np.min(shp[c[i]])]+[np.nanquantile(shp[c[i]],lims[t]) for t in range(len(lims))]+[np.max(shp[c[i]])])
   return shp,pd.DataFrame(vals)

chi,vals = make_cats_multi(chi,lims)

####################################################################################
# Check in pollution areas in Chicago # reshape
c = ['pWhite','pBlack','pLatino','income_per']

c = ['income_per','incpc_w','incpc_bl','incpc_lat']
cat = ['cat_no2','cat_pm25_y','cat_mdao3']

chi['pWhite'] = chi['pWhite'] /100
chi['pBlack'] = chi['pBlack'] /100
chi['pAsian'] = chi['pAsian'] /100
chi['pLatino'] = chi['pLatino'] /100

c = ['med_income','svi_pecentile','total_pop','pWhite','pBlack','pLatino','pAsian']
titles2 = ['Median Income ($)','SVI','Population','White (%)','Black (%)','Hispanic + Latino (%)','Asian (%)']
cat = ['cat_no2','cat_pm25','cat_mdao3']

titles = [r'NO$_2$',r'PM$_{2.5}$',r'O$_3$']
vmin = [-.05,-.05,-.05,-.05,-.05,-.05,-.05]
vmax=[160000,105,12000,105,105,105,105]
alphabet = np.array(list(string.ascii_lowercase)[:len(cat)*(len(c)+1)]).reshape([len(cat),len(c)+1])
alphabet = alphabet.T

def check_pollution_dist_over_shp(shp,figname,save=False,show=True):
#for i in range(1):
   fig,ax = plt.subplots(len(cat),len(c)+1,figsize=(16,7))
   ax=ax.T
   for k in range(1,len(c)+1):
      for j in range(len(cat)):
         i = k-1
         hi = shp[shp[cat[j]] == 'high'].reset_index(drop=True)
         lo = shp[shp[cat[j]] == 'low'].reset_index(drop=True)
         hilo = pd.DataFrame([hi[c[i]],lo[c[i]]]).T
         hilo.columns = ['High '+titles[j],'Low '+titles[j]]
         sns.boxplot(data=hilo,ax=ax[k][j], palette="Set1", width = 0.5)
         ax[k][j].set_ylim(vmin[i],vmax[i])
         sns.despine(offset=10, trim=True,ax=ax[k][j])
         ax[k][j].text(0.02, 0.9,'('+alphabet[k][j]+')',transform=ax[k][j].transAxes,fontsize = 10,c='k')
         #
         if j == 0: ax[k][j].set_title(titles2[i])
         #if i == 0: ax[i][j].set_ylabel(tiltes[i])
#
   for j in range(len(cat)):
      hi = shp[shp[cat[j]] == 'high'].reset_index(drop=True)
      lo = shp[shp[cat[j]] == 'low'].reset_index(drop=True)
      chi_union.plot(facecolor='gray',ax=ax[0][j])
      hi.plot(ax=ax[0][j],facecolor='#FF0033');
      lo.plot(ax=ax[0][j],facecolor='#336699');
      ax[0][j].axis('off')
      ax[0][j].text(0.02, 0.95,'('+alphabet[0][j]+')',transform=ax[0][j].transAxes,fontsize = 10,c='k')
      ax[0][j].text(0.02, 0.5,titles[j],transform=ax[0][j].transAxes,fontsize = 12,c='k',rotation=90)
      #ax[-1][j].set_title(alphabet[-1][j],ha='left')
   #
   plt.tight_layout()
   if save: plt.savefig(figname,dpi=300)
   if show: plt.show()


#chi = chi.to_crs(orig_crs)
#msa = msa.to_crs(orig_crs)

check_pollution_dist_over_shp(chi,'exposure_chi_soc.png',save=True,show=True)

####################################################################################
# Check in pollution areas in Chicago # reshape
c = ['pWhite','pBlack','pLatino','income_per']

c = ['income_per','incpc_w','incpc_bl','incpc_lat']
cat = ['cat_no2','cat_pm25_y','cat_mdao3']

#chi['pWhite'] = chi['pWhite'] *100
#chi['pBlack'] = chi['pBlack'] *100
#chi['pAsian'] = chi['pAsian'] *100
#chi['pLatino'] = chi['pLatino'] *100

c = ['med_income','svi_pecentile','total_pop','pWhite','pBlack','pLatino','pAsian']
titles = ['Median Income ($)','SVI','Population','White (%)','Black (%)','Hispanic + Latino (%)','Asian (%)']

titles = ['Income','SVI','Population','White','Black','Hispanic+Latino','Asian']
cat = ['cat_no2','cat_pm25','cat_mdao3']


def check_pollution_dist_over_shp(shp,alphabet,titles2,titles,vmin,vmax,cat,c,figname,save=False,show=True):
#for i in range(1):
   fig,ax = plt.subplots(len(cat),len(c)+1,figsize=(12,8.5))
   ax=ax.T
   for k in range(1,len(c)+1):
      for j in range(len(cat)):
         i = k-1
         hi = shp[shp[cat[j]] == 'high'].reset_index(drop=True)
         lo = shp[shp[cat[j]] == 'low'].reset_index(drop=True)
         hilo = pd.DataFrame([hi[c[i]],lo[c[i]]]).T
         hilo.columns = ['High '+titles[j],'Low '+titles[j]]
         sns.boxplot(data=hilo,ax=ax[k][j], palette="Set1", width = 0.5)
         ax[k][j].set_ylim(vmin[i],vmax[i])
         sns.despine(offset=10, trim=True,ax=ax[k][j])
         ax[k][j].text(0.02, 0.9,'('+alphabet[k][j]+')',transform=ax[k][j].transAxes,fontsize = 10,c='k')
         if j == 0: ax[k][j].set_title(titles2[i])
         #if i == 0: ax[i][j].set_ylabel(tiltes[i])
#
   for j in range(len(cat)):
      hi = shp[shp[cat[j]] == 'high'].reset_index(drop=True)
      lo = shp[shp[cat[j]] == 'low'].reset_index(drop=True)
      chi_union.plot(facecolor='gray',ax=ax[0][j])
      hi.plot(ax=ax[0][j],facecolor='#FF0033');
      lo.plot(ax=ax[0][j],facecolor='#336699');
      ax[0][j].axis('off')
      #ax[0][j].set_title(cat[j])
      ax[0][j].text(0.02, 0.95,'('+alphabet[0][j]+')',transform=ax[0][j].transAxes,fontsize = 10,c='k')
      ax[0][j].text(0.02, 0.1,titles3[j],transform=ax[0][j].transAxes,fontsize = 12,c='k',rotation=90)
   #
   plt.tight_layout()
   if save: plt.savefig(figname,dpi=300)
   if show: plt.show()


#chi = chi.to_crs(orig_crs)
#msa = msa.to_crs(orig_crs)
#titles = ['\nWhite','\nBlack','\nHispanic+\nLatino','\nAsian']
#titles3 = ['White','Black','Hispanic+\nLatino','Asian']
titles3 = ['Income','SVI','Population']
titles= ['\nIncome','\nSVI','\nPopulation']

c = ['no2','pm25','mdao3']
ci = ['pWhite','pBlack','pLatino','pAsian']
ci = ['med_income','svi_pecentile','total_pop']
cat = ['cat_'+ci[i] for i in range(len(ci))]

alphabet = np.array(list(string.ascii_lowercase)[:len(cat)*(len(c)+1)]).reshape([len(cat),len(c)+1])
alphabet = alphabet.T

titles2 = [r'NO$_2$',r'PM$_{2.5}$',r'O$_3$']
vmin = [10,7,35,-.05,-.05,-.05]
vmax = [25,15,50]
check_pollution_dist_over_shp(chi,alphabet,titles2,titles,vmin,vmax,cat,c,'exposure_chi_soccat_second.png',save=True,show=True)

chi['pWhite'] =  chi.White/chi.total_pop
chi['pBlack'] =  chi.Black/chi.total_pop
chi['pLatino'] =  chi.Latino/chi.total_pop
chi['pAsian'] =  chi.Asian/chi.total_pop


c = ['pWhite','pBlack','pLatino','pAsian']
cat = ['cat_asthRT_no2','cat_asthRT_pm25','cat_mortRT_no2','cat_mortRT_pm25']
titles = [r'Asthma NO$_2$',r'Asthma PM$_{2.5}$',r'Mort  NO$_2$',r'Mortality PM$_{2.5}$']
titles2 = ['White','Black','Hispanic+Latino','Asian']
titles3 = ['Asthma \n'+r'NO$_2$','Asthma \n'+r'PM$_{2.5}$','Mortality \n'+r'NO$_2$','Mortality \n'+r'PM$_{2.5}$']
titles=titles3
vmin = [0]*len(cat)
vmax= [1]*len(cat)
alphabet = np.array(list(string.ascii_lowercase)[:len(cat)*(len(c)+1)]).reshape([len(cat),len(c)+1])
alphabet = alphabet.T
check_pollution_dist_over_shp(chi,alphabet,titles2,titles,vmin,vmax,cat,c,'exposure_chi_soccat_health.png',save=True,show=True)




########################################################################################################################################################################
# Make pollution for multiple cats

def check_pollution_dist_over_shp_multicat(shp,alphabet,titles2,titles,vmins,vmaxs,cat,c,palette,vals,figname,save=False,show=True):
   #for i in range(1):
   # titles2 == along the top
   # titles = title on the botto of plot
   # titles3 = title on the side
#for z in range(1):
   col,row = len(cat),len(c)
   fig,ax = plt.subplots(len(c),len(cat),figsize=(2.7*col,2.2*row))
   ax=ax.T
   hilos = []
   for k in range(len(cat)):
      for j in range(len(c)):
         i = k
         vmin,vmax = vmins[j],vmaxs[j]
         #print('here')
         #
         catk = cat[i]
         print(catk)
         qlabel = np.array(vals[vals[0]==catk.split('mcat_')[1]]).ravel()[1:]*100
         print(qlabel)
         h4 = shp[shp[catk] == 0].reset_index(drop=True)
         h3 = shp[shp[catk] == 1].reset_index(drop=True)
         h2 = shp[shp[catk] == 2].reset_index(drop=True)
         h1 = shp[shp[catk] == 3].reset_index(drop=True)
         h0 = shp[shp[catk] == 4].reset_index(drop=True)
         hilo = pd.DataFrame([h4[c[j]],h3[c[j]],h2[c[j]],h1[c[j]],h0[c[j]]]).T
         #print('here1')
         hilo.columns = ['Q1','Q2','Q3','Q4','Q5']
         #print(hilo)
         if j == len(c)-1: hilo.columns = ['Q1\n    '+str(round(qlabel[0],1)),'Q2\n    '+str(round(qlabel[1],1)),'Q3\n     '+str(round(qlabel[2],1)),'Q4\n      '+str(round(qlabel[3],1)),'Q5']#\n >'+str(round(qlabel[3],1))]
         #print('here2')
         sns.boxplot(data=hilo,ax=ax[k][j], palette=palette, width = 0.5,linewidth=1.5,fliersize=1,flierprops={"marker": "o"})
         ax[k][j].set_ylim(vmin,vmax)
         sns.despine(offset=10, trim=True,ax=ax[k][j])
         #print('here3')
         ax[k][j].text(-0.03, 0.95,'('+alphabet[k][j]+')',transform=ax[k][j].transAxes,fontsize = 10,c='k')
         #print('here4')
         if j == 0: ax[k][j].set_title(titles2[i])
         #print('here5')
         #ax[0][j].text(0.02, 0.1,titles3[j],transform=ax[0][j].transAxes,fontsize = 12,c='k',rotation=90)
         if i == 0: ax[i][j].set_ylabel(titles[j])
         #print('here3')
         hilos.append(hilo)
   plt.tight_layout()
   if save: plt.savefig(figname,dpi=300)
   if show: plt.show()
   return hilos

c = ['pWhite','pBlack','pLatino','pAsian']
cat = ['mcat_asthRT_no2','mcat_asthRT_pm25','mcat_mortRT_no2','mcat_mortRT_pm25','mcat_mortRT_no2','mcat_mortRT_pm25']
#cat = ['mcat_asth_no2','mcat_asth_pm25','mcat_mort_no2','mcat_mort_pm25']
cat = ['mcat_asthRT_no2','mcat_asthRT_pm25','mcat_mortRT_no2','mcat_mortRT_pm25','mcat_mortRT_o3']


titles2 = [r'Asthma NO$_2$',r'Asthma PM$_{2.5}$',r'Mort  NO$_2$',r'Mortality PM$_{2.5}$']
titles = ['NH White','Black','Hispanic+Latino','Asian']
titles2 = ['Asthma Rate \n'+r'NO$_2$','Asthma Rate\n'+r'PM$_{2.5}$','Mortality Rate\n'+r'NO$_2$','Mortality Rate\n'+r'PM$_{2.5}$','Mortality Rate\n'+r'O$_3$']

cat = ['mcat_mortRT_no2','mcat_mortRT_pm25','mcat_mortRT_o3']
titles = ['NH White','Black','Hispanic+Latino','Asian']
titles2 = ['Mortality Rate\n'+r'NO$_2$','Mortality Rate\n'+r'PM$_{2.5}$','Mortality Rate\n'+r'O$_3$']

#titles=titles3
vmin = [0]*len(cat)
vmax= [1]*len(cat)
alphabet = np.array(list(string.ascii_lowercase)[:len(cat)*(len(c))]).reshape([len(cat),len(c)])
#alphabet = np.array(list(string.ascii_lowercase)+['aa','bb','cc','dd'][:len(cat)*(len(c)+1)]).reshape([len(cat),len(c)+1])
alphabet = alphabet
palette="BuPu"
check_pollution_dist_over_shp_multicat(chi,alphabet,titles2,titles,vmin,vmax,cat,c,palette,vals,'exposure_chi_soccat_healthmort_RT.png',save=True,show=True)



titles2 = [r'Asthma NO$_2$',r'Asthma PM$_{2.5}$',r'Mort  NO$_2$',r'Mortality PM$_{2.5}$']
titles = ['NH White','Black','Hispanic+Latino','Asian']
titles2 = ['Asthma Rate \n'+r'NO$_2$','Asthma Rate\n'+r'PM$_{2.5}$','Mortality Rate\n'+r'NO$_2$','Mortality Rate\n'+r'PM$_{2.5}$','Mortality Rate\n'+r'O$_3$']

cat = ['mcat_pWhite','mcat_pBlack','mcat_pLatino','mcat_pAsian']
c = ['no2','pm25','o3']
titles2 = ['NH White','Black','Hispanic+Latino','Asian']
titles = [r'NO$_2$',r'PM$_{2.5}$',r'O$_3$']

#titles=titles3
vmin = [5,5,15]
vmax= [30,20,45]
alphabet = np.array(list(string.ascii_lowercase)[:len(cat)*(len(c))]).reshape([len(cat),len(c)])
#alphabet = np.array(list(string.ascii_lowercase)+['aa','bb','cc','dd'][:len(cat)*(len(c)+1)]).reshape([len(cat),len(c)+1])
alphabet = alphabet
palette="BuPu"
hilos = check_pollution_dist_over_shp_multicat(chi,alphabet,titles2,titles,vmin,vmax,cat,c,palette,vals,'exposure_chi_soccat_healthmort_RT.png',save=True,show=True)





mdao3warm
c = ['pWhite','pBlack','pLatino','pAsian']
cat = ['mcat_no2','mcat_pm25','mcat_mdao3']
titles = [r'NO$_2$',r'PM$_{2.5}$',r'O$_3$']
titles2 = ['White','Black','Hispanic+Latino','Asian']

cat = ['pWhite','pBlack','pLatino','pAsian']
c = ['mcat_no2','mcat_pm25','mcat_mdao3']
titles2 = [r'NO$_2$',r'PM$_{2.5}$',r'O$_3$']
titles = ['White','Black','Hispanic+Latino','Asian']
titles3 = ['White','Black','Hispanic+Latino','Asian']


#titles3 = ['\n'+r'NO$_2$','Asthma \n'+r'PM$_{2.5}$','Mortality \n'+r'NO$_2$','Mortality \n'+r'PM$_{2.5}$']
titles3 = titles
#titles=titles3
vmin = [-0.05]*len(cat)
vmax= [1.05]*len(cat)
alphabet = np.array(list(string.ascii_lowercase)[:len(c)*(len(cat))]).reshape([len(c),len(cat)])
alphabet = alphabet.T
palette="RdPu"
check_pollution_dist_over_shp_multicat(chi,alphabet,titles,titles2,vmin,vmax,cat,c,palette,vals,'exposure_chi_soccat_pollutants_90-10.png',save=True,show=True)





c = ['pWhite','pBlack','pLatino','pAsian']
cat = ['mcat_Under19ER_no2_rt','mcat_Under19ER_pm25_rt']

titles = [r'Asthma ER NO$_2$',r'Asthma ER PM$_{2.5}$']
titles2 = ['White','Black','Hispanic+Latino','Asian']
titles3 = ['Asthma ER RT \n'+r'NO$_2$','Asthma ER RT \n'+r'PM$_{2.5}$']
#titles=titles3
vmin = [0]*len(c)
vmax= [1]*len(c)
alphabet = np.array(list(string.ascii_lowercase)[:len(cat)*(len(c)+1)]).reshape([len(cat),len(c)+1])
#alphabet = np.array(list(string.ascii_lowercase)+['aa','bb','cc','dd'][:len(cat)*(len(c)+1)]).reshape([len(cat),len(c)+1])
alphabet = alphabet.T
palette="BuPu"
check_pollution_dist_over_shp_multicat(chi,alphabet,titles2,titles,vmin,vmax,cat,c,palette,vals,'exposure_chi_soccat_health_AsthmaERRT_90-10.png',save=True,show=True)



########################################################################################################################################################################
# Make pollution for multiple cats -- just averages

import matplotlib.ticker as mtick

def make_barplot(shp,alphabet,cat,c,palette,vals,clabels,figname,titles,save=False,show=True):
#for zzz in range(1):
   # titles2 == along the top
   # titles = title on the botto of plot
   # titles3 = title on the side
   col,row = len(cat),len(c)
   fig,axs = plt.subplots(1,len(cat),figsize=(len(cat)*2.6,3))
   #for x in range(1):
   import matplotlib
   cmap = matplotlib.cm.get_cmap(palette)
   fl = np.linspace(0.1,0.9,len(c))
   colors = [cmap(fl[i]) for i in range(len(fl))]+['darkgray']
   for k in range(len(cat)):
      ax = axs[k]
      hilos = []
      for j in range(len(c)):
         print('here')
         catk = cat[k]
         qlabel = np.array(vals[vals[0]==catk.split('mcat_')[1]]).ravel()[1:]
         print(qlabel)
         vi = np.array(vals[vals[0]==catk.split('mcat_')[1]]).ravel()[0]
         qlabel = list(qlabel)+np.max(chi[vi])
         h4 = shp[shp[catk] == 0].reset_index(drop=True)
         h3 = shp[shp[catk] == 1].reset_index(drop=True)
         h2 = shp[shp[catk] == 2].reset_index(drop=True)
         h1 = shp[shp[catk] == 3].reset_index(drop=True)
         h0 = shp[shp[catk] == 4].reset_index(drop=True)
         hilo = pd.DataFrame([h4[c[j]],h3[c[j]],h2[c[j]],h1[c[j]],h0[c[j]]]).T
         #print('here1')
         hilo.columns = ['Q1','Q2','Q3','Q4','Q5']
         hilos.append(hilo.mean())      
      #print(hilo)
      hilos = pd.DataFrame(np.array(hilos))
      hilos.columns = hilo.columns
      hilos.columns = ['Q1\n '+str(round(qlabel[0],1)),'Q2\n '+str(round(qlabel[1],1)),'Q3\n '+str(round(qlabel[2],1)),'Q4\n '+str(round(qlabel[3],1)),'Q5\n '+str(round(qlabel[4],1))]#\n >'+str(round(qlabel[3],1))]
      hilos.index = c
      hilos = hilos.T
      hilos['Other'] = 1-hilos.sum(axis=1)
      ci = c+['Other']
      for i in range(len(ci)):
         if i == 0: 
            sns.barplot(hilos.index,hilos[ci[i]],label=clabels[i],color=colors[i],ax=ax,zorder=2.5,edgecolor='k',linewidth=0.1)
            ax.set_ylabel('$%$ of Census Tracts')
         else: 
            ax.bar(hilos.index,hilos[ci[i]],bottom  = np.sum(hilos[ci[0:i]],axis=1),label=clabels[i],color=colors[i],zorder=2.5,edgecolor='k',linewidth=0.1)
            ax.set_ylabel('')
      ax.set_title(titles[k])
      ax.set_ylim([0,1])
      ax.grid(zorder=-1,axis='y')
      ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
      #
   handles, labels = ax.get_legend_handles_labels()
   [sns.despine(offset=7, trim=True,ax=ax) for ax in axs]
   fig.legend(handles, labels, loc = (0.15, 0.02),ncol=5)
   plt.tight_layout()
   fig.subplots_adjust(bottom=0.3)
   if save: plt.savefig(figname,dpi=300)
   if show: plt.show()
   return hilos

########################################################################################################################################################################


# o3
#

c = ['pWhite','pBlack','pLatino','pAsian']
clabels = ['White','Black','Hispanic/Latino','Asian','Other']
cat = ['mcat_no2','mcat_pm25','mcat_mdao3']
titles = [r'NO$_2$',r'PM$_{2.5}$',r'O$_3$']
alphabet = np.array(list(string.ascii_lowercase)[:len(c)*(len(cat))]).reshape([len(c),len(cat)])
alphabet = alphabet.T
palette="RdPu_r"
make_barplot(chi,alphabet,cat,c,palette,vals,clabels,'exposure_stack_barplot.png',titles,save=True,show=True)
#make_barplot(msa,alphabet,cat,c,palette,vals,clabels,'msa_exposure_stack_barplot.png',titles,save=True,show=True)


c = ['pWhite','pBlack','pLatino','pAsian']
clabels = ['White','Black','Hispanic/Latino','Asian','Other']
cat = ['mcat_med_income','mcat_med_income']#','mcat_pm25','mcat_mdao3']
titles = ['Income','Income']#,r'PM$_{2.5}$',r'O$_3$']
alphabet = np.array(list(string.ascii_lowercase)[:len(c)*(len(cat))]).reshape([len(c),len(cat)])
alphabet = alphabet.T
palette="RdPu_r"
make_barplot(chi,alphabet,cat,c,palette,vals,clabels,'income_baseplot.png',titles,save=True,show=True)
#make_barplot(msa,alphabet,cat,c,palette,vals,clabels,'msa_exposure_stack_barplot.png',titles,save=True,show=True)






c = ['pWhite','pBlack','pLatino','pAsian']
cat = ['mcat_mortRT_no2','mcat_mortRT_pm25','mcat_mortRT_o3']
titles = [r'Mortality: NO$_2$',r'Mortality: PM$_{2.5}$',r'Mortality: O$_3$']
alphabet = np.array(list(string.ascii_lowercase)[:len(c)*(len(cat))]).reshape([len(c),len(cat)])
alphabet = alphabet.T
palette="RdPu_r"
make_barplot(chi,alphabet,cat,c,palette,vals,clabels,'mortality_stack_barplot.png',titles,save=True,show=True)
#make_barplot(msa,alphabet,cat,c,palette,vals,clabels,'msa_mortality_stack_barplot.png',titles,save=True,show=True)


c = ['pWhite','pBlack','pLatino','pAsian']
cat = ['mcat_Under19ER_no2_rt','mcat_Under19ER_pm25_rt','mcat_Under19ER_o3_rt']
titles = [r'Pedatric Asthma ER: NO$_2$',r'Pedatric Asthma ER: PM$_{2.5}$',r'Pedatric Asthma ER: MDAO$_3$']
alphabet = np.array(list(string.ascii_lowercase)[:len(c)*(len(cat))]).reshape([len(c),len(cat)])
alphabet = alphabet.T
palette="RdPu_r"
make_barplot(chi,alphabet,cat,c,palette,vals,clabels,'pedEReasthma_stack_barplot.png',titles,save=True,show=True)

cat = ['mcat_asthRT_no2','mcat_asthRT_pm25','mcat_asthRT_o3']
titles = [r'Adult Asthma: NO$_2$',r'Adult Asthma: PM$_{2.5}$',r'Adult Asthma: MDAO$_3$']
palette="RdPu_r"
make_barplot(chi,alphabet,cat,c,palette,vals,clabels,'adultasthma_stack_barplot.png',titles,save=True,show=True)


########################################################################################################################################################################

########################################################################################################################################################################
########################################################################################################################################################################
###################

def make_barplot_BMR_Exposure(shp,alphabet,mi,mi2,cat,c,palette,vals,poll,ylim,clabels,figname,save=False,show=True,extend=False):
   # titles2 == along the top
   # titles = title on the botto of plot
   # titles3 = title on the side
   hilos_out = []
   col,row = len(cat),len(c)
   fig,axs = plt.subplots(1,len(cat),figsize=(len(cat)*2.6,3))
   import matplotlib
   cmap = matplotlib.cm.get_cmap(palette)
   fl = np.linspace(0.1,0.9,len(c))
   colors = [cmap(fl[i]) for i in range(len(fl))]+['darkgray']
   for k in range(len(cat)):
      ax = axs[k]
      hilos = []
      catk = cat[k]
      print(catk)
      qlabel = np.array(vals[vals[0]==catk.split('mcat_')[1]]).ravel()[1:]
      h4 = shp[shp[catk] == 0].reset_index(drop=True)
      h3 = shp[shp[catk] == 1].reset_index(drop=True)
      h2 = shp[shp[catk] == 2].reset_index(drop=True)
      h1 = shp[shp[catk] == 3].reset_index(drop=True)
      h0 = shp[shp[catk] == 4].reset_index(drop=True)
      m = mi2#'popweight-mortall'; print(m)
      li = [h4,h3,h2,h1,h0]
      bmr = [(h4[m].mean()-shp[m].mean())/shp[m].mean()*100 for h4 in li]
      m = poll[k]
      print(m)
      exp = [(h4[m].mean()-shp[m].mean())/shp[m].mean()*100 for h4 in li]
      m = mi+poll[k]; print(m)
      mortrtall = [(h4[m].mean()-shp[m].mean())/shp[m].mean()*100 for h4 in li]
      m = 'asthRT_'+poll[k]; print(m)
      asthrtall = [(h4[m].mean()-shp[m].mean())/shp[m].mean()*100 for h4 in li]
      m = 'Under19ER_'+poll[k]+'_rt'; print(m)
      pedasthrtall = [(h4[m].mean()-shp[m].mean())/shp[m].mean()*100 for h4 in li]
      #
      ped_air = [(h4['asthma_age_adj_rate'].mean()-shp['asthma_age_adj_rate'].mean())/shp['asthma_age_adj_rate'].mean()*100 for h4 in li]
      adult_air = [(h4['AdAsthRt'].mean()-shp['AdAsthRt'].mean())/shp['AdAsthRt'].mean()*100 for h4 in li]
      #
      hilos = pd.DataFrame([bmr,mortrtall,ped_air,adult_air,exp,asthrtall,pedasthrtall])
      hilos.index = ['BMR','MortRate','PedAIR','AdulatAIR','Exposure','AsthmaRate','PediatricAsthmaER']
      hilos.columns = ['Q1','Q2','Q3','Q4','Q5']
      hilos = hilos.T
      hilos_out.append(hilos)
      ci = ['BMR','Exposure']
      for i in range(len(ci)):
         if i == 0: 
            _=sns.barplot(hilos.index,hilos[ci[i]],label=ci[i],color=colors[i+1],ax=ax,zorder=2.5,edgecolor='k',linewidth=0.1)
            _=ax.set_ylabel('$%$ of Census Tracts')
         else: 
            _=ax.bar(hilos.index,hilos[ci[i]],bottom  = 0,label=ci[i],color=colors[i+2],zorder=2.5,edgecolor='k',linewidth=0.1,alpha=0.6)
            _=ax.set_ylabel('')
      _ = ax.plot(hilos.index,hilos['MortRate'],label='Mortality Rate',marker='o',markersize=5,color='k',zorder=2.5,linestyle=None,linewidth =0)
      print((extend == True) & (poll[k] != 'o3'))
      if (extend == True) & (poll[k] != 'o3'):
         _ = ax.plot(hilos.index,hilos['AsthmaRate'],label='Adult Asthma Rate',marker='^',markersize=5,color='k',zorder=2.5,linestyle=None,linewidth =0)
         _ = ax.plot(hilos.index,hilos['PediatricAsthmaER'],label='Pediatric Asthma ER Rate',marker='s',markersize=5,color='k',zorder=2.5,linestyle=None,linewidth =0)
      #
      #
      _ = ax.set_title(titles[k])
      _ = ax.set_ylim(ylim)
      _=ax.grid(zorder=-1,axis='y')
      if k == 0:
         _=ax.set_ylabel('% Difference from\nDomain Average')
      #ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
      #
      handles, labels = axs[0].get_legend_handles_labels()
      _=sns.despine(offset=7, trim=True,ax=ax)
      print(labels)
      _=fig.legend(handles, clabels, loc = (0.1, 0.02),ncol=5,fontsize=8)
      plt.tight_layout()
      fig.subplots_adjust(bottom=0.3)
   if save: plt.savefig(figname,dpi=300)
   if show: plt.show()
   return hilos_out

##########################################################################

c = ['pWhite','pBlack','pLatino','pAsian']

poll = ['no2','pm25','o3']
cat = ['mcat_no2','mcat_pm25','mcat_mdao3']
titles = [r'NO$_2$',r'PM$_{2.5}$',r'O$_3$']
alphabet = np.array(list(string.ascii_lowercase)[:len(c)*(len(cat))]).reshape([len(c),len(cat)])
alphabet = alphabet.T
palette="BuPu"
ylim=[-25,25]
mi='mortRT_'
mi2 ='popweight-mortall'
clabels = ['Mortality Rate', 'BMR', 'Exposure']
hilos_out = make_barplot_BMR_Exposure(chi,alphabet,mi,mi2,cat,c,palette,vals,poll,ylim,clabels,'mcat_exposure_BMR.png',save=True,show=True,extend=True)

hilos_out[0]['Poll'] = 'no2'
hilos_out[1]['Poll'] = 'pm25'
hilos_out[2]['Poll'] = 'o3'

hilos_out[0].append(hilos_out[1]).append(hilos_out[2]).to_csv('quantile_differences_by_mcatpoll.csv')

#make_barplot_BMR_Exposure(msa,alphabet,mi,mi2,cat,c,palette,vals,poll,ylim,clabels,'msa_mcat_exposure_BMR.png',save=True,show=True,extend=False


ylim=[-100,100]
cat = ['mcat_mortRT_no2','mcat_mortRT_pm25','mcat_mortRT_o3']
make_barplot_BMR_Exposure(shp,alphabet,mi,mi2,cat,c,palette,vals,poll,ylim,clabels,'mcat_mortRT_BMR.png',save=True,show=True)
#make_barplot_BMR_Exposure(msa,alphabet,mi,mi2,cat,c,palette,vals,poll,ylim,clabels,'msa_mcat_mortRT_BMR.png',save=True,show=True)


ylim=[-60,60]
poll = ['no2','pm25','o3']
cat = ['mcat_asthRT_no2','mcat_asthRT_pm25','mcat_asthRT_o3']
mi='asthRT_'
mi2 ='AdAsthRt'
clabels = ['Asthma Rate', 'IAR', 'Exposure']

make_barplot_BMR_Exposure(shp,alphabet,mi,mi2,cat,c,palette,vals,poll,ylim,clabels,'mcatASTHMA_'+mi+mi2+'.png',save=True,show=True)


ylim=[-60,60]
poll = ['no2','pm25','o3']
cat = ['mcat_Under19ER_no2_rt','mcat_Under19ER_pm25_rt','mcat_Under19ER_o3_rt']
#cat = ['mcat_asthRT_no2','mcat_asthRT_pm25','mcat_asthRT_o3']
mi='asthRT_'
mi2 ='asthma_age_adj_rate'
clabels = ['Asthma ER Rate', 'ERVisits', 'Exposure']

make_barplot_BMR_Exposure(shp,alphabet,mi,mi2,cat,c,palette,vals,poll,ylim,clabels,'mcatASTHMA_'+mi+mi2+'.png',save=True,show=True)




#########################################################################################################################################################################

def pop_weighted_avg(group,weight,tobeweighted):
     return np.average(group[tobeweighted], weights=np.sum(group[weight],axis=1))


def make_exposure_diff(chi,popgroup,Bpopgroup,Hpopgroup,Apopgroup,Wpopgroup,weights):
   exposures=[]
   for p in poll:
      avg = pop_weighted_avg(chi,popgroup,p)
      b = pop_weighted_avg(chi,Bpopgroup,p)
      w = pop_weighted_avg(chi,Wpopgroup,p)
      h = pop_weighted_avg(chi,Hpopgroup,p)
      a = pop_weighted_avg(chi,Apopgroup,p)
      #
      exposures.append([(b-avg)/avg*100,(w-avg)/avg*100,(h-avg)/avg*100,(a-avg)/avg*100])
   exposures = pd.DataFrame(exposures)
   exposures.columns = ['Black','White','Hispanic','Asian']
   exposures.index = poll
   return exposures

pre=['','W','B','H','A']
popgroup = ['pop_30-35', 'pop_35-44', 'pop_45-54', 'pop_55-64', 'pop_65-74', 'pop_75-84', 'pop_85+']

Bpopgroup = ['B'+popgroup[i] for i in range(len(popgroup))]
Wpopgroup = ['W'+popgroup[i] for i in range(len(popgroup))]
Hpopgroup = ['H'+popgroup[i] for i in range(len(popgroup))]
Apopgroup = ['A'+popgroup[i] for i in range(len(popgroup))]


def make_mortall_weighted(chi,popgroup,mortgroup,pre):
#for i in range(1):
   A = chi[popgroup]; B = chi[mortgroup]
   #A.columns = B.columns
   chi[pre] = np.sum(np.array(A)*np.array(B),axis=1)/np.sum(A,axis=1)
   return chi


# make_mortall_weighted(chi,popgroup,mortgroup,pre):
#make race-weighted mortality rates
chi = make_mortall_weighted(chi,Bpopgroup,mortgroup,'Bpopweight-mortall')
chi = make_mortall_weighted(chi,Hpopgroup,mortgroup,'Hpopweight-mortall')
chi = make_mortall_weighted(chi,Apopgroup,mortgroup,'Apopweight-mortall')
chi = make_mortall_weighted(chi,Wpopgroup,mortgroup,'Wpopweight-mortall')


chi = 
pop_weighted_avg(chi,Bpopgroup,['asthma_age_adj_rate'])
chi = make_mortall_weighted(chi,Hpopgroup,['asthma_age_adj_rate'],'H_asth_popweight')
chi = make_mortall_weighted(chi,Apopgroup,['asthma_age_adj_rate'],'A_asth_popweight')
chi = make_mortall_weighted(chi,Wpopgroup,['asthma_age_adj_rate'],'W_asth_popweight')

np.nanmax(chi[['Bpopweight-mortall','Hpopweight-mortall','Apopweight-mortall','Wpopweight-mortall']],axis=0)


poll = ['no2', 'pm25', 'o3']
def make_exposure_diff(chi,popgroup,Bpopgroup,Hpopgroup,Apopgroup,Wpopgroup,chihi,chilo,poll):
   exposures=[]
   for p in poll:
      avg = pop_weighted_avg(chi,popgroup,p)
      b = pop_weighted_avg(chi,Bpopgroup,p)
      w = pop_weighted_avg(chi,Wpopgroup,p)
      h = pop_weighted_avg(chi,Hpopgroup,p)
      a = pop_weighted_avg(chi,Apopgroup,p)
      hi_inc = pop_weighted_avg(chihi,popgroup,p)
      lo_inc = pop_weighted_avg(chilo,popgroup,p)
      #
      exposures.append([(b-avg)/avg*100,(w-avg)/avg*100,(h-avg)/avg*100,(a-avg)/avg*100,(hi_inc-avg)/avg*100,(lo_inc-avg)/avg*100])
   exposures = pd.DataFrame(exposures)
   exposures.columns = ['Black','White','Hispanic','Asian','High Income','Low Income']
   exposures.index = poll
   return exposures

mortgroup = ['mort_25-34', 'mort_35-44', 'mort_45-54', 'mort_55-64', 'mort_65-74', 'mort_75-84', 'mort_85-99']
#mortgroup = 'AdAsthRt'


def make_mortality_diff(chi,mortgroup,popgroup,Bpopgroup,Wpopgroup,Hpopgroup,Apopgroup,chihi,chilo,here=False):
#for i in range(1):
   if here == True:
      chimortgroup = np.ma.masked_array(chi[mortgroup],np.isnan(chi[mortgroup]))
      avg = np.average(chimortgroup,weights=np.sum(chi[popgroup],axis=1))
      b = np.average(chimortgroup,weights=np.sum(chi[Bpopgroup],axis=1))
      w= np.average(chimortgroup,weights=np.sum(chi[Wpopgroup],axis=1))
      h = np.average(chimortgroup,weights=np.sum(chi[Hpopgroup],axis=1))
      a = np.average(chimortgroup,weights=np.sum(chi[Apopgroup],axis=1))
      chihii = np.ma.masked_array(chihi[mortgroup],np.isnan(chihi[mortgroup]))
      chiloo = np.ma.masked_array(chilo[mortgroup],np.isnan(chilo[mortgroup]))
      hi_inc = np.average(chihii,weights=np.sum(chihi[popgroup],axis=1))
      lo_inc = np.average(chiloo,weights=np.sum(chilo[popgroup],axis=1))
   else:
      chimortgroup = np.ma.masked_array(chi[mortgroup],np.isnan(chi[mortgroup]))
      avg = np.average(chimortgroup,weights=chi[popgroup])
      b = np.average(chimortgroup,weights=chi[Bpopgroup])
      w= np.average(chimortgroup,weights=chi[Wpopgroup])
      h = np.average(chimortgroup,weights=chi[Hpopgroup])
      a = np.average(chimortgroup,weights=chi[Apopgroup])
      hi_inc = np.average(chihi[mortgroup],weights=chihi[popgroup])
      lo_inc = np.average(chilo[mortgroup],weights=chilo[popgroup])
   #avg = (w+a+h)/3
   #exposures = [np.mean((b-avg)/avg)*100,np.mean((w-avg)/avg)*100,np.mean((h-avg)/avg)*100,np.mean((a-avg)/avg)*100]
   exposures = [(b-avg)/avg*100,(w-avg)/avg*100,(h-avg)/avg*100,(a-avg)/avg*100,(hi_inc-avg)/avg*100,(lo_inc-avg)/avg*100]
   #exposures = [b,w,h,a,hi_inc,lo_inc]
   #
   exposures = pd.DataFrame(exposures).T
   exposures.columns = ['Black','White','Hispanic','Asian','High\nIncome','Low\nIncome']
   return exposures

chilo = chi[chi.mcat_med_income < 2].reset_index(drop=True)
chihi = chi[chi.mcat_med_income > 3].reset_index(drop=True)

polls = ['no2','pm25','o3']
exposure_diff = make_exposure_diff(chi,popgroup,Bpopgroup,Hpopgroup,Apopgroup,Wpopgroup,chihi,chilo,polls).T
mortality_diffs = make_mortality_diff(chi,mortgroup,popgroup,Bpopgroup,Wpopgroup,Hpopgroup,Apopgroup,chihi,chilo).T
asthma_diffs = make_mortality_diff(chi,'AdAsthRt',popgroup,Bpopgroup,Wpopgroup,Hpopgroup,Apopgroup,chihi,chilo,True).T
asthmaER_diffs = make_mortality_diff(chi,'asthma_age_adj_rate','Under19','under18bla','under18whi','under18his','under18asi',chihi,chilo).T
health_diff = pd.DataFrame([mortality_diffs[0],asthma_diffs[0],asthmaER_diffs[0]]).T
health_diff.columns = 'Mortality','Adult Asthma','Pediatric Asthma ER'

mortality_o3 = make_mortality_diff(chi,'mortRT_o3',popgroup,Bpopgroup,Wpopgroup,Hpopgroup,Apopgroup,chihi,chilo,True).T
mortality_pm25 = make_mortality_diff(chi,'mortRT_pm25',popgroup,Bpopgroup,Wpopgroup,Hpopgroup,Apopgroup,chihi,chilo,True).T
mortality_no2 = make_mortality_diff(chi,'mortRT_no2',popgroup,Bpopgroup,Wpopgroup,Hpopgroup,Apopgroup,chihi,chilo,True).T

mort_diffs = pd.DataFrame([mortality_no2[0],mortality_pm25[0],mortality_o3[0]]).T
mort_diffs.columns = 'NO2-Attributable Mortality','PM25-Attributable Mortality','O3-Attributable Mortality'

asthma_no2 = make_mortality_diff(chi,'asthRT_no2',popgroup,Bpopgroup,Wpopgroup,Hpopgroup,Apopgroup,chihi,chilo,True).T
asthma_pm25 = make_mortality_diff(chi,'asthRT_pm25',popgroup,Bpopgroup,Wpopgroup,Hpopgroup,Apopgroup,chihi,chilo,True).T
asthma_o3 = make_mortality_diff(chi,'asthRT_o3',popgroup,Bpopgroup,Wpopgroup,Hpopgroup,Apopgroup,chihi,chilo,True).T
asthmaER_no2 = make_mortality_diff(chi,'Under19ER_no2_rt',popgroup,Bpopgroup,Wpopgroup,Hpopgroup,Apopgroup,chihi,chilo,True).T
asthmaER_pm25 = make_mortality_diff(chi,'Under19ER_pm25_rt',popgroup,Bpopgroup,Wpopgroup,Hpopgroup,Apopgroup,chihi,chilo,True).T
asthmaER_o3 = make_mortality_diff(chi,'Under19ER_o3_rt',popgroup,Bpopgroup,Wpopgroup,Hpopgroup,Apopgroup,chihi,chilo,True).T

asthma_diffs = pd.DataFrame([asthma_no2[0],asthma_pm25[0],asthma_o3[0],asthmaER_no2[0],asthmaER_pm25[0],asthmaER_o3[0]]+[mortality_no2[0],mortality_pm25[0],mortality_o3[0]]).T
asthma_diffs.columns = 'NO2-Attributable Asthma','PM25-Attributable Asthma','O3-Attributable Asthma','NO2-Attributable ER','PM25-Attributable ER','O3-Attributable ER','NO2-Attributable Mortality','PM25-Attributable Mortality','O3-Attributable Mortality'


# allpopgroup = ['pop_25-34', 'pop_35-44', 'pop_45-54', 'pop_55-64', 'pop_65-74', 'pop_75-84', 'pop_85+']

# for i in range(len(allpopgroup)):
#    chi['SYN'+allpopgroup[i]] = chi[[Bpopgroup[i],Wpopgroup[i],Hpopgroup[i],Apopgroup[i]]].sum(axis=1)

# synpopgroup = ['SYN'+allpopgroup[i] for i in range(len(allpopgroup))]

# def makehealthoutcome_weighted(chi,popgroup,healthgroup,pre):
#    A = np.array(chi[popgroup]); B = np.array(chi[healthgroup]); C = np.array(chi[synpopgroup])
#    if A.shape == B.shape:
#       chi[pre+healthgroup[0]] = np.sum(A*B,axis=1)/np.sum(C,axis=1)
#    else:
#       chi[pre+healthgroup] = np.sum(A.T*B,axis=0)/np.sum(C.T,axis=0)
#    return chi

# chi = makehealthoutcome_weighted(chi,Bpopgroup,'asthma_age_adj_rate','B')
# chi = makehealthoutcome_weighted(chi,Apopgroup,'asthma_age_adj_rate','A')
# chi = makehealthoutcome_weighted(chi,Hpopgroup,'asthma_age_adj_rate','H')
# chi = makehealthoutcome_weighted(chi,Wpopgroup,'asthma_age_adj_rate','W')

# chi = makehealthoutcome_weighted(chi,Wpopgroup,mortgroup,'W')
# chi = makehealthoutcome_weighted(chi,Apopgroup,mortgroup,'A')
# chi = makehealthoutcome_weighted(chi,Hpopgroup,mortgroup,'H')
# chi = makehealthoutcome_weighted(chi,Bpopgroup,mortgroup,'B')


# chi = makehealthoutcome_weighted(chi,Wpopgroup,'med_income','W')
# chi = makehealthoutcome_weighted(chi,Apopgroup,'med_income','A')
# chi = makehealthoutcome_weighted(chi,Hpopgroup,'med_income','H')
# chi = makehealthoutcome_weighted(chi,Bpopgroup,'med_income','B')

# healthgroup = 'med_income'
# np.array(chi[['B'+healthgroup,'W'+healthgroup,'H'+healthgroup,'A'+healthgroup]])


# [np.mean(chi['Basthma_age_adj_rate']/chi[Bpopgroup].sum(axis=1)),np.nanmean(chi['Wasthma_age_adj_rate']/chi[Wpopgroup].sum(axis=1)),
# np.nanmean(chi['Hasthma_age_adj_rate']/chi[Hpopgroup].sum(axis=1)),np.nanmean(chi['Aasthma_age_adj_rate']/chi[Apopgroup].sum(axis=1))]

# #
#exposure_diff_top20 = make_exposure_diff(chi[chi.mcat_no2 == '4'],popgroup,Bpopgroup,Hpopgroup,Apopgroup,Wpopgroup)
#mortality_diffs_top20 = make_exposure_diff(chi[chi.mcat_no2 == '4'],['popweight-mortall'], ['Bpopweight-mortall'], ['Hpopweight-mortall'],
##       ['Apopweight-mortall'], ['Wpopweight-mortall'])

#exposure_diff_bot20 = make_exposure_diff(chi[chi.mcat_no2 == '1'],popgroup,Bpopgroup,Hpopgroup,Apopgroup,Wpopgroup)
#mortality_diffs_bot20 = make_exposure_diff(chi[chi.mcat_no2 == '1'],['popweight-mortall'], ['Bpopweight-mortall'], ['Hpopweight-mortall'],
#       ['Apopweight-mortall'], ['Wpopweight-mortall'])



palette = 'RdPu'
import matplotlib
cmap = matplotlib.cm.get_cmap(palette)
fl = np.linspace(0.1,0.9,10)
colors = [cmap(fl[i]) for i in range(len(fl))]+['darkgray']

fig,axs = plt.subplots(3,sharey=True,sharex=True,figsize=(4,5))
ax=axs[0]
ax.barh(mortality_diffs.index,mortality_diffs[0],color = 'white',alpha=1,edgecolor='k',hatch='//',zorder=2.5,height=0.6,label='Baseline Mortality Rate')
ax.barh(exposure_diff.index,exposure_diff['no2'],color = 'k',edgecolor='k',zorder=2.5,label='Exposure')
#ax.scatter(mortrt_diffs['no2'],mortrt_diffs.index,color = 'k',edgecolor='k',zorder=3,label='Mortality Rate')
handles, labels = ax.get_legend_handles_labels()

ax.barh(mortality_diffs.index,mortality_diffs[0],color = 'white',alpha=1,edgecolor=colors[9],hatch='//',zorder=2.5,height=0.6)
ax.barh(exposure_diff.index,exposure_diff['no2'],color = colors[9],edgecolor=colors[9],zorder=2.5)
#ax.scatter(mortrt_diffs['no2'],mortrt_diffs.index,color = 'k',edgecolor='k',zorder=3,label='Mortality Rate')
ax.set_title(r'NO$_2$')

ax=axs[1]
ax.barh(mortality_diffs.index,mortality_diffs[0],color = 'white',alpha=1,edgecolor=colors[6],hatch='//',zorder=2.5,height=0.6)
ax.barh(exposure_diff.index,exposure_diff['pm25'],color = colors[6],edgecolor=colors[6],zorder=2.5)
#ax.scatter(mortrt_diffs['no2'],mortrt_diffs.index,color = 'k',edgecolor='k',zorder=3,label='Mortality Rate',s=5)
ax.set_title(r'PM$_{2.5}$')

ax=axs[2]
ax.barh(mortality_diffs.index,mortality_diffs[0],color = 'white',alpha=1,edgecolor=colors[3],hatch='//',label='Mortality Rate',zorder=2.5,height=0.6)
ax.barh(exposure_diff.index,exposure_diff['o3'],color = colors[3],edgecolor=colors[3],label='Exposure',zorder=2.5)
#ax.scatter(mortrt_diffs['no2'],mortrt_diffs.index,color = 'k',edgecolor='k',zorder=3,label='Mortality Rate')

ax.set_xlabel('Population-Weighted Differences\nRelative to Chicago (%)')
ax.set_title(r'O$_3$')

_ = [axs[i].grid(zorder=-1,axis='x') for i in range(len(axs))]
plt.tight_layout()


fig.subplots_adjust(bottom=0.2)
_=fig.legend(handles, labels, loc = (0.01, 0.02),ncol=5,fontsize=8)

#plt.savefig('AMR_population_weighted_differences_relativetoChicago.png',dpi=300)

plt.show()


##################################################################################################################

palette = 'RdPu'
import matplotlib
cmap = matplotlib.cm.get_cmap(palette)
fl = np.linspace(0.1,0.9,10)
colors = [cmap(fl[i]) for i in range(len(fl))]+['darkgray']

cmap = matplotlib.cm.get_cmap('BuPu_r')
fl = np.linspace(0.1,0.9,10)
colors2 = [cmap(fl[i]) for i in range(len(fl))]+['darkgray']


fig,ax = plt.subplots(3,1,sharey=True,figsize=(5,4))
ax=ax.T.ravel()
datas = [[exposure_diff['no2'],exposure_diff['pm25'],exposure_diff['o3']],[health_diff['Mortality'],health_diff['Adult Asthma'],health_diff['Pediatric Asthma ER']]]

cs = [[colors[9],colors[6],colors[3]],[colors2[9],colors2[6],colors2[3]]]

titles = [r'NO$_2$',r'PM$_{2.5}$',r'O$_3$','Mortality','Asthma Incidence','Pedatric Asthma ER']

for i in range(len(datas[0])):
   ax[i].barh(exposure_diff.index,datas[1][i],color = 'white',alpha=1,edgecolor=cs[0][i],hatch='//',zorder=2.5,height=0.6)
   ax[i].barh(exposure_diff.index,datas[0][i],alpha=1,color=cs[0][i],zorder=2.5)
   ax[i].set_title(titles[i])

_ = [ax[i].grid(zorder=-1,axis='x') for i in range(len(ax))]
plt.tight_layout()


#fig.subplots_adjust(bottom=0.2)
#_=fig.legend(handles, labels, loc = (0.01, 0.02),ncol=5,fontsize=8)

#plt.savefig('AMR_population_weighted_differences_relativetoChicago.png',dpi=300)

plt.show()




fig,ax = plt.subplots(1,2,sharey=True,figsize=(5,4))
ax=ax.T.ravel()
datas = [exposure_diff['no2'],exposure_diff['pm25'],exposure_diff['o3'],health_diff['Mortality'],health_diff['Adult Asthma'],health_diff['Pediatric Asthma ER']]
cs = [colors[9],colors[6],colors[3],colors2[9],colors2[6],colors2[3]]
xlim = [[-5,5]]*3+[[-50,50]]*3

titles = [r'NO$_2$',r'PM$_{2.5}$',r'O$_3$','Mortality','Asthma Incidence','Pedatric Asthma ER']

for j in range(len(datas)):
   for i in range(3):
      ax[j].barh(exposure_diff.index,datas[i],alpha=1,color=cs[i],edgecolor='k',linewidth=0.5,zorder=2.5)
   else: 
      ax[i].barh(exposure_diff.index,datas[i],color=cs[i],edgecolor='k',linewidth=0.5,zorder=2.5)
   #ax
   #ax[i].barh(exposure_diff.index,datas[0][i],alpha=1,color=cs[0][i],zorder=2.5)
   ax[i].set_title(titles[i])
   ax[i].set_xlim(xlim[i])

_ = [ax[i].grid(zorder=-1,axis='x') for i in range(len(ax))]
plt.tight_layout()

#fig.subplots_adjust(bottom=0.2)
#_=fig.legend(handles, labels, loc = (0.01, 0.02),ncol=5,fontsize=8)

#plt.savefig('AMR_population_weighted_differences_relativetoChicago.png',dpi=300)

plt.show()

asthma_diffs = asthma_diffs.reset_index()

fig,ax = plt.subplots(1,1,sharey=True,figsize=(7,7))
asthma_diffs.plot(x='index',y=['NO2-Attributable Asthma', 'PM25-Attributable Asthma',
       'O3-Attributable Asthma', 'NO2-Attributable ER', 'PM25-Attributable ER',
       'O3-Attributable ER', 'NO2-Attributable Mortality',
       'PM25-Attributable Mortality', 'O3-Attributable Mortality'],
        kind='barh',cmap='PuBuGn',ax=ax,edgecolor='k',linewidth=0.5,zorder=2.5,legend=False)

ax.set_xlim([-50,50])
_ = ax.grid(zorder=-1,axis='x')

handles, labels = ax.get_legend_handles_labels()
fig.subplots_adjust(bottom=0.18)
fig.legend(handles,labels, loc = (0.05, 0.03),ncol=3,fontsize=9)

plt.savefig('allhealthoutcomes_relativetochicago.png',dpi=300)
plt.show()
chi.to_


import matplotlib.ticker as mtick

exposure_diff = exposure_diff.reset_index() 
health_diff = health_diff.reset_index()
asthma_diffs = asthma_diffs.reset_index()

fig,ax = plt.subplots(1,3,sharey=True,figsize=(8,3.5))
exposure_diff.plot(x='index',y=['no2', 'pm25', 'o3'],kind='barh',cmap='RdPu',ax=ax[0],edgecolor='k',linewidth=0.5,zorder=2.5,legend=False)
health_diff.plot(x='index',y=['Mortality','Adult Asthma','Pediatric Asthma ER'],kind='barh',cmap='BuPu',ax=ax[1],edgecolor='k',linewidth=0.5,zorder=2.5,legend=False)
asthma_diffs.plot(x='index',y=['NO2-Attributable Mortality', 'PM25-Attributable Mortality',
       'O3-Attributable Mortality'],
                  kind='barh',cmap='PuBuGn',ax=ax[2],edgecolor='k',linewidth=0.5,zorder=2.5,legend=False)

ax[0].set_xlim([-5,5]);ax[1].set_xlim([-50,50])

_ = [ax[i].grid(zorder=-1,axis='x') for i in range(len(ax))]

ax[0].set_title('(a) Pollution \nExposure')
ax[1].set_title('(b) Baseline \nHealth Incidence')
ax[2].set_title('(c) Pollution-Attributable \nHealth Outcomes')

for i in range(len(ax)):
   ax[i].xaxis.set_major_formatter(mtick.PercentFormatter())
   ax[i].tick_params(axis='x', labelsize=8)


handles, labels = ax[0].get_legend_handles_labels()
labels=[r'NO$_2$',r'PM$_{2.5}$',r'O$_3$']

handles1, labels1 = ax[1].get_legend_handles_labels()
handles2, labels2 = ax[2].get_legend_handles_labels()
labels2=[r'NO$_2$-Attr. Mortality Rt',r'PM$_{2.5}$-Attr. Mortality Rt',r'O$_3$-Attr. Mortality Rt']


la = [labels,labels1,labels2]
ha = [handles,handles1,handles2]
dists = [0.19,0.41,0.68]
fig.subplots_adjust(bottom=0.3)

for i in range(len(ha)):
   _=fig.legend(ha[i],la[i], loc = (dists[i], 0.03),ncol=1,fontsize=9)

#_=fig.legend(handles1, labels1, loc = (0.5, 0.02),ncol=3,fontsize=8)
plt.savefig('baselinehealthstuff.png',dpi=300)
plt.show()






