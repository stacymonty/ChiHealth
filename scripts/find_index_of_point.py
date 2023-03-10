# find index

def find_index(stn_lon, stn_lat, wrf_lon, wrf_lat):
   # stn_lon/lat -- 1-d array or list of lat and lon 
   # wrf_lon/lat -- 2-d array of lat and lon
   # returns the indices of where the stn_lon/lat lies within the array of lon,lat data
   xx=[];yy=[]
   for i in range(len(stn_lat)):
      abslat = np.abs(wrf_lat-stn_lat[i])
      abslon= np.abs(wrf_lon-stn_lon[i])
      c = np.maximum(abslon,abslat)
      latlon_idx = np.argmin(c)
      x, y = np.where(c == np.min(c))
      #add indices of nearest wrf point station
      xx.append(x)
      yy.append(y)
   #
   xx=[xx[i][0] for i in range(len(xx))];yy=[yy[i][0] for i in range(len(yy))]
   #return indices list
   return xx, yy
