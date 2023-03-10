def mask_given_shapefile(lon,lat,shapefile):
   '''
   Make sure to import these packages:
   import geopandas
   from shapely.ops import unary_union 
   
   Make a mask given a shapefile:
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
