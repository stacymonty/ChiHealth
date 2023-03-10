#############################
# crop to relevant areas
#############################

def crop_return_shp(fin,target_crs,gdf_in):
   # have to: import geopandas as gpd
   # Warning -- will not work if you're cropping points?
   # fin = shapefile to crop to (a geopandas dataframe)
   # gdf_in = geopandas dataframe with your data 
   # target_crs = the crs you want to crop to (recommended to be the fin.crs or gdf_in.crs)
   # gdf_out = data cropped to shapefile
   # gdf_out_union = outside bounds of the gdf you just cropped
   fin = fin.to_crs(target_crs)
   gdf_in = gdf_in.to_crs(target_crs)
   fin_union = gpd.GeoSeries(unary_union(fin.geometry))
   gdf_out = gpd.clip(gdf_in,fin_union).reset_index(drop=True)
   # removing anything that was clipped weird
   gdf_out[‘dt’] = [(gdf_out.geometry[i].type != ‘Point’)  & (gdf_out.geometry[i].type != ‘MultiLineString’)  
                & (gdf_out.geometry[i].type != ‘LineString’) for i in range(len(fin))]
   gdf_out=gdf_out[fin.dt].reset_index(drop=True)
   gdf_out = gdf_out.reset_index(drop=True)
   gdf_out_union = gpd.GeoDataFrame({‘geometry’:gdf_out},crs=target_crs)
   #fin = fin[fin.geometry.type == ‘Polygon’].reset_index(drop=True)
   return gdf_out,gdf_out_union
