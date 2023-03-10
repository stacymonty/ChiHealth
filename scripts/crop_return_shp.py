#############################
# crop to relevant areas
#############################

def crop_return_shp(fin,target_crs,gdf_in):
   fin = fin.to_crs(target_crs)
   fin_union = gpd.GeoSeries(unary_union(fin.geometry))
   fin = gpd.clip(gdf_in,fin_union).reset_index(drop=True)
   fin[‘dt’] = [(fin.geometry[i].type != ‘Point’)  & (fin.geometry[i].type != ‘MultiLineString’)  & (fin.geometry[i].type != ‘LineString’) for i in range(len(fin))]
   fin=fin[fin.dt].reset_index(drop=True)
   fin = fin.reset_index(drop=True)
   fin_union = gpd.GeoDataFrame({‘geometry’:fin_union},crs=target_crs)
   #fin = fin[fin.geometry.type == ‘Polygon’].reset_index(drop=True)
   return fin,fin_union
