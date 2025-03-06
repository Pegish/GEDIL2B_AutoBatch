# GEDIL2B_AutoBatch

# Import 
  import os
  import h5py
  import numpy as np
  import pandas as pd
  import geopandas as gp
  from shapely.geometry import Point
  from shapely.geometry import box
  from shapely.geometry import MultiPolygon
  import geoviews as gv
  from geoviews import opts, tile_sources as gvts
  import holoviews as hv
  gv.extension('bokeh', 'matplotlib')
  import shapely
  import warnings
  from shapely.errors import ShapelyDeprecationWarning
  warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning) 

# Verify the change and directory
  current_directory = os.getcwd()
  print("Current Directory:", current_directory)
  new_directory = "D:/Python/NewFilter-GEDIL2BGS2022"
  os.chdir(new_directory)
  
  inDir = os.getcwd()
  print("Updated Directory:", inDir)
  output_folder = os.path.join(new_directory, "outputpolygonWithinFOREST")
  os.makedirs(output_folder, exist_ok=True)
  print("Output Folder:", output_folder)

# Boundary of Study area (ROI)
  shapefile_path = 'D:/Python/NewFilter-GEDIL2BGS2022/forestproject.shp'

# Step 1: Read the shapefile using geopandas
  forest = gp.read_file(shapefile_path)


# Step 2: Initialize lists to store lengths
  lengths_variable_selection = []
  lengths_after_filtering = []
  lengths_after_clipping = []
# Step 3: List GEDI L2A .h5 files in the inDir
  gediFiles = [g for g in os.listdir(inDir) if g.startswith('processed_GEDI02_B_2022270205255_O21479_03_T00950_02_003_01_V002') and g.endswith('.h5')]  

# Step 4: Explore metadata and beam names
  for L2B_file in gediFiles:
      L2B_path = os.path.join(inDir, L2B_file)
      #print(f"Processing file: {L2A_path}")
  
      with h5py.File(L2B_path, 'r') as gediL2B:
  
          print("Keys in the file:", list(gediL2B.keys()))
          print("METADATA:", list(gediL2B['METADATA']))
          for g in gediL2B['METADATA']['DatasetIdentification'].attrs:
              print(g)
          #print("Purpose:", gediL2A['METADATA']['DatasetIdentification'].attrs['purpose'])
  
          beamNames = [g for g in gediL2B.keys() if g.startswith('BEAM')]
          #print("Beam Names:", beamNames)

  # Step 5: Initialize lists to store data
          gediL2B_objs = []
          gediL2B.visit(gediL2B_objs.append)
          gediSDS = [o for o in gediL2B_objs if isinstance(gediL2B[o], h5py.Dataset)]
          
          shotNum, cover, coverZ, fhdNormal, pai, paiz, pavdz, dem, zElevation, zHigh, zLat, zLon, rh100, quality, degrade, sensitivity, SolarElev, Solarazi, surface, beamI, selectedAlgorithm = ([] for _ in     
          range(21))
          coverZ = [[] for _ in range(30)]
          paiZ = [[] for _ in range(30)]
          pavdZ = [[] for _ in range(30)]
          
  # Step 6: Read data from each beam and append to lists
          for beamName in beamNames:
             [shotNum.append(str(h)) for h in gediL2B[[g for g in gediSDS if g.endswith('/shot_number') and beamName in g][0]][()]]
             [dem.append(h) for h in gediL2B[[g for g in gediSDS if g.endswith('/digital_elevation_model') and beamName in g][0]][()]]
             [cover.append(h) for h in gediL2B[[g for g in gediSDS if g.endswith('/cover') and beamName in g][0]][()]]
             [fhdNormal.append(h) for h in gediL2B[[g for g in gediSDS if g.endswith('/fhd_normal') and beamName in g][0]][()]]
             [pai.append(h) for h in gediL2B[[g for g in gediSDS if g.endswith('/pai') and beamName in g][0]][()]]
             cover_data = gediL2B[f'{beamName}/cover_z'][()]  # Retrieve all rh data for this beam
             for i in range(30):  # Loop over indices
                 [coverZ[i].append(h[i]) for h in cover_data]
             pai_data = gediL2B[f'{beamName}/pai_z'][()]  # Retrieve all rh data for this beam
             for i in range(30):  # Loop over indices
                 [paiZ[i].append(h[i]) for h in pai_data]  
             pavd_data = gediL2B[f'{beamName}/pavd_z'][()]  # Retrieve all rh data for this beam
             for i in range(30):  # Loop over indices
                 [pavdZ[i].append(h[i]) for h in pavd_data]               
             [zElevation.append(h) for h in gediL2B[[g for g in gediSDS if g.endswith('/elev_lowestmode') and beamName in g][0]][()]]  
             [zHigh.append(h) for h in gediL2B[[g for g in gediSDS if g.endswith('/elev_highestreturn') and beamName in g][0]][()]]  
             [zLat.append(h) for h in gediL2B[[g for g in gediSDS if g.endswith('/lat_lowestmode') and beamName in g][0]][()]]  
             [zLon.append(h) for h in gediL2B[[g for g in gediSDS if g.endswith('/lon_lowestmode') and beamName in g][0]][()]]  
             [rh100.append(h) for h in gediL2B[[g for g in gediSDS if g.endswith('/rh100') and beamName in g][0]][()]]  
             [quality.append(h) for h in gediL2B[[g for g in gediSDS if g.endswith('/l2b_quality_flag') and beamName in g][0]][()]]  
             [degrade.append(h) for h in gediL2B[[g for g in gediSDS if g.endswith('/degrade_flag') and beamName in g][0]][()]]  
             [sensitivity.append(h) for h in gediL2B[[g for g in gediSDS if g.endswith('/sensitivity') and beamName in g][0]][()]]  
             [SolarElev.append(h) for h in gediL2B[[g for g in gediSDS if g.endswith('/solar_elevation') and beamName in g][0]][()]]
             [Solarazi.append(h) for h in gediL2B[[g for g in gediSDS if g.endswith('/solar_azimuth') and beamName in g][0]][()]]
             [surface.append(h) for h in gediL2B[[g for g in gediSDS if g.endswith('/surface_flag') and beamName in g][0]][()]]
             [beamI.append(h) for h in [beamName] * len(gediL2B[[g for g in gediSDS if g.endswith('/shot_number') and beamName in g][0]][()])]
             [selectedAlgorithm.append(h) for h in gediL2B[[g for g in gediSDS if g.endswith('/selected_l2a_algorithm') and beamName in g][0]][()]]
          
   # Step 7: Store the data in a DataFrame
          dataframe = pd.DataFrame({
              'Shot Number': shotNum,
              'Beam': beamI,
              'Latitude': zLat,
              'Longitude': zLon,
              'Digital Elevation Model': dem,
              'Cover': cover,
              'fhd_normal': fhdNormal,
              'PAI': pai,
              'RH100': rh100,
              'Ground Elevation (m)': zElevation,
              'Canopy Elevation (m)': zHigh,
              **{f'PAI_z{i}': paiZ[i] for i in range(30)}, 
              **{f'PAVD_z{i}': pavdZ[i] for i in range(30)}, 
              **{f'cover_z{i}': coverZ[i] for i in range(30)}, 
              'L2B Quality Flag': quality,
              'Degrade Flag': degrade,
              'Sensitivity': sensitivity,
              'Solar Elevation': SolarElev,
              'Solar Azimuth': Solarazi,
              'Surface Flag': surface,
              'Selected L2A Algorithm': selectedAlgorithm
          })
        lengths_variable_selection.append(len(dataframe))
        # Read RH data from HDF5 file and append to the list
        del shotNum, cover, coverZ, fhdNormal, pai, paiz, pavdz, dem, zElevation, zHigh, zLat, zLon, rh100, quality, degrade, sensitivity, SolarElev, Solarazi, surface, beamI, selectedAlgorithm
        print(len(dataframe))
        
  # Step 8: Import GeoJSON as GeoDataFrame
        ROI = gp.GeoDataFrame.from_file('ROI.geojson')  
        ROI ['geometry'][0]  # Plot GeoDataFrame      
        ROI.envelope[0].bounds
        minLon, minLat, maxLon, maxLat = ROI.envelope[0].bounds  # Define the min/max lat/lon from the bounds of Redwood NP    
        dataframe = dataframe.where(dataframe['Latitude'] > minLat)
        dataframe = dataframe.where(dataframe['Latitude'] < maxLat)
        dataframe = dataframe.where(dataframe['Longitude'] > minLon)
        dataframe = dataframe.where(dataframe['Longitude'] < maxLon)
        dataframe = dataframe.dropna()
        print(len(dataframe))
        print("dataframe after variable selection:", len(dataframe))
        
  # Step 9: Apply filters to remove specific variables
        dataframe = dataframe.where(dataframe['L2B Quality Flag'].ne(0))
        dataframe = dataframe.where(dataframe['Degrade Flag'] < 1) 
        dataframe = dataframe.where(dataframe['Sensitivity'] > 0.9)
        dataframe = dataframe.where((dataframe['Digital Elevation Model']-dataframe['Ground Elevation (m)']) < abs(2))
        #dataframe = dataframe.where(dataframe['Solar Elevation'] < 0)
        dataframe = dataframe.where(dataframe['Surface Flag'] == 1)
        dataframe = dataframe.where(dataframe['RH100'] < 3600) 
        dataframe = dataframe.dropna()
        lengths_after_filtering.append(len(dataframe))
        print(len(dataframe))
        #print(dataframe.dtypes)
        print("dataframe after filtering:", len(dataframe))
        
  # Step 10: Convert the DataFrame to GeoDataFrame
        geometry = [Point(xy) for xy in zip(dataframe['Longitude'], dataframe['Latitude'])]
        crs = {'init': 'epsg:4326'}
        gdf = gp.GeoDataFrame(dataframe, crs=crs, geometry=geometry)

        # Create circles around each point with a radius of 12.5 meters
        radius = 0.0001126126  # Radius in meters
        gdf['geometry'] = gdf['geometry'].buffer(radius)
        
        gdf_forest = gp.sjoin(gdf, forest, how="inner", op="within")
        lengths_after_clipping.append(len(gdf_forest))
        
  # Step 11: Save the GeoDataFrame as a shapefile
        out_shapefile = os.path.join(output_folder, L2B_file.replace('.h5', '_polygonsWithinFOREST.shp'))
        gdf_forest.to_file(out_shapefile)
        # Save the GeoDataFrame with polygons within the forest as GeoJSON
        out_geojson = os.path.join(output_folder, L2B_file.replace('.h5', '_Withinforest.json'))
        gdf_forest.to_file(out_geojson, driver='GeoJSON')

        # Clean up
        del dataframe, gdf


# Step 12: After processing all files, create a dataframe from the lists
  lengths_df = pd.DataFrame({
      'Length after variable selection': lengths_variable_selection,
      'Length after filtering': lengths_after_filtering,
      'Length after clipping' : lengths_after_clipping
  })

# Step 13: Save the dataframe to a CSV file
  lengths_df.to_csv('dataframe_lengthsWithinFOREST.csv', index=False)
