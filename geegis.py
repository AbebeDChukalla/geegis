# -*- coding: utf-8 -*-
"""
Created on Thu Nov 01 11:40:22 2018

@author: bec
"""
import ee
import urllib
import zipfile
import numpy as np
ee.Initialize()

def convert_to_yearly(image_collection, start_end, reducer = ee.Reducer.sum()):

    
    band_name = ee.Image(image_collection.first()).bandNames().getInfo()[0]
    
    def _fixTIMERES(year):

        year = ee.Number(year)

        VARnew = image_collection.filterDate(ee.Date.fromYMD(year, 1, 1), ee.Date.fromYMD(year.add(1), 1, 1))
        
        VARyear = VARnew.reduce(reducer)
      
        return VARyear.copyProperties(VARnew.first(), VARnew.first().propertyNames())

    years = ee.List.sequence(start_end[0], start_end[1], 1)
            
    VARyearly = ee.ImageCollection(years.map(_fixTIMERES))
    
    VARyearly = VARyearly.select([0], [band_name])
    
    return VARyearly
    
def downloadImage(image, output_fh, shape, scale):
    
    region = ee.Geometry(shape.geometry().bounds(1).getInfo()).toGeoJSONString()
    
    params = {
              'name':'test',
              'crs': 'EPSG:4326',
              'scale': scale,
              'region': region
             }
    
    url = image.getDownloadURL(params)

    succes = False
    while not succes:
        try:
            print("start download")
            urllib.urlretrieve(url, output_fh)
            zip_ref = zipfile.ZipFile(output_fh, 'r')
            zip_ref.extractall(output_fh[:-4])
            zip_ref.close()
            succes = True
        except:
            pass
        
def createTS(image_collection, geometry, scale, reducer_names = ['mean'], copy_props = ['system:time_start']):
    
    reducers_dict = {'mean':    ee.Reducer.mean(), 
                     'p95':     ee.Reducer.percentile([95]),
                     'product': ee.Reducer.product(),
                     'sum':     ee.Reducer.sum()}
    
    reducer = reducers_dict[reducer_names[0]]
    for name in reducer_names[1:]:
        reducer = reducer.combine(reducer2 = reducers_dict[name], sharedInputs = True)
    
    NDV = -9999

    band_name = ee.Image(image_collection.first()).bandNames().getInfo()[0]
    
    def _createTS(VARimg):

        # Create dictionary with stats
        dict1 = VARimg.reduceRegion(reducer = reducer, geometry = geometry, maxPixels = 1e15, scale = scale)

        # Adjust key-name incase of using 1 reducer
        if len(reducer_names) == 1:
            dict1 = dict1.rename([band_name], [band_name + '_' + reducer_names[0]])

        # Set no-data-value
        dict1 = dict1.map(lambda key, val: ee.List([val, NDV]).reduce(ee.Reducer.firstNonNull()))

        # Create feature from dictionary
        VARft = ee.Feature(None, dict1)
    
        # Copy the requested properties
        if copy_props == 'all':
            VARft = VARft.copyProperties(VARimg, VARimg.propertyNames())
        elif isinstance(copy_props, list):
            VARft = VARft.copyProperties(VARimg, copy_props)

        return VARft.copyProperties(VARimg, copy_props)

    VARts = image_collection.map(_createTS).set('NDV', NDV)       
    
    properties = [s for s in VARts.first().propertyNames().getInfo()]
    
    results = dict()
    
    for name in properties:
        
        vals = np.array(VARts.aggregate_array(name).getInfo())
        vals = np.where(vals == NDV, NDV, vals)

        results[name] = vals
    
    return results