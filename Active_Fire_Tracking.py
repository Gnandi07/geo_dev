import ee
import streamlit as st
import geemap.foliumap as geemap
import pandas as pd
import numpy as np
import altair as alt
from datetime import datetime, date, timedelta
import requests
import json
import time

scale = 50

def setTime(image):
  return ee.ImageCollection(ee.List(image.get('images'))).mosaic().set('system:time_start', ee.Date(image.get('date')).millis())

def formatTime(image):
  return image.set('date', image.date().format('yyyy-MM-dd'))

def getDailyMosaic(img_coll):
  img_coll = img_coll.map(algorithm = formatTime)

  daily = ee.ImageCollection(ee.Join.saveAll('images').apply(**{
      'primary': img_coll,
      'secondary': img_coll,
      'condition': ee.Filter.And(
          ee.Filter.equals(**{'leftField': 'date','rightField': 'date'}),
          ee.Filter.equals(**{'leftField': 'SPACECRAFT_NAME','rightField': 'SPACECRAFT_NAME'}),
          ee.Filter.equals(**{'leftField': 'SENSING_ORBIT_NUMBER','rightField': 'SENSING_ORBIT_NUMBER'})
          )
    })
  ).map(algorithm = setTime)
  daily = daily.sort('date',bool(1))
  #daily = daily.first()
  return daily

def getCover(image):
  pix = ee.Image(1).clip(AOI)
  totPixels = ee.Number(pix.reduceRegion(**{
      'reducer': ee.Reducer.count(),
      'scale': scale,
      #'maxPixels': 9999999,
      #'bestEffort': 'true',
      'geometry': AOI
  }).values().get(0))

  cldPixels = ee.Number(ee.Image(image.select('CloudMask')).reduceRegion(**{
      'reducer': ee.Reducer.sum(),
      'scale': scale,
      #'maxPixels': 9999999,
      #'bestEffort': 'true',
      'geometry': AOI
  }).values().get(0))

  actPixels = ee.Number(image.reduceRegion(**{
      'reducer': ee.Reducer.count(),
      'scale': scale,
      #'maxPixels': 9999999,
      #'bestEffort': 'true',
      'geometry': AOI
  }).values().get(0))

  percAoiCover = actPixels.divide(totPixels).multiply(100)
  percCldCover = cldPixels.divide(totPixels).multiply(100)
  image = image.set('actPixels', actPixels)
  image = image.set('cldPixels', cldPixels)
  image = image.set('totPixels', totPixels)
  image = image.set('AoiCover', percAoiCover)
  image = image.set('perCloud', percCldCover)
  return image

MAX_CLOUD_PROBABILITY = 85;

def maskClouds(img):
  clouds = ee.Image(img.get('cloud_mask')).select('probability')
  isCloud = clouds.gt(MAX_CLOUD_PROBABILITY)
  isCloud = isCloud.rename('CloudMask')
  img = img.addBands(isCloud)
  return img

def calc_geomArea(feat):
  pix = ee.Image(1).clip(feat)
  totPixels = ee.Number(pix.reduceRegion(**{
    'reducer': ee.Reducer.count(),
    'scale': 10,
    'geometry': feat,
    }).values().get(0))
  return totPixels.multiply(100).divide(1000000)

def findLakes1(sat_img):
  ndwi = sat_img.normalizedDifference(['B3', 'B8'])
  lakes = ndwi.gt(0)
  lakes = lakes.updateMask(lakes)
  return lakes

def findLakes2(sat_img):
  MNDWI_threshold=0.42; #testing shows recommended 0.42 for Sentinel-2 and Landsat 8. For the scene in article [1] it was 0.8.
  NDWI_threshold=0.4; #testing shows recommended 0.4 for Sentinel-2 and Landsat 8. For the scene in article [1] it was 0.5.
  ndvi = sat_img.normalizedDifference(['B8', 'B4'])
  mndwi = sat_img.normalizedDifference(['B3', 'B11'])
  ndwi = sat_img.normalizedDifference(['B3', 'B8'])
  ndwi_leaves = sat_img.normalizedDifference(['B8', 'B11'])
  aweish = sat_img.expression(
    'B2 + (2.5 * B3) - (1.5 * (B8 + B11)) - (0.25 * B12)',
    {'B2': sat_img.select('B2'),
     'B3': sat_img.select('B3'),
     'B8': sat_img.select('B8'),
     'B11': sat_img.select('B11'),
     'B12': sat_img.select('B12')
    })
  aweinsh = sat_img.expression(
    '4 * (B3 - B11) - ((0.25 * B8) + (2.75 * B11))',
    {'B3': sat_img.select('B3'),
     'B8': sat_img.select('B8'),
     'B11': sat_img.select('B11')
    })
  dbsi = sat_img.normalizedDifference(['B11', 'B3'])
  dbsi = dbsi.subtract(ndvi)
  wii = (sat_img.select('B8').pow(2)).divide(sat_img.select('B4'))
  wri = sat_img.expression(
    '(B3 + B4)/(B8 + B11)',
    {'B3': sat_img.select('B3'),
     'B4': sat_img.select('B4'),
     'B8': sat_img.select('B8'),
     'B11': sat_img.select('B11')
    })
  puwi = sat_img.expression(
    '(5.83 * B3) - (6.57 * B4) - (30.32 * B8) + 2.25',
    {'B3': sat_img.select('B3'),
     'B4': sat_img.select('B4'),
     'B8': sat_img.select('B8')
    })
  uwi = sat_img.expression(
    '(B3 - (1.1 * B4) - (5.2 * B8) + 0.4)',
    {'B3': sat_img.select('B3'),
     'B4': sat_img.select('B4'),
     'B8': sat_img.select('B8')
    })
  uwi2 = sat_img.expression(
    'B3 - (1.1 * B4) - (5.2 * B8)',
    {'B3': sat_img.select('B3'),
     'B4': sat_img.select('B4'),
     'B8': sat_img.select('B8')
    })
  uwi = uwi.divide(uwi2.abs())
  usi = sat_img.expression(
    '(0.25 * (B3/B4)) - (0.57 * (B8/B3)) - (0.83 * (B2/B3)) + 1',
    {'B2': sat_img.select('B2'),
     'B3': sat_img.select('B3'),
     'B4': sat_img.select('B4'),
     'B8': sat_img.select('B8')
    })
  lakes = (mndwi.gt(MNDWI_threshold)).Or(
    ndwi.gt(NDWI_threshold)).Or(aweinsh.gt(0.1879)).Or(aweish.gt(0.1112)).Or(
      ndvi.lt(-0.2)).Or(ndwi_leaves.gt(1))
  lakes = lakes.And(((aweinsh.gt(-0.03)).And(dbsi.lte(0))).Not())
  lakes = lakes.updateMask(lakes)
  lakes = lakes.rename('lakes')
  return lakes

def findLakeArea(img):
  plakes = findLakes2(img)
  lakes = findLakes1(img)
  lakes = (plakes.unmask()).Or(lakes.unmask())
  # int_lakes = lakes.updateMask(lakes)
  return lakes

addressLat = -10.9282
addressLon = 22.2698

tabWidth = 1000
tabHeight = 800

api_key = '66e7817752589937542751ydwf45f5c'
address = 'Los Angeles'

locationMode = st.sidebar.radio(
    "Select location input mode",
    ["Address", "Coordinates"],
    index=None,
)

if locationMode=="Address":
  address = st.sidebar.text_input("Location", value = "Los Angeles")
  fwdGeoCodeUrl = 'https://geocode.maps.co/search?q={address}&api_key={api_key}'.format(api_key=api_key, address = address)
  response = requests.get(fwdGeoCodeUrl)
  data = response.json()
  # addressResponse = json.loads(str(response))
  # formattedReponse = json.dumps(data, indent=2)
  # print(formattedReponse)
  print('Selected Address: ', address)
  addressLat = float(data[0]['lat'])
  addressLon = float(data[0]['lon'])
  dispName = data[0]['display_name']
  st.sidebar.write("Centroid Address: ", dispName)
elif locationMode=="Coordinates":
  sidebarCol1, sidebarCol2 = st.sidebar.columns((2, 2), gap='medium', vertical_alignment="center")
  addressLon = st.sidebar.number_input("Lon", min_value = -180.0, max_value = 180.0, value = addressLon)
  addressLat = st.sidebar.number_input("Lat", min_value = -90.0, max_value = 90.0, value = addressLat)
else:
  addressLon = st.sidebar.number_input("Lon", min_value = -180.0, max_value = 180.0, value = addressLon)
  addressLat = st.sidebar.number_input("Lat", min_value = -90.0, max_value = 90.0, value = addressLat)

radius = st.sidebar.slider("Seach radius (km)", 1, 1, 10)
maxPerCloud = st.sidebar.slider("Pick max permissible clouds cover", 0, 75, 30)
today = date.today()
earliestDate = today + timedelta(days=-180)
print('Earliest date: ', earliestDate)
print('Today: ', today)
refDate = st.sidebar.date_input(
    "Reference Date for Burn Scar Detection", value=datetime(2024,8,11), min_value=earliestDate, max_value=today, format="YYYY-MM-DD",
)
timeWindow = st.sidebar.number_input("Pick time window for change detection", min_value = 5, max_value = 180, value = 7)
buttonState = st.sidebar.button("Request Data", type="primary")

print('Search Radius: ', radius)
print('Max Cloud Cover: ', maxPerCloud)
print("Reference date: ", refDate)
print("Time window: ", timeWindow)

tab1Col1, tab1Col2 = st.columns((4, 1), gap='medium', vertical_alignment="center")

filteredResults2 = 0

with tab1Col1:
  Map = geemap.Map(center=[addressLat, addressLon], zoom=12)
  Map.add_basemap("Satellite")

  if buttonState and locationMode:    
      AOI = ee.Geometry.Point([addressLon, addressLat]).buffer(1000*radius).bounds()
      Map.centerObject(AOI)

      if radius and buttonState:
          # Map.addLayer(AOI, {'color': 'blue'}, 'Search Area')
          # Map.to_streamlit(height=tabHeight)

          # stop_date = date.today()
          stop_date = ee.Date(str(refDate))
          # stop_date = ee.Date(str(stop_date))
          start_date = stop_date.advance(-1*timeWindow,'day')
          minAoiCover = 95
          # scale = 50

          Sen2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED").filter(ee.Filter.date(start_date, stop_date)).filterBounds(AOI)
          s2Clouds = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY').filter(ee.Filter.date(start_date, stop_date)).filterBounds(AOI)

          s2SrWithCloudMask = ee.Join.saveFirst('cloud_mask').apply(**{
            'primary': Sen2,
            'secondary': s2Clouds,
            'condition': ee.Filter.equals(**{'leftField': 'system:index', 'rightField': 'system:index'})
          })

          datasetCloudMasked = ee.ImageCollection(s2SrWithCloudMask).map(algorithm = maskClouds)
          datasetCloudMasked = getDailyMosaic(datasetCloudMasked)
          datasetCloudMasked = datasetCloudMasked.map(algorithm = getCover)

          queriedOutput = datasetCloudMasked.filter(ee.Filter.gte('AoiCover', minAoiCover))
          print('Results with AOI >= ', minAoiCover,'%: ',queriedOutput.size().getInfo())
          queriedOutput = queriedOutput.filter(ee.Filter.lt('perCloud', maxPerCloud))
          filteredResults2 = queriedOutput.size().getInfo()
          print('Results with Cloud cover < ', maxPerCloud,'%: ',filteredResults2)


          if filteredResults2 >= 2:
            waterRes = queriedOutput.min()
            ndwi = waterRes.normalizedDifference(['B3', 'B8'])
            water = ndwi.gt(-0.1)

            # queriedOutput = queriedOutput.distinct('system:index')

            maxResults = 10
            filteredDataset = queriedOutput.limit(maxResults,'system:time_start', False)
            print('Number of cloud filtered & sorted results: ', filteredDataset.size().getInfo())

            prevImg = filteredDataset.sort('system:time_start', True).first()
            # prev_waterMask = findLakeArea(prevImg)
            prevImg1 = prevImg.updateMask(water.eq(0))
            print('previous date: ', ee.Date(prevImg.get('system:time_start')).format('yyyy-MM-dd').getInfo())

            BAI_1 = prevImg1.expression('(1-sqrt((B06*B07*B8A)/B04))*((B12-B8A)/(sqrt(B12+B8A))+1)',
              {'B04': prevImg.select('B4').multiply(0.0001),
              'B06': prevImg.select('B6').multiply(0.0001),
              'B07': prevImg.select('B7').multiply(0.0001),
              'B8A': prevImg.select('B8A').multiply(0.0001),
              'B12': prevImg.select('B12').multiply(0.0001)
              })
            
            fire_prev = prevImg1.select('B12').divide(prevImg1.select('B11'))
            fire_prev_mask = fire_prev.gt(1.5)
            fire_prev_mask = fire_prev_mask.updateMask(fire_prev_mask)
            fire_1 = fire_prev_mask.reduceToVectors(**{
              'geometry': AOI,
              'scale': 20,
              # 'maxPixels': 9999999,
              # 'bestEffort': True,
              'geometryType': 'centroid',
              'labelProperty': 'PrevFires',
            })

            latestImg = filteredDataset.first()
            # latest_waterMask = findLakeArea(latestImg)
            latestImg1 = latestImg.updateMask(water.eq(0))
            print('latest date: ', ee.Date(latestImg.get('system:time_start')).format('yyyy-MM-dd').getInfo())

            BAI_2 = latestImg1.expression('(1-sqrt((B06*B07*B8A)/B04))*((B12-B8A)/(sqrt(B12+B8A))+1)',
              {'B04': latestImg1.select('B4').multiply(0.0001),
              'B06': latestImg1.select('B6').multiply(0.0001),
              'B07': latestImg1.select('B7').multiply(0.0001),
              'B8A': latestImg1.select('B8A').multiply(0.0001),
              'B12': latestImg1.select('B12').multiply(0.0001)
              })
            
            BAI_diff = BAI_2.subtract(BAI_1)

            fire_latest = latestImg1.select('B12').divide(latestImg1.select('B11'))
            fire_latest_mask = fire_latest.gt(1.5)
            fire_latest_mask = fire_latest_mask.updateMask(fire_latest_mask)
            fire_2 = fire_latest_mask.reduceToVectors(**{
              'geometry': AOI,
              'scale': 20,
              # 'maxPixels': 9999999,
              # 'bestEffort': True,
              'geometryType': 'centroid',
              'labelProperty': 'PrevFires',
            })

            burnedAreaMask = BAI_diff.gt(0.2);
            burnedAreaMask = burnedAreaMask.updateMask(burnedAreaMask)
            burnScars = burnedAreaMask.reduceToVectors(**{
              'geometry': AOI,
              'scale': 20,
              # 'maxPixels': 9999999,
              # 'bestEffort': true,
              'geometryType': 'polygon',
              'labelProperty': 'burnt_shape',
            })

            burnArea = calc_geomArea(burnScars).getInfo()

            visualization = {
              'min': 0,
              'max': 3500,
              'bands': ['B4', 'B3', 'B2']
            }

            swirVisualization = {
              'min': 1000,
              'max': 5500,
              'bands': ['B12', 'B8A', 'B4']
            }

            Map.addLayer(latestImg.clip(AOI), swirVisualization, 'SWIR');
            # Map.addLayer(latestImg.clip(AOI), visualization, 'VIS');
            print('Past fire count: ', fire_1.size().getInfo())
            # Map.addLayer(fire_1, {'color':'yellow'}, 'FireCenters1')
            # Map.addLayer(water.updateMask(water),{'min': 0, 'max': 1, 'palette': ['white', 'blue']}, 'Prev Water')
            # Map.addLayer(latest_waterMask.updateMask(latest_waterMask),{'min': 0, 'max': 1, 'palette': ['white', 'blue']}, 'Water')
            print('Latest fire count: ', fire_2.size().getInfo())
            print('-------------------------------------------------------------')
            Map.addLayer(burnScars, {'color': ' brown'}, 'Burn Scars')
            Map.addLayer(fire_2, {'color':'orange'}, 'Fire Locations')
          Map.to_streamlit(height=tabHeight)
          # print(datetime.fromtimestamp(dateList[0]).strftime('%Y-%m-%d'))
  else:
    Map.to_streamlit(height=tabHeight)

with tab1Col2:
  if radius and buttonState:
    if filteredResults2 >= 2:
      st.metric(label="Active Fires", value=fire_2.size().getInfo(),
                      delta=(fire_2.size().getInfo() - fire_1.size().getInfo()), delta_color="inverse", help="Actively burning fires")
      st.metric(label="Fresh Burned Area [km2]", value="{:.2f}".format(burnArea), help="Newly burned area")
