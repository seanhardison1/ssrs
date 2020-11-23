var dictionary = ee.Image.pixelLonLat().reduceRegion({
  reducer: ee.Reducer.toCollection(['longitude', 'latitude']), 
  geometry: station_locs2, 
  scale: 500
});
// print(dictionary)
var points = ee.FeatureCollection(dictionary.get('features'))
.map(function(feature) {
  var lon = feature.get('longitude');
  var lat = feature.get('latitude');
  return ee.Feature(ee.Geometry.Point([lon, lat]), {
    'featureID': ee.Number(lon).multiply(10000).round().format('%6.0f')
    .cat('_')
    .cat(ee.Number(lat).multiply(10000).round().format('%6.0f'))
  });
});
print('points', points);
// Map.addLayer(points, {palette: ["red"]})
//print('stations', station_locs)

var dataset = ee.ImageCollection("NASA/OCEANDATA/MODIS-Terra/L3SMI")
.filterDate('2010-02-26', '2020-07-30')
.select('Rrs_645');

// print('dataset', dataset)


Map.addLayer(dataset, {bands: 'Rrs_645', min: 0, max: 0.011, palette: ['blue', 'green', 'red'], opacity: 0.2});
Map.addLayer(station_locs, {});

var triplets = dataset.map(function(image) {
  return image.reduceRegions({
    collection: points, 
    reducer: ee.Reducer.first().setOutputs(image.bandNames()), 
    scale: 500,
  }).map(function(feature) {
    return feature.set({
      'imageID': image.id(),
      'timeMillis': image.get('system:time_start')
    });
  });
}).flatten();
//print("triplets",triplets) 

var format = function(table, rowId, colId, rowProperty, colProperty) {
  var rows = table.distinct(rowId); 
  //print('rows',rows)
  var joined = ee.Join.saveAll('matches').apply({
    primary: rows, 
    secondary: table, 
    condition: ee.Filter.equals({
      leftField: rowId, 
      rightField: rowId
    })
  });
  return joined.map(function(row) {
    var values = ee.List(row.get('matches'))
    .map(function(feature) {
      feature = ee.Feature(feature);
      return [feature.get(colId), feature.get(colProperty)];
    }).flatten();
    return row.select([rowId, rowProperty]).set(ee.Dictionary(values));
  });
};

var results = format(triplets, 'imageID', 'featureID', 'timeMillis', 'Rrs_645');
// print('results',results)

Export.table.toDrive({
  collection: results, 
  description: '2010-2020', 
  fileNamePrefix: 'modis_rrs_645_2010_2020', 
  fileFormat: 'KML'
});

