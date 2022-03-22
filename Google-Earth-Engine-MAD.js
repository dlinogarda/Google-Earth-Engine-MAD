var L8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR"),
    point = 
    /* color: #d63000 */
    /* shown: false */
    ee.Geometry.Polygon(
        [[[119.7549582444061, 22.91017694961182],
          [120.1097823105682, 22.829195314437055],
          [120.47559270485532, 22.748165467752404],
          [120.72587499245297, 23.679588057011674],
          [119.97193822487485, 23.83418946226891],
          [119.7549582444061, 22.91017694961182]]], null, false);


var Geometry = /* color: #d63000 */ee.Geometry.Polygon(
        [[[119.7549582444061,22.91017694961182],
          [120.1097823105682,22.829195314437055],
          [120.47559270485532,22.748165467752404],
          [120.72587499245297,23.679588057011674],
          [119.97193822487485,23.83418946226891],
          [119.7549582444061,22.91017694961182]]],null,false);

//call image and select band respectively
var image1 = L8.filterDate('2018-01-01', '2018-05-31')
       .filter(ee.Filter.and(ee.Filter.eq('WRS_PATH', 118),         
                      ee.Filter.eq('WRS_ROW', 44)))
       .filter(ee.Filter.lt('CLOUD_COVER', 40))
       //.map(maskL8sr)
       .select("B4", "B3", "B2", "B5")
       .first()
       .clip(point);
       
print(image1)

var image2 =ee.Image('LANDSAT/LC08/C01/T1_SR/LC08_118044_20171011')
       .select("B4", "B3", "B2", "B5").clip(point)
print(image2)  

//register image
var image2 = image2.register(image1,60)

//mean centered operations
var Scale = image1.projection().nominalScale()
var bandNames = image1.bandNames()

var meanDict1 =image1.reduceRegion({
              reducer: ee.Reducer.mean(),
              geometry: Geometry,
              scale   : Scale,
              maxPixels: 1e9
});

var means1 = ee.Image.constant(meanDict1.values(bandNames));
var centered1 = image1.subtract(means1);

var meanDict2 =image2.reduceRegion({
              reducer: ee.Reducer.mean(),
              geometry: Geometry,
              scale   : Scale,
              maxPixels: 1e9
});
var means2 = ee.Image.constant(meanDict2.values(bandNames));
var centered2 = image2.subtract(means2);

//convert image to Array
var array1 = centered1.toArray()
var arrays1 = array1.toArray(1) // array image1
print(arrays1)


var array2 = centered2.toArray()
var arrays2 = array2.toArray(1) //array image2
print(arrays2)

//covariance array
var covar1 = array1.reduceRegion({
                reducer: ee.Reducer.centeredCovariance(),
                geometry: Geometry,
                scale   : Scale,
                maxPixels: 1e9
                });
                
var covA = ee.Array(covar1.get('array'));
print(covA,'covariance array ref')

var covar2 = array2.reduceRegion({
                reducer: ee.Reducer.centeredCovariance(),
                geometry: Geometry,
                scale   : Scale,
                maxPixels: 1e9
                });
                
var covB = ee.Array(covar2.get('array'));
print(covB,'covariance array trg')

var covAA = covA.matrixMultiply(covA.matrixInverse())
print(covAA, 'covariance xx')

var covAB = covA.matrixMultiply(covB.matrixInverse())
print(covAB, 'covariance xy')

var covBA = covB.matrixMultiply(covA.matrixInverse())
print(covBA, 'covariance yx')

var covBB = covB.matrixMultiply(covB.matrixInverse())
print(covBB, 'covariance yy')

//CCA

var A = (covAA.matrixInverse()).matrixMultiply(covAB).matrixMultiply(covBB.matrixInverse()).matrixMultiply(covBA)
print(A,'CCA A')

var B = (covBB.matrixInverse()).matrixMultiply(covBA).matrixMultiply(covAA.matrixInverse()).matrixMultiply(covAB)
print(B,'CCA B')

//eigens solutions

var eigenA = A.eigen()
print(eigenA,'eigenA')

  //eigen value and eigen vector solutions reference
var eigvalA = eigenA.slice(1,0,1)
print(eigvalA,'eigen values A')

var eigvecA = eigenA.slice(1,1)
print(eigvecA,'eigen vector A')


var eigenB = B.eigen()
print(eigenB,'eigenB')

  //eigen value and eigen vector solutions reference
var eigvalB = eigenB.slice(1,0,1)
print(eigvalB,'eigen values B')

var eigvecB = eigenB.slice(1,1)
print(eigvecB,'eigen vector B')

  //MAD Variates
var U = ee.Image(eigvecA).matrixMultiply(arrays1)
print(U)
var V = ee.Image(eigvecB).matrixMultiply(arrays2)
print(V)
var MAD = U.subtract(V)
print(MAD)

var nBands = bandNames.length().divide(4)
print(nBands)
var bNames = centered1.bandNames()
print(bNames)
var bNames1 = bNames.slice(nBands)
var bNames2 = bNames.slice(nBands)

var MADs = MAD.arrayProject([0]).arrayFlatten([bNames]).float()
print(MADs)

Map.addLayer(image1,{ min :0, max:5000, gamma: 1},'ref')
Map.addLayer(image2,{ min :0, max:5000, gamma: 1},'trg')
Map.addLayer(MADs,{bands:['B4'], min :0, max:255},'MAD')

Map.setCenter(120.1331, 23.065,12)

  //Normalized MAD
var stdevmad =MADs.reduceRegion({
              reducer: ee.Reducer.stdDev(),
              geometry: Geometry,
              scale   : Scale,
              maxPixels: 1e9
});

var stdevsmad = ee.Image.constant(stdevmad.values(bandNames));
print(stdevsmad,'STD MAD')

var N1 = (MADs.select("B4").pow(2)).divide(stdevsmad.select('constant_0').pow(2))
var N2 = (MADs.select("B3").pow(2)).divide(stdevsmad.select('constant_1').pow(2))
var N3 = (MADs.select("B2").pow(2)).divide(stdevsmad.select('constant_2').pow(2))
var N4 = (MADs.select("B5").pow(2)).divide(stdevsmad.select('constant_3').pow(2))


var NMAD = N1.add(N2).add(N3).add(N4).rename('NMAD');
Map.addLayer(NMAD, {min :0, max:255},'NMAD')


// define function to display histogram

var showHistogram = function (image,name,band) {
    var options = {
                   title    : 'Histogram of ' + name,
                   fontSize : 20,
                   hAxis    : {title:'DN'},
                   vAxis    : {title:'Count of DN'},
                   series   :{ 0: {color: 'blue'},}
    };


 
  var histogram = ui.Chart.image.histogram(image.select(band),point,1024,2048)
    .setSeriesNames([band])
    .setOptions(options);
  print(histogram);
};

//show histogram ----> showHistogram = function (image,name,band)
showHistogram(NMAD,'NMAD','NMAD');


// apply threshold if NMAD value less then or equal 0 == 0 else will be 1
 var thres = NMAD.lte(2).rename('thres')// white (1) = no change - black (0) = change
 
 var tNMAD = NMAD.addBands(thres)
 print(tNMAD,'Threshold NMAD')
 
 Map.addLayer(tNMAD, {bands :['thres'],min :0, max:1, gamma: 0.75},' thresholded NMAD')
 
 
// Apply threshold as PIF
// Select the land/water mask.
var PIF = tNMAD.select('thres');

// Create a binary mask.
//all the pixels that do not have the value of 1 in the PIF band (those that are change data) get a value of 0 in the resulting image.
var mask = PIF.eq(1);

// Update the composite mask with the water mask.
var maskeref = image1.updateMask(mask);
var masketrg = image2.updateMask(mask);

Map.addLayer(maskeref,{ min :0, max:5000, gamma: 1},'ref PIF')
Map.addLayer(masketrg,{ min :0, max:5000, gamma: 1},'trg PIF')


// Radiometric Normalization 
  //Slope and Intercept Definition
  //Slope
var stdevPIFref =maskeref.reduceRegion({
              reducer: ee.Reducer.stdDev(),
              geometry: Geometry,
              scale   : Scale,
              maxPixels: 1e9
});

var stdevsPIFref = ee.Image.constant(stdevPIFref.values(bandNames));

var stdevPIFtrg =masketrg.reduceRegion({
              reducer: ee.Reducer.stdDev(),
              geometry: Geometry,
              scale   : Scale,
              maxPixels: 1e9
});

var stdevsPIFtrg = ee.Image.constant(stdevPIFtrg.values(bandNames));
print(stdevsPIFref,'STD pif Reference')
print(stdevsPIFtrg,'STD pif Target')

var slope = stdevsPIFref.divide(stdevsPIFtrg)
print(slope,'slope')

    //Mean
var meanref =maskeref.reduceRegion({
              reducer: ee.Reducer.mean(),
              geometry: Geometry,
              scale   : Scale,
              maxPixels: 1e9
});

var meansref = ee.Image.constant(meanref.values(bandNames));

var meantrg =masketrg.reduceRegion({
              reducer: ee.Reducer.mean(),
              geometry: Geometry,
              scale   : Scale,
              maxPixels: 1e9
});

var meanstrg = ee.Image.constant(meantrg.values(bandNames));

print(meansref,'mean pif Reference')
print(meanstrg,'mean pif Target')


var intercept = meansref.subtract(slope.multiply(meanstrg))
print(intercept,'intercept')


  //Normalization step (linier regresion)
var Norm = (image2.multiply(slope)).add(intercept)  
Map.addLayer(Norm,{ min :0, max:5000, gamma: 1},'Normalize image')