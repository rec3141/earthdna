/**
 * This class is an example of a custom geocoder. It provides geocoding by searching inside a LabelCollection.
 * @alias LabelCollectionGeocoder
 * @constructor
 */

function LabelCollectionGeocoder() {
}

/**
 * The function called to geocode using this geocoder service.
 *
 * @param {String} input The query to be sent to the geocoder service
 * @returns {Promise<GeocoderResult[]>}
 */

// I could output bounding boxes for all of the tree taxa...

LabelCollectionGeocoder.prototype.geocode = function (input) {
    
    var searchtext = input;
    var searchlist = [];

	if (input.length < 2) return null;
//    console.log(viewer.scene.primitives);
    var plen = viewer.scene.primitives.length;
    
    for (var i = 0; i < plen; ++i) {
      var gcLC = viewer.scene.primitives.get(i); 
      if (gcLC.constructor.name == "LabelCollection" || gcLC.constructor.name == "x") {
        var len = gcLC.length;
        for (var j = 0; j < len; ++j) {
          var k = gcLC.get(j);
          if ( k.text.toLowerCase().indexOf( searchtext.toLowerCase() ) > -1 ) {
            searchlist.push(k);
          }
        }
      }
    }
    
// if the searchlist is empty, look through the pointprimitives
	if (searchlist.length < 1) {
		for (var i = 0; i < plen; ++i) {
		  var gcPP = viewer.scene.primitives.get(i); 
		  if (gcPP.constructor.name == "PointPrimitiveCollection" || gcPP.constructor.name == "A") {
			var len = gcPP.length;
			for (var j = 0; j < len; ++j) {
			  var k = gcPP.get(j);
			  if ( k.id.taxon.toLowerCase().indexOf( searchtext.toLowerCase() ) > -1 ) {
				searchlist.push(k);
			  }
			}
		  }	
		}
    }
    
    //console.log(searchlist);

// UNIMPLEMENTED
// if the searchlist is still empty, do a MySQL query for synonyms in the full SILVA ref

//	if (searchlist.length < 1) return null;

// carry on and parse the searchlist
    return Cesium.loadText("")
        .then(function() {
            var bboxDegrees;
            var heightmin = 10000;
            var heightmax = 10000;
            var texto = "unknown";

            return searchlist.map(function (resultObject) {
                var lonlat = Cesium.Ellipsoid.WGS84.cartesianToCartographic(resultObject.position);
				if (typeof resultObject.distanceDisplayCondition !== "undefined") {
                	heightmin = resultObject.distanceDisplayCondition.near;
                	heightmax = resultObject.distanceDisplayCondition.far;
				};
                var horizdeg = Math.sqrt(.5*6371000*(heightmax+heightmin)/2)/111000;
                var nwlat = lonlat.latitude + Math.PI/180*horizdeg/2; if (nwlat > Math.PI/2) nwlat=(1-(nwlat/Math.PI/2) % 1) * Math.PI/2;
                var nwlon = lonlat.longitude + Math.PI/360*horizdeg; if (nwlon > Math.PI) nwlon=(nwlon/Math.PI - 1) % 1 * Math.PI;
                var swlat = lonlat.latitude - Math.PI/180*horizdeg/2; if (swlat < -Math.PI/2) swlat=(1+(swlat/Math.PI/2) % 1) * Math.PI/2;
                var swlon = lonlat.longitude - Math.PI/360*horizdeg; if (swlon < -Math.PI) swlon=(swlon/Math.PI + 1) % 1 * Math.PI;
                var carto = [
                        new Cesium.Cartographic(swlon, swlat, heightmin),
                        new Cesium.Cartographic(nwlon, nwlat, heightmax)
                            ];
                var recto = Cesium.Rectangle.fromCartographicArray(carto);
				if (typeof resultObject.text !== "undefined") { texto = resultObject.text }
				else if (typeof resultObject.id.taxon !== "undefined") { texto = resultObject.id.taxon }

				var returnObject = {
                    displayName: texto,
                    destination: recto
                };
//                console.log(returnObject);
                return returnObject;
            });
        });
};
