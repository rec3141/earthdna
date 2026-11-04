/**
 *  	THIS EVENTUALLY NEEDS TO BE BIOMDataSource
 *		USING THE BIOM FORMAT https://github.com/biocore/biom-format/blob/master/examples/rich_sparse_otu_table.biom
 *  
 * 		RIGHT NOW I'M JUST USING THIS TO LOAD RAW NR DATA QUICKLY
 *		
 * @alias NRPrimDataSource
 * @constructor
 *
 * @param {String} [name] The name of this data source.  If undefined, a name
 *                        will be derived from the url.
 *
 * @example
 * var dataSource = new Cesium.NRPrimDataSource();
 * dataSource.loadUrl('sample.json');
 * viewer.dataSources.add(dataSource);
 * 
 *
 * more comments available at the Cesium Sandbox from which this was copied
 * https://cesiumjs.org/Cesium/Apps/Sandcastle/index.html?src=Custom%20DataSource.html&label=DataSources
 */
// console.log(1)

function NRPrimDataSource(name) {
    //All public configuration is defined as ES5 properties
    //These are just the "private" variables and their defaults.
    this._name = name;
    this._changed = new Cesium.Event();
    this._error = new Cesium.Event();
    this._isLoading = false;
    this._loading = new Cesium.Event();
    this._entityCollection = new Cesium.PointPrimitiveCollection();
//     this._seriesNames = [];
//     this._seriesToDisplay = undefined;
//     this._abundanceScale = 1;
//     this._entityCluster = new Cesium.EntityCluster();
}
// console.log(2)

Object.defineProperties(NRPrimDataSource.prototype, {
    name : { get : function() { return this._name; } },
    clock : { value : undefined, writable : false },
    entities : { get : function() { return this._entityCollection; } },
    isLoading : { get : function() { return this._isLoading; } },
    changedEvent : { get : function() { return this._changed; } },
    errorEvent : { get : function() { return this._error; } },
    loadingEvent : { get : function() { return this._loading; } },
//     seriesNames : { get : function() { return this._seriesNames; } },

//     seriesToDisplay : {
//         get : function() {
//             return this._seriesToDisplay;
//         },
//         set : function(value) {
//             this._seriesToDisplay = value;
//             
//             //Iterate over all entities and set their show property
//             //to true only if they are part of the current series.
//             var collection = this._entityCollection;
//             var entities = collection.values;
//             collection.suspendEvents();
//             for (var i = 0; i < entities.length; i++) {
//                 var entity = entities[i];
//                 entity.show = value === entity.seriesName;
//             }
//             collection.resumeEvents();
//         }
//     },

//     abundanceScale : {
//         get : function() { return this._abundanceScale; },
//         set : function(value) { if (value > 0) { throw new Cesium.DeveloperError('value must be greater than 0'); } this._abundanceScale = value; } 
//     },
    show : {
        get : function() { return this._entityCollection; },
        set : function(value) { this._entityCollection = value; } 
    },
//     clustering : {
//         get : function() { return this._entityCluster; },
//         set : function(value) { if (!Cesium.defined(value)) { throw new Cesium.DeveloperError('value must be defined.'); } this._entityCluster = value; } 
//     }
});

// console.log(3)

NRPrimDataSource.prototype.loadUrl = function(url) {
//     if (!Cesium.defined(url)) { throw new Cesium.DeveloperError('url is required.'); }

    var name = Cesium.getFilenameFromUri(url);
// console.log("3a");

    if (this._name !== name) { this._name = name; this._changed.raiseEvent(this); }
// console.log("3b");

    //Use 'when' to load the URL into a json object
    //and then process it with the `load` function.
    var that = this;
    return Cesium.when(Cesium.loadJson(url), function(json) {
	// console.log("3c");

        return that.load(json, url);
    }).otherwise(function(error) {
        //Otherwise will catch any errors or exceptions that occur 
        //during the promise processing. When this happens, 
        //we raise the error event and reject the promise. 
        this._setLoading(false);
        that._error.raiseEvent(that, error);
		// console.log("3d");
        return Cesium.when.reject(error);
    });

};
 
// makesequence = function(start, count) {
//   return Array.apply(0, Array(count))
// 	.map(function (element, index) { 
// 	  return index + start;  
//   });
// };

/**
 * Loads the provided data, replacing any existing data.
 * @param {Object} data The object to be processed.
 */

// console.log(4);

NRPrimDataSource.prototype.load = function(data) {
    //>>includeStart('debug', pragmas.debug);
//     if (!Cesium.defined(data)) {
//         throw new Cesium.DeveloperError('data is required.');
//     }
    //>>includeEnd('debug');
// console.log("4a");

    //Clear out any data that might already exist.
    this._setLoading(true);
//     this._seriesNames.length = 0;
//     this._seriesToDisplay = undefined;
// console.log("4b");

//     var abundanceScale = this.abundanceScale;
    var entities = this._entityCollection;
// console.log("4c");

    //It's a good idea to suspend events when making changes to a 
    //large amount of entities.  This will cause events to be batched up
    //into the minimal amount of function calls and all take place at the
    //end of processing (when resumeEvents is called).
//     entities.suspendEvents();
//    entities.removeAll();

    //REC JSON is an array of series, where each series itself is an
    //array of two items, the first containing the series name, the second
    //containing the heirarchical taxon level that defines fonts/viewability, and the third
    //being an array of repeating name, latitude, longitude, abundance  values.
    //
    //Here's a more visual example.
    //[["series1", taxlevel, color, ["name", latitude, longitude, abundance, ... ]
    // ["series2", taxlevel, color, ["name", latitude, longitude, abundance, ... ]]

    // Loop over each series
    for (var x = 0; x < data.length; x++) {
        var series = data[x];
        var seriesName = series[0];
        var taxlevel = series[1];
        var color = Cesium.Color.fromCssColorString(series[2]);
        var coordinates = series[3];
		
        //Add the name of the series to our list of possible values.
//         this._seriesNames.push(seriesName);

        //Make the first series the visible one by default
//         var show = x === 0;
//         if (show) {
//             this._seriesToDisplay = seriesName;
//         }

		//get range of sizes in order to properly scale dots
		var abundances = [];
		for (var i = 0; i < coordinates.length; i += 4) {			
			abundances.push(coordinates[i + 3]);
		}
		var maxabund = Math.max(...abundances);
		var minabund = Math.min(...abundances);
		var constantsize = 0;
		if (minabund == maxabund) { constantsize=1; }
		var maxradius = 20; //px
		var minradius = 4; //px

		var taxlevels = ["domain","major_clade","superkingdom","kingdom","subkingdom","superphylum","phylum","subphylum","infraphylum","superclass","class","subclass","infraclass","superorder","order","suborder","superfamily","family","subfamily","genus", "species","strain","otu"];
		var taxscalar = {};
		for (var i = 0; i < taxlevels.length; i += 1) {
			var range = {};
			range.max = 5e3 - (5e3 - 1e3)/taxlevels.length * i;
			range.max = range.max * range.max;
			range.min = range.max/2;
			if (i > 18) range.min = 1;
			taxscalar[taxlevels[i]] = range;
		}
			// console.log(minabund);
			// console.log(maxabund);

        //Now loop over each coordinate in the series and create
        // our entities from the data.
        for (var i = 0; i < coordinates.length; i += 4) {
        	var taxon = coordinates[i];
            var latitude = coordinates[i + 1];
            var longitude = coordinates[i + 2];
            var abundance = coordinates[i + 3];
            var normabundance = maxabund;
            var normsize = 3*maxabund;
            if (constantsize==0) {
                    normsize = minradius + maxradius*Math.sqrt(abundance/maxabund);
//                 normabundance = Math.sqrt( (1 + abundance - minabund)/(maxabund - minabund) );
//                 normsize = minradius + (maxradius-minradius)*normabundance;
            }

            //Ignore lines of zero abundance. may want to remove this? hahahah
//            if(abundance === 0) {
                //continue;
//            }

//            var color = Cesium.Color.fromHsl((0.6 - (normabundance * 0.5)), 1.0, 0.5);

			entities.add({
              id : {taxon: taxon, abundance: abundance},
//               abundance : abundance,
 			  position : Cesium.Cartesian3.fromDegrees(longitude, latitude),
 			  pixelSize : normsize,
 			  color : color,
 			  translucencyByDistance : new Cesium.NearFarScalar(taxscalar[taxlevel].min, 1.0, taxscalar[taxlevel].max, 0.0)
			});

// 			var entity = new Cesium.Entity({
//               id : seriesName + ' index ' + i.toString(),
// 			  name : taxon,
// 			  position : Cesium.Cartesian3.fromDegrees(longitude, latitude),
// 			  point : {
// 				pixelSize : normsize,
// 				color : Cesium.Color.BLACK,
// // 				outlineColor : Cesium.Color.WHITE,
// // 				outlineWidth : 2
// 				translucencyByDistance : new Cesium.NearFarScalar(taxscalar[taxlevel].min, 1.0, taxscalar[taxlevel].max, 0.0)
// 			  },
// //              description : '<h1>' + taxon + '</h1><h1>' + abundance + '</h1>',
// 			  label : {
// 				text : taxon,
// 				font : '12pt sans-serif',
// //				style: Cesium.LabelStyle.FILL,
// //				outlineWidth : 2,
// //				verticalOrigin : Cesium.VerticalOrigin.BOTTOM,
// //				pixelOffset : new Cesium.Cartesian2(0, -9),
// 				translucencyByDistance : new Cesium.NearFarScalar(taxscalar[taxlevel].min, 1.0, taxscalar[taxlevel].max/5, 0.0)
// 			  },
//               show : show,
// //              seriesName : seriesName //Custom property to indicate series name
// 			});
// 
//             //Add the entity to the collection.
//             entities.add(entity);

        }
        
        // console.log(entities);

    }
// console.log(6);

    //Once all data is processed, call resumeEvents and raise the changed event.
     entities.resumeEvents();
     this._changed.raiseEvent(this);
     this._setLoading(false);
};
// console.log(7);

NRPrimDataSource.prototype._setLoading = function(isLoading) {
    if (this._isLoading !== isLoading) {
        this._isLoading = isLoading;
        this._loading.raiseEvent(this, isLoading);
    }
};
// console.log(8);


