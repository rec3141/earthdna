// first load the phenotype database
var phenotypes = {};
Cesium.loadText('data/ProTraits_positive.txt').then(function(textData) {
    var datafrom = "ProTraits: ";
	var lines = textData.split("\n");
	var len = lines.length;
	for (var i = 0; i < len; i++) {
		var line = lines[i].split("\t");
		var dataindex = datafrom + line[2];
		//665942	Desulfovibrio sp. 6_1_46AFAA	metabolism=sulfatereducer	0.805
		if( !Cesium.defined(phenotypes[dataindex]) ) {
			phenotypes[dataindex] = { taxid: [line[0]], taxon: [line[1]], value: [line[3]] };	
		} else {
			phenotypes[dataindex].taxid.push(line[0]);
			phenotypes[dataindex].taxon.push(line[1]);
			phenotypes[dataindex].value.push(line[3]);
		}
	}
}).otherwise(function(error) {
	console.log("failed to load phenotypes " + error);
});


// this is an event listener which is called when the geocoder is finished
// if a phenotype was selected, it plots the selected phenotype
// and make a toggle switch to turn it off and on

function plotPhenotype(){

    var pheno = viewer.geocoder.viewModel.searchText;
    var phenolist = phenotypes[pheno];
    console.log(phenolist);
    if (typeof phenolist == "undefined") return null;

    var taxaIndex = makeTaxaIndex();
	var phenoPPC = viewer.scene.primitives.add(new Cesium.PointPrimitiveCollection());
    var color = Cesium.Color.fromRandom( {minimumRed : 0.25, minimumGreen : 0.25, minimumBlue : 0.25, alpha : 1} );
    var len = phenolist.taxon.length;
    for (var i=0; i<len; i++) {
        var matched = taxaIndex[phenolist.taxon[i]] || taxaIndex[phenolist.taxon[i].split(" ")[0]] || undefined;
        if (typeof matched == "undefined") continue;
        var position = matched.position;

        phenoPPC.add(
        {
            description : 
                '<h2>Taxon: ' + phenolist.taxon[i] + '</h2>' +
                '<h2>Taxid: ' + phenolist.taxid[i] + '</h2>' +
                '<h2>Support: ' + phenolist.value[i]*100 + '%</h2>',
            position : matched.position,
            color : color,
            outlineColor : Cesium.Color.WHITE,
            pixelSize : 25 * Math.sqrt(phenolist.value[i]),
            id : {taxon: phenolist.taxon[i], abundance: phenolist.value[i]},
            translucencyByDistance : new Cesium.NearFarScalar(4e4, 1.0, 1e8, 0.0),
            show : true
        }
        );

        
    };

    console.log("plotphenotype3");

    function togglePhenotype() {
        var len = phenoPPC.length;
            for (var i = 0; i < len; ++i) {
              var l = phenoPPC.get(i);
              l.show = !l.show;
            };
    };

    console.log("plotphenotype4");

    Sandcastle.addToolbarButton('Toggle '+pheno, function(){
        togglePhenotype();
    });


viewer.camera.zoomOut(1000000.0);

};





/**
 * This class is an example of a custom geocoder. It provides geocoding by searching inside a LabelCollection.
 * @alias LabelCollectionGeocoder
 * @constructor
 */

function PhenotypeGeocoder() {
}

/**
 * The function called to geocode using this geocoder service.
 *
 * @param {String} input The query to be sent to the geocoder service
 * @returns {Promise<GeocoderResult[]>}
 */

// I could output bounding boxes for all of the tree taxa...

PhenotypeGeocoder.prototype.geocode = function (input) {
    
    var searchtext = input;
    var searchlist = [];
    var returnObject = {};
	if (input.length < 2) return null;
    for (pheno in phenotypes) {
          if ( pheno.toLowerCase().indexOf( searchtext.toLowerCase() ) > -1 ) {
            searchlist.push(pheno);
	}
    }
    
// console.log(searchlist);

// carry on and parse the searchlist
    return Cesium.loadText("")
        .then(function() {
            return searchlist.map(function (resultObject) {

				returnObject = {
                    displayName: resultObject,
                    destination: viewer.camera.computeViewRectangle()
                };
                return returnObject;
            });
        });
};
