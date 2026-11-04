//// SETUP CESIUM VIEWER

/// transparent base layer as single pixel
var transparentBaseLayer = new Cesium.SingleTileImageryProvider({
    url : "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNgYAAAAAMAASsJTYQAAAAASUVORK5CYII="
});

/// we don't have terrain for now
var noTerrain = new Cesium.EllipsoidTerrainProvider({});

/// setup the viewer
var viewer = new Cesium.Viewer('cesiumContainer', {
	baseLayerPicker : false,
    skyBox : false,
    homeButton : false,
    skyAtmosphere : false,
    navigationHelpButton: true,
    animation: false,
	timeline: false,
    imageryProvider : transparentBaseLayer,
	terrainProvider : noTerrain,
    geocoder: [new PhenotypeGeocoder(), new LabelCollectionGeocoder()],
    contextOptions : {
        webgl: {
            alpha: true
        }
    }
});


//// MAPS

viewer.scene.imageryLayers.removeAll();
//viewer.scene.fxaa = false;
viewer.scene.backgroundColor = Cesium.Color.BLACK;
viewer.scene.globe.baseColor = Cesium.Color.ROYALBLUE; //DEEPSKYBLUE; //STEELBLUE; //DEEPSKYBLUE
viewer.scene.globe.show = true;

var bg_map = viewer.imageryLayers.addImageryProvider(new Cesium.SingleTileImageryProvider({
    url : 'data/map.mds-big.render.png',
    rectangle : Cesium.Rectangle.fromDegrees(-180, -90, 180, 90)
}));

function toggleMap() {
	// console.log(bg_map);
	bg_map.show=!bg_map.show;
}

Sandcastle.addToolbarButton('Toggle Map', function(){
    toggleMap();
});


//// where to start

viewer.camera.flyTo({
    destination : Cesium.Cartesian3.fromDegrees(0, 0, 1.1e7),
    duration : 0.1
});

//// Add event listener to rotate the globe
var spinRate=0.002;
function changeRate(scene, time){
	if (this.viewer.scene.mode == 3) {
	this.viewer.scene.camera.rotate(Cesium.Cartesian3.UNIT_Z, -spinRate);
	}
};
this.viewer.scene.postRender.addEventListener(changeRate);

// some code to remove default cesium logo
// Cesium.CreditDisplay.removeDefaultCredit();


//// GEOCODER
var searchbox = document.getElementsByClassName("cesium-geocoder-input")[0];
searchbox.setAttribute('placeholder', 'Enter a name, accession, or phenotype ...');
searchbox.style = "{width:300px}";

var suggestions = document.getElementsByClassName("search-results")[0];
suggestions.style = "{max-height:500px}";

// cesium-geocoder-input cesium-geocoder-input-wide
// 
// add listener to geocoder viewmodel
viewer.geocoder.viewModel.complete.addEventListener(plotPhenotype);
viewer.geocoder.viewModel.keepExpanded = true;
viewer.geocoder.viewModel.flightDuration = 0.5; //this can't be 0 or the plotPhenotype call doesn't work
// console.log(viewer.geocoder);


//// POINT INFORMATION / PICKER

viewer.infoBox.frame.removeAttribute('sandbox');
viewer.infoBox.frame.style = "{max-height:600px; height=100%;}";


function makeBox(picklist) {

    // console.log(picklist);
    var entity = new Cesium.Entity(); 
    var deschtml = '<div width=40% height=100%>';
    var taxa = {};
    //add onto object with taxa = {"taxa1" : [1,3,5,7,9], "taxa2" : [2,4,6,8] }
    for(var j=0; j<picklist.length; j++) {
        var pick = picklist[j];
        var taxon = picklist[j].id.taxon || picklist[j].id || undefined;
        if(Cesium.defined(taxa[taxon])) {
            taxa[taxon].push(j);
        } else {
            taxa[taxon] = [];
            taxa[taxon].push(j);
        }
    };

    //  		console.log(taxa);	 
    //  		console.log(Object.keys(taxa));
    
    for(var k in taxa) {
        var desc; var abund = "";
        var arr = taxa[k];
        for(var i=0; i<arr.length; i++) {
            var j = arr[i];
            var pick = picklist[j];
            var sample = pick.collection.name || "value";
            desc = pick.id.taxon || pick.id || pick.text || pick.description || "unknown";
            var pos = Cesium.Cartographic.fromCartesian(pick.primitive.position);
            // console.log(pos);
            desc = desc;// + ' ' + Cesium.Math.toDegrees(pos.latitude) + ' ' + Cesium.Math.toDegrees(pos.longitude);
            var whichabund = pick.id.abundance || "";
            var cattext = '<li>' + sample + ': ' + pick.id.abundance + '</li>';
            if(whichabund !== "") {
                abund = abund.concat(cattext);
            };
            // console.log(abund);
    };
    
    //     	on link click make a new entity with a new infobox that has the ncbi site in an iframe 
    var target = desc.split(".");
    var link = "taxonomy/?term=";
    if (target.length>1) { link="nuccore?term=" }
    target=target[0];

    deschtml = deschtml + '<p>taxon: <a target="_blank" href="https://www.ncbi.nlm.nih.gov/' + link + target + '">' + desc + '</a><ul>' + abund + '</ul></p>';
    }


    deschtml = deschtml + '</div>';
    entity.description = deschtml;
    entity.name = pick.name || "Taxonomy and Abundance";
    //    	console.log(entity);

    viewer.selectedEntity = entity;

}

// If the mouse is over the entity, add it to the infoBox
function doOnHover(movement) {
 	var picklist = viewer.scene.drillPick(movement.endPosition);
 	if (picklist.length>0 && typeof picklist[0].id !== "undefined") {
 	    makeBox(picklist);
 	};
};

var hover_handler = new Cesium.ScreenSpaceEventHandler(viewer.scene.canvas);
hover_handler.setInputAction(doOnHover, Cesium.ScreenSpaceEventType.MOUSE_MOVE);


// If the mouse clicks on the entity, turn off the hover handler
// If the mouse clicks anywhere, stop the globe rotation

function doOnSC(movement) {
// need to add handler for when there are multiple taxa under one click
 	var picklist = viewer.scene.drillPick(movement.position);
 	if (picklist.length>0 && typeof picklist[0].id !== "undefined") { 
 	    makeBox(picklist);
        hover_handler.removeInputAction(Cesium.ScreenSpaceEventType.MOUSE_MOVE);    
 	} else {
        hover_handler.setInputAction(doOnHover,Cesium.ScreenSpaceEventType.MOUSE_MOVE);
 	};
    
	var remove_rotate = viewer.scene.postRender.removeEventListener(changeRate);  
};

var sc_handler = new Cesium.ScreenSpaceEventHandler(viewer.scene.canvas);
sc_handler.setInputAction(doOnSC, Cesium.ScreenSpaceEventType.LEFT_CLICK);

// If the mouse double clicks on the globe, start the rotation

var dc_handler = new Cesium.ScreenSpaceEventHandler(viewer.scene.canvas);
dc_handler.setInputAction(function (movement) {
		var remove = viewer.scene.postRender.removeEventListener(changeRate); 
		if(!remove) {
			viewer.scene.postRender.addEventListener(changeRate);
 		}
}, Cesium.ScreenSpaceEventType.LEFT_DOUBLE_CLICK);

function toggleInfoBox(){
    var infobox = document.getElementsByClassName("cesium-infoBox")[0];
    // console.log(infobox);
    if (infobox.style.getPropertyValue('visibility') == 'hidden') { infobox.style.setProperty('visibility','visible') }
    else { infobox.style.setProperty('visibility','hidden') }
};

Sandcastle.addToolbarButton('Toggle InfoBox', function(){
    toggleInfoBox();
});

//// ADD LABELS


var tax_labels = viewer.scene.primitives.add(new Cesium.LabelCollection() );
Cesium.loadText('tax_labels_2.js').then(function(textData) {
	eval(textData);
// 	console.log(tax_labels);

	var len = tax_labels.length;

	for (var i=0; i < len; i++) {
		var lab =  tax_labels.get(i);
		lab.style = Cesium.LabelStyle.FILL_AND_OUTLINE;
		lab.outlineColor = Cesium.Color.BLACK;
		lab.outlineWidth = 3;
		lab.fillColor = Cesium.Color.LIGHTGREY;
		lab.scale = 0.5;
	};

}).otherwise(function(error) {
	console.log("failed tax_labels");
});


function toggleLabels() {
	// console.log(tax_labels);
//	tax_labels.show=!tax_labels.show;
	var len = tax_labels.length;
	for (var i = 0; i < len; ++i) {
	  var l = tax_labels.get(i);
	  l.show = !l.show;
	}
}

function darkLabelColors() {
	var len = tax_labels.length;
	for (var i = 0; i < len; ++i) {
	  var l = tax_labels.get(i);
	  l.fillColor = Cesium.Color.BLACK;
	};

};

function lightLabelColors() {
	var len = tax_labels.length;
	for (var i = 0; i < len; ++i) {
	  var l = tax_labels.get(i);
	  l.fillColor = Cesium.Color.LIGHTGREY;
	};

};

function toggleLabelColors() {};

function increaseLabelSize() {};

function decreaseLabelSize() {};

Sandcastle.addToolbarButton('Toggle Labels', function(){
    toggleLabels();
});


//// ADD DATA POINTS


var tax_points = viewer.scene.primitives.add(new Cesium.PointPrimitiveCollection() );
Sandcastle.addToolbarButton('Toggle SEED', function(){toggleSEED()});

var nr_points = viewer.scene.primitives.add(new Cesium.PrimitiveCollection() );
Sandcastle.addToolbarButton('Toggle NR', function(){toggleNR()});



function toggleSEED() {
	var len = tax_points.length;
	if(len > 0) {
		for (var i = 0; i < len; ++i) {
		  var l = tax_points.get(i);
		  l.show = !l.show;
		};
	} else {
		loadReference("tax_points.js");
	}
};


function toggleNR() {
	var len = nr_points.length;

	if(len > 0) {
		for(var i=0; i < 10; i++) {
			var NRdata = nr_points.get(i);
// 			NRdata.show = !NRdata.show;
			
		for (var j=0; j < NRdata.length; ++j) {
		  var l = NRdata.get(j);
		  l.show = !l.show;
		}
		};
	
	} else {
		for(var i=0; i < 10; i++) {
			var NRdata = new NRPrimDataSource();
			NRdata.loadUrl("./data/nr_all_" + i + ".json");
			nr_points.add(NRdata.entities);
		};
	}
}



