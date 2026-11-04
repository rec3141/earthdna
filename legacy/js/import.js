
// okay next level shit here
// these need to process the input into visualizations that are EASY e.g. circles
// should be FAST because it's all client-side, don't need to upload the file
// should take input from... 
// 1) BLAST hits tabular output
// 2) .biom files (QIIME)
// 3) mothur .report files
// 4) mothur .tax.summary files
// 4) kraken .summary files
// 5) anything that krona can take -- use krona as translator

// start with #4 because that's the easiest
// results should have been classified with silva alignment
// some kind of kmer short-string classifier to automatically parse titles would be sweet

if (window.File && window.FileReader && window.FileList && window.Blob) {
  // Great success! All the File APIs are supported.
} else {
  alert('The File APIs are not fully supported in this browser.');
}

  function handleFileSelect(evt) {
    evt.stopPropagation();
    evt.preventDefault();

    var files = evt.dataTransfer.files; // FileList object.
    var f=files[0];
    
    // console.log(f);
    draggedfileReader(f);
    
    // files is a FileList of File objects. List some properties.
//    for (var i = 0, f; f = files[i]; i++) {

//only using the first file if multiple are selected


			// show selected file info
// 		  output.push('<li><strong>', escape(f.name), '</strong> (', f.type || 'n/a', ') - ',
// 					  f.size, ' bytes, last modified: ',
// 					  f.lastModifiedDate ? f.lastModifiedDate.toLocaleDateString() : 'n/a',
// 					  '</li>');
// 			document.getElementById('list').innerHTML = '<ul>' + output.join('') + '</ul>';

		// Only process text files.
		// this fails for some local files
// 		  if (!f.type.match('text.*')) {
// 			continue;
// 		  }

//     }	//end for each file
  }; // end handleFileSelect

  function handleDragOver(evt) {
    evt.stopPropagation();
    evt.preventDefault();
    evt.dataTransfer.dropEffect = 'copy'; // Explicitly show this is a copy.
  };



function clickFile(serverFile) {

        Cesium.loadText(serverFile, {
          'X-Custom-Header' : 'some value'
        }).then(function(text) {
//             console.log(text);
            importTaxSummaryFile(text);
        }).otherwise(function(error) {
            // an error occurred
        });
        
    };


function loadReference(serverFile) {
    // console.log(serverFile);
    Cesium.loadText(serverFile, {
    }).then(function(text) {
        eval(text);
    }).otherwise(function(error) {
        console.log("can't import "+serverFile);
    });

}

function importReference(fileObject) {
    
}


function draggedfileReader(f) {
    
		// console.log(f);

 		var reader = new FileReader();
		reader.onload = function(progressEvent){

        importTaxSummaryFile(progressEvent.target);
    
		};		 //end onload

	  reader.readAsText(f);
}




function makeTaxaIndex() {
    var index = {};
    var tax_len = tax_labels.length; // tax_labels is a global variable
    for(var j=0; j < tax_len; j++) {
        var k = tax_labels.get(j);
        // this could be dangerous, maybe do exact matching
        // if ( k.text.toLowerCase().indexOf( taxon.toLowerCase() ) > -1 ) {
        index[k.text] = k; // this will overwrite any duplicate taxanames
    }
    return index;						
}

function importTaxSummaryFile(fileObject) {
		// Entire file
    
        var inputtext = fileObject;
        if(Cesium.defined(fileObject.result)) { inputtext = fileObject.result };
        
		// Will parse this into a nice clickable table...
		// linking each taxon to it's latlon
		// and each table to its visibility	
		// console.log("start adding table");

    	// make taxon index

	    var taxaIndex = makeTaxaIndex();

        // logic to turn different file types into tables
        // default table has samples are in [1,] and taxa in [,3]
        
        
        // IMPORT FROM KRONA INDEX.HTML FILES

        // look for a krona file
        var kronamatch = inputtext.match(/<krona/);
        if (Cesium.defined(kronamatch)) {
        
            reg = new RegExp(/<node[\s\S]*?<\/count>/g); 

            var stack = [];
            var header = ["taxlevel","rankID","taxon","daughterlevels","total"];
            stack.push(header.join('\t'));
            var element;
            while((element = reg.exec(inputtext)) !== null) {
                var name = element[0].match(/name=\"([\s\S]+?)\">/);
                var value = element[0].match(/<count><val>(\d+)/);
                var arry = ["na","na",name[1],"na",value[1]];
                stack.push(arry.join('\t'));
            }
            inputtext = stack.join('\n');
 
        }
        
        
        // IMPORT FROM BIOM FILES
        
        var biommatch = inputtext.match(/biom-format/);
        
        if (Cesium.defined(biommatch)) {
            var jsontext = JSON.parse(inputtext);
            var biom = new Biom(jsontext); // "creates new biom object"
            biom.matrix_type = 'dense';
            console.log(biom);
//                 [Array(95), Array(95), ... ]

            var stack = [];
            var header = ["taxlevel","rankID","taxon","daughterlevels","total"];
            var biomsamples = [];
            for(var sample = 0; sample < biom.columns.length; sample++) {
                biomsamples.push(biom.columns[sample].id);
            }
            header = header.concat(biomsamples);
            stack.push(header.join('\t'));
            // console.log(biomsamples);
            
            // need to make Root
            var biomsum = [];
            for(var sample = 0; sample < biom.columns.length; sample++) {
                var sum = 0;
                for(var taxa = 0; taxa < biom.rows.length; taxa++) {
                    sum = sum + biom.data[taxa][sample];
                }
                biomsum.push(sum);                
            }
            
            var taxasum = 0;
            var arry = [depth,"na","root",7,taxasum];
            arry = arry.concat(biomsum);
            stack.push(arry.join('\t'));
            
            for(var taxa = 0; taxa < biom.rows.length; taxa++) {
                var name = biomsamples[taxa];
                var values = biom.data[taxa];
                var depth = biom.rows[taxa].metadata.taxonomy.length - 1;
                while(depth>0) {
                    var biomtaxa = biom.rows[taxa].metadata.taxonomy[depth];
                    var taxasplit = biomtaxa.split('__');
                    biomtaxa = taxasplit[1];
                    if(biomtaxa == "") {depth--} else {depth = 0};
                }
                var taxasum = 0;
                var arry = [depth,"na",biomtaxa,7-depth,taxasum];
                arry = arry.concat(values);
                stack.push(arry.join('\t'));

            }

//             biomtaxa = biom.getMetadata({dimension: 'rows', attribute: 'taxonomy'});


            inputtext = stack.join('\n');
 
        }

        // console.log(inputtext);
        
        
        // TRANSPOSE 
		// now we're going to transpose the table to make it easier to select samples
		var hlines = inputtext.split(/[\r\n]+/);
		transtable = new Array(hlines.length);
		for(var hline = 0; hline < hlines.length; hline++){
			var hcols = hlines[hline].split('\t');
			transtable[hline] = new Array(hcols.length);
			for(var hcol = 0; hcol < hcols.length; hcol++){
				transtable[hline][hcol] = hcols[hcol];
			}
		}
	    // console.log(transtable);
	    
		// transtable has been transposed so that samples are in [,1] and taxa are in [3,]
		
		
		// for now let's keep javascript and html separate for clarity
		// want to make a pointprimitive collection for each sample
		// later want to subdivide it by taxa, for now just use genus
		// scale so that abundance is proportional to sqrt area
		// tax_points.add( {position : Cesium.Cartesian3.fromDegrees(-100.909687068042,14.2350191313423), color : Cesium.Color.BLUE, pixelSize : 6, id : "AB000662.TchCifer", translucencyByDistance : new Cesium.NearFarScalar(4e6, 1.0, 2e7, 0.0)});
		// tax_points.add( {position : Cesium.Cartesian3.fromDegrees(46.4881233339125,2.23382593968493), color : Cesium.Color.BLUE, pixelSize : 6, id : "AY232255.StpTend2", translucencyByDistance : new Cesium.NearFarScalar(4e6, 1.0, 2e7, 0.0)});
		// PrimitiveCollection contains PrimitivePointCollections contains PrimitivePoints
		// PC_list -> PPC_member -> PP
		// .add(primitive) returns the primitive itself
		// not sure why nothing is plotting at the moment, probably a mixup in naming, should heirarchicalize the variable names
		// 
			var PC_list = window.PC_added;
			if(Cesium.defined(PC_list)) { PC_list.destroy(); }
			PC_list = new Cesium.PrimitiveCollection();
			PC_list.show = true;

			starthighlight = 4; //zero-based
	
				
// 			console.log(tax_len);
			var samples = [];
// 			var lines = this.result.split(/[\r\n]+/);
// 			var lines_length = lines.length;
			var max_taxa = transtable.length; // 1000; // max taxa to include
			for(var taxcol = 0; taxcol < max_taxa; taxcol++){
// 				document.getElementById('loadingOverlay').innerHTML = '<h1>' + Math.floor(line/lines.length*100) + '%</h1>';
// 				var cols = lines[line].split('\t');
				var taxon = "unknown";
				var matches = []; // get match from tax_labels, for now don't worry if there's no match
// 				var cols_length = cols.length // max samples to include
// 				cols_length = 10;
//				console.log(transtable[0]);
				for(var samplerow = 0; samplerow < transtable[0].length; samplerow++){
					if (taxcol == 0) {
						var PPC_member = new Cesium.PointPrimitiveCollection(); //one per sample
						samples.push(transtable[taxcol][samplerow]);
						PPC_member.name = transtable[taxcol][samplerow];
						PPC_member.samplecolors = Cesium.Color.fromRandom( {minimumRed : 0.25, minimumGreen : 0.25, minimumBlue : 0.25, alpha : 1} );
						PPC_member.show = false;
							if (samplerow == starthighlight) {PPC_member.show = true;}
						PC_list.add(PPC_member); // everyone gets a PPC even the junk columns, this maintains the mapping
					} else {


					/// match with the taxonomy in tax_labels
						if (samplerow == 2) {
							taxon = transtable[taxcol][samplerow];
								if ( taxaIndex[taxon]) {
									matches.push(taxaIndex[taxon]);
								}
                            // console.log(taxon);
							if (matches.length < 1) { break } // just take first exact match
						}
						
						if (samplerow < starthighlight) { continue }

// 						console.log(taxon);
// 						console.log(matches);
// 						if (matches.length > 1) {alert("duplicate matches detected"); console.log(matches)}
						firstmatch = matches[0]; // just take first exact match //but some taxon names are not unique!
						// this needs to be switched to NCBI ids
						// this adds the points for each taxa found

						var PPC_member = PC_list.get(samplerow);
						var pixelSize = 55 * Math.sqrt( transtable[taxcol][samplerow]/(transtable[1][samplerow]+1 ) );
						// oh this assumes that transtable[1] is 'root'
						// console.log(pixelSize);
						if (transtable[taxcol][samplerow] > 0) pixelSize = pixelSize+5;
						
						PPC_member.add(
							{
								description : 
									'<h2>Sample: ' + PPC_member.name + '</h2>' +
									'<h2>Taxon: ' + taxon + '</h2>' +
									'<h2>Count: ' + transtable[taxcol][samplerow] + '</h2>',
								position : firstmatch.position,
								color : PPC_member.samplecolors,
								outlineColor : Cesium.Color.BLACK,
								pixelSize : pixelSize,
								id : {taxon: firstmatch.text, abundance: transtable[taxcol][samplerow]},
//									scaleByDistance : new Cesium.NearFarScalar(1, 100, 2e7, 0.0),
								translucencyByDistance : new Cesium.NearFarScalar(4e4, 1.0, 1e8, 0.0),
								show : PPC_member.show
							}
						);
// 							console.log(PPC_member);
					}
				} // end for mrow
				// console.log(PC_list);
			} //end for mcol

// 			console.log(PC_save);
// 			console.log(PC_list[PC_save[5]]);
//			var PC_here = PC_list.get(5);
//			PC_here.show = true;
			
 			// console.log("done processing collections");

			window.PC_added = viewer.scene.primitives.add(PC_list);
// 			for(var l=0; l < PC_list.length; ++l) {
// 				var PC_added = viewer.scene.primitives.add(PC_list.get(l));
// 			}
			

//			output HTML table
// 			console.log(transtable);
			var output = [];
			output.push('<span id="close" style="color:white;" onclick="' + "this.parentNode.style.setProperty('display','none');" + ' return false;">[close]</span><br>');  //this.parentNode.parentNode.removeChild(this.parentNode)
			output.push('<table id="taxatable" style="border-spacing: 10px;border-collapse: separate;">');
			for(var hcol = 0; hcol < transtable[0].length; hcol++){
				output.push('<tr>');
				for(var hrow = 0; hrow < transtable.length; hrow++){
					var hclass = 'class="highlighted"';
					var hstyle = "";
					if (hrow == 0) {
						var PC_member = PC_list.get(hcol);
						var PC_color = PC_member.samplecolors;
						PC_color = Cesium.Color.fromAlpha(PC_color,1).toCssColorString();
						var hstyle = 'style="background-color: ' + PC_color  + '; color: black;"';
						// take away the highlight to reveal the inline color
						if(hcol == starthighlight) {
							hclass = '';
						}
					}
						output.push(
						'<td ' + hclass + ' ' + hstyle + '>' + escape(transtable[hrow][hcol]) + '</td>'
					);
				}
				output.push('</tr>');
			}
			output.push('</table>');
			
		  	document.getElementById('output').innerHTML = output.join('');
		  	document.getElementById('output').style.setProperty('display','inline-block');

 			// console.log("adding clickability");

//			adding clickability
			var tbl = document.getElementById("taxatable");
			if (tbl != null) {
 				for (var trow = 0; trow < tbl.rows.length; trow++) { 
					var tcol = 0;
					if (trow < starthighlight) { continue }	
					tbl.rows[trow].cells[tcol].onclick = function (funrow) {
						return function () {
							this.classList.toggle("highlighted");
							var PPC_member = PC_list.get(funrow);
							var color = Cesium.Color.fromRandom( {minimumRed : 0.25, minimumGreen : 0.25, minimumBlue : 0.25, alpha : 1} );
							console.log(color);
							PC_member.samplecolors = color;
							this.style.backgroundColor = Cesium.Color.fromAlpha(color,1).toCssColorString();
							for (var point = 0; point < PPC_member.length; point++) {
								var PPC_point = PPC_member.get(point);
								PPC_point.show = !PPC_point.show;
								PPC_point.color = color;	
							}
						};
					}(trow);				}
			}
// console.log("added clickability");
}






  // Setup the dnd listeners.
  var dropZone = document;//.getElementById('body');
  dropZone.addEventListener('dragover', handleDragOver, false);
  dropZone.addEventListener('drop', handleFileSelect, false);
  



