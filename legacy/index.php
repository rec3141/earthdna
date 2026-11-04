<head>
<script src="./../../Build/Cesium/Cesium.js"></script>
<script src="./../Sandcastle/Sandcastle-header.js"></script>
<style>
    @import url(./../../Build/Cesium/Widgets/widgets.css);
</style>

<style>
td {
	background-color: black;
	color: white;
}
td.highlighted {
    background-color: black !important;
    color: white !important;
}
#output, #cesiumContainer, #about {
	float: left;
	display: inline-block;
	vertical-align: top;
	overflow:auto;
}
#maincontainer {
    width:100%;
    height: 100%;
}
a {
	color:white
}

p {
    color:skyblue;
    padding: 5 5px;
}

ul,li {
	color:skyblue;
}

h1,h2,h3 {
    color:lightblue;
    padding: 5 5px;
}
</style>

<script>
  (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

  ga('create', 'UA-97387016-1', 'auto');
  ga('send', 'pageview');

</script>

<script src="require.js"></script>
<script src="md5.js"></script>
<script src="NRPrimDataSource.js"></script>
<script src="import.js"></script>
<script src="PhenotypeGeocoder.js"></script>
<script src="LabelCollectionGeocoder.js"></script>
<script src="biom.js"></script>
<script type="text/javascript">var Biom = require('biojs-io-biom').Biom;</script>

</head>

<body style="overflow:hidden;position:relative;">
<!-- <div id="output_holder" style="z-index:1;position:relative;width:40%;height:100%;background-color:transparent;"> ~~> -->
	<div id="output" style="z-index:1;position:absolute;left:20%;width:15%;height:100%;background-color:black;display:none;"></div>
<!-- </div> -->

<div id="maincontainer">
<div id="cesiumContainer" class="fullSize" style="width:80%;height:100%;float:right"></div>
<div id="toolbar" name="toolber" style="z-index:1;position:absolute;left:40%;bottom:60px;overflow:auto;"></div>
<!-- <div id="toolbar" name="toolber" style="z-index:1;position:absolute;left:40%;"></div>   -->
<!-- position:absolute;left:20%; -->

<script src="sedna.js"></script>

<div id="about" style="width:20%;height:100%;background-color:black;float:left;overflow-y:auto;">

<h1 style="text-align:center">Welcome to the <strong>SEDNA WorldTree</strong></h1>

<h2>About</h2>
<p>The WorldTree is based on the <a href="https://www.arb-silva.de/documentation/release-123/">SILVA SSU rRNA Database (Release 123)</a>. The locations of reference sequences were determined with MDS using the <a href="https://www.rdocumentation.org/packages/smacof">smacof</a> package in R and interpolated with <a href="https://github.com/DSC-SPIDAL/damds">DA-MDS</a> using phylogenetic distances calculated in <a href="http://www.mothur.org">mothur</a>. Elevation is scaled to the number of sequences in the SILVA reference tree at each taxonomic level. The globe is visualized with <a href="http://cesiumjs.org">Cesium</a>. Please direct questions and issues to the <a href="https://github.com/rec3141/sedna">Github repo</a>.</p>

<p>Click the globe once to stop it from rotating. Click twice to resume.</p>

<h2>View data as a community network on <a href="./../sedna_comm/">SEDNA Community View</a></h2>

<h2>Sample Display</h2>

<h3>Available Datasets:</h3>
<ul>
    <li>SILVA v123 reference datasets
    <ul>
    <li><a href="#" onclick='toggleSEED();'>SSU SEED (6,902 sequences)</a></li>
    <li><a href="#" onclick='toggleNR();'>SSU NR (597,607 sequences)</a></li>
    </ul>
    </li>

    <li>16S rRNA gene (Bacteria and Archaea)
    <ul>
    <li><a href="#" onclick='clickFile("data/chukchi16S.trim.contigs.good.unique.good.filter.unique.precluster.pick.full.wang.tax.txt")'>Chukchi Sea seawater (BPU)</a></li>
    <li><a href="#" onclick='clickFile("data/barents16S.trim.contigs.good.unique.good.filter.unique.precluster.pick.full.wang.tax.txt")'>Barents Sea seawater (BPU)</a></li> 
    <li><a href="#" onclick='clickFile("data/barrowiceald16S.trim.contigs.unique.good.filter.unique.seed_v119.wang.tax.txt")'>Barrow sea ice (ALD)</a></li> 
    <li><a href="#" onclick='clickFile("data/cmi_16S.trim.contigs.unique.good.filter.unique.seed_v119.wang.tax.txt")'>CMI oil in ice tank (KBD)</a></li> 
    <li><a href="#" onclick='clickFile("data/daneborg16S.trim.contigs.unique.good.filter.unique.seed_v119.wang.tax.txt")'>Daneborg sea ice (REC)</a></li> 
    </ul>
    </li>
    
    <li>18S rRNA gene (Eukaryotes)</li>
    <ul>
    <li><a href="#" onclick='clickFile("data/chukchi18S.trim.contigs.good.unique.good.filter.unique.precluster.pick.full.wang.tax.txt")'>Chukchi Sea seawater (BPU)</a></li> 
    <li><a href="#" onclick='clickFile("data/barents18S.trim.contigs.good.unique.good.filter.unique.precluster.pick.full.wang.tax.txt")'>Barents Sea seawater (BPU)</a></li> 
    <li><a href="#" onclick='clickFile("data/barrowice18S.trim.contigs.good.unique.good.filter.unique.precluster.full.wang.tax.txt")'>Barrow sea ice (BH)</a></li> 
    <li><a href="#" onclick='clickFile("data/barrowsed18S.trim.contigs.good.unique.good.filter.unique.precluster.full.wang.tax.txt")'>Barrow sediment (BH)</a></li> 
    <li><a href="#" onclick='clickFile("data/daneborg18S.trim.contigs.unique.good.filter.unique.seed_v119.wang.tax.txt")'>Daneborg sea ice (REC)</a></li> 
    <li><a href="#" onclick='clickFile("data/svalbardice18S.trim.contigs.good.unique.good.filter.unique.precluster.full.wang.tax.txt")'>Svalbard sea ice (REC)</a></li> 
    <li><a href="#" onclick='clickFile("data/cmi_18S.trim.contigs.unique.good.filter.unique.manual.wang.tax.txt")'>CMI oil in ice tank (KBD)</a></li> 

    </ul>
    </li>
</ul>

<p>To add your own samples, drag and drop a tax.summary file from mothur anywhere on the map. Ideally this file was processed using the same SILVA reference alignment as the WorldTree. This file is processed on your computer and is not uploaded. On the resulting page, click the names of the samples you want to compare. Relative abundance is shown scaled to the area of the points.</p>

<div id="download">
<p>To add the results from a krona analysis, drag the krona graph .html file onto the map, or paste the URL of the krona graph .html file below:</p>
<form method="post" action="index.php">
<input name="krona_url" value="http://reric.org/minion/Group1-Jackie-Jay/index.html" size="50" />
<input name="submit" type="submit" />
</form>


<?php

    // maximum execution time in seconds
//     set_time_limit (24 * 60 * 60);

if (isset($_POST['submit'])) {

    // folder to save downloaded files to. must end with slash
    $destination_folder = 'downloads/';

    $url = $_POST['krona_url'];
    $newfname = $destination_folder . md5($url);

    print_r($url);
    print_r($newfname);    

    $download = file_get_contents($url);
    preg_match_all('#<node.*?</count>#s', $download, $matches);
    $stack = array();
    foreach ($matches[0] as $element) {
           preg_match('/name=\"(?P<name>.+)\">/',$element,$name);
           preg_match('/<count><val>(?P<value>\d+)/',$element,$value);
           array_push($stack, implode("\t", array("na","na",$name["name"], "na",$value["value"])) . "\n");
           
    }
    unset($value);
    file_put_contents($newfname, $stack);

    echo '<script type="text/javascript">clickFile("' . $newfname . '");</script>';
}
?>
</div>

<h2>More Information</h2>
<p>Created by Eric Collins in the <a href="http://www.cryomics.org">Cryomics Lab</a> at the <a href="http://www.uaf.edu">University of Alaska Fairbanks</a>, with funding from the <a href="http://oceanexplorer.noaa.gov/">NOAA Office of Ocean Exploration and Research</a> for the project <a href="http://oceanexplorer.noaa.gov/explorations/15arctic-microbes/welcome.html">Mapping the Uncharted Diversity of Arctic Marine Microbes</a>.</p>

</div>



</div>

</body>
