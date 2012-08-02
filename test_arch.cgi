#!/usr/bin/env perl

use strict;
use warnings;
use CGI qw(:standard);
use JSON::XS;
use CGI::Carp qw/fatalsToBrowser/;
my $cgi = CGI->new;
print $cgi->header;




print <<ENDHTML;
<html>
  <head>
    <title>archpic test page</title>
    <script type="text/javascript" src="../javascript/jquery.js"></script>
    <style type="text/css">


           html, body { height: 100%; width: 100%; padding: 0; margin: 0; }
           note { width: 20%; padding: 0; margin-left: 80%; float:left; }
    </style>
            <!-- JSON support for IE (needed to use JS API) -->
        <script type="text/javascript" src="../javascript/js/min/json2.min.js"></script>
        

        
        <!-- Cytoscape Web JS API (needed to reference org.cytoscapeweb.Visualization) -->
            <script type="text/javascript" src="../javascript/jquery.js"></script>
	    
	            
  </head>
  <body>

<h2> This is a test page for the archpic.cgi domain plotter script </h2>

		<script type="text/javascript">
            window.onload=function() {
			var data;
var json = '{"12345":{"length":"1205","names":[["hs","ENSP00000449791"]],"structures":{"strong":[["11111","801-830,940-1164","4.03e-56","A Superfamily","rgb(84,242,87)"],["11111","80-83,94-116","4.03e-56","A Superfamily","rgb(84,242,87)"]]}}}';
     var data;

if ( document.implementation.hasFeature("http://www.w3.org/TR/SVG11/feature#BasicStructure", "1.1")) {
	\$.ajax({url: 'archpic.cgi', data: {doms:json} ,async: true ,type: 'POST', dataType:'xml', success: function(reply,text) { \$('#image_holder').append(\$(reply).find('svg')); }});
}else{
	\$.ajax({url: 'archpic.cgi', data: {doms:json,png:1} ,async: true ,type: 'POST',dataType:'text', success: function(reply,text) { \$('#image_holder').append('<img src="data:image/png;base64,' + reply + '" />');}});
}
                    
 
            
};

        </script>
        
        <div id="image_holder" style="overflow-x:scroll" >
		
		</div>

       
	

</body>
	</html>






ENDHTML
