#!/usr/bin/env perl

use warnings;
use strict;

=head1 NAME

B<archpic.cgi> - Display the architecture for a specified SUPERFAMILY protein.

=head1 DESCRIPTION

Outputs an SVG rendering of the given proteins structual Weaker hits are included with their e-values specified as 'hanging' blocks. This a hack on an original script by Matt Oates

=head1 AUTHOR
B<Matt Oates> - I<Matt.Oates@bristol.ac.uk>
B<Owen Rackham> - I<Owen.Rackham@bristol.ac.uk>

=head1 NOTICE

B<Matt Oates> (Jan 2012) First features added.

=head1 LICENSE AND COPYRIGHT

B<Copyright 2012 Owen Rackham/Matt Oates>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

=head1 FUNCTIONS

=over 4

=cut

use POSIX qw/ceil floor/;
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
use Data::Dumper;
use DBI;
use JSON::XS;
use URI::Escape;
use MIME::Base64 qw(encode_base64);

#Deal with the CGI parameters here
my $cgi = CGI->new;

#Specify a json string
my $json_string;
$json_string = (defined $cgi->param('doms'))?$cgi->param('doms'):undef;
#$json_string = '{"cath101":{"length":"1205","names":[["hs","ENSP00000449791"]],"structures":{"strong":[["11111","801-830,940-1164","4.03e-56","DEFAULT","rgb(84,242,87)"],["11111","80-83,94-116","4.03e-56","gfsgdf","rgb(84,242,87)"]]}}}';


#What colouring scheme should we use (deprecated)
my $colouring = (defined $cgi->param('colouring'))?$cgi->param('colouring'):undef;

#Don't include structure!
my $nostruct = (defined $cgi->param('nostruct'))?$cgi->param('nostruct'):undef;

#What should we render?
my $draw_binding = (defined $cgi->param('binding'))?$cgi->param('binding'):0;
my $draw_consensus = (defined $cgi->param('consensus'))?$cgi->param('consensus'):0;
my $draw_conflicts = (defined $cgi->param('conflicts'))?$cgi->param('conflicts'):0;

#Is this a for-print rendering or for the web
my $forprint = (defined $cgi->param('forprint'))?1:0;

#Do you want textual labels for a key?
my $draw_labels = (defined $cgi->param('labels'))?$cgi->param('labels'):1;
my $supfam_labels = (defined $cgi->param('supfam'))?$cgi->param('supfam'):0;

#If we want SUPERFAMILY labels then we dont want to show weak hits!
my $include_weak_hits = (defined $cgi->param('weak'))?$cgi->param('weak'):0;
$include_weak_hits = ($supfam_labels)?0:$include_weak_hits;

#Do you want callout boxes around those labels, looke more web 2.0!
my $draw_callouts = (defined $cgi->param('callouts'))?$cgi->param('callouts'):1;
$draw_labels = 1 if ($draw_callouts); #Callouts imply labels but not vice versa

#Do you want a key at the end of the figure
my $draw_key = (defined $cgi->param('key'))?$cgi->param('key'):1;

#Should we have a protein co-ordinate ruler included
my $draw_ruler = (defined $cgi->param('ruler'))?$cgi->param('ruler'):1;

#Force conversion to PNG, currently based on browser support for SVG
my $force_png = (defined $cgi->param('png'))?$cgi->param('png'):undef;


#Global variables
my %sf_details;
my %predictor_details;
my %protein_details;
my $output = '';

my $current_record = 0;
my $ruler_width = 0;


#Formatting config variables
my $min_width = 400; #Deal with short proteins having their SF labels cut short
my $weak_below = 15; #How far below do you want weaker sf hits, negative is above
my $weak_note = ' (weak support)'; #What should be appended to the structure popup title if its a weaker hit
my $xpad = 5; #left indentation
my $ypad = 5; #top indentation
my $record_height = 70; #how much spacing for a single protein record
my $record_spacing = 10;
my $dy = $record_height/2; #with the amino acid chain in the middle to begin with, this is a cursor for the vertical position reused in each draw function
my $other_height = 10; #Height of othered blocks around the amino line
my $structure_height = 8; #Height of each structural block around the amino line


#JSON hack -- this should be in a sub

%protein_details = %{decode_json($json_string)};
my $label = 1;
foreach my $id (keys %protein_details){
	foreach my $hit (keys %{$protein_details{$id}{'structures'}}){
		foreach my $array (@{$protein_details{$id}{'structures'}{$hit}}){
			my @array = @{$array};
			unless(exists($sf_details{$array[0]})){
				$sf_details{$array[0]}{'colour'} = $array[4];
				$sf_details{$array[0]}{'label'} = $label++;
				$sf_details{$array[0]}{'description'} = $array[3];
			}
		}
	}
}

=item B<error>

	Print out a CGI error message with HTTP status 204 No Content
=cut
sub error  {
	my ($error) = @_;
	print $cgi->header('type'=>'text/html', 'status'=>'204 No Content');
	print "<html><body><h1>$error</h1></body></html>";
	exit;
}

=item B<range(> {from => $from, to => $to, by => $by, with => \@with} B<)> or alternatively B<range($from, $to, $by, \@with)>

	Function to create a range of values, either from the series specified or that series sliced over an array.
	Returns all the sliced values in range with the specified list if given.
	Returns a list of integers or reals depending on the 'by' clause if no list is given to slice with.

	Using the default values you can do something like:
		my @list = range({to=>100}); which will act in the same way as 0..100
=cut
sub range {
	my %range;
	my ($from, $to, $by, $with) = @_;
	
	#Allow either positional...
	if (ref $from eq 'HASH') {
		%range = %{$from};
	}
	#or hash based parameters
	else {
		%range = ( from => $from, to => $to, by => $by, with => $with );
	}
	
	#Set some defaults for ease of use, e.b.  range({to=>100}) or range(undef,100)
	$range{'from'} = 0 unless defined $range{'from'};
	$range{'by'} = 1 unless defined $range{'by'};
	
	#Make sure the required from and by parameters are set, should be because of defaults, someone might have passed in undef!
	defined $range{'from'} or die 'Required parameter "from" not defined.';
	defined $range{'by'} or die 'Required parameter "by" not defined.';
	
	#Do some checks if we are slicing with an array
	if (defined $range{'with'}) {
		
		#Make sure the reference is really an array! die otherwise
		ref $range{'with'} eq 'ARRAY' or die 'Parameter "with" does not appear to be an array ref to slice with.';
		
		#If we are slicing to a given value make sure its in range of the array index we are slicing with
		if (defined $range{'to'}) {
			'to' < scalar @{$range{'with'}} or die 'Required parameter "to" is beyond the dimensions of the list you are slicing "with".';
		}
		
		#Set the target of the range to the length of the array we are slicing with
		else {
			$range{'to'} = scalar @{$range{'with'}} - 1;
		}
	}
	
	#Make sure the required to parameter is set
	defined $range{'to'} or die 'Required parameter "to" not defined.';
	
	#Return value
	my @list = ();
	
	#If we are slicing on the array do that
	if (defined $range{'with'}) {
		warn "Non integer 'by' term used 'with' array slice." unless $range{'by'} =~ /D/;
		for (my $i = $range{'from'}; $i <= $range{'to'}; $i += $range{'by'}) {push @list, $range{'with'}[$i]; }
	} 
	
	#Otherwise create a numerical range
	else {
		for (my $i = $range{'from'}; $i <= $range{'to'}; $i += $range{'by'}) {push @list, $i; }
	}
	
	return @list;
}

=item B<hsv2rgb($Hue, $Saturation, $Value)>

	Function to convert HSV colour space values to RGB colour space.
	Returns RGB value as [R,G,B]
=cut
sub hsv2rgb {
	my ($Hue,$Saturation,$Value) = @_;
	my ($Red,$Green,$Blue) = (0,0,0);
	
	#Check the input and warn if it's a bit wrong
	warn "Invalid Hue component of HSV colour passed, with value: $Hue." unless ($Hue >= 0.0 and $Hue <= 360.0);
	warn "Invalid Saturation component of HSV colour passed, with value: $Saturation." unless($Saturation >= 0.0 and $Saturation <= 1.0);
	warn "Invalid Value component of HSV colour passed, with value: $Value." unless ($Value >= 0.0 and $Value <= 1.0);
	
	#If colour has no saturation return greyscale RGB
	if ($Saturation == 0) {
		$Red = $Green = $Blue = $Value;
		return [$Red, $Green, $Blue];
	}
	
	#Partition the Hue into the 5 different colour chroma and then map each of these to RGB based on the colour theory
	$Hue /= 60.0;
	my $Chroma = floor($Hue) % 6; 
	my $H_d = $Hue - $Chroma; 
	
	#RGB cube components
	my ($I,$J,$K) = ( $Value * ( 1 - $Saturation ),
	                            $Value * ( 1 - $Saturation * $H_d ),
				    $Value * ( 1 - $Saturation * ( 1 - $H_d ) )
                                    );
	
	#Map components to RGB values per chroma
	if ($Chroma == 0) { ($Red,$Green,$Blue) = ($Value,$K,$I); }
	elsif ($Chroma == 1) { ($Red,$Green,$Blue) = ($J,$Value,$I); }
	elsif ($Chroma == 2) { ($Red,$Green,$Blue) = ($I,$Value,$K); }
	elsif ($Chroma == 3) { ($Red,$Green,$Blue) = ($I,$J,$Value); }
	elsif ($Chroma == 4) { ($Red,$Green,$Blue) = ($K,$I,$Value); }
	else{ ($Red,$Green,$Blue) = ($Value,$I,$J); }
	
	#Return the RGB value in the integer range [0,255] rather than real [0,1]
	return [floor($Red * 255),floor($Green * 255),floor($Blue * 255)];
}

=item B<colourset($num_colours, $method)>

	Function to grab a set of well spaced colours from an HSV colour wheel.
		$num_colours - The number of colour values to produce, must be greater than 0 but no bigger than 360
		$method - The method for selecting colours over HSV colour space, either 'equal_spacing' or for around 10 colours 'chroma_bisection' is better.
	Returns an array of RGB values of the form ([R,G,B]) and undef on $num_colours out of bounds
=cut
sub colourset {
	my ($num_colours,$method) = @_;
	if ($num_colours <= 0 or $num_colours > 360) {
		warn "Number of colours requested out of bounds.";
		return undef;
	}
	$method = 'chroma_bisection' unless $method;
	
	#Internal sub to randomly shuffle an array
	sub fisher_yates_shuffle {
		my ($array) = @_;
		my $current;
		for ($current = @$array; --$current; ) {
			my $selected = int rand ($current+1);
			next if $current == $selected;
			#Reverse the slice between current position and the randomly selected
			@$array[$current,$selected] = @$array[$selected,$current];
		}
		return $array;
	}
	
	#Colours to return
	my %colours;
	
	#Default Hue Saturation and Value, saturation of 0.65 gives a more pastel feel!
	my ($Hue, $Saturation, $Value) = (0.0,0.65,0.95);
	
	#The interval to space colours around the wheel if equal
	my $hsv_interval = 360 / $num_colours;
	
	#Array of degrees for reuse to create ranged arrays with a given interval
	my @degrees = 1..360;
	
	#Iteratively bisect each chroma segment so that the first 6 colours are well spaced perceptually.
	#However after 12 colours we will have increasing pairs that are more confused as 
	#they are increasingly close to each other compared to the rest of the colours!
	#To get around this problem of selecting closely around a single bisection, we jump around the 
	#chroma randomly sampling.
	if ($method eq 'chroma_bisection') {
		#The current cycle of chroma bisection
		my $hsv_cycle = 1;
		#Number of colours selected by bisecting chroma so far
		my $colours_selected = 0;
		
		#While we still have colours to select
		while ($colours_selected != $num_colours) {
			#Work out the size of interval to use this cycle around the wheel
			$hsv_interval = 60 / $hsv_cycle;
			#Get all the Hues for this cycle that haven't already been examined and are on the line of bisection
			my @Hues = grep { not exists $colours{$_%360} } range({ from => 1, to => 360, by => $hsv_interval});
			#Shuffle so that we don't take from around the same chroma all the time, only perceptually worthwhile after 12th colour
			fisher_yates_shuffle(\@Hues) if $hsv_cycle > 2;
			
			#While we still have hues to select from in this cycle
			while (@Hues) {
				#Finish if we have enough colours
				last if $colours_selected == $num_colours;
				#Consume a Hue from this cycle
				$Hue = shift @Hues;
				#360 should be 0 for red
				$Hue %= 360;
				#Store this Hue and mark selection
				$colours{$Hue} = hsv2rgb($Hue,$Saturation,$Value) ;
				$colours_selected++;
			}
			$hsv_cycle++;
		}
	}
	
	#Just space colours even distances apart over the HSV colour wheel.
	#You have slightly odd/garish colours coming out, but you dont get uneven perceptual distance
	#between pairs of colours. This scales far better despite the horrible colours.
	elsif ($method eq 'equal_spacing') {	
		foreach $Hue (1..$num_colours) {
			$Hue = ($Hue * $hsv_interval) % 360;
			$colours{$Hue} = hsv2rgb($Hue,$Saturation,$Value) ;
		}
	}
	
	#Otherwise return nothing and warn the programmer
	else {
		warn "Colourset method not known, use either 'equal_spacing' or for fewer colours 'chroma_bisection'";
		return undef;
	}
	
	#Shuffle final colours so that even if we do use chroma_bisection closer colours will hopefully not be sequential
	@_ = values %colours;
	fisher_yates_shuffle(\@_);
	return @_;
}




=item B<header>

	Output the SVG header with some JavaScript for performing mouse popup information
=cut
sub header {
	my ($width, $height) = @_;

	return <<EOF
<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"
"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">	
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:xlink="http://www.w3.org/1999/xlink"
     width="$width" height="$height"
     viewBox="0 0 $width $height"
     preserveAspectRatio="xMidYMid meet"
     onload="init(evt)">

     <script type="text/ecmascript"><![CDATA[
      var svg_document, svg_root;
      var tooltop, tip_box, tip_title, tip_text;
      var hit_range, hit_quality;

      function init(evt) {
	 svg_document = evt.target.ownerDocument;
	 svg_root = document.documentElement;

	 tooltip = svg_document.getElementById('tooltip');
	 tip_box = svg_document.getElementById('tip_box');
	 tip_text = svg_document.getElementById('tip_text');
	 tip_title = svg_document.getElementById('tip_title');
	 tip_desc = svg_document.getElementById('tip_desc');
	 hit_range = svg_document.getElementById('hit_range');
	 hit_quality = svg_document.getElementById('hit_quality');

      };

      function show_tip(evt) {
	  var x = evt.clientX + 1;
	  var y = evt.clientY + 2;
	  var hit = evt.target;

	  var title_text = hit.getElementsByTagName('name').item(0);
	  if (! title_text ) {
	     return false;
	  }

	  title_text = title_text.firstChild.nodeValue;
	  tip_title.firstChild.nodeValue = title_text;
	  tip_title.setAttributeNS(null, 'display', 'inline' );

	  var description_text = hit.getElementsByTagName('desc').item(0);
	  if (description_text) {
	     description_text = description_text.firstChild.nodeValue;
	     tip_desc.firstChild.nodeValue = description_text;
	     tip_desc.setAttributeNS(null, 'display', 'inline' );
	  }
	  else {
	     tip_desc.setAttributeNS(null, 'display', 'none' );
	  }

	  var range_text = hit.getElementsByTagName('range').item(0);
	  if (range_text) {
	     range_text = range_text.firstChild.nodeValue;
	     hit_range.firstChild.nodeValue = range_text;
	     hit_range.setAttributeNS(null, 'display', 'inline' );
	  }
	  else {
	     hit_range.setAttributeNS(null, 'display', 'none' );
	  }

	  var quality_text = hit.getElementsByTagName('quality').item(0);
	  if (quality_text) {
	     quality_text = quality_text.firstChild.nodeValue;
	     hit_quality.firstChild.nodeValue = quality_text;
	     hit_quality.setAttributeNS(null, 'display', 'inline' );
	  }
	  else {
	     hit_quality.setAttributeNS(null, 'display', 'none' );
	  }

	  tip_title.firstChild.nodeValue = title_text;
	  tip_desc.firstChild.nodeValue = description_text;
	  hit_range.firstChild.nodeValue = range_text;
	  hit_quality.firstChild.nodeValue = quality_text;
	  
	  var box = tip_text.getBBox();
	  tip_box.setAttributeNS(null, 'width', Number(box.width) + 10);
	  tip_box.setAttributeNS(null, 'height', Number(box.height) + 10);
	  
	  tooltip.setAttributeNS(null, 'transform', 'translate(' + x + ',' + y + ')');
	  tooltip.setAttributeNS(null, 'visibility', 'visible');
      };


      function hide_tip(evt) {
	  tooltip.setAttributeNS(null, 'visibility', 'hidden');
      };

   ]]></script>

  <defs>
    <linearGradient id="other"
		    x1="0%" y1="0%"
		    x2="0%" y2="100%"
		    spreadMethod="pad">
      <stop offset="0%"   stop-color="#cccccc" stop-opacity="1"/>
      <stop offset="100%" stop-color="#666666" stop-opacity="1"/>
    </linearGradient>
  </defs>
  <g transform="translate($xpad,$ypad)">
EOF
;
}

=item B<footer>

	Print out the SVG footer containing the tooltip popup XML
=cut
sub footer {
	my $result = <<EOF
   </g>
EOF
	;
	$result .= <<EOF
   <g id='tooltip' opacity='0.8' visibility='hidden' pointer-events='none'>
      <rect id='tip_box' x='0' y='5' width='88' height='20' rx='2' ry='2' fill='white' stroke='black'/>
      <text id='tip_text' x='5' y='20' font-family='Arial' font-size='10'>
        <tspan id='tip_title' x='5' font-weight='bold' text-decoration="underline"><![CDATA[]]></tspan>
	    <tspan id='hit_range' x='5' dy='15' fill='red'><![CDATA[]]></tspan>
        <tspan id='hit_quality' x='5' dy='10' fill='green'><![CDATA[]]></tspan>
	    <tspan id='tip_desc' x='5' dy='10' fill='blue'><![CDATA[]]></tspan>
      </text>
   </g>
EOF
	unless $force_png or $forprint;
	$result .= <<EOF
</svg>
EOF
;
	return $result;
}

=item B<callout($x,$y,$label)>

	Returns the string for drawing a callout at a given position filled with $label.
=cut
sub callout {
	my ($x,$y,$label) = @_;
	my $width = length $label;
	return "<path style=\"fill:#ffffff;fill-opacity:0.7;stroke:#000000;stroke-width:0.5;stroke-linejoin:round;stroke-opacity:1\" d=\"M $x,$y L "
					.($x+2.5).",".($y+2.5)
					." L ".($x+(5*$width)).",".($y+2.5)
					." L ".($x+(5*$width)).",".($y+17)
					." L ".($x-(5*$width)).",".($y+17)
					." L ".($x-(5*$width)).",".($y+2.5)
					." L ".($x-2.5).",".($y+2.5)
				." Z\" />";
}



=item B<draw_consensus($protein)>

	Draw all the other hits for a protein
=cut
sub draw_consensus {
	my ($protein) = @_;
    return unless defined $protein_details{$protein}{'consensus'};
	my $result = '';
	my $middle = $dy + ($record_height / 2);
    my $y = $middle + (2 * $weak_below) + $structure_height + $other_height + 5; 
	
    my $num_predictors = scalar grep {$predictor_details{$_}{type} eq 'other' } keys %predictor_details;

	#Draw the consensus strip
	foreach my $dx (0..@{$protein_details{$protein}{'consensus'}}-1) {
        next if $protein_details{$protein}{'consensus'}[$dx] == 0;
        next if $draw_conflicts and $protein_details{$protein}{'conflict'}[$dx] > 0;
        my $value = 0;
        $value = $num_predictors / $protein_details{$protein}{'consensus'}[$dx] if $protein_details{$protein}{'consensus'}[$dx] != 0;
        my ($r,$g,$b) = ($value == 0)?(255,255,255):@{hsv2rgb(120.0,0.8,$value)};
		$result .= <<EOF
	<rect x="$dx" y="$y" width="1" height="$other_height" style="fill: rgb($r,$g,$b); stroke:none; stroke-width:0.0"></rect>
EOF
		;

	}
    $result .= "<text x=\"0\" y=\"".($y-2)."\" font-size=\"10\">Consensus</text>";
	return $result;
}

=item B<draw_conflicts($protein)>

	Draw all the other hits for a protein
=cut
sub draw_conflicts {
	my ($protein) = @_;
    return unless defined $protein_details{$protein}{'conflict'};
	my $result = '';
	my $middle = $dy + ($record_height / 2) + $other_height + 5;
    my $y = $middle + (2*$weak_below) + $structure_height; 
	
    my $num_predictors = scalar grep {$predictor_details{$_}{type} eq 'other' } keys %predictor_details;

	#Draw the consensus strip
	foreach my $dx (0..@{$protein_details{$protein}{'conflict'}}-1) {
        next if $protein_details{$protein}{'conflict'}[$dx] == 0;
        my $value = 0;
        $value = $num_predictors / $protein_details{$protein}{'conflict'}[$dx];
        my ($r,$g,$b) = ($value == 0)?(255,255,255):@{hsv2rgb(5.0,0.8,$value)};
		$result .= <<EOF
	<rect x="$dx" y="$y" width="1" height="$other_height" style="fill: rgb($r,$g,$b); stroke:none; stroke-width:0.0"></rect>
EOF
		;

	}
	$result .= "<rect x=\"0\" y=\"".($y)."\" width=\"$protein_details{$protein}{length}\" height=\"".($other_height)."\" style=\"fill: none; stroke:black; stroke-width:0.5;\"></rect>";
	return $result;
}

=item B<draw_structures($protein, $is_weak)>

	Draw all the structural hits for a protein, below the amino line if they are weak
=cut
sub draw_structures {
	my ($protein, $is_weak) = @_;
	my @structures;
	my $result = '';
	
	if ($is_weak) {
		@structures = @{$protein_details{$protein}{'structures'}{'weak'}};
	} else {
		@structures = @{$protein_details{$protein}{'structures'}{'strong'}};
	}
	
	foreach my $structure (@structures) {
	my @ranges = map {[split /-/, $_]} split /,/, $structure->[1];
	my $range = $structure->[1];
	my ($start, $end) = ($ranges[0][0],$ranges[0][1]);
	my $evalue = $structure->[2];
	my $link = "http://supfam.org/SUPERFAMILY/cgi-bin/scop.cgi?sunid=".$structure->[0];
	my $name = $structure->[3];
	$name .= ($is_weak)?$weak_note:'';
	my $fill = $sf_details{$structure->[0]}{'colour'};
	$start--;
	my $width = $end - $start;
	my $middle = $dy + ($record_height / 2);
	
	#Horizontally place the text box/callout
	my $tx = $start + ($end - $start)/2;
	
	#Vertically place the structure box
	my $y = $middle;
	$y += ($is_weak)?$weak_below:0; #Deal with weaker hits below the line stronger around the amino line
	$y -= $structure_height/2;
	
	my $drop_middle = $y + ($structure_height/2);
	
	#Vertically place the structure text box/callout
	my $ty = $y + $structure_height + 12;
	
	#Setup the callout arrow boxes
	my $label = $sf_details{$structure->[0]}{'label'};
	my $callout = ($draw_callouts)?callout($tx,$ty - 14,$label):'';
	$label = ($draw_labels)?"<text x=\"$tx\" y=\"$ty\" text-anchor=\"middle\" style=\"font-size:10px\">".$label."</text>":'';
				
	$result .= <<EOF
	<g>
EOF
;
	#Draw all of the blocks of structure for this hit, usually this is one but it can be more for inserted domains
	foreach my $part (@ranges) {
		($start, $end) = @{$part};
		$start--;
		$width = $end - $start;
		
		#Print drop down lines if this is a weaker hit
		if ($is_weak) {
			$result .= <<EOF
		<line x1="$start" y1="$middle" x2="$start" y2="$drop_middle" 
		   style="fill: none; stroke-dasharray:1,1; stroke: #000; stroke-width:1" />
		<line x1="$end" y1="$middle" x2="$end" y2="$drop_middle" 
		   style="fill: none; stroke-dasharray:1,1; stroke: #000; stroke-width:1" />
EOF
			;
		}
		
		#Print the structure box
		$result .= <<EOF
		<a xlink:href="$link" xlink:show="new" target="_blank">

			<rect x="$start" y="$y" width="$width" height="$structure_height" style="fill:$fill; opacity: 0.85;" rx="2" ry="2" onmouseover="show_tip(evt)" onmouseout="hide_tip(evt)">
				<name>$name</name>
				<desc>Length: $width</desc>
				<range>Range: $range</range>
				<quality>E-value: $evalue</quality>
			</rect>
			$callout
			$label
		</a>
EOF
		;
	}
	
	#Draw any connectors for inserted domain hits if we need to
	my @connector_ranges = map {split /-/, $_} split /,/, $structure->[1];
	shift @connector_ranges;
	pop @connector_ranges;
	while (my @connector = splice(@connector_ranges, 0, 2)) {
		my ($x1, $x2) = @connector;
		if ($is_weak) {
			$result .= <<EOF
			<line x1="$x1" y1="$drop_middle" x2="$x2" y2="$drop_middle" style="stroke:$fill; opacity: 0.85; stroke-width:2" />
EOF
			;
		}
		else {
			$result .= <<EOF
			<line x1="$x1" y1="$middle" x2="$x2" y2="$middle" style="stroke:$fill; opacity: 0.85; stroke-width:2" />
EOF
			;
		}
	}
	
	
	$result .= <<EOF
	</g>
EOF
	;

	}
	return $result;
}

=item B<draw_protein($protein)>

	Draws the amino line and labels for a protein
=cut
sub draw_protein {
	my ($protein) = @_;
	my $length = $protein_details{$protein}{'length'};
	$ruler_width = $length if $length > $ruler_width;
	my $names = $protein_details{$protein}{'names'};
	my @names;
	my $result = '';

	$dy = ($current_record * $record_height) + ($current_record * $record_spacing);
	my $y = $dy + 10;
	my $top = $dy;
	my $middle = $dy + ($record_height / 2);
	my $bottom = $dy + $record_height;

	#Draw the protein labels
	foreach my $seqid (@$names) {
		push @names, '<a xlink:show="new" target="_blank" xlink:href="http://supfam.org/SUPERFAMILY/cgi-bin/gene.cgi?genome='.$seqid->[0].';seqid='.$seqid->[1].'">'.$seqid->[1].'</a>';
	}
	my $name = join ', ', @names;
	$result .= "<text x=\"0\" y=\"$y\" font-family=\"Arial\" font-size=\"10\">$name</text>";
	
	#Draw protein amino line
	$y = $middle -1;
	$result .= <<EOF
	  <rect x="0" y="$y" width="$length" height="2" rx="0" ry="0"
	     style="fill: #333;" />
EOF
	;


    	$result .= draw_structures($protein);
	    $result .= draw_structures($protein, 'weaker') if $include_weak_hits;
    
	return $result;
}

=item B<amino_ruler>

	Draw a ruler showing to the longest length amino line drawn by draw_protein (this is stored in the $ruler_width global variable).
	50 amino acid intervals are demarked with the last tick being the max length value.
=cut
sub amino_ruler {
	my $num_ticks = floor($ruler_width/50);
	my $tick_interval = floor( ($ruler_width/$num_ticks) / 10 ) * 10;
	

	my $result = <<EOF
		  <line x1="0" y1="$dy" x2="$ruler_width" y2="$dy" style="stroke: #333; stroke-width: 1;" />
EOF
		;
	foreach my $tick ( range(0,$ruler_width,$tick_interval) ) {
		my $y = $dy + 3;
		$result .= "	  <line x1=\"$tick\" y1=\"$dy\" x2=\"$tick\" y2=\"$y\" style=\"stroke: #333; stroke-width: 0.5;\" />";
		$y += 10;
		$result .= "	  <text x=\"$tick\" y=\"$y\" text-anchor=\"middle\" style=\"font-size:10px\">$tick</text>";
	}

	if ($ruler_width % $tick_interval > $tick_interval/2) {
		my $y = $dy + 3;
		$result .= "	  <line x1=\"$ruler_width\" y1=\"$dy\" x2=\"$ruler_width\" y2=\"$y\" style=\"stroke: #333; stroke-width: 0.5;\" />";
		$y += 10;
		$result .= "	  <text x=\"$ruler_width\" y=\"$y\" text-anchor=\"middle\" style=\"font-size:10px\">$ruler_width</text>";
	}
	return $result;
}

=item B<sf_key>

	Draw the key for the superfamilies found in the input set of proteins indexed with both colours and numbers.
=cut
sub sf_key {
	my $top = $dy;
	my $text_dy;
	my $result = '';
    my $dx = 0; 
    
	my $is_weak_hits = 0;
    unless ($nostruct) {
	foreach my $protein (keys %protein_details) {
		$is_weak_hits = 1 if defined $protein_details{$protein}{'structures'}{'weak'};
		last if $is_weak_hits;
	}
    }
   if(undef){
	$result = "<text x=\"$dx\" y=\"$dy\" font-size=\"12\">Key:</text>";
	$dy += 5;
	$dx += 5;
	$text_dy = $dy + 8;
    
    unless ($nostruct) {
        #Annotate predicted structure block
        $result .= "<rect x=\"$dx\" y=\"$dy\" width=\"10\" height=\"10\" fill=\"#ccc\" ry=\"2\" rx=\"2\" />";
        $result .= "<text x=\"".($dx+12)."\" y=\"$text_dy\" font-size=\"10\">Predicted Structure</text>";
        $dy+=12;
        $text_dy = $dy + 8;

        #Annotate that weak hits are dropped down	
        if ($is_weak_hits) {
            my $y = $dy + 10;protein_details
            $result .= "<line x1=\"".($dx+5)."\" y1=\"$dy\" x2=\"".($dx+5)."\" y2=\"$y\" style=\"fill:none;stroke-dasharray:1,1;stroke:#000;stroke-width:1\"/>";
            $result .= "<text x=\"".($dx+12)."\" y=\"$text_dy\" font-size=\"10\">Weak Structural Hit</text>";
            $dy+=12;
            $text_dy = $dy + 8;	
        }
    }
	

	$dy = $top;
    $dx = 150;
    }


    #Draw the superfamily labels key
    unless ($nostruct) {
        $result .= "<text x=\"$dx\" y=\"$dy\" font-size=\"12\">Superfamilies:</text>";
        $dy += 5;
        $dx += 5;
        $text_dy = $dy + 8;

        foreach my $sf (sort {$sf_details{$a}{'label'} <=> $sf_details{$b}{'label'}} keys %sf_details) {
            my $colour = $sf_details{$sf}{'colour'};
            my $name = $sf_details{$sf}{'description'};
            my $label = $sf_details{$sf}{'label'};
            my $key = ($draw_labels and not $supfam_labels)?"[$label] $name":"$name";
            $result .= "<rect x=\"$dx\" y=\"$dy\" width=\"10\" height=\"10\" fill=\"$colour\" stroke=\"black\" stroke-width=\"0.5\" />";
            $result .= "<text x=\"".($dx+12)."\" y=\"$text_dy\" font-size=\"10\">$key</text>";
            $dy+=12;
            $text_dy = $dy + 8;
        }
        if (scalar keys %sf_details == 0) {
            $result .= "<text x=\"$dx\" y=\"$text_dy\" font-size=\"10\">N/A No Hits</text>";
        }
    }
	#Makesure we have a height that covers the Key section
	$dy = $top + 50 if $dy < ($top + 50);
	return $result;
}

#Get protein and sf details from the database

#Foreach protein draw the entry and then iterate the vertical cursor to the next record
foreach my $protein (keys %protein_details) {
	$output .= draw_protein($protein);
	$current_record++;
	$dy += $record_height;
} #Protein foreach

#Print out amino ruler
if ($draw_ruler) {
	$dy += $record_spacing; #after some vertical white space
	$output .= amino_ruler();
	$dy += 20;
}

#Print out the figure key
$dy += $record_spacing; #after some vertical white space
$output .= sf_key() if $draw_key;

#Print out the SVG footer containing the tooltip popup markup/design
$output .= footer();

#Print out SVG and ECMA script popup header
my $width = $ruler_width + ($xpad*3);
$width = $min_width if $width < $min_width and $draw_labels;
my $height = $dy+$record_spacing;
$output = header($width,$height) . $output;

if ($force_png) {
	use File::Temp qw/tempfile/;
	my ($fh,$filename) = tempfile( 'architectureXXXXXX', UNLINK => 1, TMPDIR => 1, SUFFIX => '.svg' );
	print $fh $output;
	print $cgi->header("image/png");
	open (PNG, '-|', "convert -size ${width}x${height} svg:$filename png:fd:1") or die "Couldn't convert to PNG.";
	my $encode = encode_base64(join('',<PNG>));
	print $encode;
	close(PNG) or die "Failed to close pipe.";
}
else {
	print $cgi->header("image/svg+xml");
	print $output;
}

=back
=cut

1;
