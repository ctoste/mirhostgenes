############################
## version 0.1.3, 22.12.2011
## works from Ensembl version 56 onwards (oligo table was moved form core to funcgen database)

use lib $ENV{ENS} || $ENV{PERL5LIB};
use IO::File;
use strict;
use Getopt::Std;
use warnings;

my $version_="0.2.3";

my $max_mm=0;
my $min_probe_algn=24;
my $prop_probes = 0.8;   ## proportion of probes per probe set that have to map the transcript; for affy this means 9 or more (9/11=0.81).
my $array_type = "HG-U133_Plus_2,PrimeView";
my $species = "human";
my $tx_file = "host_tx.txt";
my $out_file = "array_features.txt";
## ensembl
my $user = "anonymous";
my $host = "ensembldb.ensembl.org";
my $pass = "";
my $ensembl_version="none";
my $ensembl_genome_version="none";
my $ensembl_species;
my $tx_id_column_name = "tx_id";
my $tx_id_idx = -1;

my %option=();
#getopts("i:o:f:",\%option);
getopts("i:a:p:s:hH:m:M:P:U:v",\%option);
if( defined( $option{ h } ) ){
  ## print help and exit.
  print( "\nget_array_features_for_host_tx.pl version ".$version_.".\n" );
  print( "-a: the type of Affymetrix GeneChip. Defaults to HG-U133_Plus_2,PrimeView. Multiple , separated values are allowed (no white spaces).\n" );
##  print( "-e (required): Ensembl version. The function ensures internally whether the submitted Ensembl version matches the Ensembl Perl API and database version.\n" );
  print( "-h print this help and exit.\n" );
  print( "-H: the hostname of the Ensembl database; defaults to ensembldb.ensembl.org.\n" );
  print( "-m: the maximal number of allowed mismatches in the alignment; defaults to 0.\n" );
  print( "-M: the required length of the probe alignment; default to 24 meaning that all nucleotides of a 25nt long probe have to be aligned\n" );
  print( "-P: the password to access the Ensembl database.\n" );
  print( "-U: the username to access the Ensembl database.\n" );
  print( "-i: the input file in which the transcript IDs can be found. Defaults to host_tx.txt.\n" );
  print( "-p: the minimal required proportion of probes in a probe set to target an exon of the transcript. Defaults to 0.8, which means for a 11 probe probe set that 9 or more probes have to map to the transcript.\n" );
  print( "-s: the species (e.g. -s mouse). Defaults to human.\n" );
  print( "-v: print messages to the terminal.\n\n" );
  exit 0;
}
if( defined( $option{ i } ) ){
  $tx_file= "host_tx.txt";
}
if( defined( $option{ s } ) ){
  $species = $option{ s };
}
if( defined( $option{ a } ) ){
  $array_type=$option{ a };
}
if( defined( $option{ U } ) ){
  $user=$option{ U };
}
if( defined( $option{ H } ) ){
  $host=$option{ H };
}
if( defined( $option{ m } ) ){
  $max_mm = $option{ m };
}
if( defined( $option{ M } ) ){
  $min_probe_algn = $option{ M };
}
if( defined( $option{ P } ) ){
  $pass=$option{ P };
}
if( defined( $option{ p } ) ){
  $prop_probes=$option{ p };
}
# if( defined( $option{ e } ) ){
#   $ensembl_version=$option{ e };
# }
# else{
#   die "the Ensembl version has to be specified with the -e parameter!";
# }

## which databases to query?
my $do_core=1;
my $do_vega=1;
my $do_otherfeatures=1;
my $do_cdna=0;
##
my $probe_feature_adaptor;
my $probe_set_adaptor;
my $slice_adaptor;
my $array_adaptor;
## other adaptors:
my $core_tx_adaptor;
my $vega_tx_adaptor;
my $otherfeatures_tx_adaptor;
my $cdna_tx_adaptor;
my %TX2PS;  ## mapping of tx_id;chip name to probe set id.
my %Arrays;  ## the array adaptors...

## what are we going to do?
## 1) process the input file line by line
## 2) check if for the tx_id;chip_type combination there is already something in the TX2PS, if not:
## 3) guess from the tx id whether we have to query the core (ENS), otherfeatures or vega (OTT) database.
## 4) fetch all exons for that tx.
## 5) loop through the exons...

## setting the connection to the Ensembl database up...
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::ApiVersion;
my $registry = 'Bio::EnsEMBL::Registry';
$ensembl_version = "".software_version()."";
# if( $ensembl_version ne "".software_version().""){
#   die( "Ensembl versions do not match! Version $ensembl_version was submitted, the API has version ".software_version()."!\n" );
# }
## check what databases we've got available to query.
$registry->load_registry_from_db(-host => $host, -user => $user,
				 -pass => $pass, -verbose => "0" );
## important! set this, otherwise we loose the connection after a certain time
$registry->set_reconnect_when_lost();
evaluateDatabaseAvailability();
## initialize the array adaptors.
initializeArrayAdaptors();

## run the stuff...
open( IN, "< $tx_file" ) or die "can't open file $tx_file!\n";
my $firstline=1;
my $tx_counter=0;
my @cells;
while( <IN> ){
  chomp;
  if( /^#/ ){
  }else{
    ## read the first line.
    if( $firstline==1 ){
      $firstline=0;
      @cells = split( /\t/,$_ );
      ## search for the column with the tx_id:
      for( my $idx=0; $idx < scalar( @cells ); $idx++ ){
	if( $cells[ $idx ] eq $tx_id_column_name ){
	  $tx_id_idx = $idx;
	}
      }
      if( $tx_id_idx < 0 ){
	die( "No column named $tx_id_column_name found in the input file!" );
      }
    }else{
      ## OK, now we're really starting to work!
      @cells = split( /\t/,$_ );
      $tx_counter++;
      processTx( $cells[ $tx_id_idx ] );
    }
  }
}
close( IN );
## now we can process the %TX2PS array.
open( OUT, ">$out_file" ) or die "can't open file $out_file for writing!\n";
## write the settings.
writeSettingsHeader();
print OUT "tx_id\tprobeset_id\tarray_id\tprobes_in_tx\n";
foreach my $key (keys %TX2PS){
  ## TX2PS entries are arrays!!! (which can also empty!!!)
  foreach my $allvals (@{$TX2PS{ $key }}){
    my @vals = split /:/,$allvals;
    print OUT "$key\t".$vals[ 1 ]."\t".$vals[ 0 ]."\t".$vals[2]."\n";
  }
  ##my @vals = split /:/,$TX2PS{ $key };
  ##print OUT "$key\t".$vals[ 1 ]."\t".$vals[ 0 ]."\t".$vals[2]."\n";
}
close( OUT );

##*************************************************************************
##
## Subs
##
##*************************************************************************

sub writeSettingsHeader{
  print OUT "# ensembl_version=".$ensembl_version."\n";
  print OUT "# species=".$ensembl_species."\n";
  print OUT "# array_type=".$array_type."\n";
  print OUT "# prop_probes=".$prop_probes."\n";
  print OUT "# max_mm=".$max_mm."\n";
  print OUT "# min_probe_algn=".$min_probe_algn."\n";
}

##
##
## That's the thing.
## check for the submitted transcript which probe sets would eventually target it.
## Note that $TX2PS is a hash of ARRAYS!!!
sub processTx{
  my $id = $_[ 0 ];
  my $transcript;
  my $strand;
  my %Probes=();
  ## check if we do have already probe sets for this transcript -> skip!
  if( exists( $TX2PS{ $id } ) ){
    if( defined( $option{ v } ) ){
      print "Have already transcript $id\n";
    }
  }else{
    if( $id =~ /^ENS*/ & !( $id =~ /(EST)/ )){
      ## processing the core database.
      if( defined( $option{ v } ) ){
	print( "core: ".$id."(".$tx_counter." of xxx)\n" );
      }
      if( !defined( $core_tx_adaptor ) ){
	## initialize adaptor if it was never used before...
	$core_tx_adaptor = $registry->get_adaptor( $species, "core", "transcript" );
      }
      $transcript = $core_tx_adaptor->fetch_by_stable_id( $id );
    }elsif( $id =~ /^OTT*/ ){
      if( defined( $option{ v } ) ){
	print( "vega: ".$id."(".$tx_counter." of xxx)\n" );
      }
      if( !defined( $vega_tx_adaptor ) ){
	$vega_tx_adaptor = $registry->get_adaptor( $species, "vega", "transcript" );
      }
      ## processing the vega database.
      $transcript = $vega_tx_adaptor->fetch_by_stable_id( $id );
    }else{
      ## assume that's an otherfeatures
      if( defined( $option{ v } ) ){
	print( "otherfeatures: ".$id."(".$tx_counter." of xxx)\n" );
      }
      if( !defined( $otherfeatures_tx_adaptor ) ){
	$otherfeatures_tx_adaptor = $registry->get_adaptor( $species, "otherfeatures", "transcript" );
      }
      ## processing the otherfeatures database.
      $transcript = $otherfeatures_tx_adaptor->fetch_by_stable_id( $id );
    }
    if( defined( $transcript ) ){
      $transcript = $transcript->transform( 'chromosome' );
      $strand = $transcript->strand;
      ## now get all exons
      my @exons = @{ $transcript->get_all_Exons() };
      foreach my $exon (@exons){
	## transform to chromosome:
	$exon = $exon->transform( 'chromosome' );
	my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $exon->seq_region_name, $exon->start, $exon->end );
	foreach my $key (keys %Arrays){
	  my $array = $Arrays{ $key };
	  my @affy_features=@{ $probe_feature_adaptor->fetch_all_by_Slice_Array( $slice, $array ) };
	  if( scalar( @affy_features ) > 0 ){
	    foreach my $feature (@affy_features){
	      ## check for strandedness!
	      ## AFFY_UTR are "antisense target", AFFY_ST "sense target".
	      if( ( $array->class eq "AFFY_UTR" and ( $feature->strand == $strand ) ) or ( $array->class eq "AFFY_ST" and ( $feature->strand != $strand ) ) ){
		## that's fine, now check if the length of alignment fits and the number of mismatches.
		if( $feature->mismatchcount <= $max_mm and ( $feature->end - $feature->start ) >= $min_probe_algn ){
		  ## now add that one.
		  my $psid = $array->name.":".$feature->probe->probeset->name;
		  if( exists( $Probes{ $psid } ) ){
		    my $dummy = $Probes{ $psid };
		    $dummy++;
		    $Probes{ $psid } = $dummy;
		  }else{
		    $Probes{ $psid } = 1;
		  }
		}
	      }
	    }
	  }
	}
      }
      my $ps_array=[];
      ## done with the exons... now check if we could add some of the probe sets...
      foreach my $key (keys %Probes){
	my @vals = split /:/,$key;
	my $probe_set = $probe_set_adaptor->fetch_by_array_probeset_name($vals[ 0 ], $vals[ 1 ]);
	my $probe_count = $Probes{ $key };
	if( defined( $option{ v } ) ){
	  print( " - ".$key." probe count: ".$probe_count." of ".$probe_set->size()."..." );
	}
	if( ( $probe_count / $probe_set->size() ) > $prop_probes ){
	  ##$TX2PS{ $id } = $key.":".$probe_count;
	  push( @$ps_array, $key.":".$probe_count );
	  if( defined( $option{ v } ) ){
	    print( "added!" );
	  }
	}
	if( defined( $option{ v } ) ){
	  print( "\n" );
	}
      }
      ## add also an empty array; this avoids multiple queries for the same tx.
      $TX2PS{ $id } = $ps_array;
    }
  }
}

##
##
## Initialize the database connections etc.
sub evaluateDatabaseAvailability{
  ## here we check if the databases are available.
  my @dbadaptors = @{ $registry->get_all_DBAdaptors( $species ) };
  if( scalar @dbadaptors == 0 ){
    die "Error: no databases found for species $species!\nYou might consider to use install the appropriate databases or query a different host.\n";
  }
  ## check if we've got the core database for the species.
  if( $do_core == 1 ){
    my $testslice = $registry->get_adaptor( $species, "core", "slice" );
    if( !defined( $testslice ) ){
      $do_core=0;
      print "WARNING: no core database available! Will exclude this database.\n";
    }else{
      $ensembl_species = $testslice->db->species;
      if( defined( $option{ v } ) ){
	print "core database available.\n";
      }
    }
  }
  if( $do_cdna == 1 ){
    my $testslice = $registry->get_adaptor( $species, "cdna", "slice" );
    if( !defined( $testslice ) ){
      $do_cdna=0;
      print "WARNING: no cdna database available! Will exclude this database.\n";
    }else{
      $ensembl_species = $testslice->db->species;
      if( defined( $option{ v } ) ){
	print "cdna database available.\n";
      }
    }
  }
  if( $do_otherfeatures == 1 ){
    my $testslice = $registry->get_adaptor( $species, "otherfeatures", "slice" );
    if( !defined( $testslice ) ){
      $do_otherfeatures=0;
      print "WARNING: no otherfeatures database available! Will exclude this database.\n";
    }else{
      $ensembl_species = $testslice->db->species;
      if( defined( $option{ v } ) ){
	print "otherfeatures database available.\n";
      }
    }
  }
  if( $do_vega == 1 ){
    my $testslice = $registry->get_adaptor( $species, "vega", "slice" );
    if( !defined( $testslice ) ){
      $do_vega=0;
      print "WARNING: no vega database available! Will exclude this database.\n";
    }else{
      $ensembl_species = $testslice->db->species;
      if( defined( $option{ v } ) ){
	print "vega database available.\n";
      }
    }
  }
  ## the most important thing: the funcgen database
  $probe_feature_adaptor = $registry->get_adaptor( $species, "funcgen", "ProbeFeature" );
  if( !defined( $probe_feature_adaptor ) ){
    die( "The funcgen database is not available!" );
  }
  $array_adaptor = $registry->get_adaptor( $species, "funcgen", "Array" );
  $probe_set_adaptor = $registry->get_adaptor( $species, "funcgen", "ProbeSet" );
  $slice_adaptor = $registry->get_adaptor( $species, "core", "slice" );

  ## now check whether we've got any database left at all.
  if( $do_core==0 and $do_cdna==0 and $do_otherfeatures==0 and $do_vega==0 ){
    die( "No database available to query! Have to stop.\n" );
  }
}

sub initializeArrayAdaptors{
  ## split the array string and add adaptors to the hash, if present...
  my @arrays = split /,/,$array_type;
  foreach my $array (@arrays){
    my $tmp_array = $array_adaptor->fetch_by_name_vendor( $array, "AFFY" );
    if( defined( $tmp_array ) ){
      $Arrays{ $array } = $tmp_array;
    }
  }
  ## test how many we've got...
  my $hash_size = keys %Arrays;
  if( $hash_size==0 ){
    die( "None of the specified arrays ($array_type) are available in Ensembl!" );
  }
}





