use lib $ENV{ENS} || $ENV{PERL5LIB};
use IO::File;
use strict;
use Getopt::Std;
use warnings;
my $user = "anonuser";
my $host = "manny.i-med.ac.at";
my $pass = "";
my $species = "human";
my $max_mm=0;
my $min_probe_algn=24;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::ApiVersion;
##my $tx_id="CCDS10374.2";
##my $ps_id="1554015_a_at";
my $tx_id="CCDS12689.1";
my $ps_id="38269_at";
my $array_id="HG-U133_Plus_2";
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(-host => $host, -user => $user,
				 -pass => $pass, -verbose => "0" );
## important! set this, otherwise we loose the connection after a certain time
$registry->set_reconnect_when_lost();

my $probe_feature_adaptor = $registry->get_adaptor( $species, "funcgen", "ProbeFeature" );
my $probe_set_adaptor = $registry->get_adaptor( $species, "funcgen", "ProbeSet" );
my $otherfeatures_tx_adaptor = $registry->get_adaptor( $species, "otherfeatures", "transcript" );
my $slice_adaptor = $registry->get_adaptor( $species, "core", "slice" );
my $array_adaptor = $registry->get_adaptor( $species, "funcgen", "Array" );

print "species: ".$slice_adaptor->db->species."\n";

my $array = $array_adaptor->fetch_by_name_vendor( $array_id, "AFFY" );

my $transcript = $otherfeatures_tx_adaptor->fetch_by_stable_id( $tx_id );

$transcript = $transcript->transform( 'chromosome' );
my $strand = $transcript->strand;
my @exons = @{ $transcript->get_all_Exons() };
foreach my $exon (@exons){
  print "Exon: ".$exon->stable_id."\n";
  ## transform to chromosome:
  $exon = $exon->transform( 'chromosome' );
  my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $exon->seq_region_name, $exon->start, $exon->end );
  my @affy_features=@{ $probe_feature_adaptor->fetch_all_by_Slice_Array( $slice, $array ) };
  if( scalar( @affy_features ) > 0 ){
    foreach my $feature (@affy_features){
      print "Feature ".$feature->probe->probeset->name." mm: ".$feature->mismatchcount."\n";
      ## check for strandedness!
      ## AFFY_UTR are "antisense target", AFFY_ST "sense target".
      if( ( $array->class eq "AFFY_UTR" and ( $feature->strand == $strand ) ) or ( $array->class eq "AFFY_ST" and ( $feature->strand != $strand ) ) ){
	## that's fine, now check if the length of alignment fits and the number of mismatches.
	if( $feature->mismatchcount <= $max_mm and ( $feature->end - $feature->start ) >= $min_probe_algn ){
	  ## now add that one.
	  print "Got it!\n";
	}
      }
    }
  }
}


my $probe_set = $probe_set_adaptor->fetch_by_array_probeset_name( $array_id, $ps_id);
print "probe set: ".$probe_set->name." size: ".$probe_set->size."\n";
my @probes = @{$probe_set->get_all_Probes()};
print "No of probes: ".scalar( @probes )."\n";

