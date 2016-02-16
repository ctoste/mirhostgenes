#!/usr/bin/perl

# script to define possible miRNA host genes using the genomic alignments of the (pre-)miRNAs.
# The script reads the genomic alignments of the pre-miRNA sequences from the e.g. hsa.gff file from the
# mirbase and searches various databases for genes that have either an exon or intron at the specific position
# (considering the strand orientation).
# The databases in which genes are searched are: Ensembl core, otherfeatures, cdna and vega database (requires
# that these databases are either installed locally, or that the host for the databases is the primary ensembl
# host). Note, that the script verifies the availability of these databases automatically.
use lib $ENV{ENS} || $ENV{PERL5LIB};
use IO::File;
##use DBI;
use Getopt::Std;
use strict;
use warnings;
my $script_version="0.3.0";
## ensembl
my $user = "anonymous";
my $host = "ensembldb.ensembl.org";
my $pass = "";
my $ensembl_version="none";
my $ensembl_genome_version="none";
my $ensembl_database="core";
## mirbase
my $mirbase_version="none";
my $mirbase_date="none";
my $mirbase_genome_build_id="none";
my $species="none";
my $ensembl_species;
my $gff_file;
my $gff_version="none";
my %option=();
# d: string specifying the database(s)
# g: gff file
# h help
# H: host
# U: user
# P: pass
# v: verbosity
getopts( "d:g:hH:P:U:v", \%option );

if( $option{ h } ){
  ## print help and exit.
  print( "\ndefine_mirna_host_genes.pl version ".$script_version.".\n" );
  print( "Defines potential microRNA host genes based on the genomic alignment of the pre-miRNAs taken from a submitted gff file.\n\n" );
  print( "usage: perl define_mirna_host_genes.pl -dg:hH:P:U:\n" );
  print( "-d (optional): a character string specifying the database(s) to query. Defaults to 'core' but also 'cdna', 'otherfeatures' and 'vega' are allowed. For multiple databases the database names have to be submitted separated by a , with no white spaces." );
  print( "-g (required): the gff file from which the genomic alignments of the pre-miRNAs should be taken. This should be one of the gff3 files provided by the mirbase.\n" );
  print( "-H (optional): the hostname of the Ensembl database; defaults to ensembldb.ensembl.org.\n" );
  print( "-h (optional): print this help message.\n" );
  print( "-P (optional): the password to access the Ensembl database.\n" );
  print( "-U (optional): the username to access the Ensembl database.\n" );
  print( "-v (optional): be verbose.\n" );
  print( "\n\nThe script will generate the following tables:\n" );
  print( "- pre_mirna.txt: contains the chromosomal alignment of the pre-miRNAs.\n" );
  print( "- mat_mirna.txt: contains the chromosomal alignment of the mature miRNAs along with the information from which pre-miRNA they derive.\n" );
  print( "- host_tx.txt: predicted host transcript (defined by its transcript id) for a pre-miRNA, one line per pre-miRNA. Columns in_intron and in_exon specify in which intron or which exon the pre-miRNA is localte (introns and exons numbered from 5' to 3' of the transcript). The pre-miRNA is partially exonic if both the in_intron and in_exon is not 0. is_outside is 1 if a miRNA is partially located outside of the transcript.  Note: the same transcript id might be in more than one row of this table, if the transcript hosts more than one pre-miRNA.\n" );
  print( "- host_gene.txt: table with (unique) host genes; provides gene_id, gene_name, gene_biotype and eventually NCBI EntrezGene identifiers.\n" );
  print( "- metadata.txt: some additional informations: mirbase version, Ensembl version etc.\n" );
  print( "- chromosome.txt: the information for the chromosomes.\n");
  print( "\n" );
  exit 0;
}
if( !defined( $option{ g } ) ){
  die "No input gff file provided!\n"
}else{
  $gff_file = $option{ g };
}
if( defined( $option{ U } ) ){
  $user=$option{ U };
}
if( defined( $option{ H } ) ){
  $host=$option{ H };
}
if( defined( $option{ P } ) ){
  $pass=$option{ P };
}
## checking for databases to query...
my $databasestring="core";
if( defined( $option{ d } ) ){
  $databasestring=$option{d};
}

## let's start...
## which databases to query?
my $do_core=0;
my $do_vega=0;
my $do_otherfeatures=0;
my $do_cdna=0;
## now checking which database we want to query...
evaluateDatabasestring();

## some definition of the elements in the gff3:
my $idx_seq_name = 0;
my $idx_seq_start = 3;
my $idx_seq_end = 4;
my $idx_seq_strand = 6;
my $idx_descr = 8;
###########################
## setting the connection to the Ensembl database up...
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::ApiVersion;
my $registry = 'Bio::EnsEMBL::Registry';
$ensembl_version="".software_version()."";
## will store the pre-mirna into an array and process it later; so we can query
## more than one database.
my @MIR2START = ();          # Start coord
my @MIR2ID = ();             # miRNA ID
my @MIR2END = ();            # end coordinate
my @MIR2STRAND = ();         # strand
my @MIR2CHROMOSOME = ();     # chromosome
my $pre_mirna_pk=0;          ## that'll be the primary key of the pre-miRNA (i.e. the alignment of it...)
my %PK2PREMIR;

## read some informations from the GFF file header.
parseGffHeader();
if( defined( $option{ v } ) ){
  printGffHeaderInfo();
}

## we've set the stage: the MIR2 arrays contain the information we need, now we can
## check what databases we've got available to query.
$registry->load_registry_from_db(-host => $host, -user => $user,
				 -pass => $pass, -verbose => "0" );
## important! set this, otherwise we loose the connection after a certain time
$registry->set_reconnect_when_lost();
evaluateDatabaseAvailability();

## preparing the output files:
open(pre_mirna , ">pre_mirna.txt");
print pre_mirna "pre_mirna_algn_id\tpre_mirna_id\tpre_mirna_name\tseq_name\tseq_strand\tpre_mirna_seq_start\tpre_mirna_seq_end\tpre_mirna_confidence\tpre_mirna_read_count\tpre_mirna_experiment_count\n";
open(mat_mirna , ">mat_mirna.txt");
print mat_mirna "mat_mirna_id\tmat_mirna_name\tpre_mirna_algn_id\tmat_mirna_seq_start\tmat_mirna_seq_end\tmat_mirna_confidence\tmat_mirna_read_count\tmat_mirna_experiment_count\n";

### that's the host_tx table:
## pre_mirna_id: the ID of the pre-miRNA.
## tx_id: the transcript id
## tx_biotype: the biotype of the transcript
## in_intron: in which intron is the miRNA? 0 if not intronic.
## in_exon: in which exon is the miRNA? 0 if not exonic.
## Note: if both in_intron and in_exon are true it means it is partially exonic.
## exon_id: the id of the exon (if it's in an exon).
## gene_id: the gene id for the transcript.
open( host_tx, ">host_tx.txt" );
print host_tx "pre_mirna_algn_id\ttx_id\ttx_biotype\tin_intron\tin_exon\tis_outside\texon_id\tgene_id\tdist_next_exon_start\tdist_last_exon_end\n";

### that's the host_gene table:
## gene_id: the gene id
## gene_biotype: the gene biotype
## gene_name: the name of the gene
## entrezid: the entrezgene ids (if any).
## database: the database in which the gene was defined.
open( host_gene, ">host_gene.txt" );
print host_gene "gene_id\tgene_biotype\tgene_name\tentrezid\tdatabase\tsource\n";

### that's the chromosome information table:
open( chrom_file, ">chromosome.txt" );
print chrom_file "seq_name\tseq_length\tis_circular\n";
my %done_chromosomes=();

## get all genomic alignments for the miRNAs from the GFF.
processGff();

my $current_database;
my $gene_adaptor;
my $slice_adaptor;
if( $do_core==1 ){
  ## query the core database to search for potential miRNA host genes.
  $current_database="core";
  $gene_adaptor = $registry->get_adaptor( $species, $current_database, "gene" );
  $slice_adaptor = $registry->get_adaptor( $species, $current_database, "slice" );
  lookForHostGenes();
}
if( $do_cdna==1 ){
  $current_database="cdna";
  $gene_adaptor = $registry->get_adaptor( $species, $current_database, "gene" );
  $slice_adaptor = $registry->get_adaptor( $species, $current_database, "slice" );
  lookForHostGenes();
}
if( $do_otherfeatures==1 ){
  $current_database="otherfeatures";
  $gene_adaptor = $registry->get_adaptor( $species, $current_database, "gene" );
  $slice_adaptor = $registry->get_adaptor( $species, $current_database, "slice" );
  lookForHostGenes();
}
if( $do_vega==1 ){
  $current_database="vega";
  $gene_adaptor = $registry->get_adaptor( $species, $current_database, "gene" );
  $slice_adaptor = $registry->get_adaptor( $species, $current_database, "slice" );
  lookForHostGenes();
}


## at last writing some informations.
open( info_file, ">metadata.txt" );
print info_file "name\tvalue\n";
print info_file "Db type\tMirhostDb\n";
print info_file "Supporting package\tmirhostgenes\n";
print info_file "Db created by\tmirhostgenes package from Bioconductor\n";
print info_file "Creation time\t".localtime()."\n";
print info_file "script_version\t$script_version\n";
print info_file "mirbase_version\t$mirbase_version\n";
print info_file "ensembl_version\t$ensembl_version\n";
print info_file "genome_build\t$mirbase_genome_build_id\n";
print info_file "mirbase_date\t$mirbase_date\n";
print info_file "Organism\t$ensembl_species\n";
print info_file "did_core\t$do_core\n";
print info_file "did_cdna\t$do_cdna\n";
print info_file "did_otherfeatures\t$do_otherfeatures\n";
print info_file "did_vega\t$do_vega\n";
print info_file "DBSCHEMAVERSION\t1.0\n";
close( info_file );

close( pre_mirna );
close( mat_mirna );
close( host_tx );
close( host_gene );
close( chrom_file );
$registry -> disconnect_all();


#######################################################################################################
##
##                     SUBS
##
##


## processing the GFF file:
## reading line by line, if it's a "miRNA_primary_transcript" (i.e. the pre-miRNA/hairpin), we're looking
## for genes in Ensembl, if it's a "miRNA" we just export that info.
## writing output tables:
## pre_mirna: pre_mirna_id, pre_mirna_name, seq_name, seq_strand, seq_start, seq_end
## mat_mirna: mat_mirna_id, mat_mirna_name, pre_mirna_id, seq_name, seq_strand, seq_start, seq_end
sub processGff{
  my @cells;
  open( IN, "< $gff_file" ) or die "can't open file $gff_file!\n";
  while( <IN> ){
    chomp;
    if( /^#/ ){
      ## a comment...
    }
    else{
      @cells = split( /\t/,$_ );
      if( $cells[ 2 ] eq "miRNA_primary_transcript" ){
	$pre_mirna_pk++;   ## that's what I call cool autoincremental primary key ;) (just kidding)
	processPreMirna( @cells );
      }
      if( $cells[ 2 ] eq "miRNA" ){
	processMatMirna( @cells );
      }
    }
  }
  close(IN);
}

## process a line for a pre-miRNA (pri-miRA)
sub processPreMirna{
  my @cells=@_;
  my $seq_name = $cells[ $idx_seq_name ];
  $seq_name =~ s/chr//g;
  my $seq_start = $cells[ $idx_seq_start ];
  my $seq_end = $cells[ $idx_seq_end ];
  my $seq_strand = strand2int( $cells[ $idx_seq_strand ] );
  $cells[ $idx_descr ] =~ /ID=(.*);Alias.+/;
  my $id = $1;
  ## be aware! the id could also be something like MI0004127_2!!!s
  my @tmp = split( /_/,$id );
  $id = $tmp[0];
  $cells[ $idx_descr ] =~ /Name=(.*)/;
  my $name = $1;
  ## saving values.
  print pre_mirna "$pre_mirna_pk\t$id\t$name\t$seq_name\t$seq_strand\t$seq_start\t$seq_end\t0\t0\t0\n";
  ## adding these to the arrays for later annotation.
  push( @MIR2ID, $pre_mirna_pk );
  push( @MIR2CHROMOSOME, $seq_name );
  push( @MIR2STRAND, $seq_strand );
  push( @MIR2START, $seq_start );
  push( @MIR2END, $seq_end );
  ## that's to cross-check that mature miRNAs and pre-miRNAs really map correctly!
  $PK2PREMIR{ $pre_mirna_pk } = $id;
}


## process a line for a mature miRNA
sub processMatMirna{
  my @cells=@_;
  my $seq_name = $cells[ $idx_seq_name ];
  $seq_name =~ s/chr//g;
  my $seq_start = $cells[ $idx_seq_start ];
  my $seq_end = $cells[ $idx_seq_end ];
  my $seq_strand = strand2int( $cells[ $idx_seq_strand ] );
  $cells[ $idx_descr ] =~ /ID=(.*);Alias.+/;
  my $id = $1;
  ## be aware! the id could also be something like MI0004127_2!!!s
  my @tmp = split( /_/,$id );
  $id = $tmp[0];
  $cells[ $idx_descr ] =~ /Name=(.*);/;
  my $name = $1;
  $cells[ $idx_descr ] =~ /Derives_from=(.*)/;
  my $pre_mirna_id = $1;
  if( $PK2PREMIR{ $pre_mirna_pk } ne $pre_mirna_id ){
    die( "pre-miRNA id does not match for mature miRNA $id: $pre_mirna_id != ".$PK2PREMIR{ $pre_mirna_pk } );
  }
  ## saving values.
  ## Note that we assume that after each pre-miRNA its mature miRNAs are listed in the file!!!
  print mat_mirna "$id\t$name\t$pre_mirna_pk\t$seq_start\t$seq_end\t0\t0\t0\n";
}


## now that's the function that looks for host genes based on the chromosomal coordinates stored in the MIR2 arrays
## we're using the pre-defined slice_adapter
sub lookForHostGenes{
  if( defined( $option{ v } ) ){
    print "looking for host genes.\n";
  }
  my $seq_name;
  my $seq_strand;
  my $seq_start;
  my $seq_end;
  my $pre_mirna_id;
  my $tx_id;
  my $tx_biotype;
  my %GENEBIOTYPE=();
  my %GENENAME=();
  my %GENEENTREZID=();
  my %GENESOURCE=();
  for( my $i=0; $i < scalar @MIR2ID; $i++ ){
    $seq_name=$MIR2CHROMOSOME[ $i ];
    $seq_strand=$MIR2STRAND[ $i ];
    $seq_start=$MIR2START[ $i ];
    $seq_end=$MIR2END[ $i ];
    $pre_mirna_id=$MIR2ID[ $i ];
    if( defined( $option{ v } ) ){
      print "$PK2PREMIR{ $pre_mirna_id }; $i of ".scalar(@MIR2ID)." ";
    }
    my $slice=$slice_adaptor->fetch_by_region('chromosome', $seq_name, $seq_start, $seq_end );
    if( exists( $done_chromosomes{ $seq_name } ) ){
    }else{
      $done_chromosomes{ $seq_name } = "done";
      my $chr_slice=$slice_adaptor->fetch_by_region('chromosome', $seq_name );
      my $name = $chr_slice->seq_region_name;
      my $length = $chr_slice->length;
      my $is_circular = $chr_slice->is_circular;
      print chrom_file "$name\t$length\t$is_circular\n";
    }
    ## check for genome version if i=0;
    if( $i==0 ){
      $ensembl_genome_version=$slice->coord_system()->version();
      ##print "$ensembl_genome_version $mirbase_genome_build_id\n";
      ## first check if versions match per se
      if( $mirbase_genome_build_id ne $ensembl_genome_version ){
	## maybe one of the two has some subversion?
	my @subs = split( /\./,$mirbase_genome_build_id );
	if( $subs[ 0 ] ne $ensembl_genome_version ){
	  die "Incompatible genome build versions: Ensembl has '$ensembl_genome_version', mirbase has '$mirbase_genome_build_id'\n";
	}
      }
    }
    ## now get all transcripts:
    my @txs = @{$slice -> get_all_Transcripts( 1 )};
    foreach my $tx (@txs){
      $tx = $tx->transform( 'chromosome' );
      if( $tx->strand == $seq_strand ){
	$tx_biotype="";
	my $is_intronic=0;
	my $is_exonic=0;
	my $exon_id="";
	my $is_outside=0;      ## that's if a miRNA is partially exonic, but not all of the miRNA is within the exon/gene.
	##$tx_id = $tx->stable_id;
	$tx_id = $tx->display_id();
	$tx_biotype = $tx->biotype;
	## get exon stuff
	my $exon_number=0;
	my $current_exon_start=0;
	my $current_exon_end=0;
	my $last_exon_end=0;
	my $last_exon_start=0;
	my $dist_next_exon_start="";
	my $dist_last_exon_end="";
	## note: the loop below doesn't check if a pre-miRNA matches two exons, but rather assumes that a pre-miRNA is not big
	## enough to span a complete intron!
	## idea is:
	## if start is >= exon_start and <= exon_end it's at least partially exonic.
	##    if in the above situation also the end is within the exon -> completely exonic.
	## if end is >= exon_start and <= exon_end it's at least partially exonic.
	## to test for intronic is tricky:
	## if on + strand: start > last_exon_end & end < exon_start: is in intron exon_number -1.
	## if on - strand: start > exon_end & end < last_exon_start: is in intron exon_number -1.
	my @exons = @{$tx->get_all_Exons};
	foreach my $exon (@exons){
	  $exon_number++;
	  $current_exon_start = $exon->start;
	  $current_exon_end = $exon->end;
	  ## check: is the start within the exon?
	  if( $seq_start >= $current_exon_start and $seq_start <= $current_exon_end ){
	    $is_exonic = $exon_number;  ## keep record of the exon.
	    $exon_id = $exon->stable_id();
	    ## it's at least parially exonic.
	    ## now check if the end is also within the exon.
	    if( $seq_end <= $current_exon_end ){
	      ## ok, so it's only within the exon!
	      last;
	    }
	    if( $seq_strand > 0 ){
	      ## now, if that was the last exon, then it is parially outside of the transcript
	      ## for transcripts on the + strand.
	      if( $exon_number == scalar( @exons ) ){
		$is_outside=1;
	      }else{
		## so, there is also another exon.
		## record also the intron; thus it's a partially exonic
		$is_intronic = $exon_number;
	      }
	    }else{
	      ## if the transcript is on - strand, and only the start is within the exon:
	      ## if it's the first exon: the end is outside:
	      if( $exon_number == 1 ){
		$is_outside=1;
	      }else{
		## otherwise it's in the exon and the previous intron.
		$is_intronic = $exon_number -1;
	      }
	    }
	    last;
	  }
	  ## check: is the end within the exon?, the start wasn't otherwise we wouldn't be here.
	  if( $seq_end >= $current_exon_start and $seq_end <= $current_exon_end ){
	    ## keep track of the exon:
	    $is_exonic = $exon_number;
	    $exon_id = $exon->stable_id();
	    ## however, the miRNA could still be outside of the transcript (for + strand)
	    ## and first exon.
	    if( $seq_strand > 0 ){
	      if( $exon_number==1 ){
		$is_outside=1;
	      }else{
		## otherwise the start of the miRNA was in the last intron.
		$is_intronic=$exon_number-1;
	      }
	    }else{
	      ## if it's the last exon then the start of the miRNA is outside the transcript.
	      if( $exon_number == scalar( @exons ) ){
		$is_outside=1;
	      }else{
		## otherwise it is in the intron following that exon.
		$is_intronic=$exon_number;
	      }
	    }
	    last;
	  }
	  ## check if it's completely intronic: neither start nor end is in the exon.
	  if( $seq_strand > 0 ){
	    if( $seq_start > $last_exon_end and $seq_end < $current_exon_start ){
	      ## is in the intron between the last exon and the current one.
	      $is_intronic = $exon_number-1;
	      ## distance to exons is given 5-3' of transcript; distance is in nt (i.e. nt between the miRNA and exon)
	      $dist_next_exon_start=$current_exon_start-$seq_end-1;
	      $dist_last_exon_end=$seq_start-$last_exon_end-1;
	      last;
	    }
	  }
	  if( $seq_strand < 0 ){
	    if( $seq_start > $current_exon_end and $seq_end < $last_exon_start ){
	      ## is in the intron between the last exon and the current one.
	      $is_intronic = $exon_number -1;
	      ## distance to exons is given 5-3' of transcript; distance is in nt (i.e. nt between the miRNA and exon)
	      $dist_next_exon_start=$seq_start-$current_exon_end-1;
	      $dist_last_exon_end=$last_exon_start-$seq_end-1;
	      last;
	    }
	  }
	  $last_exon_start = $current_exon_start;
	  $last_exon_end = $current_exon_end;
	}
	## get all gene stuff for the transcript.
	my $gene = $gene_adaptor->fetch_by_transcript_stable_id( $tx_id );
	my $gene_id="";
	my $gene_biotype="";
	my $gene_entrezid="";
	my $gene_name="";
	my $gene_source="";
	if( defined( $gene ) ){
	  ##$gene_id=$gene->stable_id;
	  $gene_id = $gene->display_id();
	  $gene_source = $gene->source();
	  if( exists( $GENEBIOTYPE{ $gene_id } ) ){
	    ## did already have that one...
	  }else{
	    my $gn = $gene->external_name;
	    if( !defined( $gn ) ){
	      $gn = "";
	    }
	    $GENENAME{ $gene_id }=$gn;
	    my $gbt = $gene->biotype;
	    if( !defined( $gbt ) ){
	      $gbt = "";
	    }
	    $GENEBIOTYPE{ $gene_id}=$gbt;
	    ## try to get the Entrezgene ID(s)
	    my $all_entries = $gene->get_all_DBLinks( "EntrezGene" );
	    my %entrezgene_hash=();
	    foreach my $dbe ( @{$all_entries} ){
	      $entrezgene_hash{ $dbe->primary_id } = 1;
	    }
	    my $hash_size = keys %entrezgene_hash;
	    my $entrezid = "";
	    if( $hash_size > 0 ){
	      $entrezid = join( ";", keys %entrezgene_hash );
	    }
	    ## gene biotype is OK gene name and entrezid don't work
	    ## for otherfeatures database.
	    if( $gene_source eq "refseq" ){
	      $GENEENTREZID{ $gene_id } = $gene_id;
	    }else{
	      $GENEENTREZID{ $gene_id } = $entrezid;
	    }
	    ## new stuff... to remove?
	    $GENESOURCE{ $gene_id } = $gene_source;
	  }
	}
	## write the transcript information...
	## pre_mirna_id, tx_id, tx_biotype, in_intron, in_exon, exon_id, gene_id
	print host_tx "$pre_mirna_id\t$tx_id\t$tx_biotype\t$is_intronic\t$is_exonic\t$is_outside\t$exon_id\t$gene_id\t$dist_next_exon_start\t$dist_last_exon_end\n";
	if( defined( $option{ v } ) ){
	  print ".";
	  ##print " in intron $is_intronic, in exon $is_exonic";
	}
      }
    }
    if( defined( $option{ v } ) ){
      print "\n";
    }
  }
  ## at last also processing the host_gene table...
  foreach my $key (sort keys %GENEBIOTYPE){
    print host_gene "$key\t$GENEBIOTYPE{ $key }\t$GENENAME{ $key }\t$GENEENTREZID{ $key }\t$current_database\t$GENESOURCE{ $key }\n";
  }
}

## transforms a + or - to 1 or -1
sub strand2int{
  my ($strand)=@_;
  if( $strand eq "+" ){
    return 1;
  }
  if( $strand eq "-" ){
    return -1;
  }
  return 0;
}

## read the comments from the GFF file and define the required parameters.
sub parseGffHeader{
  open( IN, "< $gff_file" ) or die "can't open file $gff_file!\n";
  while( <IN> ){
    chomp;
    if( /^#/ ){
      if( $_ =~ /gff-version (.*)/ ){
	$gff_version = "$1";
	next;
      }
      if( $_ =~ /date (.*)/ ){
	$mirbase_date = "$1";
	next;
      }
      if( $_ =~ /Chromosomal coordinates of (.*) microRNAs/ ){
      	$species = "$1";
      	next;
      }
      if( $_ =~ /miRBase (.*)/ ){
      	$mirbase_version="$1";
      	next;
      }
      if( $_ =~ /genome-build-id\:\s+(.*)/ ){
      	$mirbase_genome_build_id="$1";
##	$mirbase_genome_build_id=
      	next;
      }
    }
    else{
      last;
    }
  }
  close(IN);
  ## now we are evaluating whether we have all we need:
  if( $mirbase_date eq "none" | $species eq "none" | $mirbase_genome_build_id eq "none" | $mirbase_version eq "none" ){
    die "The GFF file $gff_file is not in the expected format!";
  }
  if( $gff_version ne "3" ){
    die "The GFF has to be in GFF version 3 format! Please submit one of the gff3 files from miRBase."
  }
}


sub printGffHeaderInfo{
  print "mirbase date: ".$mirbase_date."\n";
  print "mirbase version: ".$mirbase_version."\n";
  print "species: ".$species."\n";
  print "mirbase genome build: ".$mirbase_genome_build_id."\n";
}

sub evaluateDatabasestring{
  ## the $databasestring string has to be a , separated list of database names!!!
  my @databases = split /,/,$databasestring;
  foreach my $database (@databases){
    if( $database eq "core" ){
      $do_core = 1;
    }
    if( $database eq "vega" ){
      $do_vega = 1;
    }
    if( $database eq "cdna" ){
      $do_cdna = 1;
    }
    if( $database eq "otherfeatures" ){
      $do_otherfeatures = 1;
    }
  }
  if( $do_core==0 and $do_vega==0 and $do_cdna==0 and $do_otherfeatures==0 ){
    die "The specified databases are not valid! Allowed are core, cdna, otherfeatures or vega!";
  }
}

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
	print "core database abvailable.\n";
      }
    }
  }
  if( $do_cdna == 1 ){
    my $testslice = $registry->get_adaptor( $species, "cdna", "slice" );
    if( !defined( $testslice ) ){
      $ensembl_species = $testslice->db->species;
      $do_cdna=0;
      print "WARNING: no cdna database available! Will exclude this database.\n";
    }else{
      $ensembl_species = $testslice->db->species;
      if( defined( $option{ v } ) ){
	print "cdna database abvailable.\n";
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
	print "otherfeatures database abvailable.\n";
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
	print "vega database abvailable.\n";
      }
    }
  }
  ## now check whether we've got any database left at all.
  if( $do_core==0 and $do_cdna==0 and $do_otherfeatures==0 and $do_vega==0 ){
    die( "No database available to query! Have to stop.\n" );
  }
}


