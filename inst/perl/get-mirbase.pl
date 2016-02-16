#!/usr/bin/perl
# small script to download the latest (or the specified) miRBase database
# options:
# -v: miRBase version; optional, if not specified the latest will be downloaded.
# -f force downloading the miRBase version.
# -p: the local path where to save the mirbase; defaults to "."
use Net::FTP;
#use Net::FTP::File;
use IO::File;
use DBI;
use Getopt::Std;
use strict;
use warnings;
use FileHandle;
my $script_version_="0.2.0";
###########################################################
## settings:
my $mirbase_path_local=".";
my $mirbase_path_remote="/pub/mirbase/";
my $mirbase_host="mirbase.org";
my $username="anonymous";
my $password='';
###########################################################

my $ftp = Net::FTP -> new($mirbase_host, Timeout=>1200);
unless(defined $ftp){
    print "$@\n";
    die "Can't connect!\n";
}
$ftp -> login($username,$password) || die "Can't login $!";

my $mirbase_version;
my $force_download=0;
my %option=();
getopts("v:p:f",\%option);
if( $option{ f } ){
  $force_download=1;
}
if( defined( $option{ p } ) ){
  $mirbase_path_local = $option{ p };
}
if(  defined( $option{ v } ) ){
	$mirbase_version = $option{ v };
}else{
  ## try to get the latest miRBase.
  $ftp -> cwd( $mirbase_path_remote."CURRENT" );
  ## have to read the version information from the README file!
  my ($remote_file_content, $remote_file_handle);
  open($remote_file_handle, '>', \$remote_file_content);
  $ftp->get('README', $remote_file_handle )
    or die "get failed ", $ftp->message;
  my @lines = split /\n/,$remote_file_content;
  $lines[0] =~ /Release{1}\s(.+)$/;
  if( !defined $1 ){
    die "error: can not get mirbase version from string ".$lines[0]."!\n";
  }
  $mirbase_version=$1;
}

print "This is get-mirbase.pl version $script_version_:\ngetting miRBase version $mirbase_version from $mirbase_host and installing it into $mirbase_path_local\n";

my $localpath = $mirbase_path_local."/".$mirbase_version;
if( ! -e $localpath ){
  mkdir( $localpath );
  if( ! -e $localpath ){
    die "error: directory $localpath can not be created!\n";
  }
  downloadFiles();
}
else{
  if( $force_download==1 ){
    print "over-writing existing version $mirbase_version.\n";
    downloadFiles();
  }
  else{
    print "have already version $mirbase_version locally.\n";
  }
}

sub downloadFiles{
  $ftp -> binary;
  # download files in root dir
  $ftp -> cwd( $mirbase_path_remote.$mirbase_version );
  my $files=$ftp->ls();
  foreach my $file ( @$files ){
    $ftp -> get( $file, $localpath."/".$file );
  }
  # database_files
  $ftp -> cwd( $mirbase_path_remote.$mirbase_version."/database_files" );
  $files=$ftp->ls();
  mkdir( $localpath."/database_files/" );
  foreach my $file ( @$files ){
    $ftp -> get( $file, $localpath."/database_files/".$file );
  }
  # genomes
  $ftp -> cwd( $mirbase_path_remote.$mirbase_version."/genomes" );
  $files=$ftp->ls();
  mkdir( $localpath."/genomes/" );
  foreach my $file ( @$files ){
    $ftp -> get( $file, $localpath."/genomes/".$file );
  }
  $ftp->quit;
}

## unzip all txt and sql tables in the folder database_files
sub unzipDatabase{
  if( system( "gunzip $localpath/database_files/*.gz" ) !=0 ){
    die "Error while unzipping of database files! Last message: $!\n";
  }
}

