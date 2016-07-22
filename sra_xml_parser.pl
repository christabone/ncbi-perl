#!usr/bin/perl
#sra_xml_parser.pl
#A script to parse the XML fields of sra.xml files obtained from NCBI.
#Derived from the processing_from_ncbi.pl script used in-house at FlyBase.
#Version 0.1

use warnings;
use strict;
use Data::Dumper;
use XML::LibXML;

my $OUTFILE = "output.log";
open(STDERR, "| tee -i $OUTFILE") or die "WARNING: ERROR: Cannot open output.log\n";
open(STDOUT, "| tee -i $OUTFILE") or die "WARNING: ERROR: Cannot open output.log\n";

# Reading information from SRA files in SRA directory.
my $sra_directory = "SRA";
opendir SRADIR, $sra_directory or die "Can't open SRA directory: $!";
my @SRAfiles = grep {$_ ne '.' && $_ ne '..' } readdir SRADIR; # Read the SRA file names into an array.
closedir SRADIR;

# Subroutine to parse SRA files.
sub parse {
	my $parser = XML::LibXML->new("1.0","UTF-8");

}