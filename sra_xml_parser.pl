#!usr/bin/perl
#sra_xml_parser.pl
#A script to parse the XML fields of sra.xml files obtained from NCBI.
#Derived from the processing_from_ncbi.pl script used in-house at FlyBase.
#Version 0.1
#
#Built for SRA XSD 1-6a: http://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/sra/doc/SRA_1-6a/

use warnings;
use strict;
use Carp;
use Data::Dumper;
use XML::LibXML;

my %main_hash; # Hash to store all parsed XML data.

my $OUTFILE = "output.log";
open(STDERR, "| tee -i $OUTFILE") or die "WARNING: ERROR: Cannot open output.log\n";
open(STDOUT, "| tee -i $OUTFILE") or die "WARNING: ERROR: Cannot open output.log\n";

# Reading information from SRA files in SRA directory. We're looking for the directory "SRA" containing a list of XML files.
my $sra_directory = "SRA";
opendir SRADIR, $sra_directory or die "Can't open SRA directory: $!";
my @SRAfiles = grep {$_ ne '.' && $_ ne '..' } readdir SRADIR; # Read the SRA file names into an array.
closedir SRADIR;

foreach my $filename (@SRAfiles) {
	my $parser = XML::LibXML->new("1.0","UTF-8");
	my $dir_filename = "SRA/" . $filename; # Add the directory name.
	my $sra_xml = $parser->parse_file($dir_filename) or die "Cannot parse document: $@";
	foreach my $MainNodes ($sra_xml->findnodes('//EXPERIMENT_PACKAGE_SET/EXPERIMENT_PACKAGE')) {
		foreach my $ExperimentNode ($MainNodes->findnodes('./EXPERIMENT')) {
			my $accession = $ExperimentNode->getAttribute('accession');
        	if(!($accession)) {
            	print "WARNING: No $accession found for $alias!\n";
        	}
        	my $alias = $ExperimentNode->getAttribute('alias');
			$main_hash{$accession}{alias} = $alias; # Assigning values to a two dimensional hash. 
		}
	}
}

print Data::Dumper->Dump([\%main_hash], ['*main_hash']);