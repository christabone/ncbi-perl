#!usr/bin/perl
#sra_xml_parser.pl
#A script to parse the XML fields of sra.xml files obtained from NCBI.
#Derived from the processing_from_ncbi.pl script used in-house at FlyBase.
#
#TO USE: This script will look for the folder "SRA" in the directory where it is placed.
#Please create this folder and fill it with SRA XML files to be parsed.
#
#Built for SRA XSD 1-6a: http://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/sra/doc/SRA_1-6a/
#
#Nodes and data retrieved:
#//EXPERIMENT_PACKAGE_SET/EXPERIMENT_PACKAGE
#---->/EXPERIMENT::accession (extracted)
#---->/EXPERIMENT::alias (extracted)
#---->/SAMPLE/SAMPLE_ATTRIBUTES
#-------->/SAMPLE_ATTRIBUTE
#------------>/TAG (extracted)
#------------>/VALUE (extracted)
use warnings;
use strict;
use Carp;
use Data::Dumper;
use XML::LibXML;
use Storable;

use version; our $VERSION = qv(0.0.1);

my %main_hash; # Hash to store all parsed XML data.

my $OUTFILE = 'output_parser.log';
open(STDERR, "| tee -i $OUTFILE") or die "WARNING: ERROR: Cannot open output.log\n";
open(STDOUT, "| tee -i $OUTFILE") or die "WARNING: ERROR: Cannot open output.log\n";

# Reading information from SRA files in SRA directory. We're looking for the directory "SRA" containing a list of XML files.
my $sra_directory = 'SRA';
opendir SRADIR, $sra_directory or croak "Can't open SRA directory: $!";
my @SRAfiles = grep {$_ ne q{.} && $_ ne q{..} } readdir SRADIR; # Read the SRA file names into an array.
closedir SRADIR;

foreach my $filename (@SRAfiles) {
	my $parser = XML::LibXML->new('1.0','UTF-8');
	my $dir_filename = "SRA/$filename"; # Add the directory name.
	my $sra_xml = $parser->parse_file($dir_filename) or croak "Cannot parse document: $@";
	foreach my $main_nodes ($sra_xml->findnodes('//EXPERIMENT_PACKAGE_SET/EXPERIMENT_PACKAGE')) {
		my $accession;
		foreach my $experiment_nodes ($main_nodes->findnodes('./EXPERIMENT')) {
			$accession = $experiment_nodes->getAttribute('accession');
			my $alias = $experiment_nodes->getAttribute('alias');
        	if(!($accession)) {
            	print "WARNING: No $accession found for $alias!\n";
        	}
			$main_hash{$accession}{alias} = $alias; # Assigning values to a two dimensional hash.
		}
	    foreach my $sample_attributes_node ($main_nodes->findnodes('./SAMPLE/SAMPLE_ATTRIBUTES')) {
	    	foreach my $individual_sample ($sample_attributes_node->findnodes('./SAMPLE_ATTRIBUTE')) {
	    	my $tag = $individual_sample->findvalue('./TAG');
	    	my $value = $individual_sample->findvalue('./VALUE');
	    	$main_hash{$accession}{$tag} = $value;
	    	}
	    }
	}
}

# Print messages and provide option to save output below:
print `clear`;
print "SRA XML parser $VERSION\n";
print scalar (keys %main_hash) . " entries extracted. Save to tab-delimited file \"extracted.tsv\"? (Y/n): ";
my $answer = <>;
chomp $answer;
if ($answer eq "" || $answer eq "Y") {
	open (TAB_OUTPUT, '>:encoding(UTF-8)', 'extracted.tsv') or die "Could not create file 'extracted.tsv' $!";
	foreach my $key (keys %main_hash) { # Print the SRA accession number.
		foreach my $hash_entry (keys %{$main_hash{$key}}){
			print TAB_OUTPUT "$key\t$hash_entry\t$main_hash{$key}{$hash_entry}\n"; # Print out the contents of the hash-of-hashes.
		}
	}
	close TAB_OUTPUT;
} else {
	print "\"n\" entered.\n";
}

print "Export storable file \"extracted.results\" for use with curation script? (Y/n): ";
my $answer2 = <>;
chomp $answer2;
if ($answer2 eq "" || $answer eq "Y") {
	store \%main_hash, 'extracted.results';
} else {
	print "\"n\" entered.\n";
}

close $OUTFILE;
exit;
# print Data::Dumper->Dump([\%main_hash], ['*main_hash']);