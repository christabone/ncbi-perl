#!usr/bin/perl
#sra_curator.pl
#A script to curate fields derived from SRA entries obtained from NCBI.
#Derived from the processing_from_ncbi.pl script used in-house at FlyBase.
#
#TO USE: This script will look for the file "extracted.results" from the script sra_xml_parser.pl in the directory where it is placed.
#
use warnings;
use strict;
use Carp;
use Data::Dumper;
use Storable;
use OBO::Parser::OBOParser;
use Text::Metaphone;
use Text::Levenshtein qw(distance);

use version; our $VERSION = qv(0.0.1);

my $OUTFILE = 'output_curator.log';
open(STDERR, "| tee -i $OUTFILE") or die "WARNING: ERROR: Cannot open output.log\n";
open(STDOUT, "| tee -i $OUTFILE") or die "WARNING: ERROR: Cannot open output.log\n";

my %main_hash = %{retrieve('extracted.results')}; # Load the extracted.results file from the sra_xml_parser.pl script.
my %dv_hash;
my %bt_hash;

my $parser = OBO::Parser::OBOParser->new;

my $dv_ontology = $parser->work("fbdv.obo");
my $bt_ontology = $parser->work("fbbt.obo");

my @dv_terms = @{$dv_ontology->get_terms()}; # Store the terms in an array.
my @bt_terms = @{$bt_ontology->get_terms()}; # Store the terms in an array.

foreach my $entry (@dv_terms) { # Obtain names and ids. Store into a hash.
	my $name = $entry->name();
	my $id = $entry->id();
	$dv_hash{$id} = $name;
}

foreach my $entry (@bt_terms) { # Obtain names and ids. Store into a hash.
	my $name = $entry->name();
	my $id = $entry->id();
	$bt_hash{$id} = $name;
}

my %count_keys_hash;
my %count_values_hash;

foreach my $key (keys %main_hash) {
	foreach my $hash_entry (keys %{$main_hash{$key}}){
		$count_keys_hash{$hash_entry}++;
	}
}

my $value_total;
my $value_match;
foreach my $key (keys %main_hash) {
	foreach my $sub_key (keys %{$main_hash{$key}}) {
		$value_total++;
		my $value = %{$main_hash{$key}}->{$sub_key};
		foreach my $key2 (keys %bt_hash) {
			my $value2 = $bt_hash{$key2};
			print "$value\t$value2\n";
			my $dist_calc = distance($value, $value2);
			if ($dist_calc < 2) {
				$count_values_hash{$value}++;
				$value_match++;
				print "$value\t$value2\t$dist_calc"
			}

		}
	}
}

foreach my $key (keys %count_values_hash) {
	if ($count_values_hash{$key} > 1) {
		print $key . "\t" . $count_values_hash{$key} . "\n";
	}
}

print "Total values queries: $value_total" . "\n";
print "Total values matched: $value_match" . "\n";
print "Precent matched: " . sprintf("%.2f", (($value_match/$value_total) * 100)) . "%" . "\n";

my %transform_hash;

# print Data::Dumper->Dump([\%count_hash], ['count_hash']);
# print Data::Dumper->Dump([\%bt_hash], ['*bt_hash']);
# print Data::Dumper->Dump([\%dv_hash], ['*dv_hash']);
# print Data::Dumper->Dump([\%main_hash], ['*main_hash']);
# print Data::Dumper->Dump([\%transform_hash], ['*transform_hash']);