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
# use Text::Compare;
use Text::CSV_XS;
use Text::Metaphone;
use Text::LevenshteinXS qw(distance);

use version; our $VERSION = qv(0.0.1);

my $OUTFILE = 'output_curator.log';
open(STDERR, "| tee -i $OUTFILE") or die "WARNING: ERROR: Cannot open output.log\n";
open(STDOUT, "| tee -i $OUTFILE") or die "WARNING: ERROR: Cannot open output.log\n";

# my %store_hash = %{retrieve('extracted.results')}; # Load the extracted.results file from the sra_xml_parser.pl script.
my %dv_hash;
my %bt_hash;
my %cv_hash;

my $parser = OBO::Parser::OBOParser->new;

my $dv_ontology = $parser->work("fbdv.obo");
my $bt_ontology = $parser->work("fbbt.obo");
my $cv_ontology = $parser->work("fbcv-simple.obo");

my @dv_terms = @{$dv_ontology->get_terms()}; # Store the terms in an array.
my @bt_terms = @{$bt_ontology->get_terms()}; 
my @cv_terms = @{$cv_ontology->get_terms()}; 

my $tsv = Text::CSV_XS->new({
  sep_char => "\t",
  binary => 1,
  auto_diag => 1,
  eol => $\
});

my %store_hash;
open(my $tsv_file, "<:encoding(UTF-8)", "metadata.tsv") || croak "can't open \"metadata.tsv\": $!";
while ( my $row = $tsv->getline( $tsv_file ) )
{
    $store_hash{$row->[0]} = $row->[25]; # read into a hash.
}
close($tsv_file);

foreach my $entry (@dv_terms) { # Obtain names and ids. Store into a hash.
	my $name = $entry->name();
	my $id = $entry->id();
	$dv_hash{$id} = $name;
}

foreach my $entry (@bt_terms) {
	my $name = $entry->name();
	my $id = $entry->id();
	$bt_hash{$id} = $name;
}

foreach my $entry (@cv_terms) {
	my $name = $entry->name();
	my $id = $entry->id();
	$cv_hash{$id} = $name;
}

my %count_keys_hash;
my %count_values_hash;

# foreach my $key (keys %store_hash) {
# 	foreach my $hash_entry (keys %{$store_hash{$key}}){
# 		$count_keys_hash{$hash_entry}++;
# 	}
# }

my $value_total;
my $value_match;
my %matched_hash; # hash for successful matches
my %updated_hash; # hash for updated store_hash values.
my %cv_id_hash;
foreach my $key (keys %store_hash) {
	if ($store_hash{$key}) {
		my $query = $store_hash{$key};
		if ($matched_hash{$query}) { # if we've already matched this term.
			$updated_hash{$key} = $matched_hash{$query}; # set the value from the matched_hash.
			print "$query\t$updated_hash{$key}\t$cv_id_hash{$updated_hash{$key}}\n";
		} else {
			foreach my $entry (keys %cv_hash) {
				my $value = $cv_hash{$entry};
				my $score = distance($query, $value);
				if ($score < 2) {
					print "$query\t$value\t$entry\n";
					$matched_hash{$query} = $value;
					$cv_id_hash{$value} = $entry;
				}
			}
			foreach my $entry (keys %bt_hash) {
				my $value = $bt_hash{$entry};
				my $score = distance($query, $value);
				if ($score < 2) {
					print "$query\t$value\t$entry\n";
					$matched_hash{$query} = $value;
					$cv_id_hash{$value} = $entry;
				}
			}
			foreach my $entry (keys %dv_hash) {
				if ($dv_hash{$entry}) {
					my $value = $dv_hash{$entry};
					my $score = distance($query, $value);
					if ($score < 2) {
						print "$query\t$value\t$entry\n";
						$matched_hash{$query} = $value;
						$cv_id_hash{$value} = $entry;
					}
				}
			}
		}
	}
}
	# 	foreach my $key2 (keys %bt_hash) {
	# 		my $value2 = $bt_hash{$key2};
	# 		print "$value\t$value2\n";
	# 		my $dist_calc = distance($value, $value2);
	# 		if ($dist_calc < 2) {
	# 			$count_values_hash{$value}++;
	# 			$value_match++;
	# 			print "$value\t$value2\t$dist_calc"
	# 		}

	# 	}
	# }


# foreach my $key (keys %count_values_hash) {
# 	if ($count_values_hash{$key} > 1) {
# 		print $key . "\t" . $count_values_hash{$key} . "\n";
# 	}
# }

# print "Total values queries: $value_total" . "\n";
# print "Total values matched: $value_match" . "\n";
# print "Precent matched: " . sprintf("%.2f", (($value_match/$value_total) * 100)) . "%" . "\n";

my %transform_hash;

# print Data::Dumper->Dump([\%count_hash], ['count_hash']);
# print Data::Dumper->Dump([\%bt_hash], ['*bt_hash']);
# print Data::Dumper->Dump([\%cv_hash], ['*cv_hash']);
# print Data::Dumper->Dump([\%dv_hash], ['*dv_hash']);
# print Data::Dumper->Dump([\%store_hash], ['*store_hash']);
# print Data::Dumper->Dump([\%transform_hash], ['*transform_hash']);