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
use Text::CSV::Hashify;
use Text::CSV;
use Text::Metaphone;
use Text::Levenshtein qw(distance);

use version; our $VERSION = qv(0.0.1);

my $OUTFILE = 'output_curator.log';
open(STDERR, "| tee -i $OUTFILE") or die "WARNING: ERROR: Cannot open output.log\n";
open(STDOUT, "| tee -i $OUTFILE") or die "WARNING: ERROR: Cannot open output.log\n";

my %store_hash = %{retrieve('extracted.results')}; # Load the extracted.results file from the sra_xml_parser.pl script.

# Initialize the development, anatomy, and controlled vocabulary hashes.
my %dv_hash; # development
my %bt_hash; # anatomy
my %cv_hash; # controlled vocabulary
my %tc_hash; # cell lines
my %gn_hash; # genes
my %sn_hash; # strains

# Create a new parser object.
my $parser = OBO::Parser::OBOParser->new;

# Parse the controlled vocabularies with the OBO parser object.
my $dv_ontology = $parser->work("fbdv-simple.obo");
my $bt_ontology = $parser->work("fbbt-simple.obo");
my $cv_ontology = $parser->work("fbcv-simple.obo");

# Store all the controlled vocabulary terms in arrays.
my @dv_terms = @{$dv_ontology->get_terms()};
my @bt_terms = @{$bt_ontology->get_terms()}; 
my @cv_terms = @{$cv_ontology->get_terms()}; 

# Initialize the metadata hash.
my $hash_obj = Text::CSV::Hashify->new( {
        file        => 'metadata2.csv',
        format      => 'hoh', # hash of hashes, which is default
        key         => 'RunSRR' 
} );

my $hash_ref = $hash_obj->all; # extract the entire file into a huge hash of hashes.

# Initialize the user hash and Text::CSV object.
my %user_hash;
my $tsv_obj = Text::CSV->new ( {
	sep_char => "\t",
	auto_diag => 1,
	binary => 1
});

# Load various tsv lists into hashes using a subroutine.
&parse_for_hash('user_curation.tsv', \%user_hash); # user curated entries.
&parse_for_hash('cell_lines.tsv', \%tc_hash); # cell lines.
&parse_for_hash('genes.tsv', \%gn_hash); # genes.
&parse_for_hash('strains.tsv', \%sn_hash); # strains.

# Parse the OBO objects in the array and extract name/id information. Store into a hash.
&parse_obo(\@dv_terms, \%dv_hash); # development.
&parse_obo(\@bt_terms, \%bt_hash); # anatomy.
&parse_obo(\@cv_terms, \%cv_hash); # controlled vocabulary.

# foreach my $key (keys $hash_ref) {
# 	print $key . "\n";
# 	my @stage_words;
# 	my $stage_magic = $hash_ref->{$key}->{'Stage (in Magic)'};
# 	if ($stage_magic) {
# 		print "Stage (in Magic)\t$stage_magic\n";
# 		push @stage_words => $stage_magic;
# 	}
# 	my $stage_zhenxia = $hash_ref->{$key}->{'Stage in Zhenxia'};
# 	if ($stage_zhenxia) {
# 		print "Stage in Zhenxia\t$stage_zhenxia\n";
# 		push @stage_words => $stage_zhenxia;

# 	}
# 	my $stage_nlm = $hash_ref->{$key}->{'Selected NLM annotated stage'};
# 	if ($stage_nlm) {
# 		print "Selected NLM annotated stage\t$stage_nlm\n";
# 		push @stage_words => $stage_nlm;
# 	}
	
# }


my %matched_hash; # hash for successful matches
my %updated_hash; # hash for updated store_hash values.
my %cv_id_hash;
# foreach my $key (keys %store_hash) {
# 	if ($store_hash{$key}) {
# 		my $query = $store_hash{$key};
# 		if ($matched_hash{$query}) { # if we've already matched this term.
# 			$updated_hash{$key} = $matched_hash{$query}; # set the value from the matched_hash.
# 			print "$query\t$updated_hash{$key}\t$cv_id_hash{$updated_hash{$key}}\n";
# 		} else {
# 			foreach my $entry (keys %cv_hash) {
# 				my $value = $cv_hash{$entry};
# 				my $score = distance($query, $value);
# 				if ($score < 2) {
# 					print "$query\t$value\t$entry\n";
# 					$matched_hash{$query} = $value;
# 					$cv_id_hash{$value} = $entry;
# 				}
# 			}
# 			foreach my $entry (keys %bt_hash) {
# 				my $value = $bt_hash{$entry};
# 				my $score = distance($query, $value);
# 				if ($score < 2) {
# 					print "$query\t$value\t$entry\n";
# 					$matched_hash{$query} = $value;
# 					$cv_id_hash{$value} = $entry;
# 				}
# 			}
# 			foreach my $entry (keys %dv_hash) {
# 				if ($dv_hash{$entry}) {
# 					my $value = $dv_hash{$entry};
# 					my $score = distance($query, $value);
# 					if ($score < 2) {
# 						print "$query\t$value\t$entry\n";
# 						$matched_hash{$query} = $value;
# 						$cv_id_hash{$value} = $entry;
# 					}
# 				}
# 			}
# 		}
# 	}
# }

my %transform_hash;

# Data Dumper commands for internal testing.
# print Data::Dumper->Dump([\%count_hash], ['count_hash']);
# print Data::Dumper->Dump([\%bt_hash], ['*bt_hash']);
# print Data::Dumper->Dump([\%cv_hash], ['*cv_hash']);
# print Data::Dumper->Dump([\%dv_hash], ['*dv_hash']);
# print Data::Dumper->Dump([\%store_hash], ['*store_hash']);
# print Data::Dumper->Dump([\%gn_hash], ['*gn_hash']);
# print Dumper $hash_ref;
# print Data::Dumper->Dump([\%main_hash], ['*main_hash']);
# print Data::Dumper->Dump([\%transform_hash], ['*transform_hash']);

sub parse_for_hash {
	my $user_list = $_[0];
	my $hash_ref = $_[1];

	open(my $data, '<:encoding(utf8)', $user_list) or die "Could not open '$user_list' $!\n";
	while (my $fields = $tsv_obj->getline( $data )) {
			$hash_ref->{$fields->[0]} = $fields->[1]; # Parse the two columns of the user list into a hash.
		}
	close $data;	
}

sub parse_obo {
	my $array_ref = $_[0];
	my $hash_ref = $_[1];
	
	foreach my $entry (@$array_ref) { 
		my $name = $entry->name();
		my $id = $entry->id();
		$hash_ref->{$id} = $name;
	}
}