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

$|++;

use version; our $VERSION = qv(0.0.1);

my $OUTFILE = 'output_curator.log';
open(STDERR, "| tee -i $OUTFILE") or die "WARNING: ERROR: Cannot open output.log\n";
open(STDOUT, "| tee -i $OUTFILE") or die "WARNING: ERROR: Cannot open output.log\n";

my $profout = 'report_curator.log';
open (PROUT,">$profout") or die "WARNING: ERROR: Cannot open report output\n";
binmode(PROUT, ":utf8");

my %store_hash = %{retrieve('extracted.results')}; # Load the extracted.results file from the sra_xml_parser.pl script.

# Initialize the development, anatomy, and controlled vocabulary hashes.
# Making these global hashes to be used in subroutines.
our %dv_hash; # development
our %bt_hash; # anatomy
our %cv_hash; # controlled vocabulary
our %tc_hash; # cell lines
our %gn_hash; # genes
our %sn_hash; # strains

our %successful_match_id; # A hash to keep track of previously successful match id's.
our %successful_match_output; # Successful match outputs.

# Initialize the output hash for metadata.
my %metadata_hash_output;

my %cat_HoA; # metadata catagories. This is extracted from a TSV where the first column is the titles from the metadata
# tsv file and the second column is which controlled vocabulary or data type to search. See 'categories.tsv'.

# An array of words used to indicate an unknown value in the metadata. TODO Load this from a file.
my @empty_words = ('mix','other','unknown');

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
        file        => 'metadata.csv',
        format      => 'hoh', # hash of hashes, which is default
        key         => 'RunSRR' 
} );

my $metadata_hash_ref = $hash_obj->all; # extract the entire file into a huge hash of hashes.

my %metadata_hash = %{$metadata_hash_ref}; # dereference the hash for ease of use.
undef $metadata_hash_ref; # remove the old reference.

# Initialize the user hash and Text::CSV object.
our %user_hash;
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

# Load the categories file into a HoA.
&parse_for_hash_of_arrays('categories.tsv', \%cat_HoA); # categories for metadata.

# Parse the OBO objects in the array and extract name/id information. Store into a hash.
&parse_obo(\@dv_terms, \%dv_hash); # development.
&parse_obo(\@bt_terms, \%bt_hash); # anatomy.
&parse_obo(\@cv_terms, \%cv_hash); # controlled vocabulary.

# Main classifying algorithm. Hold on to your butts.
# Iterate through our main metadata hash.
foreach my $key (sort keys %metadata_hash) {
	foreach my $subkey (keys $metadata_hash{$key}) { # Iterate through each SRR entry.
		my $entry_to_query = $metadata_hash{$key}{$subkey}; # The actual query ($subkey is the category name)
		# If our category is not listed as "none", our entry is not blank, and it doesn't match the empty words.
		if (!('none' ~~ @{$cat_HoA{$subkey}}) && ($entry_to_query ne '') && (!($entry_to_query ~~ @empty_words))) { 
			# Our main comparison search begins here.
			# We're at the bottom level of an entry within each category within each SRR entry (e.g. an entry in the original metadata tsv file).
			# To review, $key = the SRR entry. $subkey = a column from the metadata tsv. $cat_HoA = the categories tsv file.
			foreach my $search_type (@{$cat_HoA{$subkey}}) { # Perform a search with the query for each search type specified in the categories.tsv file (e.g. cell_line, tissue, etc.)
				# print PROUT "Searching $key\t$subkey\t$search_type\t$entry_to_query\n";
				my ($search_output, $id_output) = &main_search($entry_to_query, $search_type); # The main search.
				if ($search_output ne "FAILED") {
					print PROUT "Match found $key\t$subkey\t$search_type\t$entry_to_query => $search_output\t$id_output\n";
					$successful_match_output{$search_type}{$entry_to_query} = $search_output; # Build a successful match hash to speed up future searches.
					$successful_match_id{$search_type}{$entry_to_query} = $id_output; # Successful match hash for ids too.
					$metadata_hash_output{$key}{$subkey}{$search_type} = $search_output; # Set the output in the output hash.
					my $search_type_id = $search_type . '_id';
					$metadata_hash_output{$key}{$subkey}{$search_type_id} = $id_output; # Set the output id in the output hash.
				}
			}
		}
	}
}

# TODO
# Consolidation algorithm.
# Works through the main output hash and consolidates multiple entries into a single answer.

# Data Dumper commands for internal testing.
# print Data::Dumper->Dump([\%metadata_hash], ['*metadata_hash']);
print Data::Dumper->Dump([\%metadata_hash_output], ['*metadata_hash_output']);
# print Data::Dumper->Dump([\%cat_HoA], ['*cat_HoA']);
# print Data::Dumper->Dump([\%user_hash], ['*user_hash']);

sub parse_for_hash {
	my $user_list = $_[0];
	my $hash_ref = $_[1];

	open(my $data, '<:encoding(utf8)', $user_list) or die "Could not open '$user_list' $!\n";
	while (my $fields = $tsv_obj->getline( $data )) { # Why does tsv_obj work here if we don't explicitly pass it? TODO google it.
			$hash_ref->{$fields->[0]} = $fields->[1]; # Parse the two columns of the user list into a hash.
		}
	close $data;	
}

sub parse_for_hash_of_arrays {
	my $user_list = $_[0];
	my $hash_ref = $_[1];

	open(my $data, '<:encoding(utf8)', $user_list) or die "Could not open '$user_list' $!\n";
	while (my $line = <$data>){
		chomp($line);
		my @array = split(/\t/, $line);
		foreach my $entry (@array) {
			if ($entry ne $array[0]) {
				push @{$hash_ref->{$array[0]}} => $entry;
			}
		}
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

sub main_search {
	my $entry_to_query = $_[0];
	my $search_type = $_[1];
	my ($search_output, $id_output);

	# Check for user curated entries first.
	if ($user_hash{$entry_to_query}) {
		# print PROUT "Converted $entry_to_query => $user_hash{$entry_to_query}\n";
		$entry_to_query = $user_hash{$entry_to_query};
	}

	# Check for previously curated entries.
	if ($successful_match_output{$search_type}{$entry_to_query}){
		$search_output = $successful_match_output{$search_type}{$entry_to_query};
		$id_output = $successful_match_id{$search_type}{$entry_to_query};
		# print PROUT "Previous match found for $entry_to_query.\n";
		return ($search_output, $id_output);
	}

	# Sending out the queries to the relevant subroutines.
	if ($search_type eq "sample_type") {
		($search_output, $id_output) = &sample_type_search($entry_to_query, $search_type);
	} elsif ($search_type eq "stage") {
		($search_output, $id_output) = &stage_search($entry_to_query, $search_type);
	} elsif ($search_type eq "tissue") {
		($search_output, $id_output) = &tissue_search($entry_to_query, $search_type);
	} elsif ($search_type eq "cell_line") {
		($search_output, $id_output) = &cell_line_search($entry_to_query, $search_type);
	} elsif ($search_type eq "strain") {
		($search_output, $id_output) = &strain_search($entry_to_query, $search_type);
	} elsif ($search_type eq "genotype") {
		($search_output, $id_output) = &genotype_search($entry_to_query, $search_type);
	} elsif ($search_type eq "key_genes") {
		($search_output, $id_output) = &key_genes_search($entry_to_query, $search_type);
	} elsif ($search_type eq "sex") {
		($search_output, $id_output) = &sex_search($entry_to_query, $search_type);
	}

	return ($search_output, $id_output);
}

sub sample_type_search {
	my $entry_to_query = $_[0];
	my $search_type = $_[1];
	my ($search_output, $id_output);

	$search_output = 'FAILED';
	$id_output = 'FAILED';

	return ($search_output, $id_output);
}

sub stage_search {
	my $entry_to_query = $_[0];
	my $search_type = $_[1];
	my ($search_output, $id_output);

	$search_output = 'FAILED';
	$id_output = 'FAILED';

	return ($search_output, $id_output);
}

sub tissue_search {
	my $entry_to_query = $_[0];
	my $search_type = $_[1];
	my ($search_output, $id_output);

	$search_output = 'FAILED';
	$id_output = 'FAILED';

	return ($search_output, $id_output);
}

sub cell_line_search {
	my $entry_to_query = $_[0];
	my $search_type = $_[1];
	my ($search_output, $id_output);

	foreach my $cell_line (keys %tc_hash) {
		my $id = $tc_hash{$cell_line};
		my $lc_cell_line = lc($cell_line);
		my $lc_entry_to_query = lc($entry_to_query);
		my $score = distance($lc_entry_to_query, $lc_cell_line);
		if ($score < 1) {
			$search_output = $cell_line;
			$id_output = $id;
			print "$entry_to_query\t$search_output\n";
			return ($search_output, $id_output);
		} else {
			$search_output = 'FAILED';
			$id_output = 'FAILED';
		}
	}

	return ($search_output, $id_output);
}

sub strain_search {
	my $entry_to_query = $_[0];
	my $search_type = $_[1];
	my ($search_output, $id_output);

	$search_output = 'FAILED';
	$id_output = 'FAILED';

	return ($search_output, $id_output);
}

sub genotype_search {
	my $entry_to_query = $_[0];
	my $search_type = $_[1];
	my ($search_output, $id_output);

	$search_output = 'FAILED';
	$id_output = 'FAILED';

	return ($search_output, $id_output);
}

sub key_genes_search {
	my $entry_to_query = $_[0];
	my $search_type = $_[1];
	my ($search_output, $id_output);

	$search_output = 'FAILED';
	$id_output = 'FAILED';

	return ($search_output, $id_output);
}

sub sex_search {
	my $entry_to_query = $_[0];
	my $search_type = $_[1];
	my ($search_output, $id_output);

	# This is one of the most straightfoward searches. At least we hope.
	# TODO 
	# Separate query into individual words and search for male/female or equivalent.

	$entry_to_query = lc($entry_to_query); # change to lower case.
	if ($entry_to_query eq 'male' || $entry_to_query eq 'm') {
		$search_output = 'male';
		$id_output = 'FBcv:0000333';
	} elsif ($entry_to_query eq 'female' || $entry_to_query eq 'f') {
		$search_output = 'female';
		$id_output = 'FBcv:0000334';
	} elsif ($entry_to_query eq 'mated male' || $entry_to_query eq 'mated_male' || $entry_to_query eq 'mated-male') {
		$search_output = 'mated male';
		$id_output = 'FBcv:0000729';
	} elsif ($entry_to_query eq 'mated female' || $entry_to_query eq 'mated_female' || $entry_to_query eq 'mated-female') {
		$search_output = 'mated female';
		$id_output = 'FBcv:0000727';
	} elsif ($entry_to_query eq 'virgin male' || $entry_to_query eq 'virgin_male' || $entry_to_query eq 'virgin-male') {
		$search_output = 'virgin male';
		$id_output = 'FBcv:0000728';
	} elsif ($entry_to_query eq 'virgin female' || $entry_to_query eq 'virgin_female' || $entry_to_query eq 'virgin-female') {
		$search_output = 'virgin female';
		$id_output = 'FBcv:0000726';
	} else {
		$search_output = 'FAILED';
		$id_output = 'FAILED';
	}
	return ($search_output, $id_output);
}

sub consolidate_output {
	# A subroutine to consolidate multiple results in each category into a single result.
}

sub longest_element {
	# Taken from http://stackoverflow.com/questions/4182010/the-fastest-way-execution-time-to-find-the-longest-element-in-an-list
    my $max = -1;
    my $max_ref;
    for (@_) {
        if (length > $max) {  # no temp variable, length() twice is faster
            $max = length;
            $max_ref = \$_;   # avoid any copying
        }
    }
    $$max_ref
}