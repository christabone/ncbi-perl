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
# Making these global hashes to be used in subroutines.
our %dv_hash; # development
our %bt_hash; # anatomy
our %cv_hash; # controlled vocabulary
our %tc_hash; # cell lines
our %gn_hash; # genes
our %sn_hash; # strains

# Initialize the output hash for metadata.
my %metadata_hash_output;

my %cat_HoA; # metadata catagories. This is extracted from a TSV where the first column is the titles from the metadata
# tsv file and the second column is which controlled vocabulary or data type to search. See 'categories.tsv'.

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
		my $entry_to_query = $metadata_hash{$key}{$subkey}; # The final query ($subkey is the category name)
		if (!('none' ~~ @{$cat_HoA{$subkey}})) { # If our category is not listed as "none".
			# Our main comparison search begins here.
			# We're at the bottom level of an entry within each category within each SRR entry.
			foreach my $search_type (@{$cat_HoA{$subkey}}) {
				my ($search_output, $id_output) = &main_search($entry_to_query, $search_type); # The main search.
				&structure_the_output($subkey, $search_type, $search_output, $id_output, \%metadata_hash_output); # Sort and structure the output.
			}
		}
	}
}

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

# Data Dumper commands for internal testing.
# print Data::Dumper->Dump([\%metadata_hash], ['*metadata_hash']);
# print Data::Dumper->Dump([\%cat_HoA], ['*cat_HoA']);

sub parse_for_hash {
	my $user_list = $_[0];
	my $hash_ref = $_[1];

	open(my $data, '<:encoding(utf8)', $user_list) or die "Could not open '$user_list' $!\n";
	while (my $fields = $tsv_obj->getline( $data )) {
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

	# Sending out the queries to the relevant subroutines.
	if ($search_type eq "sample_type") {
		($search_output, $id_output) = &sample_type_search;
	} elsif ($search_type eq "stage") {
		($search_output, $id_output) = &stage_search;
	} elsif ($search_type eq "tissue") {
		($search_output, $id_output) = &tissue_search;
	} elsif ($search_type eq "cell_line") {
		($search_output, $id_output) = &cell_line_search;
	} elsif ($search_type eq "strain") {
		($search_output, $id_output) = &strain_search;
	} elsif ($search_type eq "genotype") {
		($search_output, $id_output) = &genotype_search;
	} elsif ($search_type eq "key_genes") {
		($search_output, $id_output) = &key_genes_search;
	} elsif ($search_type eq "sex") {
		($search_output, $id_output) = &sex_search;
	}

	return ($search_output, $id_output);
}

sub structure_the_output {
	my $original_search_category = $_[0];
	my $search_type = $_[1];
	my $search_output = $_[2];
	my $id = $_[3];
	my $hash_ref = $_[4];

	# Sorting the data into new categories based on the search parameters and original fields.
	# This data is stored on the Google Doc (Reannotation Metadata) shared with the NCBI group.

	if ($original_search_category eq )
}

sub sample_type_search {
	my $search_output = "development test";
	my $id_output = "development id";

	return ($search_output, $id_output);
}

sub stage_search {
	my $search_output = "stage test";
	my $id_output = "stage id";

	return ($search_output, $id_output);
}

sub tissue_search {
	my $search_output = "tissue test";
	my $id_output = "tissue id";

	return ($search_output, $id_output);
}

sub cell_line_search {
	my $search_output = "cell line test";
	my $id_output = "cell line id";

	return ($search_output, $id_output);
}

sub strain_search {
	my $search_output = "strain test";
	my $id_output = "strain id";

	return ($search_output, $id_output);
}

sub genotype_search {
	my $search_output = "genotype test";
	my $id_output = "genotype id";

	return ($search_output, $id_output);
}

sub key_genes_search {
	my $search_output = "key_genes test";
	my $id_output = "key_genes id";

	return ($search_output, $id_output);
}

sub sex_search {
	my $search_output = "sex test";
	my $id_output = "sex id";

	return ($search_output, $id_output);
}