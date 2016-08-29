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
use Text::Levenshtein::XS qw(distance);

require '/users/ctabone/Programming/git/ncbi-perl/sra_curator_support.pm';

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
our %dv_hash; # development
our %bt_hash; # anatomy
our %cv_hash; # controlled vocabulary
our %tc_hash; # cell lines
our %gn_hash; # genes
our %sn_hash; # strains

our %statistics; # Tracking match statistics.

our $skip_searches = 0; # For testing sorting and other late algorithms.

our %successful_match_id; # A hash to keep track of previously successful match id's.
our %successful_match_output; # Successful match outputs.
our %failed_match_output; # A hash to keep track of failed match attempts.

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
our %dictionary_HoH;
my $tsv_obj = Text::CSV->new ( {
	sep_char => "\t",
	auto_diag => 1,
	binary => 1
});

# Load various tsv lists into hashes using a subroutine.
&parse_for_hash('cell_lines.tsv', \%tc_hash); # cell lines.
&parse_for_hash('genes.tsv', \%gn_hash); # genes.
&parse_for_hash('strains.tsv', \%sn_hash); # strains.

# Load the categories file into a HoA.
&parse_for_hash_of_arrays('categories.tsv', \%cat_HoA); # categories for metadata.

# Load the dictionary file into a HoH.
&parse_for_hash_of_hashes('curation_dictionary.tsv', \%dictionary_HoH); # user curated entries.

# Parse the OBO objects in the array and extract name/id information. Store into a hash.
&parse_obo(\@dv_terms, \%dv_hash); # development.
&parse_obo(\@bt_terms, \%bt_hash); # anatomy.
&parse_obo(\@cv_terms, \%cv_hash); # controlled vocabulary.

my @object_array;

# Storing all metadata as objects.
foreach my $key (sort keys %metadata_hash) {
	my $objectkey = sra_curator_support->new($key);
	foreach my $subkey (keys $metadata_hash{$key}) {
		$objectkey->set_old_metadata($subkey, $metadata_hash{$key}{$subkey});
	}
	push @object_array => $objectkey;
}

# Main classifying algorithm. Hold on to your butts.
# Iterate through our objects.
foreach my $object (@object_array) { # Load an object. 
	my $object_hash_ref = $object->get_entire_old_metadata(); # Pull out the old metadata for that object.
	foreach my $key (keys % {$object_hash_ref}) { # For each column title for an object.
		my $query = $object_hash_ref -> {$key}; # Get the value, this is the actual query to use. 
		if (!('none' ~~ @{$cat_HoA{$key}}) && ($query ne '') && (!($query ~~ @empty_words))) {
		# If our category is not listed as "none", our entry is not blank, and it doesn't match the empty words.
			foreach my $search_type (@{$cat_HoA{$key}}) { # Perform a search with the query for each search type specified in the categories.tsv file (e.g. cell_line, tissue, etc.)
				if (!($failed_match_output{$search_type}{$query})) { # If it's failed before, we skip the search - hugely important, speeds up everything significantly.
					my ($search_output, $id_output) = &main_search($query, $search_type, $object); # The main search.
					# Once the search is returned, we need to store it in a more complicated manner.
					# We need to sort by search_type and further by the original column header.
					# This is because we need to resolve the results further (i.e. if 3 answers come back for cell_line via 3 searches, which is best?).
					$object->set_new_results($search_type, $key, $search_output, $id_output); # Sets the search results (stores FAILED as well).
					if ($search_output ne "FAILED") {
						$successful_match_output{$search_type}{$query} = $search_output; # Build a successful match hash to speed up future searches.
						$successful_match_id{$search_type}{$query} = $id_output; # Successful match hash for ids too.
					} else {
						$failed_match_output{$search_type}{$query}++;
					}
				}
			}	
		}
	}
}

# Print close matches for curator intervention.
foreach my $object (@object_array) { # Load an object.
	my $accession = $object->get_accession();
	my $hash_ref = $object->get_all_scores(); # Returns a HoH to of score results. First keys are search types. Second keys are literal search terms. Third is compared value. Final value is a score.
	if (defined $hash_ref) {
		foreach my $type (keys %{$hash_ref}){ # For each search type.
			foreach my $query (keys %{$hash_ref}->{$type}) { # For each literal search term.
				my $has_been_matched = $object->check_if_matched($type, $query); # Find out whether this score was eventually matched.
				if (!(defined $has_been_matched)) { # If it hasn't been matched, it may have some close matches to check for curation.
					foreach my $compare (keys %{$hash_ref}->{$type}{$query}) {
						my $score = $hash_ref->{$type}{$query}{$compare}; # Grab the score.
						if ($score <= 4) {
							print PROUT "CLOSE MATCH $accession ($type, $score) FLYBASE: $compare\t QUERY: $query\n";
						}
					}
				} else {
					# Do nothing, move onto the next entry.
				}
			}
		}
	}
}

# Consolidate the results.
print "Consolidation\n";
print "-------------\n";
foreach my $object (@object_array) { # Load an object. Consolidate via names (ids are also used in the subroutine).
	&consolidate_output($object, 'sex');
	&consolidate_output($object, 'cell_line');
	&consolidate_output($object, 'stage');
	&consolidate_output($object, 'strain');
	&consolidate_output($object, 'tissue');
	&consolidate_output($object, 'sample_type');
}

# Print statistics
&statistics('sex');
&statistics('cell_line');
&statistics('stage');
&statistics('strain');
&statistics('tissue');
&statistics('sample_type');

# Output to TSV
my $tsv_out = 'final_output.tsv';
open (TSVOUT,">$tsv_out") or die "WARNING: ERROR: Cannot open tsv output\n";
binmode(TSVOUT, ":utf8");
print TSVOUT "ID\tsample_type\tsample_type_id\tsex\tsex_id\tcell_line\tcell_line_id\tstage\tstage_id\tstrain\tstrain_id\ttissue\ttissue_id\n";
foreach my $object (@object_array) { 
	my $ID = $object->get_accession();
	my ($sample_type, $sample_type_id) = $object->get_consolidated_results('sample_type');
	my ($sex, $sex_id) = $object->get_consolidated_results('sex');
	my ($cell_line, $cell_line_id) = $object->get_consolidated_results('cell_line');
	my ($stage, $stage_id) = $object->get_consolidated_results('stage');
	my ($strain, $strain_id) = $object->get_consolidated_results('strain');
	my ($tissue, $tissue_id) = $object->get_consolidated_results('tissue');
	print TSVOUT "$ID\t";
	if (defined $sample_type) { print TSVOUT "$sample_type\t"; } else { print TSVOUT "\t"; }
	if (defined $sample_type_id) { $sample_type_id =~ s/://g; print TSVOUT "$sample_type_id\t"; } else { print TSVOUT "\t"; }
	if (defined $sex) { print TSVOUT "$sex\t"; } else { print TSVOUT "\t"; }
	if (defined $sex_id) { $sex_id =~ s/://g; print TSVOUT "$sex_id\t"; } else { print TSVOUT "\t"; }
	if (defined $cell_line) { print TSVOUT "$cell_line\t"; } else { print TSVOUT "\t"; }
	if (defined $cell_line_id) { $cell_line_id =~ s/://g; print TSVOUT "$cell_line_id\t"; } else { print TSVOUT "\t"; }
	if (defined $stage) { print TSVOUT "$stage\t"; } else { print TSVOUT "\t"; }
	if (defined $stage_id) { $stage_id =~ s/://g; print TSVOUT "$stage_id\t"; } else { print TSVOUT "\t"; }
	if (defined $strain) { print TSVOUT "$strain\t"; } else { print TSVOUT "\t"; }
	if (defined $strain_id) { $strain_id =~ s/://g; print TSVOUT "$strain_id\t"; } else { print TSVOUT "\t"; }
	if (defined $tissue) { print TSVOUT "$tissue\t"; } else { print TSVOUT "\t"; }
	if (defined $tissue_id) { $tissue_id =~ s/://g; print TSVOUT "$tissue_id\t"; } else { print TSVOUT "\t"; }
	print TSVOUT "\n";
}
close TSVOUT;

# Data Dumper commands for internal testing.
# print Data::Dumper->Dump([\%metadata_hash], ['*metadata_hash']);
# print Data::Dumper->Dump([\%metadata_hash_output], ['*metadata_hash_output']);
# print Data::Dumper->Dump([\%cat_HoA], ['*cat_HoA']);
# print Data::Dumper->Dump([\%dv_hash], ['*dv_hash']);
# print Data::Dumper->Dump([\%statistics], ['*statistics']);
# print Data::Dumper->Dump([\%successful_match_output], ['*successful_match_output']);

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

sub parse_for_hash_of_hashes {
	my $user_list = $_[0];
	my $hash_ref = $_[1];

	open(my $data, '<:encoding(utf8)', $user_list) or die "Could not open '$user_list' $!\n";
	while (my $fields = $tsv_obj->getline( $data )) {
			$hash_ref->{$fields->[0]}{$fields->[1]} = $fields->[2]; # Parse the two columns of the user list into a hash.
		}
	close $data;	
}

sub parse_obo {
	my $array_ref = $_[0];
	my $hash_ref = $_[1];
	
	foreach my $entry (@$array_ref) { 
		my $name = $entry->name();
		my $id = $entry->id();
		$hash_ref->{$name} = $id;
	}
}

sub main_search {
	my $entry_to_query = $_[0]; # Used for the search, may be altered by the dictionary check.
	my $original_entry = $_[0]; # Used for the previously curated entry search. Not modified by dictionary check.
	my $search_type = $_[1];
	my $object =$_[2];
	my ($search_output, $id_output);

	# Skip all searches and return now.
	if ($skip_searches == 1){
		my ($search_output, $id_output) = "FAILED";
		return ($search_output, $id_output);
	}

	# Check for previously curated entries.
	if (exists $successful_match_output{$search_type}{$original_entry}){
		$search_output = $successful_match_output{$search_type}{$original_entry};
		$id_output = $successful_match_id{$search_type}{$original_entry};
		# print "Previous match found for $entry_to_query.\n";
		return ($search_output, $id_output);
	}

	# Check for user dictionary entries.
	if ($dictionary_HoH{$search_type}->{$entry_to_query}) {
		# print "Converted $entry_to_query => $dictionary_HoH{$search_type}->{$entry_to_query}\n";
		$entry_to_query = $dictionary_HoH{$search_type}->{$entry_to_query};
	}

	# Sending out the queries to the relevant subroutines.
	if ($search_type eq "sample_type") {
		($search_output, $id_output) = &sample_type_search($entry_to_query, $search_type, $object);
	} elsif ($search_type eq "stage") {
		($search_output, $id_output) = &stage_search($entry_to_query, $search_type, $object);
	} elsif ($search_type eq "tissue") {
		($search_output, $id_output) = &tissue_search($entry_to_query, $search_type, $object);
	} elsif ($search_type eq "cell_line") {
		($search_output, $id_output) = &cell_line_search($entry_to_query, $search_type, $object);
	} elsif ($search_type eq "strain") {
		($search_output, $id_output) = &strain_search($entry_to_query, $search_type, $object);
	} elsif ($search_type eq "genotype") {
		($search_output, $id_output) = &genotype_search($entry_to_query, $search_type, $object);
	} elsif ($search_type eq "key_genes") {
		($search_output, $id_output) = &key_genes_search($entry_to_query, $search_type, $object);
	} elsif ($search_type eq "sex") {
		($search_output, $id_output) = &sex_search($entry_to_query, $search_type, $object);
	}

	return ($search_output, $id_output);
}

sub sample_type_search {
	my $entry_to_query = $_[0];
	my $search_type = $_[1];
	my $object = $_[2];
	my $type = 'sample_type';
	my ($search_output, $id_output, $score);

	foreach my $sample_type (keys %cv_hash) {
		my $id = $cv_hash{$sample_type};
		my $lc_sample_type = lc($sample_type);
		my $lc_entry_to_query = lc($entry_to_query);
		$score = distance($lc_entry_to_query, $lc_sample_type);
		if ($score == 0) {
			$search_output = $sample_type;
			$id_output = $id;
			$object->set_score($type, $entry_to_query, $sample_type, $score);
			$object->set_match($type, $entry_to_query, $score); # Flag the match as a success.
			return ($search_output, $id_output);
		} else {
			$object->set_score($type, $entry_to_query, $sample_type, $score); # Record all the failed scores for curation purposes later.
			$search_output = 'FAILED';
			$id_output = 'FAILED';
		}
	}
	return ($search_output, $id_output);
}

sub stage_search {
	my $entry_to_query = $_[0];
	my $search_type = $_[1];
	my $object = $_[2];
	my $type = 'stage';
	my ($search_output, $id_output, $score);

	foreach my $stage (keys %dv_hash) {
		my $id = $dv_hash{$stage};
		my $lc_stage = lc($stage);
		my $lc_entry_to_query = lc($entry_to_query);
		$score = distance($lc_entry_to_query, $lc_stage);
		if ($score == 0) {
			$search_output = $stage;
			$id_output = $id;
			$object->set_score($type, $entry_to_query, $stage, $score);
			$object->set_match($type, $entry_to_query, $score); # Flag the match as a success.
			return ($search_output, $id_output);
		} else {
			$object->set_score($type, $entry_to_query, $stage, $score); # Record all the failed scores for curation purposes later.
			$search_output = 'FAILED';
			$id_output = 'FAILED';
		}
	}
	return ($search_output, $id_output);
}

sub tissue_search {
	my $entry_to_query = $_[0];
	my $search_type = $_[1];
	my $object = $_[2];
	my $type = 'tissue';
	my ($search_output, $id_output, $score);

	foreach my $tissue (keys %bt_hash) {
		my $id = $bt_hash{$tissue};
		my $lc_tissue = lc($tissue);
		my $lc_entry_to_query = lc($entry_to_query);
		$score = distance($lc_entry_to_query, $lc_tissue);
		if ($score == 0) {
			$search_output = $tissue;
			$id_output = $id;
			print PROUT "$entry_to_query\t$search_output\n";
			$object->set_score($type, $entry_to_query, $tissue, $score);
			$object->set_match($type, $entry_to_query, $score); # Flag the match as a success.
			return ($search_output, $id_output);
		} else {
			$object->set_score($type, $entry_to_query, $tissue, $score); # Record all the failed scores for curation purposes later.
			$search_output = 'FAILED';
			$id_output = 'FAILED';
		}
	}
	return ($search_output, $id_output);
}

sub cell_line_search {
	my $entry_to_query = $_[0];
	my $search_type = $_[1];
	my $object = $_[2];
	my $type = 'cell_line';
	my ($search_output, $id_output, $score);

	foreach my $cell_line (keys %tc_hash) {
		my $id = $tc_hash{$cell_line};
		my $lc_cell_line = lc($cell_line);
		my $lc_entry_to_query = lc($entry_to_query);
		$score = distance($lc_entry_to_query, $lc_cell_line);
		if ($score == 0) {
			$search_output = $cell_line;
			$id_output = $id;
			print PROUT "$entry_to_query\t$search_output\n";
			$object->set_score($type, $entry_to_query, $cell_line, $score);
			$object->set_match($type, $entry_to_query, $score); # Flag the match as a success.
			return ($search_output, $id_output);
		} else {
			$search_output = 'FAILED';
			$id_output = 'FAILED';
			$object->set_score($type, $entry_to_query, $cell_line, $score); # Record all the failed scores for curation purposes later.

		}
	}
	return ($search_output, $id_output);
}

sub strain_search {
	my $entry_to_query = $_[0];
	my $search_type = $_[1];
	my $object = $_[2];
	my $type = 'strain';
	my ($search_output, $id_output, $score);

	foreach my $strain (keys %sn_hash) {
		my $id = $sn_hash{$strain};
		my $lc_strain = lc($strain);
		my $lc_entry_to_query = lc($entry_to_query);
		$score = distance($lc_entry_to_query, $lc_strain);
		if ($score == 0) {
			$search_output = $strain;
			$id_output = $id;
			print PROUT "$entry_to_query\t$search_output\n";
			$object->set_score($type, $entry_to_query, $strain, $score); # Record the match.
			$object->set_match($type, $entry_to_query, $score); # Flag the match as a success.
			return ($search_output, $id_output);
		} else {
			$search_output = 'FAILED';
			$id_output = 'FAILED';
			$object->set_score($type, $entry_to_query, $strain, $score); # Record all the failed scores for curation purposes later.

		}
	}
	return ($search_output, $id_output);
}

sub genotype_search {
	my $entry_to_query = $_[0];
	my $search_type = $_[1];
	my $object = $_[2];
	my ($search_output, $id_output);

	$search_output = 'FAILED';
	$id_output = 'FAILED';

	return ($search_output, $id_output);
}

sub key_genes_search {
	my $entry_to_query = $_[0];
	my $search_type = $_[1];
	my $object = $_[2];
	my ($search_output, $id_output);

	$search_output = 'FAILED';
	$id_output = 'FAILED';

	return ($search_output, $id_output);
}

sub sex_search {
	my $entry_to_query = $_[0];
	my $search_type = $_[1];
	my $object = $_[2];
	my ($search_output, $id_output);

	# This is one of the most straightfoward searches. At least we hope.
	# TODO 
	# Separate query into individual words and search for male/female or equivalent.

	$entry_to_query = lc($entry_to_query); # change to lower case.
	if ($entry_to_query eq 'male' || $entry_to_query eq 'm') {
		$search_output = 'male';
		$id_output = 'FBcv0000333';
	} elsif ($entry_to_query eq 'female' || $entry_to_query eq 'f') {
		$search_output = 'female';
		$id_output = 'FBcv0000334';
	} elsif ($entry_to_query eq 'mated male' || $entry_to_query eq 'mated_male' || $entry_to_query eq 'mated-male') {
		$search_output = 'mated male';
		$id_output = 'FBcv0000729';
	} elsif ($entry_to_query eq 'mated female' || $entry_to_query eq 'mated_female' || $entry_to_query eq 'mated-female') {
		$search_output = 'mated female';
		$id_output = 'FBcv0000727';
	} elsif ($entry_to_query eq 'virgin male' || $entry_to_query eq 'virgin_male' || $entry_to_query eq 'virgin-male') {
		$search_output = 'virgin male';
		$id_output = 'FBcv0000728';
	} elsif ($entry_to_query eq 'virgin female' || $entry_to_query eq 'virgin_female' || $entry_to_query eq 'virgin-female') {
		$search_output = 'virgin female';
		$id_output = 'FBcv0000726';
	} else {
		$search_output = 'FAILED';
		$id_output = 'FAILED';
	}
	# TODO set object score for sex.
	return ($search_output, $id_output);
}

sub consolidate_output {
	# A subroutine to consolidate the output based on "majority wins" (with an exception for 'FAILED' results).
	my $object = $_[0];
	my $type = $_[1];
	my $type_id = "$type" . '_id';
	my $hash_ref = $object->get_new_results($type_id); # Get all the id results for that object.
	my $hash_ref_entries = $object->get_new_results($type); # Get all the actual entries (which correspond to the ids above).
	my %hash_entry_id_xref;
	my $number_of_keys = scalar (keys %{$hash_ref}); # Get all the keys for the id results (the keys are the column names that were searched).
	if ($number_of_keys != 0) { # If the number of keys are not empty. In other words, only compute results were the subject (e.g. sex) was searched.
		my %check_hash; # Create a hash to check the results between different between columns.
		foreach my $key (keys %{$hash_ref}) {
			my $output_id = $hash_ref->{$key}; # Grab the output_id from the search.
			my $output_entry = $hash_ref_entries->{$key};
			$check_hash{$output_id}++; # Store the id in the check_hash.
			$hash_entry_id_xref{$output_id} = $output_entry # A cross-reference has to associate id's with their respective entries.
		}
		my $check_keys = scalar (keys %check_hash);
		# Count the number of different id numbers we stored in the check_hash.
		# This is the most important part of the check, so it warrants an explanantion.
		# We're basically grabbing all the id's that have been assigned to each column from the table for a certain trait (e.g. sex).
		# If all the id's are the same, then when you put them in a hash table, they all create one key.
		# If the id's are different, then you'll get more than one key created, and you'll have a discrepency that needs to be resolved.
		my $accession = $object->get_accession();
		if ($check_keys == 1) { 
			my $entry_value = (%check_hash)[0]; # Hash slice to grab first value of hash entry.
			if ($entry_value eq 'FAILED'){
				print PROUT "SINGLE OUTPUT (FAILED) => $accession\t$type\n"; # Note, we don't store FAILED outputs as consolidated.
				$statistics{$type}{'failed'}++;

			} else {
				my $final_output_id = (keys %check_hash)[0]; # Hash slice. Only 1 key. This is the output_id.
			    $object->set_consolidated_results($type, $hash_entry_id_xref{$final_output_id}, $final_output_id); # Storing the Search type, Output, and ID output.
			    print PROUT "SINGLE OUTPUT => $accession\t$type\t$final_output_id\n";
			    $statistics{$type}{'successful'}++;
			}
		} elsif ($check_keys > 1) {
			print PROUT "ATTEMPTING CONSOLIDATION for $accession $type\n";
			my $largest = 0;
			my $second_largest = 0;
			my $largest_id = "";
			my $second_largest_id = "";
			foreach my $entry (keys %check_hash) {
				my $to_compare = $check_hash{$entry};
				if (($entry ne 'FAILED') && ($to_compare >= $largest)){ # Importantly, we ignore cases where the entry is 'FAILED'. Otherwise a double FAILED group will take priority over a single successful group.
					$second_largest_id = $largest_id;
					$largest_id = $entry;
					$second_largest = $largest;
					$largest = $to_compare;
				}
			}
			if ($largest == $second_largest) { # If there are equal amounts of the first two entries (no clear majority) then we go with a choice from Magic (if it exists).
				# This ONLY works if a column named "Magic" is used once-per-category. Having "Magic column 1" and "Magic column 2" both for "sex" will break this!
				my $is_magic = 'false';
				foreach my $key (keys %{$hash_ref}) { # Go back through the original hash_ref look where keys are column names.
					my $lc_key = lc($key);
					if (index($lc_key, 'magic') != -1) { # Check if the word 'magic' is within the column name.
						my $final_output_id = $hash_ref->{$key};
						$object->set_consolidated_results($type, $hash_entry_id_xref{$final_output_id}, $final_output_id); # Storing the Search type, Output, and ID output.
						print PROUT "CONSOLIDATED (MAGIC PRIORITY) => $accession\t$type\t$final_output_id\n";
						$statistics{$type}{'successful'}++;
						$is_magic = 'true';
					}
				}
				if ($is_magic eq 'false') { # If we never found any magic entries in the last section.
					print PROUT "NOT CONSOLIDATED => $accession\t$type\t$largest_id\t$second_largest_id\n";
					$statistics{$type}{'failed'}++;
				}
			} else {
				my $final_output_id = $largest_id;
				$object->set_consolidated_results($type, $hash_entry_id_xref{$final_output_id}, $final_output_id); # Storing the Search type, Output, and ID output.
				print PROUT "CONSOLIDATED (MAJORITY) => $accession\t$type\t$final_output_id\n";
				$statistics{$type}{'successful'}++;	
			}
		}
	}
}

sub statistics {
	my $type = $_[0];
	print "\n$type\n";
	print "-----\n";
	my $successful = $statistics{$type}{'successful'}++;
	print "Successful: $successful\n";
	my $failed = $statistics{$type}{'failed'}++;
	print "Failed: $failed\n";
	my $percentage = ($successful/($successful+$failed)) * 100;
	print "Percentage: ";
	printf("%.2f", $percentage);
	print "\n\n";
}