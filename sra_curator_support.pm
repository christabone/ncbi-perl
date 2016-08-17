# An object oriented package for the sra_curator.pl script.
package sra_curator_support;

use warnings;
use strict;

# Bless and create a new object.
sub new {
    my ($class, $name) = @_;
    my $self = {
        Input => $name,
        OldMetaData => {},
        key_genes => [{}],
        stage => [{}],
        sample_type => [{}],
        tissue => [{}],
        cell_line => [{}],
        sex => [{}],
        strain => [{}],
        genotype => [{}],
        key_genes_ids => [{}],
        stage_ids => [{}],
        sample_type_ids => [{}],
        tissue_ids => [{}],
        cell_line_ids => [{}],
        sex_ids => [{}],
        strain_ids => [{}],
        genotype_ids => [{}],
        NewCategories => {}
    };
    bless ($self, $class);
    return $self;
}

sub set_old_metadata {
	my ($self, $key, $value) = @_;
	$self->{OldMetaData}{$key} = $value;
}

sub store_new_results {
	my ($self, $search_type, $column_header, $output, $id_output) = @_;
	my $search_type_id = $search_type . '_id';
	my $storage_hash_ref = { $column_header => $output}; # Create the hash_refs to be stored in the AoH.
	my $storage_id_hash_ref = { $column_header => $id_output}; # Create the hash_refs to be stored in the AoH.
	push @{$self->{$search_type}} => $storage_hash_ref; # Store the hash_ref in the appropriate search type array.
	push @{$self->{$search_type_id}} => $storage_id_hash_ref; # Store the hash_ref in the appropriate search type array.
}

sub get_new_results {
	my ($self, $search_type) = @_;
	$self->{$search_type};
}

sub get_old_keys {
	my ($self) = @_;
	my @keys = keys %{ $self->{OldMetaData} };
}

sub get_entire_old_metadata {
	my ($self) = @_;
	$self->{OldMetaData};
}

1; 