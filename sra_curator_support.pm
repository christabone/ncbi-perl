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
        key_genes => {},
        stage => {},
        sample_type => {},
        tissue => {},
        cell_line => {},
        sex => {},
        strain => {},
        genotype => {},
        key_genes_ids => {},
        stage_ids => {},
        sample_type_ids => {},
        tissue_ids => {},
        cell_line_ids => {},
        sex_ids => {},
        strain_ids => {},
        genotype_ids => {},
        ConsolidatedOutput => {},
        MatchScores => {},
        IsMatched => {}
    };
    bless ($self, $class);
    return $self;
}

sub get_accession {
	my ($self) = @_;
	$self->{Input};
}

sub set_score {
	my ($self, $stage, $key, $compare, $value) = @_;
	$self->{MatchScores}{$stage}{$key}{$compare} = $value;
}

sub get_score {
	my ($self, $stage) = @_;
	$self->{MatchScores}{$stage};
}

sub get_all_scores {
	my ($self) = @_;
	$self->{MatchScores};
}

sub set_match {
	my ($self, $stage, $key, $value) = @_;
	$self->{IsMatched}{$key} = $value;
}

sub check_if_matched {
	my ($self, $stage, $key) = @_;
	$self->{IsMatched}{$key};
}

sub set_old_metadata {
	my ($self, $key, $value) = @_;
	$self->{OldMetaData}{$key} = $value;
}

sub set_new_results {
	my ($self, $search_type, $column_header, $output, $id_output) = @_;
	my $search_type_id = $search_type . '_id';
	$self->{$search_type}{$column_header} = $output;
	$self->{$search_type_id}{$column_header} = $id_output;
}

sub set_consolidated_results {
	my ($self, $search_type, $output, $id_output) = @_;
	my $search_type_id = $search_type . '_id';
	$self->{ConsolidatedOutput}{$search_type} = $output;
	$self->{ConsolidatedOutput}{$search_type_id} = $id_output;
}

sub get_consolidated_results {
	my ($self, $search_type, $output, $id_output) = @_;
	my $search_type_id = $search_type . '_id';
	$self->{ConsolidatedOutput}{$search_type};
	$self->{ConsolidatedOutput}{$search_type_id};
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