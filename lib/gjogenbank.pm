#
# This is a SAS component.
#
package gjogenbank;

#===============================================================================
#  Parse one or more GenBank entries in a file into perl structures.
#  All of the entries in the file as a list:
#
#      @entries = parse_genbank( )           #  \*STDIN
#     \@entries = parse_genbank( )           #  \*STDIN
#      @entries = parse_genbank( \*FH )
#     \@entries = parse_genbank( \*FH )
#      @entries = parse_genbank(  $file )
#     \@entries = parse_genbank(  $file )
#      @entries = parse_genbank( \$text )
#     \@entries = parse_genbank( \$text )
#
#  One entry per call:
#
#      $entry = parse_next_genbank( )         #  STDIN
#      $entry = parse_next_genbank( \*FH )
#      $entry = parse_next_genbank( $file )
#
#  Error or end-of-file returns undef.
#
#  Each entry is a hash with key value pairs, some of which have special
#  processing.
#
#     ACCESSION  => [ Accession, ... ]
#     COMMENT    => [ Comment_line, ... ]   # Original wrap is preserved
#     DBLINK     => [ Field, ... ]
#     DBSOURCE   => [ Field, ... ]
#     DEFINITION =>   Definition
#     KEYWORDS   => [ Keyword, ... ]
#     LOCUS      =>   Locus_id
#     ORGANISM   =>   Organism_name
#     ORIGIN     =>   Description_of_sequence_origin
#     REFERENCES => [ { Field => Value, Field => Value, ... }, ... ]
#     SEQUENCE   =>   Sequence
#     SOURCE     =>   Source_string
#     TAXONOMY   => [ Taxon, Taxon, ... ]
#     VERSION    => [ Version, Other_information ]
#
#  Data that are more derived or parsed:
#
#     date       =>   Date
#     division   =>   Division
#     ftr_list   => [ [ type, loc, quals ], ... ]
#     geometry   =>   linear | circular
#     gi         =>   Entry_gi_number
#     is_protein =>   Bool
#     key_index  => { Keyword => 1, ... }
#     locus_id   =>   Locus_id
#     mol_type   =>   Mol_type
#
#  Feature records are merged by type. Slash is removed from qualifier name.
#  Surrounding quotation marks are stripped from qualifier values.
#
#     FEATURES   => { Type => [ [ Location, { Qualifier => \@values } ],
#                               [ Location, { Qualifier => \@values } ],
#                               ...
#                             ],
#                     Type => ...
#                   }
#
#
#-------------------------------------------------------------------------------
#  Access functions to parts of structure:
#
#     @types = feature_types( $entry );
#    \@types = feature_types( $entry );
#
#  Features of a type:
#
#     @ftrs = features_of_type( $entry,  @types );
#    \@ftrs = features_of_type( $entry,  @types );
#     @ftrs = features_of_type( $entry, \@types );
#    \@ftrs = features_of_type( $entry, \@types );
#
#     WARNING: The returned features DO NOT include their respective types, so
#              this function is only useful for features of a single type.
#
#-------------------------------------------------------------------------------
#  Sequence of a feature, optionally including information on partial ends.
#  Use this form for the reasons noted below.
#
#     $seq                           = ftr_seq( $ftr,  $dna   )
#     $seq                           = ftr_seq( $ftr, \$dna   )
#     $seq                           = ftr_seq( $ftr,  $entry )
#   ( $seq, $partial_5, $partial_3 ) = ftr_seq( $ftr,  $dna   )  # boolean of > or < in location
#   ( $seq, $partial_5, $partial_3 ) = ftr_seq( $ftr, \$dna   )
#   ( $seq, $partial_5, $partial_3 ) = ftr_seq( $ftr,  $entry )
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  Deprecated method to get the sequence of a feature, optionally including
#  information on partial ends.  This interface reverses the feature and the
#  sequence args relative to all other feature access functions.  Yuk.
#
#     $seq                           = ftr_dna(  $dna, $ftr )
#     $seq                           = ftr_dna( \$dna, $ftr )
#   ( $seq, $partial_5, $partial_3 ) = ftr_dna(  $dna, $ftr )
#   ( $seq, $partial_5, $partial_3 ) = ftr_dna( \$dna, $ftr )
#
#-------------------------------------------------------------------------------
#  Feature ids in order of preference.  For non-CDS features on many genomes,
#  locus_tag is the most consistently present.  For proteins, protein_id is
#  is usually the accession number for the corresponding protein database entry.
#
#     $id  = ftr_id( $ftr,  @types )  #  The first id drawn from an ordered list
#     $id  = ftr_id( $ftr, \@types )  #      of type preferences, each defined by
#                                     #    1. a feature qualifier name,
#                                     #    2. "gi", or
#                                     #    3. "xref:$db", where $db is a database id
#                                     #  It is up to the user to handle multiword ids.
#
#  Predefined lists:
#
#     $id  = ftr_id( $ftr )           #  protein_id, locus_tag, gene or gi
#     $id  = ftr_locus_tag( $ftr )    #  locus_tag, protein_id, gene or gi
#     $id  = ftr_old_tag( $ftr )      #  old_locus_tag, locus_tag, protein_id, gene or gi
#     $id  = ftr_gene_or_id( $ftr )   #  gene, locus_tag, protein_id or gi
#     $id  = ftr_gi_or_id( $ftr )     #  gi, protein_id, locus_tag or gene
#
#  Get the feature id of a specific type, or return undef:
#
#     $gi  = ftr_gi( $ftr )           #  gi number or undef
#     $id  = ftr_xref( $ftr, $type )  #  db cross reference of $type
#     @ids = ftr_xref( $ftr, $type )  #  db cross references of $type
#
#  Feature ids in old interface:
#
#     $id  = CDS_id( $ftr )           #  protein_id, locus_tag, gene or gi
#     $id  = CDS_locus_tag( $ftr )    #  locus_tag, protein_id, gene or gi
#     $id  = CDS_gi_or_id( $ftr )     #  gi, protein_id, locus_tag or gene
#     $gi  = CDS_gi( $ftr )           #  gi number or undef
#
#-------------------------------------------------------------------------------
#  Feature location (as GenBank format location string; see conversion
#  conversion functions below).
#
#     $ftr_location = location( $ftr )      #  Returns empty string on failure.
#
#  Feature location as cbdl = [ [ contig, beg, dir, len ], ... ]
#
#     $loc                           = location_as_cbdl( $ftr, $entry )
#   ( $loc, $partial_5, $partial_3 ) = location_as_cbdl( $ftr, $entry )
#
#  Feature location as a SEED or Sapling location string
#
#     $loc                           = location_as_seed( $ftr, $entry )
#   ( $loc, $partial_5, $partial_3 ) = location_as_seed( $ftr, $entry )
#
#     $loc                           = location_as_sapling( $ftr, $entry )
#   ( $loc, $partial_5, $partial_3 ) = location_as_sapling( $ftr, $entry )
#
#  Identify features with partial 5' or 3' ends.
#
#     $partial_5_prime = partial_5_prime( $ftr )
#     $partial_3_prime = partial_3_prime( $ftr )
#
#    \%ftr_qualifiers = qualifiers( $ftr )  #  Returns empty hash reference on failure.
#
#     $gene              = gene( $ftr )
#     @gene_and_synonyms = gene( $ftr )
#
#     $product = product( $ftr )
#
#     @EC_number = EC_number( $ftr )
#    \@EC_number = EC_number( $ftr )
#
#     $pseudo    = is_pseudo( $ftr )
#
#  CDS translation table number
#
#     $trans_table = CDS_trans_table( $ftr )
#
#  CDS translation (uses the supplied translation if provided)
#
#     $translation = CDS_translation( $ftr,  $entry ) # This is the preferred form
#     $translation = CDS_translation( $ftr,  $dna )   # Assumes code table 1
#     $translation = CDS_translation( $ftr, \$dna )   # Assumes code table 1
#     $translation = CDS_translation( $ftr )          # Cannot de novo translate
#
#-------------------------------------------------------------------------------
#  Convert GenBank location to [ [ $contig, $begin, $dir, $length ], ... ]
#
#    \@cbdl = genbank_loc_2_cbdl( $loc, $contig_id )
#
#  Convert GenBank location to a SEED or Sapling location string.
#
#     $loc                           = genbank_loc_2_seed( $acc, $loc )
#   ( $loc, $partial_5, $partial_3 ) = genbank_loc_2_seed( $acc, $loc )
#
#     $loc                           = genbank_loc_2_sapling( $acc, $loc )
#   ( $loc, $partial_5, $partial_3 ) = genbank_loc_2_sapling( $acc, $loc )
#
#  Convert a [ [ contig, begin, dir, length ], ... ] location to GenBank.
#
#     $gb_location            = cbdl_2_genbank( \@cbdl )
#   ( $contig, $gb_location ) = cbdl_2_genbank( \@cbdl )
#
#===============================================================================
#  Write a GenBank entry:
#
#     write_genbank(             \%entry );
#     write_genbank( $file_name, \%entry );
#     write_genbank( \*FH,       \%entry );
#
#===============================================================================

use strict;

use Data::Dumper;

require Exporter;
our @ISA    = qw( Exporter );
our @EXPORT = qw( parse_genbank
                  parse_next_genbank
                  feature_types
                  features_of_type
                  feature_list
                  ftr_dna
                  ftr_seq
                  CDS_translation
                  CDS_trans_table

                  genbank_loc_2_seed
                  genbank_loc_2_sapling
                  genbank_loc_2_cbdl
                  cbdl_2_genbank

                  write_genbank
                );

our @EXPORT_OK = qw ( next_entry
                      location
                      location_as_cbdl
                      location_as_seed
                      location_as_sapling
                      partial_5_prime
                      partial_3_prime
                      qualifiers
                      ftr_id
                      ftr_locus_tag
                      ftr_old_tag
                      ftr_gene_or_id
                      ftr_gi_or_id
                      ftr_gi
                      ftr_xref
                      gene
                      product
                      is_pseudo
                      EC_number
                    );

#  An approximate ordering of the common qualifiers in GenBank feature table
#  entries:

my %qualifier_order;
{
    my $tmp_n = 1;
    %qualifier_order = map { $_ => $tmp_n++ }
                       qw( organism
                           sub_species
                           sex
                           mating_type
                           chromosome
                           macronuclear
                           plasmid
                           organelle
                           strain
                           sub_strain
                           cultivar
                           serotype
                           serovar
                           variety
                           isolate
                           tissue_type
                           cell_type
                           dev_stage
                           cell_line
                           rearranged
                           clone_lib
                           clone
                           sub_clone
                           host
                           lab_host
                           mol_type
                           direction
                  
                           rpt_type
                           rpt_family
                           rpt_unit_range
                           rpt_unit_seq
                  
                           operon
                           gene
                           gene_synonym
                           allele
                           locus_tag
                           old_locus_tag
                           codon_start
                           transl_table
                  
                           anticodon
                           product
                           function
                           GO_component
                           GO_function
                           GO_process
                           EC_number
                           protein_id
                           locus_peptide
                           pseudo
                  
                           note
                           db_xref
                  
                           translation
                        );
}


#  Qualifiers that do not require values:

my %valueless_qual = map { $_ => 1 }
                     qw( artificial_location
                         environmental_sample
                         focus
                         germline
                         macronuclear
                         partial
                         proviral
                         pseudo
                         rearranged
                         ribosomal_slippage
                         transgenic
                         trans_splicing
                      );

#  Qualifiers whose values are not wrapped in quotation marks:

my %unquoted_qual  = map { $_ => 1 }
                     qw( anticodon
                         citation
                         codon_start
                         compare
                         direction
                         estimated_length
                         mod_base
                         number
                         rpt_type
                         rpt_unit_range
                         tag_peptide
                         transl_except
                         transl_table
                      );


my %genbank_streams;   # used by parse_next_genbank() to track open streams

#===============================================================================
#  Read GenBank entries from one or more files and/or strings.
#
#    First parameter is optional params hash.
#    @entries = parse_genbank( )                 #  \*STDIN
#   \@entries = parse_genbank( )                 #  \*STDIN
#    @entries = parse_genbank( \*FH, ... )       #  open file
#   \@entries = parse_genbank( \*FH, ... )       #  open file
#    @entries = parse_genbank(  $gb_file, ... )  #  file name
#   \@entries = parse_genbank(  $gb_file, ... )  #  file name
#    @entries = parse_genbank( \$gb_text, ... )  #  reference to data string
#   \@entries = parse_genbank( \$gb_text, ... )  #  reference to data string
#
#===============================================================================

sub parse_genbank
{
    my $params;
    if (ref($_[0]) eq 'HASH')
    {
	$params = shift @_;
    }
    my @entries;
    my $first = 1;
    while ( @_ || $first )
    {
        my $file = shift;
        my ( $fh, $close ) = &input_filehandle( $file );
        if ( $fh )
        {
            while ( my $entry = parse_one_genbank_entry( $fh, $params ) ) { push @entries, $entry }
            close $fh if $close;
        }
        $first = 0;
    }

    wantarray ? @entries : \@entries;
}


#-------------------------------------------------------------------------------
#  Read and parse a GenBank file, one entry at a time.  Successive calls with
#  same parameter will return successive entries.  Calls to different files
#  can be interlaced.
#
#      $entry = next_entry( )         #  STDIN
#      $entry = next_entry( \*FH )
#      $entry = next_entry(  $file )
#
#      $entry = parse_next_genbank( )         #  STDIN
#      $entry = parse_next_genbank( \*FH )
#      $entry = parse_next_genbank(  $file )
#
#  Error or end-of-file returns undef.
#-------------------------------------------------------------------------------

sub parse_next_genbank { next_entry( @_ ) }

sub next_entry
{
    my $file = shift;

    my $stream = $genbank_streams{ $file || '' };
    if ( ! $stream )
    {
        $stream = [ &input_filehandle( $file ) ];   #  Value is ( $fh, $close )
        $stream->[0] or return undef;               #  Got a file handle?
        $genbank_streams{ $file || '' } = $stream;
    }

    my ( $fh, $close ) = @$stream;
    my $entry = parse_one_genbank_entry( $fh );

    if ( ! $entry ) { close $fh if $close; delete $genbank_streams{ $file || '' }; }

    $entry;
}


#-------------------------------------------------------------------------------
#  If it should be necessary to close a stream openned by parse_next_genbank()
#  before it reaches the end-of-file, this will do it.
#
#      close_next_genbank( )         # does nothing
#      close_next_genbank( \*FH )
#      close_next_genbank( $file )
#
#-------------------------------------------------------------------------------
sub close_next_genbank
{
    my $file = shift;
    my $stream = $genbank_streams{ $file || '' };
    close $stream->[0] if $stream && ref $stream eq 'ARRAY' && $stream->[1];
}


#-------------------------------------------------------------------------------
#  Parse the next GenBank format entry read from an open file handle.  This is
#  primarily intended as an internal function called through parse_genbank()
#  or parse_next_genbank().
#
#     \%entry = parse_one_genbank_entry( \*FH )
#
#  Error or end-of-file returns undef
#-------------------------------------------------------------------------------
sub parse_one_genbank_entry
{
    my $fh = shift;
    my $params = shift // {};
    local $_;

    my $state = 0;
    if ( defined( $_ = <$fh> ) ) { chomp } else { $state = -1 }

    #  0 = Looking for LOCUS
    #  1 = Header information
    #  2 = Features
    #  3 = Sequence
    # -1 = Error

    my %entry = ();
    
#          1         2         3         4         5         6         7         8
# 12345678901234567890123456789012345678901234567890123456789012345678901234567890
# LOCUS       NC_000909            1664970 bp    DNA     circular BCT 03-DEC-2005
# LOCUS       DGRINCAD_6   9696 BP DS-DNA             SYN       22-AUG-2006
#
    while ( $state == 0 )
    {
        if ( s/^LOCUS\s+// )
        {
            my @parts = split;
            $entry{ locus_id } = $entry{ LOCUS } = shift @parts;
            $entry{ date }     = pop @parts if $parts[-1] =~ m/^\d+-[A-Z][A-Z][A-Z]-\d+$/i;
            $entry{ division } = pop @parts if $parts[-1] =~ m/^\S\S\S$/;
            $entry{ geometry } = pop @parts if $parts[-1] =~ m/^(lin|circ)/i;
            $entry{ mol_type } = pop @parts if $parts[-1] =~ m/na$/i;
            $state = 1;
        }
        if ( defined( $_ = <$fh> ) ) { chomp } else { $state = -1 }
    }

    #  Reading the header requires merging continuations, then dealing
    #  with the data:

    while ( $state == 1 )
    {
        if ( /^FEATURES / )
        {
            $state = 2;
            if ( defined( $_ = <$fh> ) ) { chomp } else { $state = -1 }
        }

        elsif ( /^CONTIG / )
        {
            $state = 3;
            last;
        }

        elsif ( /^ORIGIN / )
        {
            $state = 3;
            last;
        }

        elsif ( /^REFERENCE / )
        {
            my $ref;
            ( $ref, $state, $_ ) = read_ref( $_, $fh );
            push @{ $entry{ REFERENCES } }, $ref if $ref;
            defined() or $state = -1;
        }

        elsif ( /^(.{10})  +(.*\S)\s*$/ )  # Any other keyword
        {
            my ( $tag, @value ) = ( uc( $1 ), $2 );
            $tag =~ s/^ +| +$//g;

            # Merge continuations:

            if ( defined( $_ = <$fh> ) ) { chomp } else { $state = -1 }

            while ( $state >= 0 && s/^ {12}\s*// )
            {
                s/\s+$//;
                push @value, $_;
                if ( defined( $_ = <$fh> ) ) { chomp } else { $state = -1 }
            }

            #  Special case formats

            if ( $tag eq 'COMMENT' )
            {
                push @{ $entry{ COMMENT } }, @value;
            }

            elsif ( $tag eq 'DBLINK' )
            {
                my $data = '';
                my @dbs;
                foreach ( @value )
                {
                    if ( $data && /: / )
                    {
                        push @dbs, $data;
                        $data = '';
                    }
                    $data .= $_;
                }
                push @dbs, $data  if length $data;

                $entry{ DBLINK } = \@dbs  if @dbs;
            }

            elsif ( $tag eq 'DBSOURCE' )
            {
                my $data = '';
                my @dbs;
                foreach ( @value )
                {
                    if ( $data && /: / )
                    {
                        push @dbs, $data;
                        $data = '';
                    }
                    $data .= $_;
                }
                push @dbs, $data  if length $data;

                $entry{ DBSOURCE } = \@dbs  if @dbs;
            }

            elsif ( $tag eq 'ORGANISM' )
            {
                my $org = shift @value;
                #
                #  Long genome names may split into additional lines in the value; we need
                #  to bring them into the org name. We will bring in values until we get
                #  one with a ; since that begins the taxonomy.
                #
                $entry{ ORGANISM } = $org;
                while (@value && $value[0] !~ /\.$/ && $value[0] !~ /;/)
                {
                    my $ent = shift @value;
                    $entry { ORGANISM } .= " $ent";
                }
                my $tax = @value ? join( ' ', @value ) : '';
                $tax =~ s/\s*\.$//;
                $entry{ TAXONOMY } = [ split /; */, $tax ] if $tax;
            }

            else
            {
                if ( @value > 1 )
                {
                    foreach ( @value[ 0 .. @value-2 ] ) { $_ .= ' ' if ! /-$/ }
                }
                my $value = join( '', @value );

                if ( $tag eq 'ACCESSION' )
                {
                    $entry{ ACCESSION } = [ split / +/, $value ];
                }
                elsif ( $tag eq 'VERSION' )
                {
                    if ( $value =~ s/ +GI:(\d+)// )
                    {
                        $entry{ gi } = $1;
                    }
                    $entry{ VERSION } = [ split / +/, $value ];
                }
                elsif ( $tag eq 'KEYWORDS' )
                {
                    $value =~ s/\s*\.\s*$//;
                    if ( $value )
                    {
                        my @keys = split /; */, $value;
                        $entry{ KEYWORDS  } = \@keys;
                        $entry{ key_index } = { map { $_ => 1 } @keys };
                    }
                }

                #  Generic case:

                else
                {
                    $entry{ $tag } = $value;
                }
            }

            # To know that we are at end of continuations, we must have
            # read another line.

            defined() or $state = -1;
        }

        else  # This is really a format error, but let's skip it.
        {
            if ( defined( $_ = <$fh> ) ) { chomp } else { $state = -1 }
        }
    }

    #
    #  Reading FEATURES requires merging continuations, then dealing
    #  with the data:
    #
    while ( $state == 2 && ( /^     (\S+)\s+(\S+)/ ) )
    {
        my ( $type, $loc ) = ( $1, $2 );

        #  Collect the rest of the location:

        if ( defined( $_ = <$fh> ) ) { chomp } else { $state = -1 }
        while ( $state >= 0 && /^ {15}\s*([^\/\s]\S*)\s*$/ )
        {
            $loc .= $1;
            if ( defined( $_ = <$fh> ) ) { chomp } else { $state = -1 }
        }

        #  Collect qualiiers:

        my ( $qualif, $value, %qualifs, @qualif_order );
        while ( $state == 2 && ( $_ =~ /^\s*\/\w+/ ) )
        {
            #  Qualifiers without = get an undef value (intentionally)

            ( $qualif, undef, $value ) = /^\s*\/(\w+)(=(.*))?/;
            #  Quoted strings can have value lines that start with /, so
            #  we must track quotation marks.

            my $nquote = $value ? $value =~ tr/"// : 0;

            if ( defined( $_ = <$fh> ) ) { chomp } else { $state = -1 }
            while ( $state >= 0 && /^ {15}/ && ( $nquote % 2 || ( ! /^\s*\/\w/ ) ) )
            {
                s/^ +//;
                $nquote += tr/"//;
                $value  .= ( $value =~ /\S-$/ ? '' : ' ' ) . $_;
                if ( defined( $_ = <$fh> ) ) { chomp } else { $state = -1 }
            }

            if ( $nquote % 2 )
            {
                print STDERR "Feature quotation nesting error: $type, $loc, $qualif=$value\n";
                exit;
            }

            if ( $qualif )
            {
                if ( $valueless_qual{ $qualif } )
                {
                    $value = 1;
                }
                elsif ( ! defined $value || ! length $value )
                {
                    next;
                }
                else
                {
                    $value =~ s/""/"/g  if $value =~ s/^"(.*)"$/$1/;
                    $value =~ s/ +$//;
                    $value =~ s/ +//g   if $qualif eq 'translation';
                }

                push @qualif_order, $qualif  if ! $qualifs{ $qualif };
                push @{ $qualifs{ $qualif } }, $value;
            }
        }

        $qualifs{ '_qualif_order' } = \@qualif_order if @qualif_order;
        push @{ $entry{ FEATURES }->{ $type } }, [ $loc, \%qualifs ] if ( $type && $loc );
        push @{ $entry{ ftr_list } },     [ $type, $loc, \%qualifs ] if ( $type && $loc );

        defined() or $state = -1;
    }

    $state = 3 if $state > 0;

    #  Only a few fields are allowed in the standard, but:

    while ( $state == 3 )
    {
        #  Introducer to sequence data

        if ( /^ORIGIN/ )
        {
            $entry{ ORIGIN } = $1 if /^ORIGIN \s+(\S.*\S)\s*$/;
            $state = 4;
        }

        #  2011/01/23 -- GJO
        #
        #  NCBI has stopped including BASE COUNT.  Their example is a disaster
        #  (http://www.ncbi.nlm.nih.gov/genbank/genomesubmit-Examples.html):
        #
        #  BASE COUNT  1165552 a 648314 c 647106 g1169556 t
        #  ORIGIN
        #
        elsif ( s/^BASE COUNT\s+// )
        {
            #  $entry{ BASE_COUNT } = { reverse split };
            if ( defined( $_ = <$fh> ) ) { chomp } else { $state = -1 }
        }

        elsif ( /^(.{10})  (.*)$/ )  # Any other keyword
        {
            my ( $tag, @value ) = ( uc( $1 ), $2 );
            foreach ( $tag, @value ) { s/^ +| +$//g }

            # Merge continuations:

            if ( defined( $_ = <$fh> ) ) { chomp } else { $state = -1 }

            while ( $state >= 0 && s/^ {12}\s*// )
            {
                s/\s+$//;
                push @value, $_;
                if ( defined( $_ = <$fh> ) ) { chomp } else { $state = -1 }
            }

            #  Special case formats

            if ( $tag eq 'COMMENT' )
            {
                push @{ $entry{ COMMENT } }, @value;
            }

            elsif ( $tag eq 'CONTIG' )
            {
                #  The value of CONTIG is a location string and is merged
                #  without spaces.
                $entry{ CONTIG } = join( '', @value );
            }

            else
            {
                if ( @value > 1 )
                {
                    foreach ( @value[ 0 .. @value-2 ] ) { $_ .= ' ' if ! /-$/ }
                }
                $entry{ $tag } = join( '', @value );
            }
        }

        #  End of entry without sequence data.  These can be replaced by
        #  WGS and WGS_SCAFLD, or CONTIG or possibly something else.

        elsif ( $_ eq '//' )
        {
            #  We could fetch and assemble the data
            if ( 0 and $entry{ CONTIG } and eval { require NCBI_sequence } )
            {
                $entry{ SEQUENCE } = NCBI_sequence::contig( $entry{CONTIG} );
            }
            $state = 0;
        }

        else  # This is really a format error, but let's skip it.
        {
            if ( defined( $_ = <$fh> ) ) { chomp } else { $state = -1 }
        }
    }

    #  Read the sequence:

    while ( $state == 4 )
    {
        my @sequence;

        if ( defined( $_ = <$fh> ) ) { chomp } else { $state = -1 }
        while ( $state >= 0 && s/^ *\d+ +// )
        {
	    if (!$params->{skip_contigs})
	    {
                s/[^A-Za-z]+//g;
                push @sequence, $_;
	    }
            if ( defined( $_ = <$fh> ) ) { chomp } else { $state = -1 }
        }

        $entry{ SEQUENCE } = join '', @sequence if !$params->{skip_contigs};

        $state = ( $_ eq '//' ) ? 0 : -1 if $state >= 0;
    }

    $state >= 0 ? \%entry : undef;
}


#-------------------------------------------------------------------------------
#  Parse a reference.
#-------------------------------------------------------------------------------
sub read_ref
{
    my ( $line, $fh ) = @_;
    my $state = 1;
    my %ref = ();

    if ( $line =~ /\s*(\d+)/ )
    {
        $ref{ ref_num  } = $1;
        $ref{ ref_note } = $1 if $line =~ /\s*\d+\s+\((.*)\)\s*$/;
    }

    local $_;
    if ( defined( $_ = <$fh> ) ) { chomp } else { $state = -1 }

    my ( $tag, $value );
    while ( ( $state >= 0 ) && /^  / )
    {
        if ( substr( $_, 0, 10 ) =~ /\S/ )
        {
            ( $tag, $value ) = $_ =~ /\s*(\w+)\s+(.*)$/;
        }
        elsif ( /\S/ )
        {
            s/^ +//;
            $value .= " $_" if $_;
        }

        if ( defined( $_ = <$fh> ) ) { chomp } else { $state = -1 }

        if ( ( $state < 0 ) || ( /^ {0,9}\S/ ) )
        {
            if ( $tag && $value )
            {
                $ref{ $tag } = $value;
            }
            $tag = $value = undef;
        }
    }

    ( ( keys %ref ? \%ref : undef ), $state, $_ )
}



#===============================================================================
#  Access methods for some features and feature data
#===============================================================================
#
#     @types = feature_types( $entry );
#    \@types = feature_types( $entry );
#
#-------------------------------------------------------------------------------
sub feature_types
{
    my $entry = shift;
    return wantarray ? () : [] if ! ( $entry && ref $entry eq 'HASH' );

    my $ftrs = $entry->{ FEATURES };
    return wantarray ? () : [] if ! ( $ftrs && ref $ftrs eq 'HASH' );

    my @types = sort { lc $a cmp lc $b } grep { ! /^_/ } keys %$ftrs;
    wantarray ? @types : \@types;
}


#-------------------------------------------------------------------------------
#  This is really lame for more than one type because the individual features
#  are not marked with their types.  This should be more effecient than the
#  following method for finding all the features of one type.
#
#     @ftrs = features_of_type( $entry,  @types );
#    \@ftrs = features_of_type( $entry,  @types );
#     @ftrs = features_of_type( $entry, \@types );
#    \@ftrs = features_of_type( $entry, \@types );
#
#-------------------------------------------------------------------------------
sub features_of_type
{
    my $entry = shift;
    return wantarray ? () : [] if ! ( $entry && ref $entry eq 'HASH' );

    my $ftrs = $entry->{ FEATURES };
    return wantarray ? () : [] if ! ( $ftrs && ref $ftrs eq 'HASH' );

    my @types = ( ! $_[0] )              ? sort { lc $a cmp lc $b } keys %$ftrs
              : ( ref $_[0] eq 'ARRAY' ) ? @{ $_[0] }
              :                            @_;

    my @ftrs = map { @{ $ftrs->{ $_ } || [] } } @types;
    wantarray ? @ftrs : \@ftrs;
}


#------------------------------------------------------------------------------- 
#  Get a list of a features as triples [ $type, $loc, \%quals ] in genome
#  order:
#
#      @features = feature_list( $entry )
#     \@features = feature_list( $entry )
#      @features = feature_list( $entry, $type )
#     \@features = feature_list( $entry, $type )
#
#  Uses $entry->{ ftr_list } to cache the results (before filtering by type)
#------------------------------------------------------------------------------- 
sub feature_list
{
    my $entry = shift;

    if ( ! $entry->{ ftr_list } || ref $entry->{ ftr_list } ne 'ARRAY' )
    {
        my $features = $entry->{ FEATURES };

        #  Is it a list of [ $type, $location, \@qualifiers ]?

        if ( ref $features eq 'ARRAY' )
        {
            $entry->{ ftr_list } = $features;
        }

        #  Is it a hash of ( $type => [ [ $location, \@qualifiers ], ... ] )?
        #  Build a list of features:

        elsif ( ref $features eq 'HASH' )
        {
            my @features;
            foreach my $type ( keys %$features )
            {
                push @features, map { [ $type, @$_ ] } @{ $features->{ $type } };
            }
            $entry->{ ftr_list } = sort_feature_list( \@features );
        }
    }

    if ( $_[0] )
    {
        my @features = grep { $_->[0] eq $_[0] } @{ $entry->{ ftr_list } };
        return wantarray ? @features : \@features;
    }

    wantarray ? @{ $entry->{ ftr_list } } : $entry->{ ftr_list };
}


#------------------------------------------------------------------------------- 
#  Order features by location.
#
#      @feautres = sort_feature_list(  @features )
#     \@feautres = sort_feature_list( \@features )
#
#  Handles [ $location, \%quals] and [ $type, $location, \%quals ]
#------------------------------------------------------------------------------- 
sub sort_feature_list
{
    my $by_ref = $_[0]      && ( ref $_[0]      eq 'ARRAY' )
              && $_[0]->[0] && ( ref $_[0]->[0] eq 'ARRAY' );

    my @features = map  { $_->[0] }
                   sort { $a->[1] <=> $b->[1] || $b->[2] <=> $a->[2] }
                   map  { [ $_, end_coordinates( $_->[-2] ) ] }
                   $by_ref ? @{ $_[0] } : @_;

    wantarray ? @features : \@features;
}
 

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  For a feature location, find its left and right end coordinates.
#
#      ( $left, $right ) = end_coordinates( $location )
#
#  Rather than parsing, extract a list of coordinates, and take the extremes.
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub end_coordinates
{
    local $_ = shift;
    s/([,(])[^,(]+\:[^,)]+([,)])/$1$2/g;     #  Delete locations in other contigs
    ( sort { $a <=> $b } m/(\d+)/g )[0,-1];  #  Return the lowest and highest values
}


#-------------------------------------------------------------------------------
#  Sequence of a feature.  In list context, include information on partial ends.
#  Can get extract the data from a DNA string, reference to a string, or from
#  the SEQUENCE in an entry.  However, this does not handle locations that
#  specify contigs.  Also, this adjusts CDS features to the first nucleotide
#  of the first complete codon, which is not really what we should be doing.
#
#     $seq                           = ftr_seq( $ftr,  $dna   )
#     $seq                           = ftr_seq( $ftr, \$dna   )
#     $seq                           = ftr_seq( $ftr,  $entry )
#   ( $seq, $partial_5, $partial_3 ) = ftr_seq( $ftr,  $dna   )  # boolean of > or < in location
#   ( $seq, $partial_5, $partial_3 ) = ftr_seq( $ftr, \$dna   )
#   ( $seq, $partial_5, $partial_3 ) = ftr_seq( $ftr,  $entry )
#
#  Handles both [ $location, \%quals ] and [ $type, $location, \%quals ]
#
#  Deprecated interface because args are reversed relative to other methods:
#
#     $seq                           = ftr_dna(  $dna,   $ftr )
#     $seq                           = ftr_dna( \$dna,   $ftr )
#     $seq                           = ftr_dna(  $entry, $ftr )
#   ( $seq, $partial_5, $partial_3 ) = ftr_dna(  $dna,   $ftr )
#   ( $seq, $partial_5, $partial_3 ) = ftr_dna( \$dna,   $ftr )
#   ( $seq, $partial_5, $partial_3 ) = ftr_dna(  $entry, $ftr )
#

sub ftr_dna { ftr_seq( @_[1,0] ) }

#-------------------------------------------------------------------------------

sub ftr_seq
{
    my ( $ftr, $dna ) = @_;
    return undef if ! ( $ftr && $dna );

    my $dnaR =   ref $dna eq 'SCALAR'                     ?  $dna
             :   ref $dna eq 'HASH' && $dna->{ SEQUENCE } ? \$dna->{ SEQUENCE }
             : ! ref $dna                                 ? \$dna
             :                                               undef;
    return undef if ! $dnaR;

    eval { require gjoseqlib; }
        or return undef;

    my $loc = &location( $ftr );
    $loc or return undef;

    my $loc0 = $loc;
    my $complement = ( $loc =~ s/^complement\((.*)\)$/$1/ );
    $loc =~ s/^join\((.*)\)$/$1/;
    my @spans = split /,/, $loc;
    #  For each substring, see if it needs to be complemented.  This does
    #  not occur unless pieces are drawn from multiple contigs, which we
    #  are not dealing with here anyway.
    my @cspan = map { s/^complement\((.+)\)$/$1/ } @spans;
    if ( grep { ! /^<?\d+\.\.>?\d+$/ } @spans )
    {
        print STDERR "*** Feature location parse error: $loc0\n";
        return undef;
    }

    my $partial_5 = $spans[ 0] =~ s/^<//;
    my $partial_3 = $spans[-1] =~ s/\.\.>/../;
    ( $partial_5, $partial_3 ) = ( $partial_3, $partial_5 ) if $complement;

    my @parts;
    for ( my $i = 0; $i < @spans; $i++ )
    {
        $parts[$i] = $cspan[$i] ? gjoseqlib::complement_DNA_seq( extract_span( $dnaR, $spans[$i] ) )
                                :                                extract_span( $dnaR, $spans[$i] );
    }
    my $seq = join( '', @parts );
    $seq = gjoseqlib::complement_DNA_seq( $seq ) if $complement;

    #  Sequences that run off the end can start at other than the first
    #  nucleotide of a codon.
    #  This should not really be done here, even though it makes life easier
    #  for the calling program. -- GJO, 2015-08-27

    my $qual = &qualifiers( $ftr );
    my $codon_start = $qual->{ codon_start } ? $qual->{ codon_start }->[0] : 1;
    $seq = substr( $seq, $codon_start-1 ) if $codon_start > 1;

    wantarray ? ( $seq, $partial_5, $partial_3 ) : $seq;
}


sub extract_span
{
    my ( $dnaR, $span ) = @_;
    my ( $beg, $end ) = $span =~ /^<?(\d+)\.\.>?(\d+)$/;
    ( $beg > 0 ) && ( $beg <= $end ) && ( $end <= length( $$dnaR ) ) or return '';

    substr( $$dnaR, $beg-1, $end-$beg+1 );
}


#-------------------------------------------------------------------------------
#  Get the location of a feature
#
#    $ftr_location = location( $ftr )   #  Returns empty string on failure.
#
#  Handles both [ $location, \%quals ] and [ $type, $location, \%quals ]
#-------------------------------------------------------------------------------
sub location
{
    my ( $ftr ) = @_;

    defined( $ftr ) && ( ref( $ftr ) eq 'ARRAY' ) && ( @$ftr > 1 ) ? $ftr->[-2]
                                                                   : '';
}


#-------------------------------------------------------------------------------
#  Feature location as cbdl = [ [ contig, beg, dir, len ], ... ]
#
#     $loc                           = location_as_cbdl( $ftr, $entry )
#   ( $loc, $partial_5, $partial_3 ) = location_as_cbdl( $ftr, $entry )
#
#  Feature location as a SEED or Sapling location string
#
#     $loc                           = location_as_seed( $ftr, $entry )
#   ( $loc, $partial_5, $partial_3 ) = location_as_seed( $ftr, $entry )
#
#     $loc                           = location_as_sapling( $ftr, $entry )
#   ( $loc, $partial_5, $partial_3 ) = location_as_sapling( $ftr, $entry )
#
#-------------------------------------------------------------------------------

sub location_as_cbdl
{
    genbank_loc_2_cbdl( location( $_[0] ),
                        ( $_[1]->{ ACCESSION } || [] )->[0] || $_[1]->{ LOCUS },
                      );
}

sub location_as_seed
{
    genbank_loc_2_seed( ( $_[1]->{ ACCESSION } || [] )->[0] || $_[1]->{ LOCUS },
                        location( $_[0] )
                      );
}

sub location_as_sapling
{
    genbank_loc_2_sapling( ( $_[1]->{ ACCESSION } || [] )->[0] || $_[1]->{ LOCUS },
                           location( $_[0] )
                         );
}


#-------------------------------------------------------------------------------
#  Identify features with partial 5' or 3' ends.
#
#     $partial_5_prime = partial_5_prime( $ftr )
#     $partial_3_prime = partial_3_prime( $ftr )
#
#  Handles both [ $location, \%quals ] and [ $type, $location, \%quals ]
#-------------------------------------------------------------------------------
sub partial_5_prime
{
    my $ftr = shift             or return undef;
    my $loc = &location( $ftr ) or return undef;
    my $complement = ( $loc =~ s/^complement\((.*)\)$/$1/ );
    $loc =~ s/^join\((.*)\)$/$1/;
    my @spans = split /,/, $loc;

    $complement ? $spans[-1] =~ /\.\.>/ : $spans[0] =~ /^</;
}


sub partial_3_prime
{
    my $ftr = shift             or return undef;
    my $loc = &location( $ftr ) or return undef;
    my $complement = ( $loc =~ s/^complement\((.*)\)$/$1/ );
    $loc =~ s/^join\((.*)\)$/$1/;
    my @spans = split /,/, $loc;

    $complement ? $spans[0] =~ /^</ : $spans[-1] =~ /\.\.>/;
}


#-------------------------------------------------------------------------------
#  Get the qualifier hash for a feature:
#
#    \%ftr_qualifiers = qualifiers( $ftr )
#
#  Handles both [ $location, \%quals ] and [ $type, $location, \%quals ]
#  Returns empty hash reference on failure.
#  Note that this list may include  _qualif_order => \@keys
#-------------------------------------------------------------------------------
sub qualifiers
{
    my ( $ftr ) = @_;
    my $qual;
    ( defined( $ftr )
         && ( ref( $ftr ) eq 'ARRAY' )
         && ( @$ftr > 1 )
         && defined( $qual = $ftr->[-1] )
         && ( ref( $qual ) eq 'HASH' ) )
         ? $qual
         : {};
}


#-------------------------------------------------------------------------------
#  Feature gene:
#
#   $gene          = gene( $ftr )
#   @gene_and_syns = gene( $ftr )
#
#-------------------------------------------------------------------------------
sub gene
{
    my $qual = &qualifiers( @_ );
    my %seen;
    my @gene = grep { ! $seen{ $_ }++ }
               ( ( $qual->{ gene }         ? @{ $qual->{ gene } }         : () ),
                 ( $qual->{ gene_synonym } ? @{ $qual->{ gene_synonym } } : () )
               );

    wantarray ? @gene : $gene[0];
}


#-------------------------------------------------------------------------------
#  Feature ids in order of preference.  For non-CDS features on many genomes,
#  locus_tag is the most consistently present.  For proteins, protein_id is
#  
#
#   $id  = ftr_id( $ftr,  @types )  #  the first id drawn from an ordered list
#   $id  = ftr_id( $ftr, \@types )  #      of type preferences, each defined by
#                                   #    1. a feature qualifier name,
#                                   #    2. "gi", or
#                                   #    3. "xref:$db", where $db is a database id
#                                   #  It is up to the user to handle multiword ids.
#
#  Predefined lists:
#
#   $id  = ftr_id( $ftr )           #  protein_id, locus_tag, gene or gi
#   $id  = ftr_locus_tag( $ftr )    #  locus_tag, protein_id, gene or gi
#   $id  = ftr_old_tag( $ftr )      #  old_locus_tag, locus_tag, protein_id, gene or gi
#   $id  = ftr_gene_or_id( $ftr )   #  gene, locus_tag, protein_id or gi
#   $id  = ftr_gi_or_id( $ftr )     #  gi, protein_id, locus_tag or gene
#
#  Specific feature ids or undef:
#
#   $gi  = ftr_gi( $ftr )           #  gi number or undef
#   $id  = ftr_xref( $ftr, $type )  #  db cross reference of $type
#   @ids = ftr_xref( $ftr, $type )  #  db cross references of $type
#
#  Feature ids in old interface:
#
#   $id  = CDS_id( $ftr )           #  protein_id, locus_tag, gene or gi
#   $id  = CDS_locus_tag( $ftr )    #  locus_tag, protein_id, gene or gi
#   $id  = CDS_gi_or_id( $ftr )     #  gi, protein_id, locus_tag or gene
#   $gi  = CDS_gi( $ftr )           #  gi number or undef
#
#-------------------------------------------------------------------------------
sub ftr_id
{
    my $ftr  = shift;
    my $qual = &qualifiers( $ftr );

    my @types = grep { $_ } ref( $_[0] ) eq 'ARRAY' ? @$_ : @_;
    @types = qw( protein_id locus_tag gene gi )  if ! @types;

    my $id;
    foreach ( @types )
    {
        if    ( /^gi$/i       ) { $id = ftr_gi( $ftr ) }
        elsif ( /^xref:(.+)$/ ) { $id = ftr_xref( $ftr, $1 ) }
        else                    { $id = ( $qual->{$_} || [] )->[0] }
        last if $id;
    }

    $id;
}


sub ftr_locus_tag  { ftr_id( $_[0], qw( locus_tag protein_id gene gi ) ) }
sub ftr_old_tag    { ftr_id( $_[0], qw( old_locus_tag locus_tag protein_id gene gi ) ) }
sub ftr_gene_or_id { ftr_id( $_[0], qw( gene locus_tag protein_id gi ) ) }
sub ftr_gi_or_id   { ftr_id( $_[0], qw( gi locus_tag protein_id gene ) ) }


sub ftr_gi
{
    my $qual = &qualifiers( @_ );

    my ( $id ) = map { m/^GI:(.+)$/i ? $1 : () } @{ $qual->{db_xref} || [] };

    $id;
}


sub ftr_xref
{
    my ( $ftr, $type ) = @_;
    my $qual = &qualifiers( $ftr );

    my @ids = map { m/^\Q$type\E:(.+)$/i ? $1 : () } @{ $qual->{db_xref} || [] };

    wantarray ? @ids : $ids[0];
}


sub CDS_id
{
    my $qual = &qualifiers( @_ );
    my $id;

    ( $id ) =                                 @{ $qual->{ protein_id } } if          $qual->{ protein_id };
    ( $id ) =                                 @{ $qual->{ locus_tag } }  if ! $id && $qual->{ locus_tag };
    ( $id ) =                                 @{ $qual->{ gene } }       if ! $id && $qual->{ gene };
    ( $id ) = map { m/^GI:(.+)$/i ? $1 : () } @{ $qual->{ db_xref } }    if ! $id && $qual->{ db_xref };

    $id;
}


sub CDS_gi_or_id
{
    my $qual = &qualifiers( @_ );
    my $id;

    ( $id ) = map { m/^GI:(.+)$/i ? $1 : () } @{ $qual->{ db_xref } }    if          $qual->{ db_xref };
    ( $id ) =                                 @{ $qual->{ protein_id } } if ! $id && $qual->{ protein_id };
    ( $id ) =                                 @{ $qual->{ locus_tag } }  if ! $id && $qual->{ locus_tag };
    ( $id ) =                                 @{ $qual->{ gene } }       if ! $id && $qual->{ gene };

    $id;
}


sub CDS_locus_tag
{
    my $qual = &qualifiers( @_ );
    my $id;

    ( $id ) =                                 @{ $qual->{ locus_tag } }  if          $qual->{ locus_tag };
    ( $id ) =                                 @{ $qual->{ protein_id } } if ! $id && $qual->{ protein_id };
    ( $id ) =                                 @{ $qual->{ gene } }       if ! $id && $qual->{ gene };
    ( $id ) = map { m/^GI:(.+)$/i ? $1 : () } @{ $qual->{ db_xref } }    if ! $id && $qual->{ db_xref };

    $id;
}


sub CDS_gi
{
    my $qual = &qualifiers( @_ );

    my ( $id ) = map { m/^GI:(.+)$/i ? $1 : () } @{ $qual->{ db_xref } } if $qual->{ db_xref };

    $id;
}


#-------------------------------------------------------------------------------
#  Feature product:
#
#   $product = product( $ftr )
#
#-------------------------------------------------------------------------------
sub product
{
    my $qual = &qualifiers( @_ );
    my $prod;

    ( $prod ) = @{ $qual->{ product } }  if            $qual->{ product };
    ( $prod ) = @{ $qual->{ function } } if ! $prod && $qual->{ function };
    ( $prod ) = @{ $qual->{ note } }     if ! $prod && $qual->{ note };

    $prod;
}


#-------------------------------------------------------------------------------
#  Feature is pseudo gene?
#
#   $pseudo = is_pseudo( $ftr )
#
#-------------------------------------------------------------------------------
sub is_pseudo { qualifiers( $_[0] )->{ pseudo } ? 1 : 0 }


#-------------------------------------------------------------------------------
#
#   @EC_number = EC_number( $ftr )
#  \@EC_number = EC_number( $ftr )
#
#-------------------------------------------------------------------------------
sub EC_number
{
    my $qual = &qualifiers( @_ );
    my @EC = $qual->{ EC_number } ? @{ $qual->{ EC_number } } : ();

    wantarray ? @EC : \@EC;
}


#-------------------------------------------------------------------------------
#  Find a CDS translation table number:
#
#   $trans_table = CDS_trans_table( $ftr );
#
#-------------------------------------------------------------------------------
sub CDS_trans_table
{
    ( qualifiers( $_[0] )->{ transl_table } || [] )->[0] || 1;
}


#-------------------------------------------------------------------------------
#  This is the in situ translation.  Will extract from the DNA sequence if
#  necessary.
#
#   $translation = CDS_translation( $ftr )
#   $translation = CDS_translation( $ftr,  $dna )
#   $translation = CDS_translation( $ftr, \$dna )
#   $translation = CDS_translation( $ftr,  $entry )  #  <--- preferred form
#
#  We should look for translation exceptions, but ....
#
#-------------------------------------------------------------------------------
sub CDS_translation
{
    my ( $ftr, $dna ) = @_;
    my $qual = &qualifiers( $ftr );

    return $qual->{ translation }->[0] if $qual->{ translation };

    return undef if ! $dna;

    #  Beware that ftr_seq() removes a leading partial codon based on /codon_start

    my $CDS_dna = ftr_seq( $ftr, $dna ) or return undef;

    my $transl_table = $qual->{ transl_table }->[0] || 1;

    my $start_with_met = ! partial_5_prime( $ftr );

    translate_seq_with_NCBI_code( $CDS_dna, $transl_table, $start_with_met );
}


sub translate_seq_with_NCBI_code
{
    my ( $seq, $transl_table, $start_with_met ) = @_;

    eval { require gjoseqlib; }
        or return undef;

    eval { require NCBI_genetic_code; }
        or return undef;

    $seq =~ tr/-//d;     #  remove gaps (should never happen)
    $seq =~ tr/Uu/Tt/;   #  make it DNA

    my $gc = NCBI_genetic_code::genetic_code( $transl_table );

    my $ambigs = \%gjoseqlib::DNA_letter_can_be;

    #  We can now do the codon-by-codon translation:

    my @codons = map { /[a-z]/ ? lc( $_ ) : $_ }
                 $seq =~ m/(...?)/g;  #  will try to translate last 2 nt

    my @met;
    if ( $start_with_met && ( my $codon1 = shift @codons ) )
    {
        push @met, ( $codon1 =~ /^[a-z]/ ? 'm' : 'M' );
    }

    my $pep = join( '', @met, map { gjoseqlib::translate_codon_with_user_code( $seq, $gc, $ambigs ) } @codons );

    #  If it ends with stop, and it usually will, remove it

    $pep =~ s/\*$// if $pep;

    $pep;
}


#===============================================================================
#  Utilities for locations and location strings.
#===============================================================================
#  Convert GenBank location to a SEED location string.
#
#     $loc                           = genbank_loc_2_seed( $acc, $loc )
#   ( $loc, $partial_5, $partial_3 ) = genbank_loc_2_seed( $acc, $loc )
#
#-------------------------------------------------------------------------------
sub genbank_loc_2_seed
{
    my ( $acc, $loc ) = @_;
    $acc && $loc or return undef;
    genbank_loc_2_string( $acc, $loc, 'seed' );
}


#-------------------------------------------------------------------------------
#  Convert GenBank location to a Sapling location string.
#
#     $loc                           = genbank_loc_2_sapling( $acc, $loc )
#   ( $loc, $partial_5, $partial_3 ) = genbank_loc_2_sapling( $acc, $loc )
#
#-------------------------------------------------------------------------------
sub genbank_loc_2_sapling
{
    my ( $acc, $loc ) = @_;
    $acc && $loc or return undef;
    genbank_loc_2_string( $acc, $loc, 'sapling' );
}


#-------------------------------------------------------------------------------
#  Convert GenBank location to another location format.
#  At present, only 'sapling' (D) and 'seed' are supported.
#
#     $loc                           = genbank_loc_2_string( $acc, $loc, $format )
#   ( $loc, $partial_5, $partial_3 ) = genbank_loc_2_string( $acc, $loc, $format )
#
#-------------------------------------------------------------------------------
sub genbank_loc_2_string
{
    my ( $acc, $loc, $format ) = @_;
    $acc && $loc or return undef;

    my ( $cbdl, $partial_5, $partial_3 ) = genbank_loc_2_cbdl( $loc, $acc );
    my $str = cbdl_2_string( $cbdl, $format || 'sapling' );

    wantarray ? ( $str, $partial_5, $partial_3 ) : $str;
}


#-------------------------------------------------------------------------------
#  Convert GenBank location to a list of [contig, begin, dir, len ] locations.
#  order() is treated as join().  Nesting is allowed (unlike the standard).
#
#     \@cbdl                           = genbank_loc_2_cbdl( $loc, $accession )
#   ( \@cbdl, $partial_5, $partial_3 ) = genbank_loc_2_cbdl( $loc, $accession )
#
#  Elements are:
#
#   (accession:)?<?\d+..>?\d+                  # range of sites
#   (accession:)?<?\d+^>?\d+                   # site between residues
#   (accession:)?\d+                           # single residue
#   (accession:)?complement\(element\)
#                gap\(\w+\)
#   (accession:)?join\(element,element,...\)
#   (accession:)?order\(element,element,...\)
#
#-------------------------------------------------------------------------------
#  Paterns used in the parsing.  They are in each subroutine due to a very
#  strange initialization issue.
#
#  Because $paranthetical is self-referential, it must be declared before it is
#  defined.
#
#   my $paranthetical;
#      $paranthetical    = qr/\([^()]*(?:(??{$paranthetical})[^()]*)*\)/;
#   my $contigid         = qr/[^\s:(),]+/;
#   my $complement       = qr/(?:$contigid:)?complement$paranthetical/;
#   my $complement_parts = qr/(?:($contigid):)?complement($paranthetical)/;
#   my $join             = qr/(?:$contigid:)?join$paranthetical/;
#   my $join_parts       = qr/(?:($contigid):)?join($paranthetical)/;
#   my $order            = qr/(?:$contigid:)?order$paranthetical/;
#   my $order_parts      = qr/(?:($contigid):)?order($paranthetical)/;
#   my $range            = qr/(?:$contigid:)?<?\d+\.\.>?\d+/;
#   my $range_parts      = qr/(?:($contigid):)?(<?)(\d+)\.\.(>?)(\d+)/;
#   my $site             = qr/(?:$contigid:)?<?\d+^>?\d+/;
#   my $site_parts       = qr/(?:($contigid):)?(<?)(\d+)^(>?)(\d+)/;
#   my $position         = qr/(?:$contigid:)?\d+/;
#   my $position_parts   = qr/(?:($contigid):)?(\d+)/;
#   my $element          = qr/$range|$position|$complement|$join|$order/;
#   my $elementlist      = qr/$element(?:,$element)*/;
#
#-------------------------------------------------------------------------------

sub genbank_loc_2_cbdl
{
    my ( $loc, $acc ) = @_;

    my $contigid = qr/[^\s:(),]+/;

    my $range = qr/(?:$contigid:)?<?\d+\.\.>?\d+/;
    return gb_loc_range( $loc, $acc )      if $loc =~ /^$range$/;

    #  This cannot by part of any other format except complement
    my $site = qr/(?:$contigid:)?<?\d+\^>?\d+/;
    return gb_loc_site( $loc, $acc )       if $loc =~ /^$site$/;

    my $position = qr/(?:$contigid:)?\d+/;
    return gb_loc_position( $loc, $acc )   if $loc =~ /^$position$/;

    my $paranthetical;
    {
       no warnings;
       $paranthetical = qr/\([^()]*(?:(??{$paranthetical})[^()]*)*\)/;
    }
    my $complement = qr/(?:$contigid:)?complement$paranthetical/;
    return gb_loc_complement( $loc, $acc ) if $loc =~ /^$complement$/;

    my $join = qr/(?:$contigid:)?join$paranthetical/;
    return gb_loc_join( $loc, $acc )       if $loc =~ /^$join$/;

    #  Treated as a join
    my $order = qr/(?:$contigid:)?order$paranthetical/;
    return gb_loc_order( $loc, $acc )      if $loc =~ /^$order$/;

    return ();
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  A range of positions, with optional accession number prefix, optional less
#  than first position (5' partial), begin, .., optional greater than end
#  position (3' partial), and end position.
#
#    (\S+:)?<?\d+\.\.>?\d+
#
#      \@cbdl_list                           = gb_loc_range( $loc, $acc )
#    ( \@cbdl_list, $partial_5, $partial_3 ) = gb_loc_range( $loc, $acc )
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub gb_loc_range
{
    my ( $loc, $acc ) = @_;

    my $contigid = qr/[^\s:(),]+/;
    my $range_parts = qr/(?:($contigid):)?(<?)(\d+)\.\.(>?)(\d+)/;
    my ( $acc2, $p5, $beg, $p3, $end ) = $loc =~ /^$range_parts$/;
    $beg && $end or return ();
    $acc2 ||= $acc;

    #  GenBank standard is always $beg <= $end.  We will relax that.

    my $cbdl = [ [ $acc2, $beg, (($end>=$beg)?'+':'-'), abs($end-$beg)+1 ] ];

    wantarray ? ( $cbdl, $p5, $p3 ) : $cbdl;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  A range of positions, with optional accession number prefix, optional less
#  than first position (5' partial), begin, ^, optional greater than end
#  position (3' partial), and end position.
#
#    (\S+:)?<?\d+^>?\d+
#
#      \@cbdl_list                           = gb_loc_site( $loc, $acc )
#    ( \@cbdl_list, $partial_5, $partial_3 ) = gb_loc_site( $loc, $acc )
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub gb_loc_site
{
    my ( $loc, $acc ) = @_;

    my $contigid = qr/[^\s:(),]+/;
    my $site_parts = qr/(?:($contigid):)?(<?)(\d+)\^(>?)(\d+)/;
    my ( $acc2, $p5, $beg, $p3, $end ) = $loc =~ /^$site_parts$/;
    $beg && $end or return ();
    $acc2 ||= $acc;

    #  GenBank standard is always $beg <= $end.  We will relax that.

    my $cbdl = [ [ $acc2, $beg, (($end>=$beg)?'+':'-'), abs($end-$beg)+1 ] ];

    wantarray ? ( $cbdl, $p5, $p3 ) : $cbdl;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  A singe position, with optional accession number prefix.
#
#    (\S+:)?\d+
#
#      \@cbdl_list                           = gb_loc_position( $loc, $acc )
#    ( \@cbdl_list, $partial_5, $partial_3 ) = gb_loc_position( $loc, $acc )
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub gb_loc_position
{
    my ( $loc, $acc ) = @_;

    my $contigid = qr/[^\s:(),]+/;
    my $position_parts = qr/(?:($contigid):)?(\d+)/;
    my ( $acc2, $beg ) = $loc =~ /^$position_parts$/;
    $beg or return ();

    my $cbdl = [ [ $acc2 || $acc, $beg, '+', 1 ] ];

    wantarray ? ( $cbdl, '', '' ) : $cbdl;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#      \@cbdl_list                           = gb_loc_complement( $loc, $acc )
#    ( \@cbdl_list, $partial_5, $partial_3 ) = gb_loc_complement( $loc, $acc )
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub gb_loc_complement
{
    my ( $loc, $acc ) = @_;

    my $paranthetical;
    {
       no warnings;
       $paranthetical = qr/\([^()]*(?:(??{$paranthetical})[^()]*)*\)/;
    }
    my $contigid         = qr/[^\s:(),]+/;
    my $complement_parts = qr/(?:($contigid):)?complement($paranthetical)/;

    my ( $acc2, $loc2 ) = $loc =~ /^$complement_parts$/;
    $loc2 && $loc2 =~ s/^\(// && $loc2 =~ s/\)$// or return ();
    my ( $locs, $p5, $p3 ) = genbank_loc_2_cbdl( $loc2, $acc2 || $acc );
    $locs && ref( $locs ) eq 'ARRAY' && @$locs or return ();

    my $cbdl = [ map { complement_cbdl( @$_ ) } reverse @$locs ];

    wantarray ? ( $cbdl, $p5, $p3 ) : $cbdl;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#      \@cbdl_list                           = gb_loc_join( $loc, $acc )
#    ( \@cbdl_list, $partial_5, $partial_3 ) = gb_loc_join( $loc, $acc )
#
#  There is no warning about partial sequences internal to list.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub gb_loc_join
{
    my ( $loc, $acc ) = @_;

    my $paranthetical;
    {
       no warnings;
       $paranthetical = qr/\([^()]*(?:(??{$paranthetical})[^()]*)*\)/;
    }
    my $contigid         = qr/[^\s:(),]+/;
    my $complement       = qr/(?:$contigid:)?complement$paranthetical/;
    my $join             = qr/(?:$contigid:)?join$paranthetical/;
    my $join_parts       = qr/(?:($contigid):)?join($paranthetical)/;
    my $order            = qr/(?:$contigid:)?order$paranthetical/;
    my $range            = qr/(?:$contigid:)?<?\d+\.\.>?\d+/;
    my $position         = qr/(?:$contigid:)?\d+/;
    my $element          = qr/$range|$position|$complement|$join|$order/;
    my $elementlist      = qr/$element(?:,$element)*/;

    my ( $acc2, $locs ) = $loc =~ /^$join_parts$/;
    $locs && $locs =~ s/^\(// && $locs =~ s/\)$//
          && $locs =~ /^$elementlist$/
          or return ();
    $acc2 ||= $acc;

    my @elements = map { [ genbank_loc_2_cbdl( $_, $acc2 ) ] }
                   $locs =~ m/($element)/g;
    @elements or return ();

    my $cbdl = [ map { @{ $_->[0] } } @elements ];

    wantarray ? ( $cbdl, $elements[0]->[1], $elements[-1]->[2] ) : $cbdl;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  Ordered list is treated as a join:
#
#      \@cbdl_list                           = gb_loc_order( $loc, $acc )
#    ( \@cbdl_list, $partial_5, $partial_3 ) = gb_loc_order( $loc, $acc )
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub gb_loc_order
{
    my ( $loc, $acc ) = @_;

    my $paranthetical;
    {
       no warnings;
       $paranthetical = qr/\([^()]*(?:(??{$paranthetical})[^()]*)*\)/;
    }
    my $contigid         = qr/[^\s:(),]+/;
    my $complement       = qr/(?:$contigid:)?complement$paranthetical/;
    my $join             = qr/(?:$contigid:)?join$paranthetical/;
    my $order            = qr/(?:$contigid:)?order$paranthetical/;
    my $order_parts      = qr/(?:($contigid):)?order($paranthetical)/;
    my $range            = qr/(?:$contigid:)?<?\d+\.\.>?\d+/;
    my $position         = qr/(?:$contigid:)?\d+/;
    my $element          = qr/$range|$position|$complement|$join|$order/;
    my $elementlist      = qr/$element(?:,$element)*/;

    my ( $acc2, $locs ) = $loc =~ /^$order_parts$/;
    $locs && $locs =~ s/^\(// && $locs =~ s/\)$//
          && $locs =~ /^$elementlist$/
          or return ();

    gb_loc_join( "join($locs)", $acc2 || $acc );
}


#-------------------------------------------------------------------------------
#    $cbdl = complement_cbdl(   $contig, $beg, $dir, $len   )
#    $cbdl = complement_cbdl( [ $contig, $beg, $dir, $len ] )
#-------------------------------------------------------------------------------
sub complement_cbdl
{
    defined $_[0] or return ();
    my ( $contig, $beg, $dir, $len ) = ref( $_[0] ) ? @{$_[0]} : @_;

    ( $dir =~ /^-/ ) ? [ $contig, $beg -= $len - 1, '+', $len ]
                     : [ $contig, $beg += $len - 1, '-', $len ];
}


#-------------------------------------------------------------------------------
#   $loc = cbdl_2_string( \@cbdl, $format )
#-------------------------------------------------------------------------------
sub cbdl_2_string
{
    my ( $cbdl, $format ) = @_;
    $cbdl && ( ref( $cbdl ) eq 'ARRAY' ) or return undef;
    $format = 'sapling' if ! defined $format;
    return cbdl_2_genbank( $cbdl ) if $format =~ m/genbank/i;
    join( ',', map { cbdl_part_2_string( $_, $format ) } @$cbdl );
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  Support function for formatting one contiguous part of location.
#
#   $loc_part = cbdl_part_2_string( $cbdl_part, $format )
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub cbdl_part_2_string
{
    my ( $part, $format ) = @_;
    $part && ( ref( $part ) eq 'ARRAY' ) && ( @$part == 4 ) or return ();
    my ( $contig, $beg, $dir, $len ) = @$part;
    $dir = $dir =~ /^-/ ? '-' : '+';

    if ( $format =~ m/seed/i )
    {
        my $n2 = ( $dir eq '+' ) ? $beg + $len - 1 : $beg - $len + 1;
        return join( '_', $contig, $beg, $n2 );
    }

    # Default is sapling:

    return $contig . '_' . $beg . $dir . $len;
}

#-------------------------------------------------------------------------------
#  Convert a [ [ contig, begin, dir, length ], ... ] location to GenBank.
#
#    $gb_location            = cbdl_2_genbank( \@cbdl )
#  ( $contig, $gb_location ) = cbdl_2_genbank( \@cbdl )
#-------------------------------------------------------------------------------
sub cbdl_2_genbank
{
    my ( $cbdl, $contig ) = @_;
    $cbdl && ref( $cbdl ) eq 'ARRAY' && @$cbdl or return '';
    my @cbdl = ref( $cbdl->[0] ) ? @$cbdl : ( $cbdl );
    @cbdl or return '';

    my $dir = $cbdl[0]->[2];
    @cbdl = map { complement_cbdl( $_ ) } reverse @cbdl if $dir =~ /^-/;

    $contig = $cbdl[0]->[0];
    my @gb = map { cbdl_part_2_genbank( $_, $contig ) } @cbdl;

    my $gb = ( @gb > 1 ) ? 'join(' . join( ',', @gb ) . ')' : $gb[0];

    $gb = "complement($gb)" if $dir =~ /^-/;

    return wantarray ? ( $contig, $gb ) : $gb;
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  Support function for formatting one contiguous part of location.
#
#   $loc_part = cbdl_part_2_genbank( $cbdl_part, $contig )
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub cbdl_part_2_genbank
{
    my ( $part, $contig0 ) = @_;
    $part && ( ref( $part ) eq 'ARRAY' ) && ( @$part == 4 ) or return ();
    my ( $contig, $beg, $dir, $len ) = @$part;

    my $gb;
    if ( $dir =~ /^-/ )
    {
        my $end = $beg - ( $len-1);
        $gb = "complement($end..$beg)";
    }
    else
    {
        my $end = $beg + ($len-1);
        $gb = "$beg..$end";
    }

    $gb = "$contig:$gb" if $contig0 && $contig ne $contig0;

    return $gb;
}



#===============================================================================
#  Write a GenBank entry:
#
#     write_genbank(            \%entry );
#     write_genbank( $filename, \%entry );
#     write_genbank( \*FH,      \%entry );
#
#===============================================================================
#  Recognized record types are:
#
#   + ACCESSION  => [ Accession, ... ]
#     COMMENT    => [ Comment, ... ]
#     DBLINK     => [ Field, ... ]
#     DBSOURCE   => [ Field, ... ]
#     DEFINITION =>   Definition
#     KEYWORDS   => { Key_phrase => 1, Key_phrase => 1, ... }
#   + LOCUS      =>   Locus_id            # locus_id is preferred
#   * ORGANISM   =>   Organism_name
#     ORIGIN     =>   Description_of_sequence_origin
#     REFERENCES => [ { Field => Value, Field => Value, ... }, ... ]
#   * SEQUENCE   =>   Sequence
#     SOURCE     =>   Source_string
#     TAXONOMY   => [ Taxon, Taxon, ... ]
#   + VERSION    => [ Version, Other_information ]
#
#  Data that are not in record types, but supply data to build records:
#
#     date       =>   Date
#     geometry   =>   linear | circular
#     gi         =>   Entry_gi_number
#     is_protein =>   Bool
#     locus_id   =>   Locus_id
#     mol_type   =>   Mol_type
#
#  The record types marked with a * are required.
#  At least one of the types marked with a + must also be present.
#
#  Feature records are merged by type.  Slash is removed from qualifier name.
#  Surrounding quotation marks are stripped from qualifier values.
#
#     FEATURES   => { Type => [ [ Location, { Qualifier => \@values } ],
#                               [ Location, { Qualifier => \@values } ],
#                               ...
#                             ],
#                     Type => ...
#                   }
#
#          1         2         3         4         5         6         7         8
# 12345678901234567890123456789012345678901234567890123456789012345678901234567890
# LOCUS       @<<<<<<<<<<<<<<< @########## @< DNA        @<<<<<<< @<< @<<<<<<<<<<
# LOCUS       CYC_HUMAN                105 aa            linear   PRI 15-MAR-2004
# LOCUS       NC_000909            1664970 bp    DNA     circular BCT 03-DEC-2005
#
#===============================================================================


sub write_genbank
{
    my ( $fh, $close );
    if ( $_[0] && ( ( ! ref( $_[0] ) ) || ( ref( $_[0] ) eq 'GLOB' ) ) )
    {
        ( $fh, $close ) = output_filehandle( shift );
    }

    my $entry = shift;
    $entry && ref $entry eq 'HASH' or return 0;
    $entry->{ SEQUENCE } or return 0;
    $entry->{ ORGANISM } or return 0;

    #  Prepare the data for the forms:

    my $key;

    if ( ! $entry->{ locus_id } )
    {
        ( $key ) = grep { m/locus_id/i  } keys %$entry;
        ( $key ) = grep { m/locus/i     } keys %$entry if ! $key;
        ( $key ) = grep { m/contig/i    } keys %$entry if ! $key;
        ( $key ) = grep { m/accession/i } keys %$entry if ! $key;
        ( $key ) = grep { m/version/i   } keys %$entry if ! $key;
        my $locus_id = $key ? $entry->{ $key } : 'undefined';
        $locus_id = $locus_id->[0] if ref $locus_id;
        $locus_id =~ s/\s.*$//;
        $entry->{ locus_id } = $locus_id;
    }

    if ( ! $entry->{ mol_type } )
    {
        my $mol_type;
        my $is_prot  = $entry->{ is_protein };
        $mol_type = ' ' if $is_prot;
        if ( ! $mol_type )
        {
            $mol_type = mol_type( \$entry->{ SEQUENCE } );
            $mol_type = ' ' if $mol_type =~ m/^prot/i;
            $entry->{ is_protein } = ( $mol_type =~ m/NA$/i ) ? 0 : 1;
        }
        $entry->{ mol_type } = $mol_type;
    }

    if ( ! $entry->{ DEFINITION } )
    {
        ( $key ) = grep { m/definition/i } keys %$entry;
        ( $key ) = grep { m/^def/i       } keys %$entry if ! $key;
        my $def = $key ? $entry->{ $key } : $entry->{ locus_id };
        $entry->{ DEFINITION } = $def;
    }

    if ( ! $entry->{ ACCESSION } || ! ref $entry->{ ACCESSION } )
    {
        ( $key ) = grep { m/accession/i } keys %$entry;
        ( $key ) = grep { m/^acc/i      } keys %$entry if ! $key;
        ( $key ) = grep { m/version/i   } keys %$entry if ! $key;
        ( $key ) = grep { m/^ver/i      } keys %$entry if ! $key;
          $key   = 'locus_id'                          if ! $key;
        my $acc = $entry->{ $key };
        $entry->{ ACCESSION } = ref $acc ? $acc : [ $acc ];
    }

    if ( ! $entry->{ VERSION } || ! ref $entry->{ VERSION } )
    {
        ( $key ) = grep { m/version/i } keys %$entry;
        ( $key ) = grep { m/^ver/i    } keys %$entry if ! $key;
          $key   = 'ACCESSION'                       if ! $key;
        my $ver = $entry->{ $key };
        $entry->{ VERSION } = ref $ver ? $ver : [ $ver ];
    }
    
    if ( ! grep { m/^gi:/ } @{ $entry->{ VERSION } } )
    {
        ( $key ) = grep { m/^gi$/i } keys %$entry;
        push @{ $entry->{ VERSION } }, "GI:$entry->{$key}" if $key;
    }


    #  Fill out the forms

    my $old_fh;
    $old_fh = select( $fh ) if ref $fh eq 'GLOB';

    write_locus( $entry );
    write_definition( $entry );
    write_accession( $entry );
    write_version( $entry );
    write_dblink( $entry )       if $entry->{ DBLINK };
    write_dbsource( $entry )     if $entry->{ DBSOURCE };
    write_keywords( $entry );
    write_source( $entry );
    write_organism( $entry );
    write_taxonomy( $entry );
    write_references( $entry );
    write_comment( $entry )      if $entry->{ COMMENT };
    write_features( $entry );
    #  write_base_count( $entry );  #  Just not useful, and error prone
    write_origin( $entry );
    write_sequence( $entry );
    write_end_of_entry();

    select( $old_fh ) if $old_fh;

    close( $fh ) if $fh && $close;

    return 1;
}


sub mol_type
{
    return undef if ! defined $_[0];
    local $_ = ref $_[0] ? shift : \$_[0];
    my $n_nt = $$_ =~ tr/ACGNTUacgntu//;
    my $n    = length( $$_ );
    return ( $n_nt < 0.8 * $n )                ? 'protein'
         : ( $$_ =~ tr/Tt// > $$_ =~ tr/Uu// ) ? 'DNA' : 'RNA'; 
}



sub write_locus
{
    my $entry = shift;
    my ( $locus, $mol_type, $geometry, $div, $date )
        = map { $entry->{ $_ } || ' ' }
          qw( locus_id mol_type geometry division date );

    my $length = $entry->{ length } ||= length( $entry->{ SEQUENCE } );

    my $unit = ( $mol_type =~ m/NA$/i ) ? 'bp' : 'aa';

    my $stranded = ' ';
    ( $stranded, $mol_type ) = ( $1, $2 ) if $mol_type =~ /^(\S+-)(\S+)$/;

    my $key;
    ( $key ) = grep { m/^geometry/i } keys %$entry;
    ( $key ) = grep { m/^geom/i     } keys %$entry if ! $key;
    $geometry = $entry->{ geometry } ||= $key ? $entry->{ $key }
                                          : ' ';

    ( $key ) = grep { m/^div/i } keys %$entry;
    $div = $entry->{ division } ||= $key ? $entry->{ $key } : 'UNK';

    ( $key ) = grep { m/^date/i } keys %$entry;
    ( $key ) = grep { m/date/i  } keys %$entry if ! $key;
    $date = $entry->{ date } ||= $key ? $entry->{ $key }
                                         : ( map { chomp; uc } `date '+%d-%b-%Y'` )[0];

    $^A = '';
    formline( <<'End_locus', $locus, $length, $unit, $stranded, $mol_type, $geometry, $div, $date );
LOCUS       @<<<<<<<<<<<<<<<<< @######## @< @>>@<<<<<< @<<<<<<< @<< @<<<<<<<<<<
End_locus
    print $^A;
    $^A = '';
}


sub write_definition
{
    local $_ = $_[0]->{ DEFINITION } || '';
    $_ .= '.' unless /\.$/;
    output_record( 'DEFINITION', $_ );
}


sub output_record
{
    my ( $type, $text ) = @_;

    $^A = '';

    formline( <<'End_of_line_1', $type, $text );
@<<<<<<<<<  ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
End_of_line_1

    formline( <<'End_of_continuation', $text ) if defined $text && length $text;
~~          ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
End_of_continuation

    print $^A;
    $^A = '';
}


sub write_continuation
{
    local $_ = shift;

    $^A = '';
    formline( <<'End_of_continuation', $_ ) if defined $_ && length $_;
~~          ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
End_of_continuation
    print $^A;
    $^A = '';
}


sub write_accession
{
    local $_ = join( ' ', @{ $_[0]->{ ACCESSION } || [] } );
    output_record( 'ACCESSION', $_ ) if length $_;
}


sub write_version
{
    local $_ = join( ' ', @{ $_[0]->{ VERSION } || [] } );
    $_ =~ s/ GI:/  GI:/g;
    output_record( 'VERSION', $_ ) if length $_;
}


sub write_dblink
{
    local $_ = join( "\r", @{ $_[0]->{ DBLINK } || [] } );
    output_record( 'DBLINK', $_ ) if length $_;
}


sub write_dbsource
{
    local $_ = join( "\r", @{ $_[0]->{ DBSOURCE } || [] } );
    output_record( 'DBSOURCE', $_ ) if length $_;
}


sub write_keywords
{
    local $_ = join( '; ', @{ $_[0]->{ KEYWORDS } || [] } );
    $_ .= '.' if ! /\.$/;
    output_record( 'KEYWORDS', $_ );
}


sub write_source
{
    local $_ = $_[0]->{ SOURCE } || 'unknown';
    output_record( 'SOURCE', $_ );
}


sub write_organism
{
    local $_ = $_[0]->{ ORGANISM } || 'unclasified';
    output_record( '  ORGANISM', $_ );
}


sub write_taxonomy
{
    local $_ = join( '; ', @{ $_[0]->{ TAXONOMY } || [] } );
    $_ .= '.' if ! /\.$/;
    output_record( ' ', $_ );
}


sub write_references
{
    my $n = 0;
    foreach ( @{ $_[0]->{ REFERENCES } || [] } )
    {
        write_reference( $_, ++$n );
        write_authors( $_ );
        write_consortium( $_ );
        write_title( $_ );
        write_journal( $_ );
        write_pubmed( $_ );
        write_remark( $_ );
    }
}


sub write_reference
{
    local $_ = $_[0]->{ ref_num } || $_[1] || '1';
    $_ .= "  ($_[0]->{ref_note})" if $_[0]->{ ref_note };
    output_record( 'REFERENCE', $_ );
}


sub write_authors
{
    local $_ = $_[0]->{ AUTHORS };
    output_record( '  AUTHORS', $_ ) if defined $_ && length $_;
}


sub write_consortium
{
    local $_ = $_[0]->{ CONSRTM };
    output_record( '  CONSRTM', $_ ) if defined $_ && length $_;
}


sub write_title
{
    local $_ = $_[0]->{ TITLE };
    output_record( '  TITLE', $_ ) if defined $_ && length $_;
}


sub write_journal
{
    local $_ = $_[0]->{ JOURNAL };
    output_record( '  JOURNAL', $_ ) if defined $_ && length $_;
}


sub write_pubmed
{
    local $_ = $_[0]->{ PUBMED };
    output_record( '   PUBMED', $_ ) if defined $_ && length $_;
}


sub write_remark
{
    local $_ = $_[0]->{ REMARK };
    output_record( '  REMARK', $_ ) if defined $_ && length $_;
}


sub write_comment
{
    local $_ = join( "\r", @{ $_[0]->{ COMMENT } || [] } );
    output_record( 'COMMENT', $_ ) if length $_;
}


sub write_features
{
    my $ftr_list = feature_list( $_[0] );
    if ( $ftr_list && @$ftr_list )
    {
        write_feature_header();
        foreach ( @$ftr_list ) { write_feature( $_ ) }
    }
}


sub write_feature_header
{
    print "FEATURES             Location/Qualifiers\n";
}


sub write_feature
{
    my ( $type, $location, $qualifiers ) = @{ $_[0] };

    $^A = '';

    my ( $loc, $rest ) = wrap_location( $location );

    formline( <<'End_of_line_1', $type, $loc );
     @<<<<<<<<<<<<<< ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
End_of_line_1


    while ( $rest )
    {
        ( $loc, $rest ) = wrap_location( $rest );
        formline( <<'End_of_continuation', $loc ) if defined $loc && length( $loc );
                     ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
End_of_continuation
    }

    if ( defined $qualifiers && ref $qualifiers eq 'ARRAY' )
    {
        foreach ( @$qualifiers )
        {
            add_feature_qualifier( @$_ ) if $_->[0] !~ /^_/;
        }
    }
    elsif ( defined $qualifiers && ref $qualifiers eq 'HASH' )
    {
        foreach my $key ( ordered_qualifiers( $qualifiers ) )
        {
            foreach ( @{ $qualifiers->{ $key } } ) { add_feature_qualifier( $key, $_ ) }
        }
    }

    print $^A;
    $^A = '';
}


sub wrap_location
{
    local $_ = shift;
    return ( $_, '' ) if length( $_ ) <= 58;

    return ( $1, $2 ) if m/^(.{1,57},)(.*)$/;
    return ( $1, $2 ) if m/^(.{1,56}\.\.)(.*)$/;
    return ( $1, $2 ) if m/^(.{1,57}\()(.*)$/;
    return $_;
}


sub ordered_qualifiers
{
    my $qualifiers = shift;

    if ( ! $qualifiers->{ '_qualif_order' } )
    {
        @{ $qualifiers->{ '_qualif_order' } } = sort { qualifier_priority( $a ) <=> qualifier_priority( $b )
                                                    || lc $a cmp lc $b
                                                     }
                                                grep { ! /^_/ }
                                                keys %$qualifiers;
    }

    @{ $qualifiers->{ '_qualif_order' } };
}


sub qualifier_priority { $qualifier_order{ lc $_[0] } || 999999 }


sub add_feature_qualifier
{
    my ( $key, $value ) = @_;
    $key or return;
    $value = '' if ! defined $value;

    my $qualif;
    if ( $valueless_qual{$key} )
    { 
        $qualif = "/$key";
    }
    elsif ( $unquoted_qual{$key} )
    {
        $qualif = "/key=$value";
    }
    else
    {
        $value =~ s/"/""/g;
        $qualif = qq(/$key="$value");
    }

    formline( <<'End_of_continuation', $qualif ) if defined $qualif && length $qualif;
~~                   ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
End_of_continuation
}


#  2011/01/23 -- GJO
#
#  Not required, and introduces more problems than it solves.
#
sub write_base_count
{
    return;
#
#    my $entry = shift;
#    return if $entry->{ is_protein };
#    return if $entry->{ mol_type } =~ /\S/ && $entry->{ mol_type } !~ m/NA$/i;
#
#    my ( $a, $c, $g, $t, $o );
#    my $counts = $entry->{ BASE_COUNT };
#    if ( $counts && ref $counts eq 'HASH' && keys %$counts)
#    {
#        ( $a, $c, $g, $t, $o ) = map { $counts->{ $_ } || 0 }
#                                 qw( a c g t other );
#        $a || $c || $g || $t || $o or return;
#    }
#    else
#    {
#        $entry->{ SEQUENCE } or return;
#        my $s = \$entry->{ SEQUENCE };
#        $a = $$s =~ tr/Aa//;
#        $c = $$s =~ tr/Cc//;
#        $g = $$s =~ tr/Gg//;
#        $t = $$s =~ tr/Tt//;
#        $o = length( $$s ) - ( $a + $c + $g + $t );
#        $a || $c || $g || $t || $o or return;
#
#        $entry->{ BASE_COUNT } = { a => $a,
#                                   c => $c,
#                                   g => $g,
#                                   t => $t
#                                 };
#        $entry->{ BASE_COUNT }->{ other } = $o if $o;
#    }
#
#    printf "BASE COUNT  %6d a %6d c %6d g %6d t%s\n",
#           $a, $c, $g, $t, $o ? sprintf( ' %6d other', $o ) : '';
}


sub write_origin
{
    my $entry = shift;
    my $info = $entry && defined $entry->{ ORIGIN } ? $entry->{ ORIGIN } : '';
    print "ORIGIN      $info\n";  # Use explicit print to get blank padding
}


sub write_sequence
{
    my $entry = shift;
    $entry->{ SEQUENCE } or return;
    
    my $start = 1;
    foreach ( $entry->{ SEQUENCE } =~ m/(.{1,60})/g )
    {
        print join( ' ', sprintf( '%9d', $start ), m/(.{1,10})/g ), "\n";
        $start += 60;
    }
}


sub write_end_of_entry
{
    print "//\n";
}


#===============================================================================
#  Helper function for defining an input filehandle:
#
#     filehandle is passed through
#     string is taken as file name to be opened
#     scalar reference is taken as data string to be opened
#     undef or "" defaults to STDIN
#
#      \*FH           = input_filehandle( $file );
#    ( \*FH, $close ) = input_filehandle( $file );
#
#===============================================================================
sub input_filehandle
{
    my $file = shift;

    #  Null string or undef

    if ( ! defined( $file ) || ( $file eq '' ) )
    {
        return wantarray ? ( \*STDIN, 0 ) : \*STDIN;
    }

    #  FILEHANDLE?

    if ( ref( $file ) eq "GLOB" )
    {
        return wantarray ? ( $file, 0 ) : $file;
    }

    #  File name

    if ( ! ref( $file ) || ref( $file ) eq 'SCALAR' )
    {
        ref( $file ) or -f $file or die "Could not find input file '$file'.\n";
        my $fh;
        open( $fh, "<", $file ) || die "Could not open '$file' for input.\n";
        return wantarray ? ( $fh, 1 ) : $fh;
    }

    return wantarray ? ( \*STDIN, undef ) : \*STDIN;
}


#===============================================================================
#  Helper function for defining an output filehandle:
#
#     filehandle is passed through
#     string is taken as file name to be opened
#     scalar reference is taken as data string to be opened
#     undef or "" defaults to STDOUT
#
#      \*FH           = output_filehandle( $file );
#    ( \*FH, $close ) = output_filehandle( $file );
#
#===============================================================================
sub output_filehandle
{
    my $file = shift;

    #  Null string or undef

    if ( ! defined( $file ) || ( $file eq '' ) )
    {
        return wantarray ? ( \*STDOUT, 0 ) : \*STDOUT;
    }

    #  FILEHANDLE?

    if ( ref( $file ) eq "GLOB" )
    {
        return wantarray ? ( $file, 0 ) : $file;
    }

    #  File name

    if ( ! ref( $file ) || ref( $file ) eq 'SCALAR' )
    {
        my $fh;
        open( $fh, ">", $file ) || die "Could not open '$file' for output.\n";
        return wantarray ? ( $fh, 1 ) : $fh;
    }

    return wantarray ? ( \*STDOUT, undef ) : \*STDOUT;
}


1;
