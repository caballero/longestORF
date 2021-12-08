#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $fasta_file = undef;
my $out_prot   = 'sequences.pep';
my $verbose    = undef;

GetOptions(
    'fasta|f=s'     => \$fasta_file,
    'proteins|p=s'  => \$out_prot,
    'verbose|v'     => \$verbose
);

printHelp() unless (defined $fasta_file);

my %nuc2aa = (
'TCA' => 'S', # Serine
'TCC' => 'S', # Serine
'TCG' => 'S', # Serine
'TCT' => 'S', # Serine
'TTC' => 'F', # Phenylalanine
'TTT' => 'F', # Phenylalanine
'TTA' => 'L', # Leucine
'TTG' => 'L', # Leucine
'TAC' => 'Y', # Tyrosine
'TAT' => 'Y', # Tyrosine
'TAA' => '*', # Stop
'TAG' => '*', # Stop
'TGC' => 'C', # Cysteine
'TGT' => 'C', # Cysteine
'TGA' => '*', # Stop
'TGG' => 'W', # Tryptophan
'CTA' => 'L', # Leucine
'CTC' => 'L', # Leucine
'CTG' => 'L', # Leucine
'CTT' => 'L', # Leucine
'CCA' => 'P', # Proline
'CAT' => 'H', # Histidine
'CAA' => 'Q', # Glutamine
'CAG' => 'Q', # Glutamine
'CGA' => 'R', # Arginine
'CGC' => 'R', # Arginine
'CGG' => 'R', # Arginine
'CGT' => 'R', # Arginine
'ATA' => 'I', # Isoleucine
'ATC' => 'I', # Isoleucine
'ATT' => 'I', # Isoleucine
'ATG' => 'M', # Methionine
'ACA' => 'T', # Threonine
'ACC' => 'T', # Threonine
'ACG' => 'T', # Threonine
'ACT' => 'T', # Threonine
'AAC' => 'N', # Asparagine
'AAT' => 'N', # Asparagine
'AAA' => 'K', # Lysine
'AAG' => 'K', # Lysine
'AGC' => 'S', # Serine
'AGT' => 'S', # Serine
'AGA' => 'R', # Arginine
'AGG' => 'R', # Arginine
'CCC' => 'P', # Proline
'CCG' => 'P', # Proline
'CCT' => 'P', # Proline
'CAC' => 'H', # Histidine
'GTA' => 'V', # Valine
'GTC' => 'V', # Valine
'GTG' => 'V', # Valine
'GTT' => 'V', # Valine
'GCA' => 'A', # Alanine
'GCC' => 'A', # Alanine
'GCG' => 'A', # Alanine
'GCT' => 'A', # Alanine
'GAC' => 'D', # Aspartic Acid
'GAT' => 'D', # Aspartic Acid
'GAA' => 'E', # Glutamic Acid
'GAG' => 'E', # Glutamic Acid
'GGA' => 'G', # Glycine
'GGC' => 'G', # Glycine
'GGG' => 'G', # Glycine
'GGT' => 'G'  # Glycine
);

open (my $op, ">", $out_prot) or die "cannot write $out_prot\n";

warn "Loading sequence from Fasta file\n" if ($verbose);
my $oseq = '';
my $oid = '';
open (my $fh, "<", $fasta_file) or die;
while (<$fh>) {
    chomp;
    if (/^>/) {
        s/>//;
        $oid = $_;
        $oid =~ s/\s+.*//;
    }
    else { 
        $oseq .= $_;
    }
} 
close $fh;

my $frame_p1 = $oseq;
my $frame_p2 = $oseq; $frame_p2 =~ s/^\w//; 
my $frame_p3 = $oseq; $frame_p3 =~ s/^\w\w//;
my $frame_m1 = revseq($oseq);
my $frame_m2 = revseq($oseq); $frame_m2 =~ s/^\w//;
my $frame_m3 = revseq($oseq); $frame_m3 =~ s/^\w\w//;

my $orfs_p1 = translate($frame_p1);
my $orfs_p2 = translate($frame_p2);
my $orfs_p3 = translate($frame_p3);
my $orfs_m1 = translate($frame_m1);
my $orfs_m2 = translate($frame_m2);
my $orfs_m3 = translate($frame_m3);

my $longest = getLongestORF($orfs_p1, $orfs_p2, $orfs_p3, $orfs_m1, $orfs_m2, $orfs_m3);
print $op ">$oid | FRAME +1\n$orfs_p1\n";
print $op ">$oid | FRAME +2\n$orfs_p2\n";
print $op ">$oid | FRAME +3\n$orfs_p3\n";
print $op ">$oid | FRAME -1\n$orfs_m1\n";
print $op ">$oid | FRAME -2\n$orfs_m2\n";
print $op ">$oid | FRAME -3\n$orfs_m3\n";
print $op ">$oid | LONGEST ORF\n$longest\n";

close $op;

warn "All done\n" if ($verbose);

#################
## Subroutines ##
#################

sub translate {
    my $nuc = shift @_;
    my $pep = '';
    for (my $i = 1; $i <= length($nuc) + 1; $i += 3) {
        my $codon = substr($nuc, $i - 1, 3);
        if (defined $nuc2aa{$codon}) {
            $pep .= $nuc2aa{$codon};
        }
    }
    return $pep;
}

sub revseq {
    my $seq = shift @_;
    my $rev = reverse $seq;
    $rev =~ tr/ACGTacgt/TGCAtgca/;
    return $rev;
}

sub getLongestORF {
    my $longest = '';
    foreach my $frame (@_) {
        while ($frame =~ /(M\w+\*)/g) {
            my $orf = $1;
            if (length($orf) > length($longest)) {
                $longest = $orf;
            }
        }
    }
    return $longest;
}

sub printHelp {

    print <<__HELP__

script: reversion.pl [PARAM]

Parameters
    -f | --fasta      Input file with sequence in Fasta format
    -p | --proteins   Output file for protein sequences, default: $out_prot
    -v | --verbose    Verbose mode

__HELP__
;
    exit;
}