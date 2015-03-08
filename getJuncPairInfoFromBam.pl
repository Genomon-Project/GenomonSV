#! /usr/local/bin/perl

use strict;

my $input_sam = $ARGV[0];
my $juncIDs = $ARGV[1];


my %ID2info = ();
open(IN2, $juncIDs) || die "cannot open $!";
while(<IN2>) {
    s/[\r\n]//g;
    my @F = split("\t", $_);
    $ID2info{$F[3]} = join("\t", @F);
}
close(IN);


open(IN1, $input_sam) || die "cannot open $!";
while(<IN1>) {
    s/[\r\n]//g;
    my @F = split("\t", $_);
    my @flags = reverse(split("", sprintf("%012b", $F[1])));
 
    my ($chr_current, $pos_current, $dir_current, $cigar_current, $seq_current)  = ($F[2], $F[3], ($flags[4] == 0 ? "+" : "-"), $F[5], $F[9]);
    my ($chr_pair, $pos_pair, $dir_pair) = (($F[6] eq "=" ? $F[2] : $F[6]), $F[7], ($flags[5] == 0 ? "+" : "-"));
    my $seqID = $F[0] . "/" . ($flags[6] eq "1" ? 1 : 2);
    my $mapQ = $F[4];

    # skip supplementary alignment
    next if ($flags[8] == 1 or $flags[11] == 1);

    if (exists $ID2info{$seqID}) {

        my $alignmentSize_current = 0;
        while($cigar_current =~ /(\d+)[MD]/g) {
            $alignmentSize_current = $alignmentSize_current + $1;
        }
        my $end_current = $pos_current + $alignmentSize_current - 1;


        print $ID2info{$seqID} . "\t" . $chr_current . ":" . $pos_current . "-" . $end_current . "\t" . $mapQ . "\n";
    }
}
close(IN1);


