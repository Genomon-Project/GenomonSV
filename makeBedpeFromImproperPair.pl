#! /usr/local/bin/perl

use strict;

my $input = $ARGV[0];
my $junction_dist = 800;

my $tempID = "";
my $tempPairNum = "";
my $tempChr = "";
my $tempStart = "";
my $tempEnd = "";
my $tempDir = "";
my $tempMapQ = "";

open(IN, $input) || die "cannot open $!";
while(<IN>) {
    s/[\r\n]//g;
    my @F = split("\t", $_);
    $F[0] =~ s/\/([12])//;
    my $pairNum = $1;

    if ($F[0] eq $tempID and $tempPairNum == 1 and $pairNum == 2) {

        my ($chr1, $dir1, $start1, $end1, $mapQ1, $align1) = ($tempChr, $tempDir, 0, 0, $tempMapQ, $tempChr . ":" . $tempStart . "-" . $tempEnd);
        if ($dir1 eq "+") {
            $start1 = $tempEnd;
            $end1 = $tempEnd + $junction_dist;
        } else {
            $start1 = $tempStart - $junction_dist;
            $end1 = $tempStart;
        }
        
        my ($chr2, $dir2, $start2, $end2, $mapQ2, $align2) = ($F[1], $F[4], 0, 0, $F[5], $F[1] . ":" . $F[2] . "-" . $F[3]);
        if ($dir2 eq "+") {
            $start2 = $F[3];
            $end2 = $F[3] + $junction_dist;
        } else {
            $start2 = $F[2] - $junction_dist;
            $end2 = $F[2];
        }

        if (($chr1 cmp $chr2) == -1) {
            print $chr1 . "\t" . $start1 . "\t" . $end1 . "\t" . $chr2 . "\t" . $start2 . "\t" . $end2 . "\t" . $tempID . "\t" . $mapQ1 . "," . $mapQ2 . "\t" . $dir1 . "\t" . $dir2 . "\t" . $align1 . "," . $align2 . "\n";
        } elsif (($chr1 cmp $chr2) == 1) {
            print $chr2 . "\t" . $start2 . "\t" . $end2 . "\t" . $chr1 . "\t" . $start1 . "\t" . $end1 . "\t" . $tempID . "\t" . $mapQ2 . "," . $mapQ1 . "\t" . $dir2 . "\t" . $dir1 . "\t" . $align2 . "," . $align1 . "\n"; 
        } else {
            if ($start1 <= $start2) {
                print $chr1 . "\t" . $start1 . "\t" . $end1 . "\t" . $chr2 . "\t" . $start2 . "\t" . $end2 . "\t" . $tempID . "\t" . $mapQ1 . "," . $mapQ2 . "\t" . $dir1 . "\t" . $dir2 . "\t" . $align1 . "," . $align2 . "\n";
            } else {
                print $chr2 . "\t" . $start2 . "\t" . $end2 . "\t" . $chr1 . "\t" . $start1 . "\t" . $end1 . "\t" . $tempID . "\t" . $mapQ2 . "," . $mapQ1 . "\t" . $dir2 . "\t" . $dir1 . "\t" . $align2 . "," . $align1 . "\n";
            }
        }

    }

    ($tempID, $tempPairNum, $tempChr, $tempStart, $tempEnd, $tempDir, $tempMapQ) = ($F[0], $pairNum, $F[1], $F[2], $F[3], $F[4], $F[5]);
}
close(IN);

