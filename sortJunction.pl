#! /usr/local/bin/perl

use strict;

my $input = $ARGV[0];

open(IN, $input) || die "cannot open $!";
while(<IN>) {
    s/[\r\n]//g;
    my ($chr1, $start1, $end1, $chr2, $start2, $end2, $ID, $mapQ, $dir1, $dir2, $align, $inseq, $pairPos) = split("\t", $_);

    next if ($ID =~ /_s$/);
    next if ($chr1 eq "hs37d5" or $chr2 eq "hs37d5");

    if (($chr1 cmp $chr2) == -1) {
        print $chr1 . "\t" . $start1 . "\t" . $end1 . "\t" . $chr2 . "\t" . $start2 . "\t" . $end2 . "\t" . $ID . "\t" . $mapQ . "\t" . $dir1 . "\t" . $dir2 . "\t" . $align . "\t" . $inseq . "\t" . $pairPos . "\t1" . "\n";
    } elsif (($chr1 cmp $chr2) == 1) {
        print $chr2 . "\t" . $start2 . "\t" . $end2 . "\t" . $chr1 . "\t" . $start1 . "\t" . $end1 . "\t" . $ID . "\t" . $mapQ . "\t" . $dir2 . "\t" . $dir1 . "\t" . $align . "\t" . $inseq . "\t" . $pairPos . "\t2" . "\n"; 
    } else {
        if ($start1 <= $start2) {
            print $chr1 . "\t" . $start1 . "\t" . $end1 . "\t" . $chr2 . "\t" . $start2 . "\t" . $end2 . "\t" . $ID . "\t" . $mapQ . "\t" . $dir1 . "\t" . $dir2 . "\t" . $align . "\t" . $inseq . "\t" . $pairPos . "\t1" . "\n";
        } else {
            print $chr2 . "\t" . $start2 . "\t" . $end2 . "\t" . $chr1 . "\t" . $start1 . "\t" . $end1 . "\t" . $ID . "\t" . $mapQ . "\t" . $dir2 . "\t" . $dir1 . "\t" . $align . "\t" . $inseq . "\t" . $pairPos . "\t2" . "\n";
        }
    }

}
close(IN);


