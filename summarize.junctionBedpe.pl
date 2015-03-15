#! /usr/local/bin/perl

use strict;
use List::Util qw(min max);

my $input = $ARGV[0];
my $checkMarginSize = 1000;

my %mergedBedpe = ();
open(IN, $input) || die "cannot open $!";
while(<IN>) {
    s/[\r\n]//g;
    my @F = split("\t", $_);

    my $match = 0;
    foreach my $key (keys %mergedBedpe) {

        my ($tchr1, $tstart1, $tend1, $tchr2, $tstart2, $tend2, $tdir1, $tdir2, $inseqSize) = split("\t", $key);
        my ($tids, $tmqs1, $talns1, $tinseqs, $tmqs2, $talns2, $tpinds)  = split("\t", $mergedBedpe{$key});

        # if ($tids =~ /HWI-ST1021:119:C14DVACXX:3:1102:4871:5662/) {
        if ($F[5] eq "224202785") {
            $DB::single = 1;
        }

        if ($F[0] ne $tchr1 or $F[1] > $tend1 + 1000) {

            my @talns_a1 = split(";", $talns1);
            my @talns_a2 = split(";", $talns2);

            my %count = ();
            @talns_a1 = grep {!$count{$_}++} @talns_a1;
            my %count = (); 
            @talns_a2 = grep {!$count{$_}++} @talns_a2;

            if ($#talns_a1 >= 0 or $#talns_a2 >= 0) {
                print $tchr1 . "\t" . $tstart1 . "\t" . $tend1 . "\t" . $tchr2 . "\t" . $tstart2 . "\t" . $tend2 . "\t";
                print $tids . "\t" . $tmqs1 . "\t" . $tdir1 . "\t" . $tdir2 . "\t" . $talns1 . "\t" . $tinseqs . "\t" . $tmqs2 . "\t" . $talns2 . "\t" . $tpinds . "\n";
            }
            delete $mergedBedpe{$key};
            next;
        } else {

            my $flag = 0;
            if ($F[0] eq $tchr1 and $F[3] eq $tchr2 and $F[8] eq $tdir1 and $F[9] eq $tdir2) {

                if ($F[8] eq "+") {
                    my $expectedDiffSize = ($F[2] - $tend1) + (length($F[11]) - $inseqSize);
                    if (($F[9] eq "+" and $F[5] == $tend2 - $expectedDiffSize) or ($F[9] eq "-" and $F[5] == $tend2 + $expectedDiffSize)) {
                        $flag = 1;
                    }
                } else {
                    my $expectedDiffSize = ($F[2] - $tend1) + ($inseqSize - length($F[11]));
                    if (($F[9] eq "+" and $F[5] == $tend2 + $expectedDiffSize) or ($F[9] eq "-" and $F[5] == $tend2 - $expectedDiffSize)) {
                        $flag = 1;
                    }
                }

                if ($flag == 1) {
            
                    $match = 1;
                    my $newIds = $tids . ";" . $F[6];
                    my $newMqs1 = $tmqs1 . ";" . $F[7];
                    my $newAlns1 = $talns1 . ";" . $F[10];
                    my $newInseqs = $tinseqs . ";" . $F[11];
                    my $newPinds = $tpinds . ";" . $F[14];
                    my $newMqs2 = $tmqs2 . ";" . $F[12];
                    my $newAlns2 = $talns2 . ";" . $F[13];

                    $mergedBedpe{$key} = $newIds . "\t" . $newMqs1 . "\t" . $newAlns1 . "\t" . $newInseqs . "\t" . $newMqs2 . "\t" . $newAlns2 . "\t" . $newPinds;
                }

            }
    
        }

    }

    if ($match == 0) {
        my $newKey = $F[0] . "\t" . $F[1] . "\t" . $F[2] . "\t" . $F[3] . "\t" . $F[4] . "\t" . $F[5] . "\t" . $F[8] . "\t" . $F[9] . "\t" . length($F[11]);
        $mergedBedpe{$newKey} = $F[6] . "\t" . $F[7] . "\t" . join("\t", @F[10 .. 14]);
    }

}
close(IN);


