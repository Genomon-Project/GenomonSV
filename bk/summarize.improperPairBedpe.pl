#! /usr/local/bin/perl

use strict;
use List::Util qw(min max);

my $input = $ARGV[0];

my $junction_dist = 3000;

my %mergedBedpe = ();
open(IN, $input) || die "cannot open $!";
while(<IN>) {
    s/[\r\n]//g;
    my @F = split("\t", $_);

    if ($F[6] eq "HWI-ST1021:119:C14DVACXX:1:1106:19843:147527") {
        $DB::single = 1;
    }

    my $match = 0;
    foreach my $key (sort keys %mergedBedpe) {

        my ($tchr1, $tstart1, $tend1, $tchr2, $tstart2, $tend2, $tdir1, $tdir2) = split("\t", $key);
        my ($tids, $tmqs, $talns)  = split("\t", $mergedBedpe{$key});

        # if ($tids =~ /HWI-ST1021:119:C14DVACXX:3:1102:4871:5662/) {
        # if ($F[5] eq "224202785") {
        #     $DB::single = 1;
        # }

        if ($F[0] ne $tchr1 or $F[1] > $tend1 + $junction_dist) {
            my @talns_a = split(";", $talns);
            my %count = ();
            @talns_a = grep {!$count{$_}++} @talns_a;

            if ($#talns_a >= 0) {

                print $tchr1 . "\t" . $tstart1 . "\t" . $tend1 . "\t" . $tchr2 . "\t" . $tstart2 . "\t" . $tend2 . "\t";
                print $tids . "\t" . $tmqs . "\t" . $tdir1 . "\t" . $tdir2 . "\t" . $talns . "\n";
                delete $mergedBedpe{$key};
                next;

            }

        } else {

            if ($F[0] eq $tchr1 and $F[3] eq $tchr2 and $F[8] eq $tdir1 and $F[9] eq $tdir2) {
                # condition for the overlap
                if (($F[2] > $tstart1 and $F[1] <= $tend1) and ($F[5] > $tstart2 and $F[4] <= $tend2)) { 
            
                    $match = 1;
                    my $newStart1 = max($tstart1, $F[1]);
                    my $newEnd1 = min($tend1, $F[2]);
                    my $newStart2 = max($tstart2, $F[4]);
                    my $newEnd2 = min($tend2, $F[5]); 

                    my $newKey = $tchr1 . "\t" . $newStart1 . "\t" . $newEnd1 . "\t" . $tchr2 . "\t" . $newStart2 . "\t" . $newEnd2 . "\t" . $tdir1 . "\t" . $tdir2;
                    my $newIds = $tids . ";" . $F[6];
                    my $newMqs = $tmqs . ";" . $F[7];
                    my $newAlns = $talns . ";" . $F[10];

                    delete $mergedBedpe{$key};
                    $mergedBedpe{$newKey} = $newIds . "\t" . $newMqs . "\t" . $newAlns;
                    last;

                }
            }
    
        }

    }

    if ($match == 0) {
        my $newKey = $F[0] . "\t" . $F[1] . "\t" . $F[2] . "\t" . $F[3] . "\t" . $F[4] . "\t" . $F[5] . "\t" . $F[8] . "\t" . $F[9];
        $mergedBedpe{$newKey} = $F[6] . "\t" . $F[7] . "\t" . $F[10];
    }

}
close(IN);

