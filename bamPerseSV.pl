#! /usr/local/bin/perl

use strict;

my $input = $ARGV[0];

my $junctionn_dist = 800;
my $abnormal_insert_size = 2000;
my $min_clip_size_for_primary = 20;

open(IN, $input) || die "cannot open $!";
while(<IN>) {

    s/[\r\n]//g;
    my @F = split("\t", $_);
    my @flags = reverse(split("", sprintf("%012b", $F[1])));



    my ($chr_current, $pos_current, $dir_current, $cigar_current, $seq_current)  = ($F[2], $F[3], ($flags[4] == 0 ? "+" : "-"), $F[5], $F[9]);
    my ($chr_pair, $pos_pair, $dir_pair) = (($F[6] eq "=" ? $F[2] : $F[6]), $F[7], ($flags[5] == 0 ? "+" : "-"));
    my $seqID = $F[0] . "/" . ($flags[6] eq "1" ? 1 : 2);
    my $mapQ = $F[4];

    # skip supplementary alignment
    next if ($flags[8] eq "1" or $flags[11] eq "1");

    # skip reads aligned to "hs37d5"
    next if ($chr_current eq "hs37d5" or $chr_pair eq "hs37d5");

    # skip if both directions are clipped
    next if ($cigar_current =~ /(\d+)[HS]$/ and $cigar_current =~ /^(\d+)[HS]/);
 
    # check when the junction is in + direction 
    if (join("\t", @F[11 .. $#F]) =~ /SA:Z:/ and $cigar_current =~ /(\d+)[HS]$/) {

        my $clipLen_current = $1;   
        if ($clipLen_current >= $min_clip_size_for_primary) {

            join("\t", @F[11 .. $#F]) =~ /SA:Z:(\w+),(\d+),([\-\+]),(\w+),(\d+)/;
            my ($chr_SA, $pos_SA, $dir_SA, $cigar_SA) = ($1, $2, $3, $4);

            # junction information from the side of the sequence alignment information at the current record.
            my $juncChr_currennt = $chr_current;
            my $juncDir_current = "+";
            my $alignmentSize_current = 0;
            while($cigar_current =~ /(\d+)[MD]/g) {
                $alignmentSize_current = $alignmentSize_current + $1;
            }
            my $juncPos_current = $pos_current + $alignmentSize_current - 1;

            my $juncChr_SA = $chr_SA;
            my $juncDir_SA = "";
            my $juncPos_SA = 0;
            my $juncSurplus = "---";

            # checking the expected clipping sizes and directions. 
            my $readLength_current = 0;
            while($cigar_current =~ /(\d+)[HIMS]/g) {
                $readLength_current = $readLength_current + $1;
            }
            my $expected_clipLen_supp = $readLength_current - $clipLen_current;
            my $expected_clipDir_supp = $dir_current eq $dir_SA ? "-" : "+";

            my $alignmentSize_SA = 0;
            while($cigar_SA =~ /(\d+)[MD]/g) {
                $alignmentSize_SA = $alignmentSize_SA + $1;
            }

            my $clipLen_supp = 0;

            # support read pattern B or C
            if ($dir_current eq "-" and $dir_pair eq "+" and $chr_current eq $chr_pair and $pos_current - $pos_pair >= 0 and $pos_current - $pos_pair < $abnormal_insert_size) {

                if ($dir_SA eq "+" and $expected_clipDir_supp eq "+" and $cigar_SA =~ /(\d+)[HS]$/) {
                    $clipLen_supp = $1;
                    $juncDir_SA = "+";
                    $juncPos_SA = $pos_SA + ($alignmentSize_SA - 1);

                    if ($clipLen_supp < $expected_clipLen_supp) {
                        $juncPos_SA = $juncPos_SA - ($expected_clipLen_supp - $clipLen_supp);
                    }

                    if ($clipLen_supp > $expected_clipLen_supp and $readLength_current == length($seq_current)) {
                        $juncSurplus = substr($seq_current, $readLength_current - $clipLen_current, $clipLen_supp - $expected_clipLen_supp);
                    }
   
                    print $juncChr_currennt . "\t" . ($juncPos_current - 1) . "\t" . $juncPos_current . "\t";
                    print $juncChr_SA . "\t" . ($juncPos_SA - 1) . "\t" . $juncPos_SA . "\t";
                    print $seqID . "\t" . $mapQ . "\t" . $juncDir_current . "\t" . $juncDir_SA . "\t";
                    print $chr_current . ":" . $pos_current . "-" . $juncPos_current . "," . $chr_SA . ":" . $pos_SA . "-" . $juncPos_SA . "\t" . $juncSurplus . "\t";
                    print $chr_pair . ":" . $pos_pair . "\n"; 
                }

                if ($dir_SA eq "-" and $expected_clipDir_supp eq "-" and $cigar_SA =~ /^(\d+)[HS]/) {
                    $clipLen_supp = $1;
                    $juncDir_SA = "-";
                    $juncPos_SA = $pos_SA;                    

                    if ($clipLen_supp < $expected_clipLen_supp) {
                        $juncPos_SA = $juncPos_SA + ($expected_clipLen_supp - $clipLen_supp);
                    }

                    if ($clipLen_supp > $expected_clipLen_supp and $readLength_current == length($seq_current)) {
                        $juncSurplus = substr($seq_current, $readLength_current - $clipLen_current, $clipLen_supp - $expected_clipLen_supp);
                    }

                    print $juncChr_currennt . "\t" . ($juncPos_current - 1) . "\t" . $juncPos_current . "\t";
                    print $juncChr_SA . "\t" . ($juncPos_SA - 1) . "\t" . $juncPos_SA . "\t";
                    print $seqID . "\t" . $mapQ . "\t" . $juncDir_current . "\t" . $juncDir_SA . "\t";
                    print $chr_current . ":" . $pos_current . "-" . $juncPos_current . "," . $chr_SA . ":" . $juncPos_SA . "-" . ($juncPos_SA + ($alignmentSize_SA - 1)) . "\t" . $juncSurplus . "\t";
                    print $chr_pair . ":" . $pos_pair . "\n";
                }
    
            }

            # support read pattern G or H
            if ($dir_current eq "+" and $chr_SA eq $chr_pair) {

                if ($dir_SA eq "-" and $dir_pair eq "+" and $pos_SA - $pos_pair >= 0 and $pos_SA - $pos_pair < $abnormal_insert_size and $expected_clipDir_supp eq "+" and $cigar_SA =~ /(\d+)[HS]$/) {
                    $clipLen_supp = $1;
                    $juncDir_SA = "+";
                    $juncPos_SA = $pos_SA + ($alignmentSize_SA - 1);

                    if ($clipLen_supp < $expected_clipLen_supp) {
                        $juncPos_SA = $juncPos_SA - ($expected_clipLen_supp - $clipLen_supp);
                    }

                    if ($clipLen_supp > $expected_clipLen_supp and $readLength_current == length($seq_current)) {
                        $juncSurplus = substr($seq_current, $readLength_current - $clipLen_current, $clipLen_supp - $expected_clipLen_supp);
                    }

                    print $juncChr_currennt . "\t" . ($juncPos_current - 1) . "\t" . $juncPos_current . "\t";
                    print $juncChr_SA . "\t" . ($juncPos_SA - 1) . "\t" . $juncPos_SA . "\t";
                    print $seqID . "\t" . $mapQ . "\t" . $juncDir_current . "\t" . $juncDir_SA . "\t";
                    print $chr_current . ":" . $pos_current . "-" . $juncPos_current . "," . $chr_SA . ":" . $pos_SA . "-" . $juncPos_SA . "\t" . $juncSurplus . "\t";
                    print $chr_pair . ":" . $pos_pair . "\n";

                }

                if ($dir_SA eq "+" and $dir_pair eq "-" and $pos_pair - $pos_SA >= 0 and $pos_pair - $pos_SA < $abnormal_insert_size and $expected_clipDir_supp eq "-" and $cigar_SA =~ /^(\d+)[HS]/) {

                    $clipLen_supp = $1;
                    $juncDir_SA = "-";
                    $juncPos_SA = $pos_SA;
 
                    if ($clipLen_supp < $expected_clipLen_supp) {
                        $juncPos_SA = $juncPos_SA + ($expected_clipLen_supp - $clipLen_supp);
                    }

                    if ($clipLen_supp > $expected_clipLen_supp and $readLength_current == length($seq_current)) {
                        $juncSurplus = substr($seq_current, $readLength_current - $clipLen_current, $clipLen_supp - $expected_clipLen_supp);
                    }

                    print $juncChr_currennt . "\t" . ($juncPos_current - 1) . "\t" . $juncPos_current . "\t";
                    print $juncChr_SA . "\t" . ($juncPos_SA - 1) . "\t" . $juncPos_SA . "\t";
                    print $seqID . "\t" . $mapQ . "\t" . $juncDir_current . "\t" . $juncDir_SA . "\t";
                    print $chr_current . ":" . $pos_current . "-" . $juncPos_current . "," . $chr_SA . ":" . $juncPos_SA . "-" . ($juncPos_SA + ($alignmentSize_SA - 1)) . "\t" . $juncSurplus . "\t";
                    print $chr_pair . ":" . $pos_pair . "\n";

                }
            }

        }

    }

    # check when the junction is in - direction 
    if (join("\t", @F[11 .. $#F]) =~ /SA:Z:/ and $cigar_current =~ /^(\d+)[HS]/) {

        my $clipLen_current = $1;
        if ($clipLen_current >= $min_clip_size_for_primary) {

            join("\t", @F[11 .. $#F]) =~ /SA:Z:(\w+),(\d+),([\-\+]),(\w+),(\d+)/;
            my ($chr_SA, $pos_SA, $dir_SA, $cigar_SA) = ($1, $2, $3, $4);

            # junction information from the side of the sequence alignment information at the current record.
            my $juncChr_currennt = $chr_current;
            my $juncDir_current = "-";  # discrepancy between + and - junction
            my $alignmentSize_current = 0;
            while($cigar_current =~ /(\d+)[MD]/g) {
                $alignmentSize_current = $alignmentSize_current + $1;
            }
            my $juncPos_current = $pos_current; #discrepancy between + and - junction

            my $juncChr_SA = $chr_SA;
            my $juncDir_SA = "";
            my $juncPos_SA = 0;
            my $juncSurplus = "---";

            # checking the expected clipping sizes and directions. 
            my $readLength_current = 0;
            while($cigar_current =~ /(\d+)[HIMS]/g) {
                $readLength_current = $readLength_current + $1;
            }
            my $expected_clipLen_supp = $readLength_current - $clipLen_current;
            my $expected_clipDir_supp = $dir_current eq $dir_SA ? "+" : "-"; # discrepancy between + and - junction

            my $alignmentSize_SA = 0;
            while($cigar_SA =~ /(\d+)[MD]/g) {
                $alignmentSize_SA = $alignmentSize_SA + $1;
            }

            my $clipLen_supp = 0;

            # support read pattern B or C
            if ($dir_current eq "+" and $dir_pair eq "-" and $chr_current eq $chr_pair and $pos_pair - $pos_current >= 0 and $pos_pair - $pos_current < $abnormal_insert_size) {

                if ($dir_SA eq "+" and $expected_clipDir_supp eq "+" and $cigar_SA =~ /(\d+)[HS]$/) {
                    $clipLen_supp = $1;
                    $juncDir_SA = "+";
                    $juncPos_SA = $pos_SA + ($alignmentSize_SA - 1);

                    if ($clipLen_supp < $expected_clipLen_supp) {
                        $juncPos_SA = $juncPos_SA - ($expected_clipLen_supp - $clipLen_supp);
                    }

                    if ($clipLen_supp > $expected_clipLen_supp and $readLength_current == length($seq_current)) {
                        # $juncSurplus = substr($seq_current, $readLength_current - $clipLen_current, $clipLen_supp - $expected_clipLen_supp);
                        $juncSurplus = substr($seq_current, $clipLen_current - ($clipLen_supp - $expected_clipLen_supp) - 1, $clipLen_supp - $expected_clipLen_supp);
                    }

                    print $juncChr_currennt . "\t" . ($juncPos_current - 1) . "\t" . $juncPos_current . "\t";
                    print $juncChr_SA . "\t" . ($juncPos_SA - 1) . "\t" . $juncPos_SA . "\t";
                    print $seqID . "\t" . $mapQ . "\t" . $juncDir_current . "\t" . $juncDir_SA . "\t";
                    print $chr_current . ":" . $pos_current . "-" . ($pos_current + $alignmentSize_current - 1) . "," . $chr_SA . ":" . $pos_SA . "-" . $juncPos_SA . "\t" . $juncSurplus . "\t";
                    print $chr_pair . ":" . $pos_pair . "\n";

                }
 
                if ($dir_SA eq "-" and $expected_clipDir_supp eq "-" and $cigar_SA =~ /^(\d+)[HS]/) {
                    $clipLen_supp = $1;
                    $juncDir_SA = "-";
                    $juncPos_SA = $pos_SA;

                    if ($clipLen_supp < $expected_clipLen_supp) {
                        $juncPos_SA = $juncPos_SA + ($expected_clipLen_supp - $clipLen_supp);
                    }

                    if ($clipLen_supp > $expected_clipLen_supp and $readLength_current == length($seq_current)) {
                        # $juncSurplus = substr($seq_current, $readLength_current - $clipLen_current, $clipLen_supp - $expected_clipLen_supp);
                        $juncSurplus = substr($seq_current, $clipLen_current - ($clipLen_supp - $expected_clipLen_supp) - 1, $clipLen_supp - $expected_clipLen_supp);
                    }

                    print $juncChr_currennt . "\t" . ($juncPos_current - 1) . "\t" . $juncPos_current . "\t";
                    print $juncChr_SA . "\t" . ($juncPos_SA - 1) . "\t" . $juncPos_SA . "\t";
                    print $seqID . "\t" . $mapQ . "\t" . $juncDir_current . "\t" . $juncDir_SA . "\t";
                    print $chr_current . ":" . $pos_current . "-" . ($pos_current + $alignmentSize_current - 1) . "," . $chr_SA . ":" . $juncPos_SA . "-" . ($juncPos_SA + ($alignmentSize_SA - 1)) . "\t" . $juncSurplus . "\t";
                    print $chr_pair . ":" . $pos_pair . "\n";

                }
                
            }

            # support read pattern G or H
            if ($dir_current eq "-" and $chr_SA eq $chr_pair) {

                if ($dir_SA eq "-" and $dir_pair eq "+" and $pos_SA - $pos_pair >= 0 and $pos_SA - $pos_pair < $abnormal_insert_size and $expected_clipDir_supp eq "+" and $cigar_SA =~ /(\d+)[HS]$/) {
                    $clipLen_supp = $1;
                    $juncDir_SA = "+";
                    $juncPos_SA = $pos_SA + ($alignmentSize_SA - 1);

                    if ($clipLen_supp < $expected_clipLen_supp) {
                        $juncPos_SA = $juncPos_SA - ($expected_clipLen_supp - $clipLen_supp);
                    }

                    if ($clipLen_supp > $expected_clipLen_supp and $readLength_current == length($seq_current)) {
                        # $juncSurplus = substr($seq_current, $readLength_current - $clipLen_current, $clipLen_supp - $expected_clipLen_supp);
                        $juncSurplus = substr($seq_current, $clipLen_current - ($clipLen_supp - $expected_clipLen_supp) - 1, $clipLen_supp - $expected_clipLen_supp);
                    }

                    print $juncChr_currennt . "\t" . ($juncPos_current - 1) . "\t" . $juncPos_current . "\t";
                    print $juncChr_SA . "\t" . ($juncPos_SA - 1) . "\t" . $juncPos_SA . "\t";
                    print $seqID . "\t" . $mapQ . "\t" . $juncDir_current . "\t" . $juncDir_SA . "\t";
                    print $chr_current . ":" . $pos_current . "-" . ($pos_current + $alignmentSize_current - 1) . "," . $chr_SA . ":" . $pos_SA . "-" . $juncPos_SA . "\t" . $juncSurplus . "\t";
                    print $chr_pair . ":" . $pos_pair . "\n";

                }

                if ($dir_SA eq "+" and $dir_pair eq "-" and $pos_pair - $pos_SA >= 0 and $pos_pair - $pos_SA < $abnormal_insert_size and $expected_clipDir_supp eq "-" and $cigar_SA =~ /^(\d+)[HS]/) {

                    $clipLen_supp = $1;
                    $juncDir_SA = "-";
                    $juncPos_SA = $pos_SA;

                    if ($clipLen_supp < $expected_clipLen_supp) {
                        $juncPos_SA = $juncPos_SA + ($expected_clipLen_supp - $clipLen_supp);
                    }

                    if ($clipLen_supp > $expected_clipLen_supp and $readLength_current == length($seq_current)) {
                        # $juncSurplus = substr($seq_current, $readLength_current - $clipLen_current, $clipLen_supp - $expected_clipLen_supp);
                        $juncSurplus = substr($seq_current, $clipLen_current - ($clipLen_supp - $expected_clipLen_supp) - 1, $clipLen_supp - $expected_clipLen_supp);
                    }
    
                    print $juncChr_currennt . "\t" . ($juncPos_current - 1) . "\t" . $juncPos_current . "\t";
                    print $juncChr_SA . "\t" . ($juncPos_SA - 1) . "\t" . $juncPos_SA . "\t";
                    print $seqID . "\t" . $mapQ . "\t" . $juncDir_current . "\t" . $juncDir_SA . "\t";
                    print $chr_current . ":" . $pos_current . "-" . ($pos_current + $alignmentSize_current - 1) . "," . $chr_SA . ":" . $juncPos_SA . "-" . ($juncPos_SA + ($alignmentSize_SA - 1)) . "\t" . $juncSurplus . "\t";
                    print $chr_pair . ":" . $pos_pair . "\n";

                }
            }

        }
 
    }

}
close(IN);
        
 
