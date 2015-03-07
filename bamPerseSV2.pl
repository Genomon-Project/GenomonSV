#! /usr/local/bin/perl

use strict;

my $input = $ARGV[0];
my $junctionn_dist = 800;
my $abnormal_insert_size = 2000;
my $min_clip_size_for_primary = 20;
my $min_mapping_qual = 30;
my $max_clip_size_for_primary_read_pair = 5;


open(IN, $input) || die "cannot open $!";
while(<IN>) {

    s/[\r\n]//g;
    my @F = split("\t", $_);
    my @flags = reverse(split("", sprintf("%012b", $F[1])));

    my ($chr_current, $pos_current, $dir_current, $cigar_current, $seq_current)  = ($F[2], $F[3], ($flags[4] == 0 ? "+" : "-"), $F[5]);
    my ($chr_pair, $pos_pair, $dir_pair) = (($F[6] eq "=" ? $F[2] : $F[6]), $F[7], ($flags[5] == 0 ? "+" : "-"));
    my $seqID = $F[0] . "/" . ($flags[6] eq "1" ? 1 : 2) . ($flags[8] eq "1" or $flags[11] eq "1" ? "_s" : "");
    my $mapQ = $F[4];

      
    # skip unless the both the pair are mapped. 
    next if ($flags[2] == 1 or $flags[3] == 1);

    # skip unless premaliry read
    next if ($flags[8] == 1 or $flags[11] == 1);

    # skip if below the minimum mapping quality 
    next if ($mapQ < $min_mapping_qual);

    # skip if the clipping occurs in non-junction direction
    next if ($dir_current eq "+" and $cigar_current =~ /^(\d+)[HS]/ and $1 >= $max_clip_size_for_primary_read_pair);
    next if ($dir_current eq "-" and $cigar_current =~ /(\d+)[HS]$/ and $1 >= $max_clip_size_for_primary_read_pair);

    # if paired reads are not properly aligned
    # Here, not proper means
    # (I'm not sure whether the sam flags are perfect for detecting the following patterns...)
    # 1. paired reads are mapped on different chromosomes
    # 2. paired reads are mapped on the same chromosome, but very distantly ("distant" is specified by the parameter $abnormal_insert_size).
    # 3. paired reads are mapped on the same chromosome, 
    # but breaks the rule that the left-aligned read is aligned in the forward direction and the right-aligned read is in the reverse direciton.
    my $abnormal_pair = 0;
    $abnormal_pair = 1 if ($F[6] ne "=" or abs($F[8]) > $abnormal_insert_size);
    $abnormal_pair = 1 if (abs($F[8]) > $abnormal_insert_size);
    
    # check for the above pattern 3
    if ($flags[1] == 0 and $abnormal_pair == 0) {

        if (not ($F[7] - $F[3] > 0 and $flags[4] == 0 and $flags[5] == 1) and not ($F[3] - $F[7] > 0 and $flags[4] == 1 and $flags[5] == 0)) {      
            $abnormal_pair = 1
            # print join("\t", @F) . "\n";
        }
    }

    if ($abnormal_pair == 1) {
        my $alignmentSize_current = 0;
        while($cigar_current =~ /(\d+)[MD]/g) {
            $alignmentSize_current = $alignmentSize_current + $1;
        }
        print $seqID . "\t" . $chr_current . "\t" . $pos_current . "\t" . ($pos_current + $alignmentSize_current - 1) . "\t" . $dir_current . "\t" . $mapQ . "\n";
    }


}
close(IN);
