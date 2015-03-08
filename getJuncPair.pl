#! /usr/local/bin/perl

use strict;

my $input = $ARGV[0];

my $num = 1;
open(IN, $input) || die "cannot open $!";
while(<IN>) {
    s/[\r\n]//g;
    my @F = split("\t", $_);

    my $ID = $F[6];
    my $chr = "";
    my $pos = "";

    if ($F[12] =~ /(\w+)\:(\d+)/) {
        
        $chr = $1;
        $pos = $2;
    } else {
        print STDERR join("\t", @F) . "\n";
    }

    $ID =~ s/\/1$/\/3/;
    $ID =~ s/\/2$/\/1/;
    $ID =~ s/\/3$/\/2/;

    print $chr . "\t" . ($pos - 1) . "\t" . $pos . "\t" . $ID . "\t" . $num . "\n";

    $num = $num + 1;
}
close(IN);


