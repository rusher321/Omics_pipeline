#!/usr/bin/perl -w
use strict;

my $usage = "perl $0 <input> <prefix> <decimal> <output>\n";
die $usage if @ARGV < 4;
my ($in,$pfx,$dec,$out) = @ARGV;
$dec = 1 / (10**$dec);
my $count = 0 ;
open IN,"gzip -dc $in |" || die $!;
my $head = <IN>;
while (<IN>){
	chomp;
	my @a = split;
	$count ++ if $a[3] > $dec;
}
close IN;
open OUT,"> $out" || die $!;
print OUT "$pfx\t$count\n";
close OUT;
