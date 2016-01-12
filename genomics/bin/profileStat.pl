#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $usage = "USAGE:
	perl $0 -a <taks aim> -f <profile> -p <phenotype> -g <grouping row> -d <decimal> -o <output>

DISCRIPTION:
This script is written to contain kinds of statistics frequently used.
";
my ($aim, $profile, $phenotype, $row, $dec, $output);
GetOptions(
	"a:s" => \$aim,
	"f:s" => \$profile,
	"p:s" => \$phenotype,
	"g:s" => \$row,
	"d:i" => \$dec,
	"o:s" => \$output,
);
die $usage if not defined ($profile|$aim);

if ($aim =~ /geneCount/i){
	$dec ||= 30; $dec = 1 / (10**$dec);

	my %count;
	my $openMethod = ($profile =~ /gz$/)?"gzip -dc $profile |":"$profile";
	open IN,"$openMethod" || die $!;
	chomp(my $head = <IN>);
	my @heads = split(/\t/,$head);
	while (<IN>){
		chomp;
		my @a = split;
		for (my $r=1;$r<@heads;$r++){
			$count{$r} ||= 0;
			$count{$r} ++ if $a[$r] > $dec;
		}
	}
	close IN;
	open OUT,"> $output" || die $!;
	for (my $r=1;$r<@heads;$r++){
		print OUT "$heads[$r]\t$count{$r}\n";
	}
	close OUT;
}
