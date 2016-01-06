#!/usr/bin/perl -w
use strict;
use Switch;

my $usage = "USAGE:
perl $0 <input> <output>\n";

my $in = $ARGV[0];
my $out =$ARGV[1];
my(@HEAD,$head);

open OUT,">$out" ||die $!;
open IN,$in || die $!;
$head = <IN>;
chomp($head);
print OUT "$head\tdirection\n";
my @head_names = split(/_|\t/,$head);
my ($ctr,$case1,$case2) = ($head_names[1],$head_names[3],$head_names[5]);
while (<IN>){
	chomp;
	print OUT "$_\t";
	my @a=split;
	my($N,$P,$T,$pNP,$pNT,$pPT,$zNP,$zNT,$zPT)=($ctr,$case1,$case2,"=","=","=",0,0,0);
	if($a[5]<0.05){$pNP = ">"}
	if($a[6]<0.05){$pNT = ">"}
	if($a[7]<0.05){$pPT = ">"}
	if($a[8]>0){$zNP=1}
	if($a[9]>0){$zNT=1}
	if($a[10]>0){$zPT=1}
	my($pPN,$pTN,$pTP)=($pNP,$pNT,$pPT);
	my $s = 4*$zNP+2*$zNT+$zPT;
	switch($s){
		case 0 { print OUT "$T $pTP $P $pPN $N\n"}
		case 1 { print OUT "$P $pPT $T $pTN $N\n"}
		case 2 { die $!}
		case 3 { print OUT "$P $pPN $N $pNT $T\n"}
		case 4 { print OUT "$T $pTN $N $pNP $P\n"}
		case 5 { die $!}
		case 6 { print OUT "$N $pNT $T $pTP $P\n"}
		case 7 { print OUT "$N $pNP $P $pPT $T\n"}
	}
}
close IN;
close OUT;
