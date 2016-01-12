#!/usr/bin/perl -w
use strict;
unless(@ARGV==2){
	die "a"
	}
my($input,$output)=@ARGV;
open IN,"$input"||die "read $input $!\n";
open OUT,">$output"||die $!;
my $a=<IN>;
#chomp($a);
print OUT "$a";
while(<IN>){
	chomp;
	my (@a,$header,$length);
	my $num=0;
	@a=split /\s+/,$_;
	$header=shift @a;
	$length=@a;
	for(0..$#a){
		if($a[$_]=0)
		{
			$num=$num+1;
			}
        }
		my $percent=1-$num/$length;
	
    print OUT  "$header\t$percent\n"; 
	}

