#!/usr/bin/perl -w
use warnings;
use strict;

die &usage if @ARGV !=6;

my ($in,$pfx,$qt,$l,$n,$qf) = @ARGV;

open IN,"gzip -dc $in |",or die "error\n";
open OUTT,"|gzip >$pfx.trimed.fq.gz",or die "error\n";
open OUTN,"|gzip >$pfx.cleanN.fq.gz",or die "error\n";
open OUTQ,"|gzip >$pfx.cleanQ.fq.gz",or die "error\n";
open STAT,"> $pfx.clean.stat_out",or die "error\n";

my($total,$remainN,$remainQ,$sum_bp,$max_bp,$min_bp)= (0,0,0,0,0);
while(<IN>){
	chomp;
	chomp(my $seq=<IN>);
	chomp(my $num=<IN>);
	chomp(my $quality=<IN>);
	$total++;
	my $count=0;
	my $originLength = length($seq);
	# trim
	my $Tscore = &Qstat($quality,$qt,"trim",$l);
	my $len = $originLength-$count;
	$seq=substr($seq,0,$len);
	$quality=substr($quality,0,$len);
	print OUTT "$_\n$seq\n$num\n$quality\n";
	$max_bp = ($max_bp > $len)?$max_bp:$len;
	$min_bp = (-$min_bp > -$len)?$min_bp:$len;
	$sum_bp += $len;
	# filter
	$count=$seq=~tr/N/N/;
	if($count<=$n){
		print OUTN "$_\n$seq\n$num\n$quality\n";
		$remainN ++; 
		my $Qscore = &Qstat($quality,$qf,"filter");
		if(($Qscore*2)<=length($quality)){
			print OUTQ "$_\n$seq\n$num\n$quality\n";
			$remainQ ++; 
		}
	}
}
close IN;
close OUTT;
close OUTN;
close OUTQ;
my $avg  = $sum_bp / $total;
my $rate = $remainQ / $total;
print STAT "Total\tmax\tmin\tavg\tremainAfterNfilter(Q=$qt/N=$n)\tremainAfterQfilter(Q=$qf)\trate\tSampleTAG\n";
print STAT "$total\t$max_bp\t$min_bp\t$avg\t$remainN\t$remainQ\t$rate\t$pfx\n";

close STAT;

# sub
sub Qstat {
	my ($q_l,$q_n,$method,$c_n) = @_;
	my $c = 0;
	if ($method eq "filter"){
		for(my $i=0;$i<length($q_l);$i++){
			my $q=substr($q_l,$i,1);
			$q=ord$q;
			$q=$q-33;
			if($q<=$q_n){$c++};
		}
	}elsif($method eq "trim"){
		for(my $i=length($q_l)-1;$i>=0;$i--){
			my $q=substr($q_l,$i,1);
			$q=ord$q;
			$q=$q-33;
			last if($q>=$q_n || $c>=$c_n);
			$c++;
		}
	}
	return($c);
}




