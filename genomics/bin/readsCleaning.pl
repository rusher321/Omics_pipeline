#!/usr/bin/perl -w
use warnings;
use strict;

die &usage if @ARGV !=7;

my ($fq1,$fq2,$pfx,$qt,$l,$n,$qf,$lf) = @ARGV;

open FQ1,"gzip -dc $fq1 |",or die "error\n";
#open OUTT,"|gzip >$pfx.trimed.fq.gz",or die "error\n";
#open OUTN,"|gzip >$pfx.cleanN.fq.gz",or die "error\n";
open OUT,"|gzip >$pfx.clean.fq.gz",or die "error\n";
open STAT,"> $pfx.clean.stat_out",or die "error\n";

my($total,$remainN,$remainQ,$sum_bp,$max_bp,$min_bp)= (0,0,0,0,0,10000);
my %READ;
while(<FQ1>){
	chomp;
	my @a =split;
	chomp(my $seq=<FQ1>);
	chomp(my $num=<FQ1>);
	chomp(my $quality=<FQ1>);
	$total++;
	my $count=0;
	my $originLength = length($seq);
	# trim
	my $Tscore = &Qstat($quality,$qt,"trim",$l);
	my $len = $originLength-$Tscore;
	$seq=substr($seq,0,$len);
	$quality=substr($quality,0,$len);
#	print OUTT "$_\n$seq\n$num\n$quality\n";
	# filter
	$count=$seq=~tr/N/N/;
	next if($count > $n || $len < $lf);
#	print OUTN "$_\n$seq\n$num\n$quality\n";
	$remainN ++; 
	my $Qscore = &Qstat($quality,$qf,"filter");
	next if (($Qscore*2) > length($quality));
	$READ{$a[0]}{'fqs'} = ++ ;
	$READ{$a[0]}{'seq'} = "$a[0]\s$a[1]\slength=$len\n$seq\n$num\n$quality\n";
	$READ{$a[0]}{'len'} = $len;
#	print OUT "$a[0]\s$a[1]\slength=$len\n$seq\n$num\n$quality\n";
	next if $fq1 ne $fq2;
	$remainQ ++;
	$max_bp = ($max_bp > $len)?$max_bp:$len;
	$min_bp = ($min_bp < $len)?$min_bp:$len;
	$sum_bp += $len;
}
close FQ1;
if($fq1 ne $fq2){
	open FQ2,"gzip -dc $fq2 |",or die "error\n";
	while(<FQ2>){
		chomp;
		my @a =split;
		next if $READ{$a[0]}{'fqs'} == 0;
		chomp(my $seq=<FQ2>);
		chomp(my $num=<FQ2>);
		chomp(my $quality=<FQ2>);
		my $originLength = length($seq);
		# trim
		my $Tscore = &Qstat($quality,$qt,"trim",$l);
		my $len = $originLength-$Tscore;
		$seq=substr($seq,0,$len);
		$quality=substr($quality,0,$len);
		# filter
		$count=$seq=~tr/N/N/;
		next if($count > $n || $len < $lf);
		my $Qscore = &Qstat($quality,$qf,"filter");
		next if (($Qscore*2) > length($quality));
		$remainQ ++;
		($max_bp2,$min_bp2,$sum_bp2) = &lengthAvg($len,$max_bp2,$min_bp2,$sum_bp2);
		($max_bp1,$min_bp1,$sum_bp1) = &lengthAvg($len,$max_bp1,$min_bp1,$sum_bp1);
	}
}
close OUT;

my $avg  = $sum_bp / $total;
my $rate = $remainQ / $total;
print STAT "Total\tmax\tmin\tavg\tremain(trim_limit=$l,Qt=$qt,N=$n,Qf=$qf,min=$lf)\trate\tSampleTAG\n";
print STAT "$total\t$max_bp\t$min_bp\t$avg\t$remainQ\t$rate\t$pfx\n";

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



