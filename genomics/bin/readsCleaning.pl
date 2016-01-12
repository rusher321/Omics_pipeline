#!/usr/bin/perl -w
use warnings;
use strict;
use File::Basename; 

die &usage if @ARGV !=8;

my ($fq1,$fq2,$pfx,$qt,$l,$n,$qf,$lf) = @ARGV;
if($fq1 eq $fq2){
	open FQ1,"gzip -dc $fq1 |",or die "error\n";
	open OUT1,"|gzip >$pfx.clean.fq.gz",or die "error\n";
}else{
	open FQ1,"gzip -dc $fq1 |",or die "error\n";
	open FQ2,"gzip -dc $fq2 |",or die "error\n";
	open OUT1,"|gzip >$pfx.clean.fq1.gz",or die "error\n";
	open OUT2,"|gzip >$pfx.clean.fq2.gz",or die "error\n";
}
open STAT,"> $pfx.clean.stat_out",or die "error\n";

my($total,$remainN,$remainQ,$sum_bp1,$max_bp1,$min_bp1,$sum_bp2,$max_bp2,$min_bp2)= (0,0,0,0,0,10000,0,0,10000);
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
	# filter
	$count=$seq=~tr/N/N/;
	next if($count > $n || $len < $lf);					# N number & length limit judgement
	$remainN ++; 
	# filter more
	my $Qscore = &Qstat($quality,$qf,"filter");
	next if (($Qscore*2) > length($quality));			# half Q score judgement
	# quality pass
	$READ{$a[0]}{'fqs'} = "pass";
	$READ{$a[0]}{'seq'} = "$a[0] $a[1] length=$len\n$seq\n$num\n$quality\n";
	$READ{$a[0]}{'len'} = $len;
	next unless $fq1 eq $fq2;
	# stat
	print OUT1 $READ{$a[0]}{'seq'};
	$remainQ ++;
	$max_bp1 = ($max_bp1 > $len)?$max_bp1:$len;
	$min_bp1 = ($min_bp1 < $len)?$min_bp1:$len;
	$sum_bp1 += $len;
}
close FQ1;
unless($fq1 eq $fq2){
	while(<FQ2>){
		chomp;
		my @a =split;
		next unless defined $READ{$a[0]}{'fqs'};# ne "pass";
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
		my $count=$seq=~tr/N/N/;
		next if($count > $n || $len < $lf);
		my $Qscore = &Qstat($quality,$qf,"filter");
		next if (($Qscore*2) > length($quality));
		$remainQ ++;
		print OUT1 $READ{$a[0]}{'seq'};
		print OUT2 "$a[0] $a[1] length=$len\n$seq\n$num\n$quality\n";

		($max_bp2,$min_bp2,$sum_bp2) = &lengthAvg($len,$max_bp2,$min_bp2,$sum_bp2);
		($max_bp1,$min_bp1,$sum_bp1) = &lengthAvg($READ{$a[0]}{'len'},$max_bp1,$min_bp1,$sum_bp1);
	}
}
close FQ2;
close OUT1;
close OUT2;

my $avg1 = $sum_bp1 / $total;
my $avg2 = $sum_bp2 / $total;
my $rate = $remainQ / $total;
my $tag = basename($pfx);
unless ($fq1 eq $fq2){
	print STAT "Total\tmax1\tmin1\tavg1\tmax2\tmin2\tavg2\tremain\trate\tSampleTAG(trim_limit=$l,Qt=$qt,N=$n,Qf=$qf,min=$lf)\n";
	print STAT "$total\t$max_bp1\t$min_bp1\t$avg1\t$max_bp2\t$min_bp2\t$avg2\t$remainQ\t$rate\t$tag\n";
}else{
	print STAT "Total\tmax\tmin\tavg\tremain\trate\tSampleTAG(trim_limit=$l,Qt=$qt,N=$n,Qf=$qf,min=$lf)\n";
	print STAT "$total\t$max_bp1\t$min_bp1\t$avg1\t$remainQ\t$rate\t$tag\n";
}

close STAT;

# sub
sub lengthAvg{
	my($l,$max,$min,$sum) = @_;
	$max = ($max > $l)?$max:$l;
	$min = ($min < $l)?$min:$l;
	$sum += $l;
	return($max,$min,$sum);
}

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



