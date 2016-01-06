#!/usr/bin/perl -w

my $usage =" perl $0 <input 1> <input2> > <output>\n";

print $!."\n$usage" if @ARGV != 2;

my($in1,$in2)=@ARGV;

my(%MAP);

open IN1,"<$in1" || $!."$usage";
while (<IN1>) {
	chomp;
	my @a=split;
	$MAP{$a[1]}{'in1'}=$a[0];
}
close IN1;

open IN2,"<$in2" || $!."$usage";
while (<IN2>) {
	chomp;
	my @a=split;
	$a[0] =~ /(SD-B\d+)$/ || next;
	my $tid=$1;
	my $tmap= "$a[2]\_$a[3]";
	$MAP{$tid}{'in2'} = $tmap;
	$MAP{$tid}{'extra'}=$a[3];
}
close IN2;

foreach my $i(sort keys %MAP){
	next if not defined $MAP{$i}{'in1'};
	next if not defined $MAP{$i}{'in2'};
	print "$MAP{$i}{'in1'}\t$MAP{$i}{'in2'}\t$i\t$MAP{$i}{'extra'}\n";
}


