#!/usr/bin/perl -w 
use strict;

my $usage = "Usage: perl $0 <list file> <input file> [mode] <out filename> [range](optional) [for sort(optional)]
		mode	in	<include the elements in the list>
			ex	<exclude the elements in the list>
		
		range	h	only match header
			r	only match row
			t	<default> match the whole table(both header and column)
		forSort	s	if your files are sorted,this option may help match much more faster\n";

if (@ARGV<4) {
	die "insufficient ARGV\n$usage";
}

my ($lst_f,$table_f,$mode,$dir,$range,$st)=@ARGV;
$range||="t";
$st||="F";

my (@lst,@column,$switch,$oppo);

if ($mode eq "in"){
	$oppo = "ex"
}elsif($mode eq "ex"){
	$oppo = "in"
}else{ die "ERROR:<mode>:\"$mode\" is not allowed. Only \"ex\" or \"in\" available.\n$usage"
}
$dir =~ s/\/$//;

open L,"$lst_f" || die $!;
while(<L>){
	chomp;
	my @a=split;
	push @lst, $a[0];
}
close L;
#####################################
$table_f=~ /(\/|)([^\/]+)\.(\w+)$/; 
my $pt = $2;
my $suffix = $3;
$lst_f =~ /(\/|)([^\/]+)\.(\w+)$/;
open OUT,">$dir" || die $!;
#####################################
open T,"$table_f" || die $!;

my $head=<T>;
chomp ($head);
my @col=split /\t/,$head;
for (my $i=1;$i<=$#col;$i++){
	$switch = $oppo;
	for (my $j=0;$j<=$#lst;$j++){
		if ($col[$i] eq $lst[$j]){$switch = $mode}
	}
	next if $switch eq "ex" && $range ne "r" ;
	push @column,$i;
	print OUT "\t$col[$i]";
}
print OUT "\n";

my $k=0;
while(<T>){
	chomp;
	my @a=split /\t/;
	$switch = $oppo;
	$k=0 if $st ne "s";
	for (my $i=$k;$i<=$#lst;$i++){
		if ($a[0] eq $lst[$i]){$switch = $mode;$k=$i+1;last}
	}
	next if $switch eq "ex" && $range ne "h";
	print OUT $a[0];
	for (my $i=0;$i<=$#column;$i++){
		print OUT "\t$a[$column[$i]]";
	}
	print OUT "\n";
}

close T;

