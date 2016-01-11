#!/usr/bin/perl -w
use strict;
#push @INC, "/usr/lib64/perl5/vendor_perl/5.8.8/x86_64-linux-thread-multi";
#use PerlIO::gzip;
#$|=1;
################################################################################
unless(4==@ARGV) {
    &usage;
    exit;
}
################################################################################
my($list_f,$order_f,$row,$out) = @ARGV;
my(@order,%list,@info,$i,%class,%cover);
################################################################################
#open IN, "<:gzip(autopop)", "$ARGV[0]" or die "$! $ARGV[0]!\n";
open IN,"<$list_f" || die "read $list_f $!\n";
while(<IN>) {
    chomp;
    @info=split /\//;
    $info[-1]=~/^(\S+)\.abundance(|\.gz)$/;
    my $name =$1;
#	print "\$info[-1] = $info[-1]\t\$name = $name\n";
#    $name =~ s/s3w/w/;#print $name."\n"; #custom modified by fangchao@genomics.cn #2015-04-07
    $list{$name}=$_;
}
close IN;
################################################################################
open IN,$order_f or die "read $order_f $!\n";
while(<IN>) {
    chomp;
    push(@order,$_);
}
close IN;
################################################################################
for($i=0;$i<@order;++$i) {
#    open IN, "<:gzip(autopop)", "$list{$order[$i]}" or die "$! $list{$order[$i]}\n";
	open IN, "gzip -dc $list{$order[$i]} |"  or die "$!\n";
   # open IN,$list{$order[$i]} or die "read $list{$order[$i]} $!\n";
    <IN>;
    while(<IN>) {
        chomp;
        @info=split /\t/;
        $class{$info[0]}.="\t".$info[$row];
		$cover{$info[0]}||=0;
		if ($info[$row] > 0){ $cover{$info[0]} =1}; 
    }
    close IN;
}
################################################################################
open OT,">$out" or die "write $out $!\n";
open CT,">$out.cover" || die $!;
for($i=0;$i<@order;++$i) {
    print OT "\t",$order[$i];
}
print OT "\n";
foreach $i(sort {$a<=>$b} keys %class) {
    print OT $i,$class{$i},"\n";
	print CT "$i\t$cover{$i}\n";
}
close OT;
close CT;
################################################################################
sub usage {
    print STDERR<<USAGE;

    Description\n
    This programme is to combine profiling table\n
    Usage:  perl $0 [file.list] [order] [row] [outfile]\n
    Row must be in 1(pairs),2(base_abundance),3(reads_abundance),4(depth_abundance)\n
    Author Libranjie,zhouyuanjie\@genomics.org.cn\n
    updated 2014/12/5 linyuxiang\@genomics.cn\n
USAGE
}
