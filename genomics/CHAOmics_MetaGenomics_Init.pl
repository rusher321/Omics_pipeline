#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin);
use Cwd 'abs_path';

my $cwd = abs_path;
#my($f,$f_ins,$s,$out_dir) = @ARGV;
#my $usage = "usage: perl $0 <path_file> <ins_file> <steps(1234)> <output dir>
#	path_file should contained such column:
	#sample_name\t#trim_quality\t#trim_length_cut\t#N_cutoff\t#50\%ofQ_control\t#path
	#sample1\t20\t10\t1\t15\t/path/to/fq/file
#";
sub usage {
	print <<USAGE;
usage:
	perl $0 [options]
options:
	-p|path		:[essential]sample path file
	-i|ins		:[essential]insert info file
	-s|step		:functions,default 123
					1	trim+filter
					2	remove host genomic reads
					3	soap mapping to microbiotic genomics
	-o|outdir	:output directory path. Conatins the results and scripts.
	-c|config	:set parameters for each setp, default Qt=20,l=10,N=1,Qf=15,lf=0
					Qt	Qvalue for trim
					l	bp length for trim
					N	tolerance number of N for filter
					Qf	Qvalue for filter. The reads which more than half of the bytes lower than Qf will be discarded.
	-h|help		:show help info
	-v|version	:show version and author info.
USAGE
};
my($path_f,$ins_f,$step,$out_dir,$config,%CFG,$help,$version);
GetOptions(
	"p|path:s"    => \$path_f,
	"i|ins:s"     => \$ins_f,
	"s|step:i"    => \$step,
	"o|outdir:s"  => \$out_dir,
	"c|config:s"  => \$config,
	"h|help:s"    => \$help,
	"v|version:s" => \$version,
);
die &usage if ( (!defined $path_f)||(!defined $ins_f)||(defined $help));
die &version if defined $version;

# ####################
# initialize variables
# ####################
$step    ||= "123";
$out_dir ||= $cwd; $out_dir = abs_path($out_dir);
$path_f = abs_path($path_f);
$ins_f  = abs_path($ins_f);
$config  ||= "Qt=20,l=10,N=1,Qf=15,lf=0";
foreach my $par (split(/,/,$config)){
	my @a = split(/=/,$par);
	$CFG{$a[0]} = $a[1];
}

# scripts under bin
my $bin = "$Bin/bin";
#my $s_trim   = "$bin/trimReads.pl";
#my $s_filter = "$bin/filterReads.pl";
my $s_clean  = "$bin/readsCleaning.pl";
my $s_rm     = "/ifs5/PC_MICRO_META/PRJ/MetaSystem/analysis_flow/bin/program/rmhost_v1.0.pl";
my $s_soap   = "$bin/soap2BuildAbudance.dev.pl";
# public database prefix
my $s_db     = "/nas/RD_09C/resequencing/resequencing/tmp/pub/Genome/Human/human.fa.index";
# project results directiory structure
my $dir_s = $out_dir."/script";
	my $dir_sI = $dir_s."/individual";
	my $dir_sL = $dir_s."/linear";
	my $dir_sB = $dir_s."/steps";
#my $dir_t = $out_dir."/trim";
#my $dir_f = $out_dir."/filter";
my $dir_c = $out_dir."/clean";
my $dir_r = $out_dir."/rmhost";
my $dir_sp = $out_dir."/soap";

system "mkdir -p $dir_s" unless(-d $dir_s);
	system "mkdir -p $dir_sI" unless(-d $dir_sI);
	system "mkdir -p $dir_sL" unless(-d $dir_sL);
	system "mkdir -p $dir_sB" unless(-d $dir_sB);
#system "mkdir -p $dir_f" unless(-d $dir_f or $s !~ /1/);
#system "mkdir -p $dir_t" unless(-d $dir_t or $s !~ /2/);
system "mkdir -p $dir_c" unless(-d $dir_c or $step !~ /1/);
system "mkdir -p $dir_r" unless(-d $dir_r or $step !~ /2/);
system "mkdir -p $dir_sp" unless(-d $dir_sp or $step !~ /3/);

#open B1,">$dir_s/batch.trim.sh";
#open B2,">$dir_s/batch.filter.sh";
open B1,">$dir_s/batch.clean.sh";
open B2,">$dir_s/batch.rmhost.sh";
open B3,">$dir_s/batch.soap.sh";
open C1,">$out_dir/qsub_all.sh";

open IN,"<$path_f" || die $!;
my (%SAM,@samples,$tmp_out,$tmp_outN,$tmp_outQ);
while (<IN>){
	chomp;
	my @a = split;
	my ($sam,$pfx,$path) = @a;
	$SAM{$sam}{$pfx} = $path;
}

foreach my $sam (keys %SAM){
	my @fqs = sort keys %{$SAM{$sam}};
	die "$sam got more than 2 fq files, pls check it out!" if @fqs > 2;
	my ($fq1,$fq2) = ($SAM{$sam}{$fqs[0]}, $SAM{$sam}{$fqs[1]});
###############################
	if ($step =~ /1/){
		open SH,">$dir_sI/$sam.clean.sh";
		if (@fqs eq 2){
			print SH "perl $s_clean $fq1 $fq2 $dir_c/$sam $CFG{'Qt'} $CFG{'l'} $CFG{'N'} $CFG{'Qf'} $CFG{'lf'}\n";
			($SAM{$sam}{$fqs[0]}, $SAM{$sam}{$fqs[1]}) = ("$dir_c/$sam.clean.fq1.gz","$dir_c/$sam.clean.fq2.gz");
		}else{
			print SH "perl $s_clean $fq1 $fq1 $dir_c/$sam $CFG{'Qt'} $CFG{'l'} $CFG{'N'} $CFG{'Qf'} $CFG{'lf'}\n";
			$tmp_out = "$dir_c/$sam.clean.fq.gz";
		}
		close SH;
		print B1 "sh $dir_sI/$sam.clean.sh\n";
	}
###############################
	if ($step =~ /2/){
		open SH,">$dir_sI/$sam.rmhost.sh";
		if (@fqs eq 2){
			print SH "perl $s_rm -a $SAM{$sam}{$fqs[0]} -b $SAM{$sam}{$fqs[1]} -d $s_db -m 4 -s 32 -s 30 -r 1 -v 7 -i 0.9 -t 8 -f Y -p  $dir_r/$sam -q\n";
			($SAM{$sam}{$fqs[0]}, $SAM{$sam}{$fqs[1]}) = ("$dir_r/$sam.rmhost.1.fq.gz","$dir_r/$sam.rmhost.2.fq.gz");
		}else{
			print SH "perl $s_rm -a $tmp_out -d $s_db -m 4 -s 32 -s 30 -r 1 -v 7 -i 0.9 -t 8 -f Y -p  $dir_r/$sam -q\n";
			$tmp_out = "$dir_r/$sam.rmhost.fq.gz";
		}
		close SH;
		print B2 "sh $dir_sI/$sam.rmhost.sh\n";
	}
###############################
	if ($step =~ /3/){
		open SH,">$dir_sI/$sam.soap.sh";
		if (@fqs eq 2){
			print SH "perl $s_soap -i1 $fq1 -i2 $fq2 -ins $ins_f -o $dir_s -p $sam > $dir_sp/$sam.log\n";
		}else{
			print SH "perl $s_soap -i1 $tmp_out -ins $ins_f -o $dir_s -p $sam > $dir_sp/$sam.log\n";
		}
		close SH;
		print B3 "sh $dir_sI/$sam.soap.sh\n";
	}
}
$CFG{'q'}  ||= "st.q";
$CFG{'P'}  ||= "st_ms";
$CFG{'pro'}  ||= 8;
$CFG{'vf1'} ||= "0.3G";
$CFG{'vf2'} ||= "8G";
$CFG{'vf3'} ||= "15G";
$CFG{'m'} ||= 30;
print C1 "perl /home/fangchao/bin/qsub_all.pl -N B.c -d $dir_s/qsub_1 -l vf=$CFG{'vf1'} -q $CFG{'q'} -P $CFG{'P'} -r -m $CFG{'m'} $dir_s/batch.clean.sh\n" if $step =~ /1/;
print C1 "perl /home/fangchao/bin/qsub_all.pl -N B.r -d $dir_s/qsub_2 -l vf=$CFG{'vf2'} -q $CFG{'q'} -P $CFG{'P'} -r -m $CFG{'m'} $dir_s/batch.rmhost.sh\n" if $step =~ /2/;
print C1 "perl /home/fangchao/bin/qsub_all.pl -N B.s -d $dir_s/qsub_3 -l vf=$CFG{'vf3'},p=$CFG{'pro'} -q $CFG{'q'} -P $CFG{'P'} -r -m $CFG{'m'} $dir_s/batch.soap.sh\n" if $step =~ /3/;

close B1;
close B2;
close B3;


# ####################
# SUB FUNCTION
# ####################
sub version {
	print <<VERSION;
	version:	v0.12
	update:		20160111
	author:		fangchao\@genomics.cn

VERSION
};

