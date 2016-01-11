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
	-c|config	:set parameters for each setp, default Qt=20,l=10,N=1,Qf=15
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
$config  ||= "Qt=20,l=10,N=1,Qf=15";
foreach my $par (split(/,/,$config)){
	my @a = split(/=/,$par);
	$CFG{$a[0]} = $a[1];
}

# scripts under bin
my $bin = "$Bin/bin";
#my $s_trim   = "$bin/trimReads.pl";
#my $s_filter = "$bin/filterReads.pl";
my $s_clean  = "$bin/readsCleaning.pl";
my $s_rm     = "$bin/rmhost_v1.0.pl";
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
my ($tmp_out,$tmp_outN,$tmp_outQ);
while (<IN>){
	chomp;
	my @a = split;
	my ($sam,$pfx,$path) = @a;
	open S1,">$dir_sI/$pfx.01.clean.sh";
	if ($step =~ /1/){
		print S1 "perl $s_clean $path $dir_c/$pfx $CFG{'Qt'} $CFG{'l'} $CFG{'N'} $CFG{'Qf'}\n";
		$tmp_out  = "$dir_c/$pfx.trimed.fq.gz";
		$tmp_outN = "$dir_c/$pfx.cleanN.fq.gz";
		$tmp_outQ = "$dir_c/$pfx.cleanQ.fq.gz";
	}
#	if ($s =~ /2/){
#		$tmp_out ||= $path;
#		print S2 "perl $s_filter $tmp_out $dir_f/$pfx $nc\n";
#		print S2 "perl $s_tMore $dir_f/$pfx.clean.fq.gz $dir_f/$pfx $qc \n";
#		$tmp_out = "$dir_f/$pfx.cleanN.fq.gz";
#		$tmp_out2 = "$dir_f/$pfx.cleanQ.fq.gz";
#	}
	close S1;
	open S2,">$dir_sI/$pfx.02.rmhost.sh";
	if ($step =~ /2/){
		$tmp_outQ ||= $path;
		print S2 "perl $s_rm -a $tmp_outQ -d $s_db -m 4 -s 32 -s 30 -r 1 -v 7 -i 0.9 -t 8 -f Y -p  $dir_r/$pfx -q\n";
		$tmp_out = "$dir_r/$pfx.rmhost.fq.gz";
	}
	close S2;

	if ($step =~ /3/){
		open S3,">$dir_sp/$pfx.03.soap.sh";
		$tmp_out ||= $path;
		print S3 "perl $s_soap -i1 $tmp_out -ins $ins_f -o $dir_s -p $pfx > $dir_s/$pfx.log\n";
		close S3;
	}
	print B1 "sh $dir_s/$pfx.01.clean.sh\n";
	print B2 "sh $dir_s/$pfx.02.rmhost.sh\n";
	print B3 "sh $dir_s/$pfx.03.soap.sh\n";
}
print C1 "perl /home/fangchao/bin/qsub_all.pl -N bat.2 -d $dir_s/qsub_1 -l vf=0.3G -q st.q -P st_ms -r $dir_s/batch.clean.sh\n";
print C1 "perl /home/fangchao/bin/qsub_all.pl -N bat.3 -d $dir_s/qsub_2 -l vf=8G -q st.q -P st_ms -r $dir_s/batch.rmhost.sh\n";
print C1 "perl /home/fangchao/bin/qsub_all.pl -N bat.4 -d $dir_s/qsub_3 -l vf=15G -q st.q -P st_ms -r $dir_s/batch.soap.sh\n";

close B1;
close B2;
close B3;


# ####################
# SUB FUNCTION
# ####################
sub version {
	print <<VERSION;
	version:	v0.10
	update:		20160111
	author:		fangchao\@genomics.cn

VERSION
};

