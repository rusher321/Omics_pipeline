#!/usr/bin/perl -w
use strict;

my($f,$f_ins,$s,$o) = @ARGV;
my $usage = "usage: perl $0 <path_file> <ins_file> <steps(1234)> <output dir>
	path_file should contained such column:
	#sample_name\t#trim_quality\t#trim_length_cut\t#N_cutoff\t#50\%ofQ_control\t#path
	#sample1\t20\t10\t1\t15\t/path/to/fq/file
";
die $usage if @ARGV < 4;
########################
my $s_trim   = "/lfs1/ST_PMO2015G/F13ZOOYJSY1389/metagenomics.20151225/lib/trim_dev.pl";
my $s_filter = "/lfs1/ST_PMO2015G/F13ZOOYJSY1389/metagenomics.20151225/lib/filter_dev.pl";
my $s_tMore  = "/lfs1/ST_PMO2015G/F13ZOOYJSY1389/metagenomics.20151225/lib/trim_more_dev.pl";
my $s_rm     = "/ifs5/PC_MICRO_META/PRJ/MetaSystem/analysis_flow/bin/program/rmhost_v1.0.pl";
my $s_db     = "/nas/RD_09C/resequencing/resequencing/tmp/pub/Genome/Human/human.fa.index";
my $s_soap   = "/lfs1/ST_PMO2015G/F13ZOOYJSY1389/metagenomics.20151225/lib/profiling.flow.SE.dev.pl";

my $dir_s = $o."/script";
my $dir_f = $o."/filter";
my $dir_t = $o."/trim";
my $dir_rm= $o."/rmhost";
my $dir_soap = $o."/soap";

my ($tmp_out,$tmp_out2);

system "mkdir -p $dir_s" unless(-d $dir_s);
system "mkdir -p $dir_f" unless(-d $dir_f or $s !~ /1/);
system "mkdir -p $dir_t" unless(-d $dir_t or $s !~ /2/);
system "mkdir -p $dir_rm" unless(-d $dir_rm or $s !~ /3/);
system "mkdir -p $dir_soap" unless(-d $dir_soap or $s !~ /4/);

open B12,">$o/batch.1-2nd.sh";
open B3,">$o/batch.3rd.sh";
open B4,">$o/batch.4nd.sh";
open C1,">$o/qsub_all.sh";

#print B12 "wdir=$dir_t\n";
#print B3 "wdir=$dir_rm\n";
#print B4 "wdir=$dir_soap\n";

open IN,"<$f" || die $!;
while (<IN>){
	chomp;
	my @a = split;
	my ($sam,$tq,$tl,$nc,$qc,$path,$soapPfx) = @a;
	open S12,">$dir_s/$sam.12.sh";
	if ($s =~ /1/){
		print S12 "perl $s_trim $path $dir_t/$sam $tq $tl\n";
		$tmp_out = "$dir_t/$sam.trim.fq.gz";
	}
	if ($s =~ /2/){
		$tmp_out ||= $path;
		print S12 "perl $s_filter $tmp_out $dir_f/$sam $nc\n";
		print S12 "perl $s_tMore $dir_f/$sam.clean.fq.gz $dir_f/$sam $qc \n";
		$tmp_out = "$dir_f/$sam.clean.fq.gz";
		$tmp_out2 = "$dir_f/$sam.clean_more.fq.gz";
	}
	close S12;
	open S3,">$dir_s/$sam.3.sh";
	open S3M,">$dir_s/$sam.3M.sh";
	if ($s =~ /3/){
		$tmp_out ||= $path;
		print S3 "perl $s_rm -a $tmp_out -d $s_db -m 4 -s 32 -s 30 -r 1 -v 7 -i 0.9 -t 8 -f Y -p  $dir_rm/$sam -q\n";
		print S3M "perl $s_rm -a $tmp_out2 -d $s_db -m 4 -s 32 -s 30 -r 1 -v 7 -i 0.9 -t 8 -f Y -p  $dir_rm/$sam.more -q\n";
		$tmp_out = "$dir_rm/$sam.rmhost.fq.gz";
		$tmp_out2 = "$dir_rm/$sam.more.rmhost.fq.gz";
	}
	close S3;
	close S3M;

	if ($s =~ /4/){
		open S4,">$dir_s/$sam.4.sh";
		open S4M,">$dir_s/$sam.4M.sh";
		$tmp_out ||= $path;
		print S4 "perl $s_soap -i1 $tmp_out -ins $f_ins -o $dir_soap -p $sam > $dir_soap/$sam.log\n";
		print S4M "perl $s_soap -i1 $tmp_out2 -ins $f_ins -o $dir_soap -p $sam.more > $dir_soap/$sam.more.log\n";
		close S4;
		close S4M;
	}
	print B12 "sh $dir_s/$sam.12.sh\n";
#	print B3 "qsub -wd \$wdir -q st.q -P st_ms -l vf=7G $dir_s/$sam.3.sh\n";
#	print B4 "qsub -wd \$wdir -q st.q -P st_ms -l vf=14G $dir_s/$sam.4.sh\n";
	print B3 "sh $dir_s/$sam.3.sh\nsh $dir_s/$sam.3M.sh\n";
	print B4 "sh $dir_s/$sam.4.sh\nsh $dir_s/$sam.4M.sh\n";
}
print C1 "perl /home/fangchao/bin/qsub_all.pl -N bat.2 -d $o/qsub_12 -l vf=0.2G -q st.q -P st_ms -r $o/batch.1-2nd.sh\n";
print C1 "perl /home/fangchao/bin/qsub_all.pl -N bat.3 -d $o/qsub_3 -l vf=7G -q st.q -P st_ms -r $o/batch.3rd.sh\n";
print C1 "perl /home/fangchao/bin/qsub_all.pl -N bat.4 -d $o/qsub_4 -l vf=12G -q st.q -P st_ms -r $o/batch.4nd.sh\n";

close B12;
close B3;
close B4;
