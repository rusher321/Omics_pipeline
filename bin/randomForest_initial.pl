#!/usr/bin/perl
use strict;
use FindBin qw($Bin $Script);
use lib $Bin;
use Getopt::Long;
use Cwd qw(abs_path);

my $usage = "usage: perl $0 [options]
			-s      times to repeat random selecting training set
			-phe    phenotype sheet
			-todo	phenotype to classify the samples
			-pro    profile
			-todo_id	if you have specific ids to do, set the id list here,or set as \"F\"
			-todo_marker 
			-ctr	control
			-c1     case1
			-c2     case2
			-f      fold times
			-ct     cv times
			-cs     cv Steps
			-pf     prefix
";

die $usage if @ARGV < 24;

my ($setTimes, $pheF, $todo, $proF, $do_wilcox, $todo_id, $todo_marker, $control, $case1, $case2, $fold, $cvTimes, $cvSteps, $prefix);
GetOptions(
	"s:s" => \$setTimes,
	"phe:s" => \$pheF,
	"todo:s"=> \$todo,
	"pro:s" => \$proF,
	"do_w:s" => \$do_wilcox,
	"todo_id:s"=> \$todo_id,
	"todo_marker:s"=> \$todo_marker,
	"ctr:s" => \$control,
	"c1:s" => \$case1,
	"c2:s" => \$case2,
	"f:s"   => \$fold,
	"ct:s"  => \$cvTimes,
	"cs:s"  => \$cvSteps,
	"pf:s"  => \$prefix,
);

my $Work_dir = $prefix.".batch";
my $pfx = $prefix;
$pfx =~ s/^(\S+)\/(\S+?)$/$2/;

$pheF = abs_path($pheF);
$proF = abs_path($proF);
$Work_dir = abs_path($Work_dir);

mkdir($Work_dir);

open SH,">$Work_dir/step1_setSamples.sh" || die $!;
for(my $i=1;$i<=$setTimes;$i++){
	print SH "## repeat ".$i;
	print SH "\nmkdir -p $Work_dir/Repeat_$i\n";
#	print SH "cd $Work_dir/Repeat_$i\n";
	my $repeat_pfx = $prefix;
#	$repeat_pfx =~ s/(\S+)\/(\S+?)$/$1\/$2.batch\/Repeat_$i\/$2/;
	$repeat_pfx = "$Work_dir/Repeat_$i/Repeat_$i\_$pfx";
	if (-e $todo_id){
		my $signF = abs_path($todo_id);
		print SH "Rscript $Bin/avesample.R $pheF $todo $signF $proF $case1 $case2 $control $repeat_pfx \n";
	}else{
		print SH "Rscript $Bin/avesample.R $pheF $todo $proF $proF $case1 $case2 $control $repeat_pfx 1>$repeat_pfx.dunn.log \n";
	}
	if (-e $do_wilcox){

	}
}
print SH "wait\n";
close SH;

open SH,">$Work_dir/step2_randomForest.sh" || die $!;
open WC,">$Work_dir/do_wilcox.sh"|| die $!;
for(my $i=1;$i<=$setTimes;$i++){
#	print SH "export /ifs1/ST_META/USER/xiehailiang/usr/R/ && ";
	my $wd = "$Work_dir/Repeat_$i";
	my $repeat_pfx = "$wd/Repeat_$i\_$pfx";
	my @vs;
	if ($case2 eq "F"){
		@vs = ("$case1-vs-$control");
	}else{
		@vs = ("$case1-vs-$control","$case2-vs-$control","$case2-vs-$case1");
	}
	while($#vs >-1){
		my $v = shift @vs;
		if ($do_wilcox = "T"){
			my $wilcox_marker = "$repeat_pfx.$v.wilcox.xls";
			print WC "/home/fangchao/bin/Rscript $Bin/wilcox_rf.r $repeat_pfx.$v.train_x.xls $repeat_pfx.$v.train_y.xls $wilcox_marker $v &\n";
			$todo_marker = $wilcox_marker;
		}
		print SH "/home/fangchao/bin/Rscript $Bin/randomForest_CRC.R $repeat_pfx.$v.train_x.xls $repeat_pfx.$v.train_y.xls $repeat_pfx.$v.test_x.xls $repeat_pfx.$v.test_y.xls $todo_marker 10 0.9 5 0 $repeat_pfx.$v 1>$repeat_pfx.$v.stat1 2>$repeat_pfx.$v.stat2 \n";
		#print SH "/home/fangchao/bin/Rscript $Bin/randomForest_FC.R $wd/$prefix.$v.train_x.xls $wd/$prefix.$v.train_y.xls $wd/$prefix.$v.test_x.xls $wd/$prefix.$v.test_y.xls 10 0.9 5 0 $wd/$prefix 1>$wd/$prefix.stat1 2>$wd/$prefix.stat2 \n";
	}
}
print WC "wait\n";
close WC;
close SH;

open SH,">$Work_dir/qsub_Step2.sh"|| die $!;
print SH "perl /home/fangchao/bin/qsub_all.pl -d $Work_dir/shell.sh_qsub -q st.q -l vf=0.5G -P st_ms $Work_dir/step2_randomForest.sh &";
close SH;
