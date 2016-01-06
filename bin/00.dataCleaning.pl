#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin);
use lib $Bin;
use Cwd 'abs_path';

sub usage {
	print <<USAGE
Version:	DEMO
Update :	20151120 	(since 20151120)
----------------------------------------
Description:
	This pipeline is to ease my pressure for Multiple omics analysis. For this version, I'm focus on the process of data polish.
Usage:
	perl $0 -cfg <configure> -dir <work dir> -step [options]
Options
	-help	ask for help
	-cfg*	[essential] the configure file contains various fussy settings. Be patient.
	-dir*	[essential] The shell and output date will be saved under this directory (default CHAOS)
	-step	each function will be ran step by step. You can continue by any step if needed.
			1	data polish

e.g.
	perl $0 -cfg configureA -dir ./output -step 1

USAGE
}

# Globle variables
my ($help, $cfg, $dir, $step);
GetOptions(
	"help"   => \$help,
	"cfg=s"  => \$cfg,
	"dir=s"  => \$dir,
	"step=s" => \$step,
);
&usage && exit 0 if defined $help;

# Default settings
$dir  ||= "CHAOS";
$step ||= "1";
$dir = abs_path($dir);
# load tools & software configure
my %library;

&parse_lib("$Bin/../lib/library.cfg",\%library);
my $idmap_pl = $library{'idmap_pl'};
my $replace_pl = $library{'replace_pl'};
my $trim_pl = $library{'trim_pl'};
my $sort_pl = $library{'sort_pl'};

# load config associated with the input data
my %cfg;
my ($s_name,$ion_name);
open IN,$cfg || die "fail open: $cfg\n";
while(<IN>){
	chomp;
	next if (/^#/);
	if (/(\S+)\s*=\s*(\S+)/){
		my ($cfg_name,$cfg_value) = ($1,$2);
		if ($cfg_name eq "sample_pfx"){
			$s_name = $cfg_value;
			$ion_name = "";
#			$cfg{$s_name}{'sample_pfx'} = $cfg_value;
		}elsif($cfg_name eq "ion"){
			$ion_name = $cfg_value;
			
		}elsif($ion_name ne ""){
			$cfg{$s_name}{$ion_name}{$cfg_name} = $cfg_value;
		}
	}
}
close IN;


# setup directory
my $shdir=$dir.'/Shell';
system "mkdir -p $shdir" unless(-d $shdir);
my $batdir=$shdir.'/Batch';
system "mkdir -p $batdir" unless(-d $batdir);
my $samdir=$shdir.'/Sample';
system "mkdir -p $samdir " unless(-d $samdir);
#my $resdir=$dir.'/Result';
#system "mkdir -p $resdir" unless(-d $resdir);
## setup sub dir for each samples
foreach my $sample (sort keys %cfg) {
	my $specific_Samdir = $samdir."/".$sample;
	system "mkdir -p $specific_Samdir" unless(-d $specific_Samdir);
	foreach my $ion (sort keys %{$cfg{$sample}}){
		$cfg{$sample}{$ion}{'sam_dir'} = $specific_Samdir."/".$ion;
		system "mkdir -p $cfg{$sample}{$ion}{'sam_dir'}" unless(-d $cfg{$sample}{$ion}{'sam_dir'});
	}
}

######################################
#          STEP BY STEP PART         #
#          #################         #

open FIN,">$dir/pip.work.sh";

# step 1
if ($step =~ /1/) {
	my $batch_file = $batdir."/00.pretreatment.sh";
	my $result_dir = $dir."/00.data";
	open BATCH,">$batch_file";
	foreach my $sample (sort keys %cfg) {
		my $sample_sh = $samdir."/".$sample."/$sample.sh";
		my $unify_head_cmd = "\n";
		my $unify_tail_cmd = "\n";
		open UNIFY,">$sample_sh" || die "can't access to $sample_sh";
		foreach my $ion (sort keys %{$cfg{$sample}}){
			my $metabo_csv = $cfg{$sample}{$ion}{'metabo_csv'};
			my $new_ID_lst = $cfg{$sample}{$ion}{'new_ID_lst'};
			my $raw_ID_lst = $cfg{$sample}{$ion}{'raw_ID_lst'};
			my $ID_pfx     = $cfg{$sample}{$ion}{'ID_pfx'};
			my $sam_dir    = $cfg{$sample}{$ion}{'sam_dir'};
			my $res_pfx    = $result_dir."/".$sample."_".$ion;

			my $pretreatment_sh = $sam_dir."/00.pretreatment.$sample.$ion.sh";
			open PT,">$pretreatment_sh" || die "can't access to $pretreatment_sh.$!\n";
			print PT "sh $Bin/S01_1.pretreatment.sh $metabo_csv $new_ID_lst $raw_ID_lst $res_pfx $ID_pfx $idmap_pl $replace_pl $trim_pl $sort_pl > $res_pfx.pretreatment.res.xls\n";
			close PT;
			$unify_head_cmd .= "head -1 $res_pfx.pretreatment.res.xls|xargs -n1 > $res_pfx.newID.lst\n$ion=$res_pfx.newID.lst\n";
			$unify_tail_cmd .= "perl $trim_pl $result_dir/$sample.comm.lst $res_pfx.pretreatment.res.xls in $res_pfx.pretreatment.res.comm.xls h\n";
			print UNIFY "sh $pretreatment_sh \n";
		}
		print UNIFY "$unify_head_cmd\ncomm \$pos \$neg -1 -2 > $result_dir/$sample.comm.lst\n$unify_tail_cmd";
		close UNIFY;
		print BATCH "sh $sample_sh &\n";
	}
	close BATCH;
	print FIN "sh $batch_file\n";
}
close FIN;

		



######################################
#          SUB FUNCTION PART         #
#          #################         #

#usage: $parse_cfg($cfg_file,$cfg_hash,'circulation key variable'(if need));
sub parse_lib {
	my ($cfg_file, $cfg_hash, $key) = @_;
	open IN,"<$cfg_file" || die "fail open: $cfg_file";
	while(<IN>){
		chomp;
		next if (/^#/);
		if (/(\S+)\s*=\s*(\S+)/){
			my ($cfg_name,$cfg_value) = ($1,$2);
			$cfg_hash -> {$cfg_name} = $cfg_value;
		}
	}
	close IN;
}
