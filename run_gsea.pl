#!/usr/bin/perl -w


use strict;
use warnings;

no warnings qw/uninitialized/;

use Getopt::Long qw(GetOptionsFromArray :config pass_through no_auto_abbrev bundling);
use File::Basename;

use Env::Modulecmd { unload => 'java' };
use Env::Modulecmd { load => 'java/1.8.0_45' };


my @args=@ARGV;

my $rnk='';
my $gmx='';	
my $min_set="";
my $max_set="";
GetOptionsFromArray (
    \@args,
    "rnk=s" => \$rnk,
    "gmx=s"   => \$gmx,
    "min_set=i"   => \$min_set,
    "max_set=i"   => \$max_set,
    "<>"   => \&print_usage
) or die "\n";

if (!-e "$rnk") {
	usage_gsea("Can't find rnk file <$rnk>");
	return;
}		
if (!-e "$gmx") {
	usage_gsea("Can't find gmx (gene set) file <$gmx>");
	return;
}		
if($min_set ne ''){
	if(!isnum($min_set) || ($min_set < 0)){
		usage_gsea("<$min_set> is not a valid number. Please insert an integer > 0");
		return;
	}	
}
if($max_set ne ''){
	if(!isnum($max_set) || ($max_set < 0)){
		usage_gsea("<$max_set> is not a valid number. Please insert an integer > 0");
		return;
	}	
}

my @suffix=(".txt",".gmt",".gmx",".rnk");
my $basename_gmx = basename($gmx, @suffix);
my $basename_rnk = basename($rnk, @suffix);

my $runstr="java -Xmx512m -cp ./gsea3.0/gsea-3.0.jar xtools.gsea.GseaPreranked";
$runstr.=" -gmx $gmx";
$runstr.=" -collapse false -mode Max_probe -norm meandiv -nperm 1000";
$runstr.=" -rnk $rnk";
$runstr.=" -scoring_scheme weighted";
$runstr.=" -rpt_label ${basename_rnk}_${basename_gmx}";
$runstr.=" -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp";
$runstr.=" -set_max 2000 -set_min 3 -zip_report false";
#$runstr.=" -out ./gsea_results";
$runstr.=" -gui false";
`$runstr`;




sub print_usage {

	my @a=@_;

	my $usage0="\t";
	my $usage1="\tProgram: Run_GSEA3.0";
	my $usage2 = "\tDescription: Runs gene set enrichment analysis (GSEA3.0) given a ranked list of genes (-rnk) and a gene set (-gmt).";
	my $usage3 = "\tContact: Lucia Conde <l.conde\@ucl.ac.uk>";
	my $usage4="\t";
	my $usage5="\tUsage:   run_gsea --rnk FILE.rnk --gmt FILE.gmt [options]";
	my $usage6="\t";
	my $usage7="\tRequired: --rnk FILE	Ranked list of genes. Needs a column with gene names and a column with the stat";
	my $usage8="\t          --gmt FILE		Gene sets in GMT format (https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29)";
	my $usage9="\t";
	my $usage10="\tOptions: --set_min NUM     Ignore gene sets that contain less than NUM genes [15]";
	my $usage11="\t         --set_max NUM	Ignore gene sets that contain more than NUM genes [500]";
        my $usage12="\t";

	print "$usage0\n$usage1\n$usage2\n$usage3";
	print "$usage4\n$usage5\n$usage6\n$usage7\n";
	print "$usage8\n$usage9\n$usage10\n";
	print "$usage11\n$usage12\n$usage12\n";

	die "ERROR: Invalid argument '@a'. Please check above the usage of this script\n";

}

sub usage_gsea {
	my $error=shift;
	die qq(
	Program: Run_GSEA3.0
	Description: Runs gene set enrichment analysis (GSEA3.0) given a ranked list of genes (-rnk) and a gene set (-gmt)
	Contact: Lucia Conde <l.conde\@ucl.ac.uk>";

	Usage:   run_gsea --rnk FILE.rnk --gmt FILE.gmt [options]";

	Required: --rnk FILE	Ranked list of genes. Needs a column with gene names and a column with the stat";
	          --gmt FILE		Gene sets in GMT format (https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29)";
	Options: --set_max NUM	Ignore gene sets that contain more than NUM genes [500]";
	         --set_min NUM	Ignore gene sets that contain less than NUM genes [15]";

	ERROR: $error

	);
}
