#!/usr/bin/perl -w


use strict;

my $gsea="/home/regmond/Scratch/software/gsea-3.0.jar";

use Env::Modulecmd { unload => 'mpi' };
use Env::Modulecmd { unload => 'compilers' };
use Env::Modulecmd { load => 'r/recommended' };

use Env::Modulecmd { unload => 'java' };
use Env::Modulecmd { load => 'java/1.8.0_45' };


my $min_set=15;
my $max_set=500;


my %header=();
my $data;
my @geneset=();
open(INF, "GSEA_matrix.txt");
while(<INF>){
	chomp $_;
	my @a=split("\t", $_);
	if($a[0] eq 'Gene_set'){
		for (my $i=1;$i<scalar @a; $i++){
			$header{$i}=$a[$i];
		}
	}else{
		push(@geneset, $a[0]);
		for (my $i=1;$i<scalar @a; $i++){
			push(@{$data->{$header{$i}}},$a[$i]);
		}
	}
}
close(INF);

#- I'm running gsea for each combination of rank+gmx, only because James sent the GSEA
#- 	matrix with numbers for each combination.
#- But I'm also running one per rank, with all the gene_sets together, and this is what I will use 
#-	downstream for plotting and for pvalues

foreach my $i(keys %header){
	my $rnk=$header{$i};
	
	my @data=@{$data->{$rnk}};

	($rnk =~/LUAD/) && next;

	my @genesets=();
	my $all_genesets="cat ";
	for (my $j=0;$j<scalar @data;$j++){

####		($data[$j] eq 'x') && next;

		$all_genesets.="data/${geneset[$j]}.gmt ";

		my $runstr="java -Xmx512m -cp $gsea xtools.gsea.GseaPreranked";
	        $runstr.=" -gmx data/${geneset[$j]}.gmt ";
	        $runstr.=" -collapse false -mode Max_probe -norm meandiv -nperm 1000";
	        $runstr.=" -rnk data/${rnk}.rnk";
	        $runstr.=" -scoring_scheme weighted";
	        #$runstr.=" -rpt_label ${rnk}_${geneset[$j]}";
	        $runstr.=" -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp";
	        $runstr.=" -set_max $max_set -set_min $min_set -zip_report false";
	        $runstr.=" -out ./$data[$j]_${rnk}_${geneset[$j]}";
	        $runstr.=" -gui false";
#	        print $runstr."\n\n";
#		`$runstr`;
	}
	$all_genesets.=" > all.gmt";
	`$all_genesets`;

	print "run_gsea.pl --rnk data/${rnk}.rnk --gmx all.gmt --out all --min_set $min_set --max_set $max_set\n";
	system("run_gsea.pl --rnk data/${rnk}.rnk --gmx all.gmt --out all --min_set $min_set --max_set $max_set");
	system("rm all.gmt");
	
}


my %genes=();
my %ranks=();
my %leadingedgegenes=();
open(OUTF, ">gsea_table.txt");
print OUTF "RANK\tGENESET\tGENESET_ORIGINAL_SIZE\tGENESET_SIZE_IN_DATA\tLEADING_EDGE_GENES\tRATIO\tRATIO_ORIGINAL_SIZE\tFDR_p\tNES\n";
opendir(DIR2, "./");
my @dirs2=grep { /^all_.*_all/ } readdir (DIR2);
foreach my $dir2(@dirs2){
	opendir(DIR, "./$dir2");
	my $rnk;
	if($dir2 =~/^all_(.*)_all/){
		$rnk=$1;
	}
	$ranks{$rnk}=1;
	my @dirs= grep { /^my_analysis.GseaPreranked/ } readdir (DIR);
	my $dir=$dirs[0];
	#- orioginal geneset sizes
	my %originalsize=();
	open(INF, "./$dir2/$dir/gene_set_sizes.xls");
	while(<INF>){
		chomp $_;
		my @a=split("\t",$_);
		($a[0] eq 'NAME') && next;
		$originalsize{$a[0]}=$a[1];
	}
	close(INF);
	opendir(D, "./$dir2/$dir");
	my @files= grep { /.xls$/ && /^gsea_report_for_na/} readdir (D);
	foreach my $file(@files){
		open(INF, "./$dir2/$dir/$file") || die "$! ./$dir2/$dir/$file\n";
		while(<INF>){
			chomp $_;
			my @a=split("\t",$_);
			($a[0] eq 'NAME') && next;
			
			my $count=0;
			open(INF2, "./$dir2/$dir/$a[0].xls");
			while(<INF2>){
				chomp $_;
				my @a=split("\t",$_);
				($a[0] eq 'NAME') && next;
				
				$genes{$a[1]}=1;
				
				if($a[7] eq 'Yes'){
					$leadingedgegenes{$rnk}{$a[1]}=1;
					$count++;
				}
			}
			close(INF2);
			my $ratio=sprintf("%.2f", $count/$a[3]);
			my $ratiooriginal=sprintf("%.2f", $count/$originalsize{$a[0]});
			print OUTF $rnk."\t".$a[0]."\t".$originalsize{$a[0]}."\t".$a[3]."\t".$count."\t".$ratio."\t".$ratiooriginal."\t".$a[7]."\t".$a[5]."\n";
		}
		close(INF);
	}
	closedir(D);
	closedir(DIR);
}
closedir(DIR2);
close(OUTF);

`R CMD BATCH plot.R`;

my @ranks=sort (keys %ranks);
my @genes=sort (keys %genes);
open(OUTF, ">circlize_table.txt");
print OUTF "Gene\t".join("\t", @ranks)."\n";
foreach my $gene(@genes){
	my @data=();
	foreach my $rank(@ranks){
		if($leadingedgegenes{$rank}{$gene}){
			push(@data, "1");
		}else{
			push(@data, "0");
		}
	}
	print OUTF $gene."\t".join("\t", @data)."\n";
}
close(OUTF);

 
