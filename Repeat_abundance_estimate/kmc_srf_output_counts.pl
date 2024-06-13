#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum);

my $infile;my $bp = 0;
my @speci=('Ana','Bia','Bipectin','Elegans','Erecta','Eug','Ficusphi','Kikkawai','Mel','Rho','San','Sechelli','Serrata','Sim','Suzukii','Takakash','Yakuba');
my $outfile="/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/kmc_srf/srf_counts_abundance.txt";
open(OUT, ">$outfile");



for (my $i=0; $i<@speci; $i+=1){
$bp = 0;
my $file1 = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/kmc_srf/".$speci[$i]."_srf.abun";
open(FILE1,"<", "$file1")||die"$!";
	while(my $count = <FILE1>){
my @line=split("\t", $count);
$bp+=$line[1];
    }
	print OUT "$speci[$i]\t","$bp\n";
}

