#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum);

#download coding sequences, amino acid sequences, and genome annotations for 16 studied species. Unzip and save the cds_from_genomic.fna for each species in "$species_name/ncbi_dataset/data/". Download orthologous information from OrthoDB, save the files odb11v0_genes.tab and odb11v0_genes.tab

my $infile;
my %speci=(ananassae =>'7217_0',erecta=>'7220_0',melanogaster=>'7227_0',sechellia=>'7238_0',simulans=>'7240_0',yakuba=>'7245_0',serrata=>'7274_0',suzukii=>'28584_0',eugracilis=>'29029_0',takahashii=>'29030_0',elegans=>'30023_0',ficusphila=>'30025_0',kikkawai=>'30033_0',bipectinata=>'42026_0',biarmipes=>'125945_0',rhopaloa=>'1041015_0');

my %orth;
$infile = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/odb11v0_OG2genes.tab";
open(FILE,"<", "$infile")||die"$!";
	while(my $count = <FILE>){
  	chomp($count);
		my @b = split("\t", $count);
		$orth{$b[1]}=$b[0];
		}
close FILE;


foreach my $id (keys %speci) {
my %spec_gene_id;

my $outfile=$id."_orth_coding_id.txt";
open(OUT, ">$outfile");
my $i=0; my %list; my $Iso;
my $file1 = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/odb11v0_genes.tab";
open(FILE1,"<", "$file1")||die"$!";
	while(my $count = <FILE1>){
  	chomp($count);
my @a=split("\t", $count);

if ($a[1] eq $speci{$id}){
	if (exists $orth{$a[0]}){
print OUT "$orth{$a[0]}\t","$a[0]\t","$a[2]\n";
}
}
}
close FILE1;
}


my @speci=('melanogaster','ananassae','erecta','sechellia','simulans','yakuba','serrata','suzukii','eugracilis','takahashii','elegans','ficusphila','kikkawai','bipectinata','biarmipes','rhopaloa');
for (my $j=0; $j<@speci; $j+=1){
my %refIso;
$infile = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/".$speci[$j]."_orth_coding_id.txt";
open(FILE1,"<", "$infile")||die"$!";
	while(my $count = <FILE1>){
  	chomp($count);
		my @b = split("\t", $count);
my $rna=$b[2];
	$refIso{$rna}=$b[0];
	}
close FILE1;


my $outfile="/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/16species_orthologs/".$speci[$j]."_cds_orth.fa";
open(OUT, ">$outfile");

my $i=0; my %list; my $Iso;
my $file1 = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/".$speci[$j]."/ncbi_dataset/data/cds_from_genomic.fna";
open(FILE1,"<", "$file1")||die"$!";
	while(my $count = <FILE1>){
if($count eq "\n"){
	}else{
  	chomp($count);
if ($count =~ />/){
my @line=split(" ", $count);
my $gene; my $protein;
foreach my $line (@line){
	if ($line =~ /gene=/){
		my @a = split("=", $line);
		my @b = split("]", $a[1]);
		$gene=$b[0];
	}
	if ($line =~ /protein_id=/){
		my @a = split("=", $line);
		my @b = split("]", $a[1]);
		$protein=$b[0];
	}
}
if (exists $refIso{$protein}){
	$Iso=$speci[$j].",".$gene.",".$protein;
	print OUT ">","$Iso\n";
$i++;
}else{
$i=0;
}
}else{
if ($i>0){
	print OUT "$count\n";
}
}
}
}
}