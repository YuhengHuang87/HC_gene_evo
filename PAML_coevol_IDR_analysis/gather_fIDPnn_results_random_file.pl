#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum);

my %lnL;
my $infile;

my $outfile="fIDP_results_gene_random_1000_percentage.txt";
open(OUT, ">$outfile");
print OUT "ID\t","D_melanogaster\t","D_ananassae\t","D_biarmipes\t","D_bipectinata\t","D_elegans\t","D_erecta\t","D_eugracilis\t","D_ficusphila\t","D_kikkawai\t","D_rhopaloa\t","D_sechellia\t","D_serrata\t","D_simulans\t","D_suzukii\t","D_takahashii\t","D_yakuba\n";



for (my $num=1; $num<=1000; $num+=1){
opendir my $dir, "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/random_set_AA_Clustal_input/AA_random_gene/sa".$num or die "Cannot open directory: $!";
my @files = readdir $dir;
closedir $dir;
print "$num\n";
my $id; my $dis_order=0; my $total=0; my $idr_per; my $run=0;
for (my $i=0; $i<@files; $i+=1){
    	if ($files[$i] =~ m/.fa/){
        my @a=split(/\./, $files[$i]);
        my @b=split("_",$a[0]);
		$id = $b[1];
        }
	if ($files[$i] =~ m/results.csv/){	
        my @idr=(); $run++;
		$infile = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/random_set_AA_Clustal_input/AA_random_gene/sa".$num."/".$files[$i];
		open(FILE1,"<", "$infile")||die"$!";
			while(my $count = <FILE1>){
		  	chomp($count);
   
        if ($count =~ m/Residue/){}else{
		if ($count =~ m/>/){
            if ($count =~ m/melanogaster/){}else{
                $idr_per=$dis_order/$total;
                push @idr, $idr_per;
            }
            $dis_order=0; $total=0;
	}else{
	my @b = split(",", $count);
    #print "@b\n";
    if ($b[1] ne "U"){
        $total++;
        if ($b[2]>0.3){
            $dis_order++;
        }
    }
    }}
    }
$idr_per=$dis_order/$total;
push @idr, $idr_per;
if (@idr ==16){
    print OUT "$id\t";
    my $j=0;
    for ($j=0; $j<(@idr-1); $j+=1){
        print OUT "$idr[$j]\t";
    }
    print OUT "$idr[$j]\n";
}
}
}
}