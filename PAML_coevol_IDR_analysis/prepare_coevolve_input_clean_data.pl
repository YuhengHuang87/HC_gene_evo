#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum);

my $infile;

my @sig_gene=();
my $num=0; my %id;
$infile = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/PAML_pairwise_output_likelihood_focal_sig.txt";
open(FILE1,"<", "$infile")||die"$!";
while(my $count = <FILE1>){
	chomp($count);
	my @b = split("\t", $count);
my $name = $b[0];
push @sig_gene, $b[0];
#$num++;
$id{$name}=1;
}

$infile = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/paml_results_gene_random_1000_pairwise_significant_omega";
open(FILE1,"<", "$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		#push @sig_gene, $count;
        #$num++;
	}
print "$num\n";

	for (my $k=0; $k<@sig_gene; $k+=1){
        my $file_name=$sig_gene[$k].".fa.codon";
        my $index=$k+1;
        print "$file_name\n";
        
        my $outfile2="/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/coevol/focal/sa".$index."/coevol_run.sub";

				open(OUT2, ">$outfile2");
				my $line=0;

                my $file1 = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/coevol_run.sub";

				open(FILE1,"<", "$file1")||die"$!";
					while(my $count = <FILE1>){
						chomp($count);
                        my @b = split(" ", $count);
                        if ($count =~ m/#/){ 
						print OUT2 "$count\n";
			}elsif($count =~ m/coevol-master/){
            #print OUT2 "$b[0] ","$b[1] ","$file_name ","$b[3] ","$b[4] ","$b[5] ","$b[6] ","$b[7] ","$file_name","_run4\n";
            print OUT2 "$b[0] ","$b[1] ","$b[2] ","$file_name","_run4\n";
            }
	}}



