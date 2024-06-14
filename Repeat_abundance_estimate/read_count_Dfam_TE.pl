#!/usr/bin/perl
#use strict;
use warnings;
use List::Util qw(sum);

my $outfile="TE_counts_no_dup_Dfam_species";
open(OUT, ">$outfile");
my $read=0;my $chr;my $file = '';
my %neu;
$file="out_sorted_Dfam_species.sam";#need to use the correct genomes located in the species's folder
open(FILE1,"<", "$file")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
      my @b = split("\t", $count);
			    if($count =~ m/@/){
}else{
		#if (($b[4] > 20)&&($b[2] =~ m/DF/)){
		if ($b[2] ne "*"){
			if (exists $neu{$b[0]}){}else{
				$neu{$b[0]}=$b[2];
			$read++;
		#}
}
}}
}
print OUT "$read\n";
