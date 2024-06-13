#!/usr/bin/perl
#use strict;
use warnings;
use List::Util qw(sum);

my $outfile="median_coverage";
open(OUT, ">$outfile");
my $posi_index=0;my $chr;my $file = '';
my %neu;
$file="GCF_018152265.1_ASM1815226v1_genomic.fna.masked";#need to use the correct genomes!! 
open(FILE1,"<", "$file")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
    if($count =~ m/>/){
      my @b = split(" ", $count);
      #print "$b[0]\n";
      my @c = split(">", $b[0]);
  	   $chr=$c[1];
       $posi_index=0;
}else{
    my @b = split("", $count);
		for (my $i=0; $i < @b; $i++) {
			$posi_index++;
      my $posi=$chr."\t".$posi_index;
			my $nuclotide=$b[$i];
			if ($nuclotide ne "N"){
			$neu{$posi}=$nuclotide;
}}}
}

my @cov;
	$file = "output_sorted_depth.txt";
	  open I, "<$file" or print "Can't open /$file\n";
		while(my $count = <I>){
	  	chomp($count);
	    my @b = split("\t", $count);
			my $site=$b[0]."\t".$b[1];
			if (exists $neu{$site}){
				push @cov, $b[-1];
			}
		}

my $size=scalar(@cov);
my $median=median(@cov);
print OUT "GCF_018152265.1_ASM1815226v1_genomic.fna.masked\t","$size\t","$median\n";



sub average {
my @array = @_; # save the array passed to this function
my $sum; # create a variable to hold the sum of the array's values
foreach (@array) { $sum += $_; } # add each element of the array
# to the sum
return $sum/@array; # divide sum by the number of elements in the
# array to find the mean
}


sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}
