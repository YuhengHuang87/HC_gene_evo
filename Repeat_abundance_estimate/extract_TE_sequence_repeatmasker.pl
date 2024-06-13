#!/usr/bin/perl
use strict;
use warnings;

my %TE_neu; my $name; my $posi;
my @TE_pacbio;

my $outfile="/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/Bia_TE_library.fasta"; #change for different species
open(OUT, ">$outfile");
my $file_name="/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/usftp21.novogene.com/raw_data/Biarmipe/ncbi_dataset/data/GCF_018148935.1/GCF_018148935.1_ASM1814893v1_genomic.fna";#change for different genomes
my $infile1=$file_name.".out";
open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split(" ", $count);
		if (($b[10] =~ /LINE/)||($b[10] =~ /LTR/)||($b[10] =~ /DNA/)||($b[10] =~ /Unknown/)){
			my $chr = $b[4]; my $left = $b[5]; my $right = $b[6];
	 $name=$b[10].":".$chr.":".$left.":".$right;
push @TE_pacbio,$name;
	for (my $i=$left; $i<=$right; $i+=1){
		$posi = $chr."\t".$i;
		$TE_neu{$posi}=$name;
}}
}

my $num_fam=scalar(@TE_pacbio);
print "$num_fam\n";

my %TE_seq; my $chr;my $posi_index=0;
my $infile2=$file_name;
open(FILE1,"<", "$infile2")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		if($count =~ m/>/){
			my @a = split(" ", $count);
			my @c = split(">", $a[0]);
			$chr = $c[1];
			 $posi_index=0;
	}else{
		my @d = split("", $count);
		for (my $i=0; $i < @d; $i++) {
			$posi_index++;
			my $site=$chr."\t".$posi_index;
			if (exists $TE_neu{$site}){
				if (exists $TE_seq{$TE_neu{$site}}){
					$TE_seq{$TE_neu{$site}}=$TE_seq{$TE_neu{$site}}.$d[$i];
				}else{
					$TE_seq{$TE_neu{$site}}=$d[$i];
				}
			}
	}}}


foreach my $key (keys %TE_seq) {
	print OUT ">","$key\n","$TE_seq{$key}\n";
}
