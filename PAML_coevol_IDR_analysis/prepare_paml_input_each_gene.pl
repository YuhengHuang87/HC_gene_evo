#!/usr/bin/perl
use strict;
use warnings;

my %lnL;
#opendir my $dir, "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/Clustal_to_PAML/" or die "Cannot open directory: $!";
opendir my $dir, "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/Clustal_PAML_random/" or die "Cannot open directory: $!";
my @files = readdir $dir;
closedir $dir;
my $num=0;
for (my $i=0; $i<@files; $i+=1){
	if ($files[$i] =~ m/.fa.codon/){

	my $speci_num=0; my $l=0;
		my @a=split(/\./, $files[$i]);
		my $id = $a[0];
		print "$id\n";
		#my $infile = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/Clustal_to_PAML/".$files[$i];
		my $infile = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/Clustal_PAML_random/".$files[$i];
		open(FILE1,"<", "$infile")||die"$!";
			while(my $count = <FILE1>){
		  	chomp($count);
				my @line=split(" ", $count);
					$l++;
					if ($l ==1){
						$speci_num = $line[0];
						last;
					}
				}

				if ($speci_num == 16){
					$num++;
					#my $outfile="/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/Clustal_to_PAML/input_files/sa".$num."/".$files[$i];
					my $outfile="/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/Clustal_PAML_random/input_files/sa".$num."/".$files[$i];
					open(OUT, ">$outfile");
					open(FILE1,"<", "$infile")||die"$!";
						while(my $count = <FILE1>){
							chomp($count);
							print OUT "$count\n";
					}

					#my $outfile1="/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/Clustal_to_PAML/input_files/sa".$num."/16_dro_updated.tree.txt";
					my $outfile1="/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/Clustal_PAML_random/input_files/sa".$num."/16_dro_updated.tree.txt";
					open(OUT1, ">$outfile1");

					my $file1 = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/Clustal_to_PAML/16_dro_updated.tree.txt";
					open(FILE1,"<", "$file1")||die"$!";
						while(my $count = <FILE1>){
							chomp($count);
							print OUT1 "$count\n";
				}

				#my $outfile2="/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/Clustal_to_PAML/input_files/sa".$num."/codeml.ctl";
				my $outfile2="/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/Clustal_PAML_random/input_files/sa".$num."/codeml.ctl";
				open(OUT2, ">$outfile2");
				print OUT2 "seqfile = ","$files[$i]\n";
				#print OUT2 "outfile = ","$id","_pairwise_output\n";
				#print OUT2 "outfile = ","$id","_pairwise_fix_omega_output\n";


				#print OUT2 "outfile = ","$id","_clean_output\n";
				#print OUT2 "outfile = ","$id","_clean_fix_omega_output\n";

				print OUT2 "outfile = ","$id","_pairwise_one_omega_output\n";


				$file1 = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/Clustal_to_PAML/codeml_update_file_names.txt";
				open(FILE1,"<", "$file1")||die"$!";
					while(my $count = <FILE1>){
						chomp($count);
						print OUT2 "$count\n";
			}
	}
	}
	}


my $infile;
opendir my $dir, "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/random_set_AA_Clustal_input" or die "Cannot open directory: $!";

my @files = readdir $dir;
closedir $dir;

my $id;
for (my $i=0; $i<@files; $i+=1){
	if ($files[$i] =~ m/_clustalo.out/){
    }elsif($files[$i] =~ m/.fa/){
		#my @a=split(".fa", $files[$i]);
		print "$files[$i]\n";
        my $num=$i-1;

	my $outfile="/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/random_set_AA_Clustal_input/AA_random_gene/sa".$num."/AA_".$files[$i];
    open(OUT, ">$outfile");
        my $line=0;
		$infile = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/random_set_AA_Clustal_input/".$files[$i];
		open(FILE1,"<", "$infile")||die"$!";
			while(my $count = <FILE1>){
		  	chomp($count);
				if ($count =~ m/>/){
                    if ($line ==0){
                        print OUT "$count\n";
                    }else{
                        print OUT "\n","$count\n";
                    }
                }else{
                    print OUT "$count";
                }
$line++;
}
}
}