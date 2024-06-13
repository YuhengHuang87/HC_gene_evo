#!/usr/bin/perl
use strict;
use warnings;

my %lnL;
my $infile;

opendir my $dir, "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/Clustal_input/" or die "Cannot open directory: $!";
my @files = readdir $dir;
closedir $dir;

for (my $i=0; $i<@files; $i+=1){
	if ($files[$i] =~ m/clustalo.out/){
	}else{
		if ($files[$i] =~ m/.fa/){
			my $outfile1="/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/Clustal_AA_input/".$files[$i];
			open(OUT1, ">$outfile1");
		my $j=0; my $nucl_seq; my $triple;
		$infile = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/Clustal_input/".$files[$i];
		open(FILE1,"<", "$infile")||die"$!";
			while(my $count = <FILE1>){
		  	chomp($count);
				if ($count =~ m/>/){
					if ($count eq ">melanogaster"){
					print OUT1 "$count\n";
					$triple='';
					$nucl_seq='';
				}else{
					my $num_nu=0;
					my @nuc = split("", $nucl_seq);
						for (my $k=0; $k<@nuc; $k+=1){
							my $site = $k+1;
							$triple=$triple.$nuc[$k];
							if ($site%3==0){
								#print "$triple\n";
								$num_nu++;
								my $aa = codon2aa($triple);
								print OUT1 "$aa";
								if ($num_nu%60==0){
								print OUT1 "\n";
								}
								$triple='';
							}
						}
					print OUT1 "\n";
					$triple='';
					$nucl_seq='';
					print OUT1 "$count\n";
				}
					}else{
						$nucl_seq=$nucl_seq.$count;
							}
				}

						my $num_nu=0;
						my @nuc = split("", $nucl_seq);
							for (my $k=0; $k<@nuc; $k+=1){
								my $site = $k+1;
								$triple=$triple.$nuc[$k];
								if ($site%3==0){
									#print "$triple\n";
									$num_nu++;
									my $aa = codon2aa($triple);
									print OUT1 "$aa";
									if ($num_nu%60==0){
									print OUT1 "\n";
									}
									$triple='';
								}
							}
						print OUT1 "\n";
				}
			}
		}

opendir my $dir, "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/All_gene_seq" or die "Cannot open directory: $!";
my @files = readdir $dir;
closedir $dir;
my $num=0; my %select_gene;
while ($num<1000) {
my $randf = $files[rand @files];
if (exists $select_gene{$randf}){}else{
	my $speci_num=0; my $l=0;
		my @a=split(/\./, $randf);
		my $id = $a[0];
		print "$id\n";
		my $infile = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/All_gene_seq/".$randf;
		open(FILE1,"<", "$infile")||die"$!";
			while(my $count = <FILE1>){
		  	chomp($count);
				if ($count =~ m/>/){
					$l++;
					}
				}
				if ($l == 16){
					$num++;
					$select_gene{$randf}=1;
					my $outfile="/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/random_set_for_Clustal/".$randf;
					open(OUT, ">$outfile");
					open(FILE1,"<", "$infile")||die"$!";
						while(my $count = <FILE1>){
							chomp($count);
							print OUT "$count\n";
						}
						
					my $outfile1="/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/random_set_AA_Clustal_input/".$randf;
					open(OUT1, ">$outfile1");
				my $j=0; my $nucl_seq; my $triple;
				open(FILE1,"<", "$infile")||die"$!";
					while(my $count = <FILE1>){
						chomp($count);
						if ($count =~ m/>/){
							if ($count eq ">melanogaster"){
							print OUT1 "$count\n";
							$triple='';
							$nucl_seq='';
						}else{
							my $num_nu=0;
							my @nuc = split("", $nucl_seq);
								for (my $k=0; $k<@nuc; $k+=1){
									my $site = $k+1;
									$triple=$triple.$nuc[$k];
									if ($site%3==0){
										#print "$triple\n";
										$num_nu++;
										my $aa = codon2aa($triple);
										print OUT1 "$aa";
										if ($num_nu%60==0){
										print OUT1 "\n";
										}
										$triple='';
									}
								}
							print OUT1 "\n";
							$triple='';
							$nucl_seq='';
							print OUT1 "$count\n";
						}
							}else{
								$nucl_seq=$nucl_seq.$count;
									}
						}

								my $num_nu=0;
								my @nuc = split("", $nucl_seq);
									for (my $k=0; $k<@nuc; $k+=1){
										my $site = $k+1;
										$triple=$triple.$nuc[$k];
										if ($site%3==0){
											#print "$triple\n";
											$num_nu++;
											my $aa = codon2aa($triple);
											print OUT1 "$aa";
											if ($num_nu%60==0){
											print OUT1 "\n";
											}
											$triple='';
										}
									}
								print OUT1 "\n";
						}
					}
				}


sub codon2aa{
		    my($codon) = @_;
		    $codon = uc $codon;
		    my(%genetic_code) = (
		    TCA => 'S',    # Serine
		    TCC => 'S',    # Serine
		    TCG => 'S',    # Serine
		    TCT => 'S',    # Serine
		    TTC => 'F',    # Phenylalanine
		    TTT => 'F',    # Phenylalanine
		    TTA => 'L',    # Leucine
		    TTG => 'L',    # Leucine
		    TAC => 'Y',    # Tyrosine
		    TAT => 'Y',    # Tyrosine
		    TAA => 'X',    # Stop
		    TAG => 'X',    # Stop
		    TGC => 'C',    # Cysteine
		    TGT => 'C',    # Cysteine
		    TGA => 'X',    # Stop
		    TGG => 'W',    # Tryptophan
		    CTA => 'L',    # Leucine
		    CTC => 'L',    # Leucine
		    CTG => 'L',    # Leucine
		    CTT => 'L',    # Leucine
		    CCA => 'P',    # Proline
		    CCC => 'P',    # Proline
		    CCG => 'P',    # Proline
		    CCT => 'P',    # Proline
		    CAC => 'H',    # Histidine
		    CAT => 'H',    # Histidine
		    CAA => 'Q',    # Glutamine
		    CAG => 'Q',    # Glutamine
		    CGA => 'R',    # Arginine
		    CGC => 'R',    # Arginine
		    CGG => 'R',    # Arginine
		    CGT => 'R',    # Arginine
		    ATA => 'I',    # Isoleucine
		    ATC => 'I',    # Isoleucine
		    ATT => 'I',    # Isoleucine
		    ATG => 'M',    # Methionine
		    ACA => 'T',    # Threonine
		    ACC => 'T',    # Threonine
		    ACG => 'T',    # Threonine
		    ACT => 'T',    # Threonine
		    AAC => 'N',    # Asparagine
		    AAT => 'N',    # Asparagine
		    AAA => 'K',    # Lysine
		    AAG => 'K',    # Lysine
		    AGC => 'S',    # Serine
		    AGT => 'S',    # Serine
		    AGA => 'R',    # Arginine
		    AGG => 'R',    # Arginine
		    GTA => 'V',    # Valine
		    GTC => 'V',    # Valine
		    GTG => 'V',    # Valine
		    GTT => 'V',    # Valine
		    GCA => 'A',    # Alanine
		    GCC => 'A',    # Alanine
		    GCG => 'A',    # Alanine
		    GCT => 'A',    # Alanine
		    GAC => 'D',    # Aspartic Acid
		    GAT => 'D',    # Aspartic Acid
		    GAA => 'E',    # Glutamic Acid
		    GAG => 'E',    # Glutamic Acid
		    GGA => 'G',    # Glycine
		    GGC => 'G',    # Glycine
		    GGG => 'G',    # Glycine
		    GGT => 'G',    # Glycine
		    );
		    if(exists $genetic_code{$codon}) {
		        return $genetic_code{$codon};
		    }else{
		        return 'X';
		    }
		}
