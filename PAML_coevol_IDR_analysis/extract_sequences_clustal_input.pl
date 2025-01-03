#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum);

my $infile;
my @id=('Su(var)2-10','Su(var)205','Su(var)3-9','Su(var)2-HP2','HP4','Su(var)3-7','SuUR','Hmt4-20','crol','blanks','CG8108','piwi','Pep','jumu','ADD1','D1','prod','XNP','rhi','cuff','del','moon','wde','HP6','Pof','Hmr','Fkbp39','CG2129','Hrb87F','l(3)neo38','Tlk','Hr83','Atf-2','Cap-G','MED17','MED26','Rif1','vtd','SA','Orc2','mei-S332','Nipped-B','NSD','e(y)3','vig','Gas41','stwl','vig2','dmt','HP1b','Arp6','Stat92E','eff','HIPP1','cav','Trf2','mu2','Kdm2','esc','Sfmbt','ph-p','Pc','Trl','Jarid2','Psc','ph-d','Caf1-55','Pcl','Sce','Spps','Su(z)12','cg','calypso','Asx','jing','Scm','pho','Gcn5','Sgf29','wda','wde','Panx','arx','nxf2','Sam-S','Nap1','jumu','HDAC1','crol','SuUR','Hmt4-20','Su(var)2-HP2','ADD1','CG8108','Dp1','MTA1-like','chm','HP5','Kdm4A','Lhr','HP4','Odj','XNP','G9a','nub','egg','AGO2','Su(var)3-3','bin3','Su(var)3-9','Hsc70-4','bon','mbo','Sec13','Alh','velo','mh','sov','Rrp6','Mcm10','eyg','Su(z)2','E(z)','phol','Pp1-87B','RpLP0','Sam-S','snf','wapl','CDK2AP1','cg','E2f1','fl(2)d','HDAC11','HDAC3','HDAC4','HDAC6','Hel25E','His2Av','Hnf4','JHDM2','JIL-1','Kdm4B','kis','l(3)neo38','lid','Lsd-1','mael','MBD-like','MEP-1','Mi-2','mod','mof','nej','NO66','Nxt1','ova','PR-Set7','simj','Sirt1','Sirt2','Sirt4','Sirt6','Sirt7','Spps','Su(z)2','Usp47','vig','ttk','Khc','Tip60');
my %focal;
for (my $i=0; $i<@id; $i+=1){
my $id = $id[$i];
$focal{$id}=1;
}

my @speci=('melanogaster','ananassae','erecta','sechellia','simulans','yakuba','serrata','suzukii','eugracilis','takahashii','elegans','ficusphila','kikkawai','bipectinata','biarmipes','rhopaloa');

my %orth_melIso;my %melIso_orth;
$infile = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/melanogaster_orth_coding_id.txt";
open(FILE1,"<", "$infile")||die"$!";
while(my $count = <FILE1>){
	chomp($count);
	my @b = split("\t", $count);
	if (exists $orth_melIso{$b[0]}){
		$orth_melIso{$b[0]}=$orth_melIso{$b[0]}."\t".$b[2];
	}else{
	$orth_melIso{$b[0]}=$b[2];
}
	$melIso_orth{$b[2]}=$b[0];
}

my $outfile1="/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/dup_gene_random";
open(OUT1, ">$outfile1");

foreach my $key (keys %melIso_orth) {
my $seq=''; my $m=0;
my %orth_ID; my $orth_Iso; my $ID;
$infile = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/16species_orthologs/melanogaster_cds_orth.fa";
open(FILE1,"<", "$infile")||die"$!";
	while(my $count = <FILE1>){
  	chomp($count);
		if($count =~ m/>/){
		my @b = split(",", $count);
		if ($b[2] eq $key){
			$m=1;
			$orth_Iso = $melIso_orth{$b[2]};
			if (exists $orth_ID{$orth_Iso}){
				$orth_ID{$orth_Iso} = $orth_ID{$orth_Iso}."\t".$b[1];
				#print OUT1 "$b[1]\t","$b[2]\t","$focal_dup{$b[2]}\n";
			}else{
				$orth_ID{$orth_Iso} = $b[1];
			}
}else{
	$m=0;
}
}else{
	if ($m==1){
		$seq=$seq.$count;
	}
}}
my $leng=length($seq);
close FILE1;
#print "$leng\n";

my $dup=0;
$infile = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/16species_orthologs/melanogaster_cds_orth.fa";
open(FILE1,"<", "$infile")||die"$!";
	while(my $count = <FILE1>){
  	chomp($count);
		if($count =~ m/>/){
		my @b = split(",", $count);
		if ($b[2] eq $key){
		$orth_Iso = $melIso_orth{$b[2]};
		$ID = $b[1];
		if (exists $orth_ID{$orth_Iso}){
			my @iso = split("\t", $orth_ID{$orth_Iso});

		for (my $k=0; $k<@iso; $k+=1){
			if ($b[1] ne $iso[$k]){
				print OUT1 "$b[1]\t","$b[2]\n";
				$dup++;
			}
}}}}}
print "$dup\n";
if ($dup ==0){

if (exists $focal{$ID}){
my $outfile="/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/Clustal_input/".$ID.".fa";
open(OUT, ">$outfile");
}else{
my $outfile="/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/All_gene_seq/".$ID.".fa";
open(OUT, ">$outfile");
}

for (my $j=0; $j<@speci; $j+=1){
	my $hit=0;
	my @spec_iso;

$infile = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/".$speci[$j]."_orth_coding_id.txt";
open(FILE1,"<", "$infile")||die"$!";
while(my $count = <FILE1>){
	chomp($count);
	my @b = split("\t", $count);
	if ($b[0] eq $orth_Iso){
		push @spec_iso, $b[2];
		$hit++;
	}
}

#print OUT2 "$speci[$j]\t","$hit\t","$spec_iso\t";
my $iso_id; my $min_dif_length=$leng;
foreach my $spec_iso (@spec_iso){
	my $check=0;my $spec_seq;
	$infile = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/16species_orthologs/".$speci[$j]."_cds_orth.fa";
	open(FILE1,"<", "$infile")||die"$!";
		while(my $count = <FILE1>){
	  	chomp($count);
			if($count =~ m/>/){
			my @b = split(",", $count);
			if ($b[2] eq $spec_iso){
				$check++;
			}else{
				$check=0;
			}
		}else{
			if ($check == 1){
				$spec_seq=$spec_seq.$count;
			}
	}}
	my $len_dif=abs(length($spec_seq)-$leng);
	if ($len_dif<$min_dif_length){
		$iso_id=$spec_iso;
		$min_dif_length=$len_dif;
	}
}
#print "$iso_id\n";
	my $check=0;
$infile = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/16species_orthologs/".$speci[$j]."_cds_orth.fa";
open(FILE1,"<", "$infile")||die"$!";
	while(my $count = <FILE1>){
  	chomp($count);
		if($count =~ m/>/){
		my @b = split(",", $count);
		if ($b[2] eq $iso_id){
			print OUT "$b[0]\n";
			#print OUT "$count\n";
			$check++;
		}else{
			$check=0;
		}
	}else{
		if ($check == 1){
			print OUT "$count\n";
		}
}}
}
}
}
