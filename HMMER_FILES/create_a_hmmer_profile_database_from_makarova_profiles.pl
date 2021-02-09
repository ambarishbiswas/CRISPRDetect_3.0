#!/usr/bin/perl -w
use strict;


#---- first open the makarova_summary_table_of_cas_genes.tab file and get the name of the profile_file containg the sequence alignement; then create a .hmm profile using the makarova alignment file
my $skip_1=1;
if($skip_1==0)
{
my @arr_info_lines=`grep -v '^#' makarova_summary_table_of_cas_genes.tab >&1`;

foreach my $line(@arr_info_lines)
	{
		chomp $line; $line=~s/\r//;
		my @arr_t1=split('\t',$line);
		my $profile=$arr_t1[0];
		
		my $aln_file=$profile.".sr";
		
		my $cas_gene_name=$arr_t1[3];
		my $cas_gene_function=$arr_t1[2];
		
		my $output_hmm_file=$cas_gene_name."-".$cas_gene_function."-".$profile.".hmm";
		
		#system("hmmbuild --amino --cpu 50 Cas_HMMs\/$output_hmm_file makarova_cas_NRM2015_profiles\/$aln_file");
		
		#last;
		
	}
#---- now cat all the hmms into one file 
system("cat Cas_HMMs/*.hmm >all_cas_profiles.hmm");

#--- then compress the file -----------
system("hmmpress all_cas_profiles.hmm");

}



#------ search all the sub sampled sequence files with the Cas_gene's hmm profile

my $skip_2=1;
if($skip_2==0)
{
#----- finally search the translated metagenomic reads against the hmm database ----
my $metagenomic_sample1="/PROJECTS/DATA2/chris/Other_files/seqs/hot_spring_water_sample1.fa";
my $translated_metagenomic_sample1="/tmp/"."hot_spring_water_sample1.fa.translated";

system("transeq -sequence $metagenomic_sample1 -outseq $translated_metagenomic_sample1 -frame 6 -clean");


my $hmmer_tabular_output_file=$translated_metagenomic_sample1.".tblout";
system("hmmsearch --tblout $hmmer_tabular_output_file -E 1e-5 --cpu 50 all_cas_profiles.hmm $translated_metagenomic_sample1 >/dev/null");
}


#--- now create a collectors curve [similar way to the repeat spacer collectors curve analysis]












