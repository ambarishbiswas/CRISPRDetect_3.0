#!/usr/bin/perl -w
use strict;


######  README #################
# This script will Predict Genes with GeneMarkS and identify potential Cas genes with hmmsearch
#################################


my $script_dir="";
my $tmp_dir="/tmp";
my $input_seq_file="";
my $all_gene_positions_file="";
my $predicted_cas_gene_position_and_sequences="";
my $hmm_seq_evalue_cutoff=0.00001;	#0.00001; # -E 1e-5 
my $hmm_dom_evalue_cutoff=0.00001;	#0.00001; # -E 1e-5 

#----------------------------------------------------
my $no_of_threads=20;


#-----------------------------------------------------



for(my $i=0;$i<=$#ARGV;$i++)
	{
						
		#--------------------------- input options ------------------------------------------------------------
				if($ARGV[$i]=~/-script_dir$/)
					{
						$script_dir=$ARGV[$i+1];					
					}
				elsif($ARGV[$i]=~/-tmp$/)
					{
						$tmp_dir=$ARGV[$i+1];					
					}	
				elsif($ARGV[$i]=~/-i$/)
					{
						$input_seq_file=$ARGV[$i+1];					
					}	
				elsif($ARGV[$i]=~/-o$/)
					{
						$all_gene_positions_file=$ARGV[$i+1];					
					}
				elsif($ARGV[$i]=~/-s$/)
					{
						$predicted_cas_gene_position_and_sequences=$ARGV[$i+1];					
					}					
				elsif($ARGV[$i]=~/-T$/)
					{
						$no_of_threads=$ARGV[$i+1];					
					}
				elsif($ARGV[$i]=~/-hmm_seq_evalue_cutoff$/)
					{
						$hmm_seq_evalue_cutoff=$ARGV[$i+1];					
					}
				elsif($ARGV[$i]=~/-hmm_dom_evalue_cutoff$/)
					{
						$hmm_dom_evalue_cutoff=$ARGV[$i+1];					
					}			
								
	}

if(($tmp_dir!~/\S/ and $input_seq_file!~/\S/) or $all_gene_positions_file!~/\S/ or $predicted_cas_gene_position_and_sequences!~/\S/)
	{
		print "Error: script_find_cas_genes.pl : please provide the input and out put file in this manner: -tmp tmp_dir -f seq_file -o all_gene_positions_file -T no_of_threads \n\n"; exit;
	}


unless(-d $tmp_dir)
	{
		system("mkdir -p $tmp_dir");
	}

#-----------------------------------------------------
chomp $input_seq_file; $input_seq_file=~s/\r//;

#print "0: $input_seq_file\n";

my $only_file_name;
if($input_seq_file=~/\//)
	{		
		#print "1: $input_seq_file\n";
		my @arr_t1=split('\/',$input_seq_file);
		$only_file_name=$arr_t1[$#arr_t1];
	}
else{
		#print "2: $input_seq_file\n";
		$only_file_name=$input_seq_file;
	}	
	

#--- move in to the tmp_dir folder ---	
chdir "$tmp_dir"; 

#--- Step 1: run GeneMarkS to predict the Genes
#print "\tRunning: gmhmmp.pl --org $script_dir/GeneMarkS_MODS --mode combined --output $only_file_name.GFF --format GFF --faa $only_file_name \n";
system("gmhmmp.pl --org $script_dir/GeneMarkS_MODS --mode combined --output $only_file_name.GFF --format GFF --faa $only_file_name");



#---- Step2: run hmmscan[very slow] or hmmsearch [fast] against the makarova hmm models 
#print "\tRunning: hmmsearch --notextw -E $hmm_seq_evalue_cutoff --domE $hmm_dom_evalue_cutoff -o /dev/null --cpu $no_of_threads --tblout $only_file_name.GFF.tbl --domtblout $only_file_name.GFF.tbl.dom $script_dir/HMMER_FILES/all_cas_profiles.hmm $only_file_name.GFF.faa \n";
system("hmmsearch --notextw -E $hmm_seq_evalue_cutoff --domE $hmm_dom_evalue_cutoff -o /dev/null --cpu $no_of_threads --tblout $only_file_name.GFF.tbl --domtblout $only_file_name.GFF.tbl.dom $script_dir/HMMER_FILES/all_cas_profiles.hmm $only_file_name.GFF.faa ");	# -E 1e-5 




#---- Step3: Process output tables load the hmm seqID and predicted Cas gene names and CRISPR (sub)types from makarova_summary_table_of_cas_genes.tab
my %hash_of_makarova_hmm_seqID_and_cas_gene;
my %hash_of_makarova_hmm_seqID_and_crispr_type;

my @arr_table_content=`grep -v '#' $script_dir/HMMER_FILES/makarova_summary_table_of_cas_genes.tab >&1`;
foreach my $row(@arr_table_content)
	{
		chomp $row; $row=~s/\r//;
		
		my @arr_t1=split('\t',$row);
		
		my $seqID=$arr_t1[0];
		my $cas_gene=$arr_t1[3];
		my $crispr_type=$arr_t1[5]; $crispr_type=~s/CAS-//g; $crispr_type="Type-$crispr_type";
		
		$hash_of_makarova_hmm_seqID_and_cas_gene{$seqID}=$cas_gene;
		$hash_of_makarova_hmm_seqID_and_crispr_type{$seqID}=$crispr_type;
		
	}


#--- now load the gene sequences in a hash
my $gene_fasta_file="$only_file_name.GFF.faa";
my %hash_of_accession_geneID_and_sequences;
my %hash_of_accession_geneID_and_position;
my %hash_of_accession_position_and_seq;
&load_protein_sequences($gene_fasta_file,\%hash_of_accession_geneID_and_sequences,\%hash_of_accession_geneID_and_position,\%hash_of_accession_position_and_seq);



 



########### now get the matched domain position and sequences 
my %hash_of_accession_geneID_start_stop_highest_score;
my %hash_of_accession_geneID_start_stop_best_matching_domainID;
my %hash_of_accession_geneID_start_stop_best_matching_domain_seq;


my %hash_of_accession_geneID_start_stop_and_source_seq;
my %hash_of_accession_geneID_start_stop_and_nucleotide_position;

my %hash_of_accession_geneID_start_stop_and_source_seq_matching_e_value;
my %hash_of_accession_geneID_start_stop_and_domain_matching_e_value;


my @arr_hmmsearch_table_rows_domains=`grep -v '^#' $only_file_name.GFF.tbl.dom >&1`;
my $seq_count=0;
foreach my $row(@arr_hmmsearch_table_rows_domains)
	{
		chomp $row; $row=~s/\r//;
		$row=~s/>//g;
		
		$seq_count++;
		my @arr_t1=split('\s+',$row);  #gene_3017|GeneMark.hmm|305_aa|-|3140655|3141572 -          cd09634              -            1.3e-11   41.8   0.0   6.5e-08   29.6   0.0   2.1   2   0   0   2   2   2   2 >NC_017933
		
		#--- get the gene info ----		
		my @arr_t2=split('\|',$arr_t1[0]);
		my $gene_id=$arr_t2[0];		
		my $gene_start=$arr_t2[4];
		my $gene_stop=$arr_t2[5];		
		#------------------------
		
		my $accession=$arr_t1[$#arr_t1];
		
		my $makarova_domainID=$arr_t1[3];		
		my $domain_score=$arr_t1[13];
		my $aligned_domain_start=$arr_t1[17];
		my $aligned_domain_stop=$arr_t1[18];		
		my $dom_position=$aligned_domain_start."-".$aligned_domain_stop;
		
		
		my $full_sequence_e_value=$arr_t1[6];	$full_sequence_e_value=0 + $full_sequence_e_value;
		my $domain_c_value=$arr_t1[11];			$domain_c_value=0 + $domain_c_value; 	#$domain_c_value=sprintf("%.30f",$domain_c_value);
		my $domain_i_value=$arr_t1[12];			$domain_i_value=0 + $domain_i_value; 	#$domain_i_value=sprintf("%.30f",$domain_i_value);
		
		#print "$full_sequence_e_value\t$domain_c_value\t$domain_i_value\n";
		
		#---- get the Cas gene name and CRISPR type
		my $cas_gene=$hash_of_makarova_hmm_seqID_and_cas_gene{$makarova_domainID}; $cas_gene=ucfirst($cas_gene);
		my $crispr_type=$hash_of_makarova_hmm_seqID_and_crispr_type{$makarova_domainID};
		my $position=$gene_start."-".$gene_stop;
		my $gene_seq=$hash_of_accession_position_and_seq{$accession}{$position};
		
		
		
		#---- now get the aligned domain sequence ---
		my $aligned_domain_seq=substr($gene_seq,$aligned_domain_start,($aligned_domain_stop-$aligned_domain_start+1));
		
		#---- check if the reported region already exist; if exist, then check the score: only keep the best scoring domain
		if(defined $hash_of_accession_geneID_start_stop_highest_score{$accession}{$gene_id})
			{
				foreach my $existing_dom_position(keys %{$hash_of_accession_geneID_start_stop_highest_score{$accession}{$gene_id}})
					{
						my ($ex_start,$ex_stop)=split('-',$existing_dom_position);
						my $ex_dom_score=$hash_of_accession_geneID_start_stop_highest_score{$accession}{$gene_id}{$existing_dom_position};
						
						#--- my overlap found
						my $overlap_found=0;
						if($aligned_domain_start>=$ex_start and $aligned_domain_start <=$ex_stop ){$overlap_found++;}
						if($aligned_domain_stop>=$ex_start and $aligned_domain_stop <=$ex_stop ){$overlap_found++;}
						
						if($ex_start >= $aligned_domain_start and $ex_start <= $aligned_domain_stop ){$overlap_found++;}
						if($ex_stop >= $aligned_domain_start and $ex_stop <= $aligned_domain_stop ){$overlap_found++;}
						
						if($overlap_found >0 and $domain_score >$ex_dom_score)
							{
								#---- delete all previous record and push the new onces ----
								delete $hash_of_accession_geneID_start_stop_highest_score{$accession}{$gene_id}{$existing_dom_position};
								delete $hash_of_accession_geneID_start_stop_best_matching_domainID{$accession}{$gene_id}{$existing_dom_position};
								delete $hash_of_accession_geneID_start_stop_best_matching_domain_seq{$accession}{$gene_id}{$existing_dom_position};
								
								delete $hash_of_accession_geneID_start_stop_and_source_seq{$accession}{$gene_id}{$existing_dom_position};
								delete $hash_of_accession_geneID_start_stop_and_nucleotide_position{$accession}{$gene_id}{$existing_dom_position};
								
								delete $hash_of_accession_geneID_start_stop_and_source_seq_matching_e_value{$accession}{$gene_id}{$existing_dom_position};
								delete $hash_of_accession_geneID_start_stop_and_domain_matching_e_value{$accession}{$gene_id}{$existing_dom_position};
								
								
								#--- push the new ones
								$hash_of_accession_geneID_start_stop_highest_score{$accession}{$gene_id}{$dom_position}=$domain_score;
								$hash_of_accession_geneID_start_stop_best_matching_domainID{$accession}{$gene_id}{$dom_position}=$makarova_domainID;
								$hash_of_accession_geneID_start_stop_best_matching_domain_seq{$accession}{$gene_id}{$dom_position}=$aligned_domain_seq;
								
								$hash_of_accession_geneID_start_stop_and_source_seq{$accession}{$gene_id}{$dom_position}=$gene_seq;
								$hash_of_accession_geneID_start_stop_and_nucleotide_position{$accession}{$gene_id}{$dom_position}=$position;
								
								#---e-values ---
								$hash_of_accession_geneID_start_stop_and_source_seq_matching_e_value{$accession}{$gene_id}{$dom_position}=$full_sequence_e_value;
								$hash_of_accession_geneID_start_stop_and_domain_matching_e_value{$accession}{$gene_id}{$dom_position}=$domain_i_value;
							}
						
						
					}
			}
		else{
				$hash_of_accession_geneID_start_stop_highest_score{$accession}{$gene_id}{$dom_position}=$domain_score;
				$hash_of_accession_geneID_start_stop_best_matching_domainID{$accession}{$gene_id}{$dom_position}=$makarova_domainID;
				$hash_of_accession_geneID_start_stop_best_matching_domain_seq{$accession}{$gene_id}{$dom_position}=$aligned_domain_seq;
				
				$hash_of_accession_geneID_start_stop_and_source_seq{$accession}{$gene_id}{$dom_position}=$gene_seq;
				$hash_of_accession_geneID_start_stop_and_nucleotide_position{$accession}{$gene_id}{$dom_position}=$position;
				
				#---e-values ---
				$hash_of_accession_geneID_start_stop_and_source_seq_matching_e_value{$accession}{$gene_id}{$dom_position}=$full_sequence_e_value;
				$hash_of_accession_geneID_start_stop_and_domain_matching_e_value{$accession}{$gene_id}{$dom_position}=$domain_i_value;
			}	


			
	}
	



#--- now print the best scoring domain hits ---
my %hash_of_already_written_accession_position_and_cas_gene;
open(APP,">>$all_gene_positions_file");
open(WR,">>$predicted_cas_gene_position_and_sequences");
foreach my $accession(sort keys %hash_of_accession_geneID_start_stop_highest_score)
	{
		foreach my $gene_id(sort keys %{$hash_of_accession_geneID_start_stop_highest_score{$accession}})
			{				
				foreach my $dom_position(sort keys %{$hash_of_accession_geneID_start_stop_highest_score{$accession}{$gene_id}})
					{
						my($aligned_domain_start,$aligned_domain_stop)=split('-',$dom_position);
						my $aligned_domain_seq=$hash_of_accession_geneID_start_stop_and_source_seq{$accession}{$gene_id}{$dom_position};
						
						my $gene_seq=$hash_of_accession_geneID_start_stop_and_source_seq{$accession}{$gene_id}{$dom_position};	
						my $position=$hash_of_accession_geneID_start_stop_and_nucleotide_position{$accession}{$gene_id}{$dom_position};					
						my($gene_start,$gene_stop)=split('-',$position);
						
						my $makarova_domainID=$hash_of_accession_geneID_start_stop_best_matching_domainID{$accession}{$gene_id}{$dom_position};
						my $crispr_type=$hash_of_makarova_hmm_seqID_and_crispr_type{$makarova_domainID};
						my $cas_gene=$hash_of_makarova_hmm_seqID_and_cas_gene{$makarova_domainID}; $cas_gene=ucfirst($cas_gene);
						
						
						#---- the cas_protein_acc is applicable for NCBI .GBK file 
						my $cas_protein_acc="NA";
						   $cas_protein_acc=$makarova_domainID;
						
						#--- get the evalues ---
						my $full_sequence_e_value=$hash_of_accession_geneID_start_stop_and_source_seq_matching_e_value{$accession}{$gene_id}{$dom_position};
						my $domain_i_value=$hash_of_accession_geneID_start_stop_and_domain_matching_e_value{$accession}{$gene_id}{$dom_position};
						
						
						
						print WR "$accession\t$cas_gene\t$gene_start\t$gene_stop\t$full_sequence_e_value\t$gene_seq\t$makarova_domainID\t$aligned_domain_start\t$aligned_domain_stop\t$domain_i_value\t$aligned_domain_seq\t$crispr_type\n";
						
						if(not $hash_of_already_written_accession_position_and_cas_gene{$accession}{$position}{$cas_gene})	#-- many a Cas gene has many domains, this will stop it from printing duplicate lines
							{
								print APP "$accession\tCRISPR\t$gene_start\t$gene_stop\t$cas_gene\t$cas_protein_acc\t$crispr_type\t$full_sequence_e_value\n";
								
								$hash_of_already_written_accession_position_and_cas_gene{$accession}{$position}{$cas_gene}=1; 
							}
					}
			}		
	}
close(WR);
close(APP);













################ subs 

sub load_protein_sequences()
	{
		my($gene_fasta_file,$hash_of_accession_geneID_and_sequences,$hash_of_accession_geneID_and_position,$hash_of_accession_position_and_seq)=@_;
		
		
		
		#----------- now open the user fasta file and extract the individual sequences, and pass them to first pilercr, then to CRISPRDetect -----
		open(RD,"$gene_fasta_file") or print "$!: $gene_fasta_file not found\n";


		my $seq_index=0;
		
		my $last_accession;		
		my $last_gene_id;
		my $last_start;
		my $last_stop;
		my $last_position;
		my $seq="";
		my $o_id="";
		while( my $line=<RD> )
			{				
				chomp $line; $line=~ s/\r//; $line=~ s/^\s+//; 
				#print "$line\n";
				
				if($line=~/^>/)
					{
						#print "matched\n";
						if($seq_index>0)
							{
								$seq=uc($seq);
								
								#----- store the id and seq
								$hash_of_accession_geneID_and_sequences->{$last_accession}->{$last_gene_id}=$seq;
								$hash_of_accession_geneID_and_position->{$last_accession}->{$last_gene_id}=$last_start."-".$last_stop;
								$hash_of_accession_position_and_seq->{$last_accession}->{$last_position}=$seq;
								
								$last_accession="";
								$last_gene_id="";
								$seq="";
							}						
						
						#---- assign the 
						$o_id=$line;	chomp $o_id;$o_id=~ s/\r//; 
						
						#--- remove all the >
						$o_id=~ s/>//g;
							
						my @arr_t1=split(' ',$o_id);  #>gene_1|GeneMark.hmm|146_aa|-|364|804	>NC_017933
						$last_accession=$arr_t1[1];
						
						my @arr_t2=split('\|',$arr_t1[0]);
						$last_gene_id=$arr_t2[0];	
						
						
						$last_start=$arr_t2[4];
						$last_stop=$arr_t2[5];
						$last_position=$last_start."-".$last_stop;
						
						
						#---- check if NCBI accession style ---

												
						$seq_index++;
					}
				else{
						#---------- get the sequence ---
						#print "matched\n";
						
						$line=uc($line);
						$seq=$seq.$line;
					}																	 								 						 
				#print WR "$line\n";
												
			 }
		close(RD); 
		#------- now write the last sequence
		$seq=uc($seq);
		
		$hash_of_accession_geneID_and_sequences->{$last_accession}->{$last_gene_id}=$seq;
		$hash_of_accession_geneID_and_position->{$last_accession}->{$last_gene_id}=$last_start."-".$last_stop;	 
		$hash_of_accession_position_and_seq->{$last_accession}->{$last_position}=$seq;
			
		return 1;
	}


exit;
