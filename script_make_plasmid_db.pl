#!/usr/bin/perl 


use strict;
use Term::ANSIColor;

use FindBin qw($Bin);
use lib "$Bin/lib";


#use lib '/PROJECTS/custom_scripts/lib1';
use Parallel::ForkManager;
my $pm = new Parallel::ForkManager(5);





#-----

my $ftp="ftp.ncbi.nlm.nih.gov/refseq/release/plasmid"; 



my $tmp_directory="/tmp/REFSEQ_PLASMIDS";
unless(-d "$tmp_directory")
	{
		system("mkdir -p $tmp_directory");
	}


my $output_directory="$Bin/REFSEQ_PLASMIDS_DB";
unless(-d "$output_directory")
	{
		system("mkdir -p $output_directory");
	}



my @organisms;
push(@organisms,$ftp);



foreach my $ftp_dir(@organisms)
	{
		my @files_to_download;
		system("rm ./index.htm*");
		 
		system("wget --no-remove-listing $ftp_dir\/");
		
		open(RD,"index.html");
		my @arr_rd=<RD>;

		foreach my $line(@arr_rd)
			{
				chomp $line; $line=~s/\r//g;
				
				
																							################# change the following line for redundant protein
				if($line=~/genomic\.fna\.gz/ )
					{				
						
						my @arr_t1=split('>|<', $line);
						my $file=$arr_t1[2];
						
						
						push(@files_to_download,$file);							
						print "\tTo be downloaded: $file\n";
						
					}
			}
		close(RD);
		
		
		foreach my $file(@files_to_download)
			{
				chomp $file; $file=~s/\r//; $file=~s/\s+//g;
				
				my $file_size=-s "$tmp_directory/$file";
				#--- skip the files already downloaded ---------------------------------
				if(-e "$tmp_directory/$file" and $file_size >10)
					{
						print "File exist. Skipping..\n";next;
					}
				#-----------------------------------------------------------------------
				
				select(undef,undef,undef,0.5);
				
				$pm->start and next;
				
				print "\nDownloading file $file..\n";
				unless(-e "$tmp_directory/$file")
					{
						system("wget -N -q $ftp_dir\/$file -O $tmp_directory\/$file");
					}
				
				
				
				my $output_file=$file;
				$output_file=~s/\.gz$//;
				
				print "\nExtracting $file..";
				if(not -e "$tmp_directory/$output_file" and $file=~/\.fna\.gz$/)
					{
						system("gunzip -c $tmp_directory\/$file >$tmp_directory\/$output_file");
					}
				
				
				print "\tdone.\n";
				
				$pm->finish;
			}
		$pm->wait_all_children;	
		
		
		#--- now cat the sequences in plasmids.fa file
		system("cat $tmp_directory/plasmid.*.fna >$output_directory/plasmids.fa");
		
		
		#--- make a blastdb
		system("makeblastdb -in $output_directory/plasmids.fa -dbtype nucl");
	
	}



