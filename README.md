CRISPRDetect Version 3.0 help:
-----------------------------
CRISPREDetect3 is a perl program developed and tested in Fedora 21 Linux operating system. <b>CRISPRDetect3</b> should run under any unix based 
operating system that has a working 'perl' executable [comes with default installations under all *nix based operating systems]. If all the
3rd party dependencies are installed and available in user/system $PATH CRISPRDetect3 should run without any issues. 

The cd-hit-est program installation in Mac operating system were reported to have issues. If you face issues installing cd-hit-est, please 
refer to the link: https://github.com/weizhongli/cdhit/issues/24   


INSTALLATION:
------------
Please make sure that the following 3rd party tools are installed in your system. The only CPAN perl package 'Parallel' is provided in the 'lib' folder 
and should work. However, if it doesn't then either execute "cpan Parallel::ForkManager" to install it or execute the following commands after making sure 
you have its dependencies (POSIX, Storable, File::Spec, File::Temp, File::Path 2.00 and Test::More 0.81_01) installed in your system :

	wget http://search.cpan.org/CPAN/authors/id/Y/YA/YANICK/Parallel-ForkManager-1.19.tar.gz
	tar xvzf Parallel-ForkManager-1.19.tar.gz
	cd Parallel-ForkManager-1.19
	perl Makefile.PL && make test && make install



CRISPRDetect3 dependencies: 
--------------------------
The following dependencies are needed by CRISPRDetect3 to run. 

	clustalw 	Download from ftp://ftp.ebi.ac.uk/pub/software/clustalw2/2.1/ 
	water 		Comes with EMBOSS:6.3.1+ tools : 	Download from ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz
	seqret 		Comes with EMBOSS:6.3.1+ tools : 	Download from ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz
	RNAfold 	Comes with Vienna RNA package  :	Download the correct version specific for your operating system from http://www.tbi.univie.ac.at/RNA/#download
	cd-hit-est 	Comes with cdhit package  	   :	Download from http://weizhongli-lab.org/cd-hit/download.php  	
	blastn 		Comes with ncbi-blast+ package :	Download from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
	makeblastdb 	Comes with ncbi-blast+ package : 	Download from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

Other dependencies [Optional]:	
-----------------------------

	gmhmmp.pl 	Developed by Georgia Institute of Technology :	Download from http://exon.gatech.edu/Genemark/license_download.cgi	
	hmmsearch 	Developed by Howard Hughes Medical Institute :	Download from http://eddylab.org/software/hmmer/hmmer-3.3.1.tar.gz 


CRISPRDetect3 requires clustalw2 to be available in the $PATH as clustalw. You can do that by <b>"cp /PATH_TO_CLUSTALW2/clustalw2 /bin/clustalw"</b> or by creating 
a symbolic link <b>"ln -s /PATH_TO_CLUSTALW2/clustalw2 /bin/clustalw"</b> where PATH_TO_CLUSTALW2 is the full path of the directory where the executable clustalw2 is
present. 



Once all the dependencies are installed, please check that they are successfully installed and available in the user/system PATH by typing the following:

	clustalw -help
	water -help
	seqret -help
	RNAfold -help
	cd-hit-est -help
	blastn -help

Optional dependencies check [Only required if using '-annotate_cas_genes 1' ]:		
	hmmsearch -h
	gmhmmp.pl



NOTE:   
----	
The 'tmp' folder in the CRISPRDetect3 installation directory should have read and write permissions. An easy way to do that is by issuing the command 'chmod -R 755 . && chmod 777 tmp' from the CRISPRDetect3 installation directory. 


Also make sure that CRISPRDetect3 installation directory is in your system/user path. You can do that by adding the full path of the CRISPRDetect_3.0 installation directory to ./bashrc file (see example below):
export PATH="/FULL_PATH_TO_CRISPRDETECT_3.0:\$PATH"


The '-annotate_cas_genes 1' requires GeneMarkS v.4.30 package. For this to work, you will require to install GeneMarkS and obtain a License from http://exon.gatech.edu/Genemark/license_download.cgi . 
Download the gm_key_64.gz file, uncompress it (e.g. gunzip gm_key_64.gz ) and copy the license it in your home directory (e.g. cp gm_key_64 ~/.gm_key).


The '-wgs 1' parameter will instruct CRISPRDetect3 to use any/all predicted cas genes (applicable for multiFASTA input files) for all the predicted CRISPRs (suitable for WGS sequences,
where all the contigs are from the same genome). If the sequences are from different genome/organism, please specify -wgs 0 [i.e. Cas genes present in the same contig will be used for the CRISPRs scoring]. 


The '-check_plasmid 1' option compares the input sequence(s) against RefSeq plasmid sequences and if significant similarity is found between them, the relevant information is shown in the output.  



Basic syntax:	
---------------------

	CRISPRDetect3 -g NZ_CP006019.gbk -o NZ_CP006019_CRISPRDetect.txt > NC_003106_CRISPRDetect.log

The above command runs CRISPRDetect3 with default paramaeter on a complete gbk file (file containing both annotation and sequence) that has cas1 or cas2 annotated.
	
	CRISPRDetect3 -f test_multifasta.fa -o test_CRISPRDetect3 -check_direction 0 -array_quality_score_cutoff 3 -T 0 > test.log
        
The above command runs CRISPRDetect3 with a lower score cutoff on a fasta file (i.e. cutoff 3 rather than 4, as cas1 and cas2 are not annotated. GenerallY appropriate for contigs/fasta files), -T 0 will use all CPUs rather than the default of 4 and will not check CRISPR direction (not recommended) ]
	
	CRISPRDetect3 -f test_multifasta.fa -o test_CRISPRDetect3 -T 0 -wgs 1 -annotate_cas_genes 1 -check_plasmid 1 > test.log
        
The above command will instruct CRISPRDetect3 to identify Cas genes [based on Cas profiles published by Makarova et. al ] and use them for scoring of all the CRISPRs, as well as to check if the input sequence belongs to a potential plasmid sequence  [ based on plasmid sequences available in RefSeq DB ]




Compulsory input parameters for CRISPRDetect3:
---------------------------------------------
 
	-f/-g		FASTA/Genbank	Specify a FASTA formatted [e.g. -f test.fa] or Genbank formatted file [e.g. -g NZ_CP006019.gbk] containing the sequence. 		
	
	Note: the default cutoff of 4 is appropriate for genbank files that have cas1 or cas2 annoated, 3 is more appropriate for fa.
	
	-o		TEXT		Specify a text file that will contain the output [e.g. -o NC_003106_CRISPRDetect] 			
						Note: CRISPRDetect3 will provide two additional output files - one containing the filtered out arrays (e.g. NC_003106_CRISPRDetect.fp) 
						and a gff annotation file (e.g. NC_003106_CRISPRDetect.gff)
						
Basic options:
-------------						
	-h/-help	HELP		shows this help text
	-q/-quiet	0 or 1		Switch off/on step by step progress reporting [default 0]	
	-T			Threads		Specify number of parallel processes CRISPRDetect3 should use  [default 4; specify '-T 0' to use all processors]		
	-tmp_dir	tmp/		This is the default directory where temporary files generated by CRISPRDetect3 and its dependencies will be stored.	


Parameters for putative CRISPR identification [optional]:
--------------------------------------------------------	
	-word_length			11	This is the default word length CRISPRDetect3 uses to find the putative CRISPRs. Any positive integer >=6 can be used.
	-minimum_word_repeatation	3	By default CRISPRDetect3 uses 3 repeating identical words to find putative CRISPRs. To find CRISPRs with 2 repeats, use -minimum_word_repeatation 2	
	-max_gap_between_crisprs	125	By default the maximum gap is set tp 125 nucleotides between the repeating identical seed words.
	-repeat_length_cutoff		17	After the intial processing, putative CRISPRs with repeat lengths less than this value will be rejected.


Filtering parameters [optional]:
-------------------------------	
	-minimum_repeat_length		23	Minimum length of repeats 
	-minimum_no_of_repeats		3	Predicted CRISPRs with number of repeats less than this value will be excluded. To include CRISPRs with only 2 repeats, use -minimum_no_of_repeats 2
	-array_quality_score_cutoff	4	Predicted CRISPRs with score less than this value will be excluded from the output file. 
									The CRISPRs with score >=0 and less than the specied value will be moved to the output.fp file [output refers to user given output filename]. Cutoff of 3 is more appropriate for fasta files.
						

Additional parameters [optional]:
--------------------------------
	-left_flank_length		500	This is the default length of the 5' (upstream) region of the CRISPRs.
	-right_flank_length		500	This is the default length of the 3' (downstream) region of the CRISPRs.		
	


Advanced options [optional]:
---------------------------	
To test different methods as specified in the literature, open the CRISPRDetect.pl program with any text editor [e.g. gedit in RHEL/Fedora/CentOS, or vi in any *nix OS, or 
notepad in Windows OS] and change the parameters in the top most section of the script. To toggle individual methods, locate the '$check_' prefix and change the value to 1 
(i.e. the method will be applied) or 0 (i.e. the method will not be applied). 
	
	Examples:
		
		Direction specific options:
		--------------------------
			$check_direction=0;			[ Default is 1, making it 0 will turn the method off.] 
				
		To change the parameter(s) of a particular method (such as check_array_degeneracy) change the nested variables under that particular method.
		
			$check_array_degeneracy=1;	 
				$array_degeneracy_score=0.41; 		[ Default: PPV (0.91) - 0.50 ]
				$permitted_mutation_per_array=0; 	[ Default 0 ]
	
	Changing to '$permitted_mutation_per_array=2;' will instruct the program to allow maximum 2 bases as permitted mutations per CRISPR array.



Newly added options [for CRISPRDetect3 or higher]:
-------------------------------------------------
	-annotate_cas_genes		1/0	 By default this option is turned off. Use 1 to enable this option [Requires GeneMarkS package installed and gmhmmp.pl file in the PATH].
	-hmm_seq_evalue_cutoff	 0.01	 Specify sequence e-value cutoff for hmmsearch. Default value is 0.00001
	-hmm_dom_evalue_cutoff	 0.01	 Specify domain e-value cutoff for hmmsearch. Default value is 0.00001
	-wgs 	1/0	 Specify if all sequences in the input sequence file belongs to the same genome [useful with option '-annotate_cas_genes 1']. Default value is 0 [Zero]. 
	-check_plasmid 	1/0	 By default this option is turned off. Use 1 to enable this option [Requires RefSeq Plasmid sequence DB. You can generate the required DB by running script_make_plasmid_db.pl from the CRISPRDetect_3.0 installation directory.
	



If you use CRISPRDetect3, then please cite:
-----------------------------------------
	Biswas, A., Staals, R. H., Morales, S. E., Fineran, P. C. & Brown, C. M. CRISPRDetect: a flexible algorithm to define CRISPR arrays. BMC Genomics 17, 356 (2016).

For version updates and bug fixes refer to https://github.com/ambarishbiswas/CRISPRDetect_3.0 or email ambarishbiswas[at]gmail[dot]com  

 	
