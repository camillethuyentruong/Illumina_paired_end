# Pipeline


# Trimmomatic: denoising
# Perl script to automate the Trimmomatic command for all files in a folder
# Make the script executable
chmod +x my_Make_Trimmomatic.sh
./my_Make_Trimmomatic.sh
# Execute the script
make_Trimmomatic.pl -d /path_to_Illmina_read_folder -p/path_to_Trimmomatic_application


# move paired reads, singe forward reads R1 (reverse read was descarded) and single reverse reads R2 (forward read was descarded) into separate folders
mv *paired.fastq ./paired
mv *R1_unpaired.fastq ./unpairedR1
mv *R2_unpaired.fastq ./unpairedR2


# QIIME: join reads from paired folder
multiple_join_paired_ends.py -i /path_to_paired_read_folder -o join_paired_ends -p join_paired_ends_parameters.txt


# Perl script to reorganize the data from the output directories of join_paired_ends.py in Qiime, moving joined and unjoined forward reads to separate folders
# Make the script executable
chmod +x make_separate_forSplitlibraries.sh
./make_separate_forSplitlibraries.sh
# Execute the script
make_separate_forSplitlibraries.pl -d /path_to_join_paired_ends_folder


# QIIME: convert to fasta format
#joined reads from joined_paired_ends
multiple_split_libraries_fastq.py -i /path_to_folder_joined_fastq  -o split_libraries_join --demultiplexing_method sampleid_by_file -p split_libraries_parameters.txt 

# unjoined forwards reads from joined_paired_ends
multiple_split_libraries_fastq.py -i /path_to_folder_un1_fastq  -o split_libraries_un1 --demultiplexing_method sampleid_by_file -p split_libraries_parameters.txt --sampleid_indicator .fastq 

# single forwards reads R1 from Trimmomatic
multiple_split_libraries_fastq.py -i /path_to_folder_unpairedR1 -o split_libraries_singleR1 --demultiplexing_method sampleid_by_file -p split_libraries_parameters.txt --sampleid_indicator _unpaired 

# single reverse reads R2 from Trimmomatic
multiple_split_libraries_fastq.py -i /path_to_folder_unpairedR2 -o split_libraries_singleR2 --demultiplexing_method sampleid_by_file -p split_libraries_parameters.txt --sampleid_indicator _unpaired 


# QIIME: eliminate reads of size < 150 bp
# for each read type fasta file
filter_short_reads.py input_file.fasta > output_file.fasta


# QIIME: eliminate potential chimeric sequences
# for each read type fasta file
identify_chimeric_seqs.py -i input_file.fasta -o output_file.fasta -m usearch61 --non_chimeras_retention intersection -r /path_to_UNITE_dynamic_ITS_reference_database.fna
filter_fasta.py -f input_file.fasta -s output_file.fasta/chimeras.txt -o output_file.nochimera.fasta -n


# ITSx: filter non-ITS and non-fungal reads
# for each read type fasta file
itsx -i /path_to_output_file.nochimera.fasta -o ITSx_output_file -t F --preserve T --save_regions ITS1 --reset T --graphical


# QIIME: create OTU table by merging all read type fasta files
pick_open_reference_otus.py -m -m sortmerna_sumaclust -i ITSx_output_file_join,ITSx_output_file_un1,ITSx_output_file_singleR1,ITSx_output_file_singleR2 -o output_folder/ -p pick_otu_parameters.txt -s 0.5 --min_otu_size 1 --suppress_align_and_tree --new_ref_set_id name_of_new_sequence_reference_database  -r /path_to_custom/UNITE_sequence_reference_database.fna


# Perl script to filter OTUs from an OTU table (in csv format) based upon a percent threshold of reads per sample [default = 0.05%]
# Make the script executable
chmod +x my_BIOM_filter_by_sample.sh
./my_BIOM_filter_by_sample.sh
# Execute the script 
BIOM_filter_by_sample.pl -t 0.0005 -c otu_table.csv
