#!/usr/bin/perl -w

#use warnings;
use strict;
use Getopt::Std;
use File::Basename;
#use Bio::SeqIO;
#use Bio::DB::Fasta;


sub reverse_complement($);
sub trim($);
sub usage;
sub isint($);

#declare the global variables
#hashes


#arrays
my @File_array;

#strings


#integers


#declare the options string and associated global variables
#d:directory,p: path to trimmomatic jar n:number-of-threads, l:leading, t:trailing, s:sliding_window, m:minimum_length c:headcut amount
my $opt_string='hd:p:n:l:t:s:m:c:';
my $prg = basename $0;
my %opt;
my $Seq_directory;
my $Trim_path;
my $Threads;
my $Lead;
my $Trail;
my $Sliding_window;
my $Min_length;
my $Headcrop;

####process the options and error check input
getopts("$opt_string", \%opt) or usage();
$opt{h} and usage();
usage() unless ($opt{d} && $opt{p});
print "Could not read directory: $opt{d}\n" and die(usage()) unless (-e $opt{d});
print "File does not exist: $opt{p}\n" and die(usage()) unless (-e $opt{p});
if($opt{n}){print "Number of threads must be an integer\n" and die(usage()) unless(isint($opt{n}));}
if($opt{c}){print "Headclip argument must be an integer\n" and die(usage()) unless(isint($opt{c}));}
if($opt{l}){print "Leading trim quality score minimum must be an integer\n" and die(usage()) unless(isint($opt{l}));}
if($opt{t}){print "Trailing trim quality score minimum must be an integer\n" and die(usage()) unless(isint($opt{t}));}
if($opt{s}){print "Sliding window parameter must follow the format (Integer):(Integer)" and die(usage()) unless($opt{s} =~ m/\d+:\d+/);}
if($opt{m}){print "Minimum sequence length must be an integer\n" and die(usage()) unless(isint($opt{m}));}


#set the input files and default parameters
$Seq_directory = $opt{d};
$Trim_path= $opt{p};
$Threads = $opt{n} ? $opt{n} : 1;
$Lead = $opt{l} ? $opt{l} : 20;
$Trail = $opt{t} ? $opt{t} : 20;
$Sliding_window = $opt{s} ? $opt{s} : "4:20";
$Min_length = $opt{m} ? $opt{m} : 200;
$Headcrop = $opt{c} ? $opt{c} : 10;

######inform the user of program settings
print STDERR "Creating trimmomatic command script for all .gz fastq files stored in the directory: $opt{d}\n";
print STDERR "Path to trimmomatic .jar file (including filename): $Trim_path\n";
print STDERR "Threads for each analysis: $Threads\n";
print STDERR "Number of bases at the head of the sequence: $Headcrop\n";
print STDERR "Minimum quality score for bases at the leading edge of the sequence: $Lead\n";
print STDERR "Minimum quality score for bases at the trailing edge of the sequence: $Trail\n";
print STDERR "Sliding_window parameter set to: $Sliding_window\n";
print STDERR "Minimum sequence length to be retained after trimming: $Min_length\n";



# open up BLAST and ASSEMBLY directories and reference file 
opendir(SEQDIR, $Seq_directory) or die "can't open directory: $Seq_directory $!\n";

# if input directory does not end with '/', add it.
#This is for the purposes of printing files later on, no need for the opening of the directory.
if ($Seq_directory !~ /\/$/){
        $Seq_directory .= "/";
}

###loop through all of the files in the input directory and grab the file names
foreach my $f (readdir (SEQDIR)){
	next if ($f =~ /^\./); # skip if this file is default '.' or '..' directory
	next if ($f !~ /\.gz$/); # skip if this file does not end with the right file suffix
	push(@File_array,$f);

}

#now sort the array of filenames
@File_array=sort @File_array;

#Initialize the output file (script) and print the header line to the file
my $Output_script= "my_Make_Trimmomatic.sh";
open(OUTPUT, ">$Output_script");
print OUTPUT "#!/bin/bash\n";
print OUTPUT "mkdir my_Trimmomatic_OUTPUT\n";
print OUTPUT "PID=\"\"\n";
if(@File_array % 2 !=0){print STDERR "Warning, files not paired correctly (there is an odd number of read files in the directory).\nFiles will be paired and commands generated where possible.";};
while(@File_array){
	#load in a pair of files from the array.
	my $f1 = shift @File_array;
	my $f2 = shift @File_array;
	my @file_array1=split(/\_/,$f1);	
	my @file_array2=split(/\_/,$f2);
	my $id1 = $file_array1[0]; 
	my $id2 = $file_array2[0];
	#make sure the files are properly paired
	if($id1 ne $id2){
		print STDERR "Warning: $id1 from file $f1 does not equal $id2 from file $f2, skipping these files\n";
		next;
	}
	#make sure the files are in the correct order
	if($f1 !~ m/_R1_/){
		print STDERR "Warning: Files $f1 and $f2 are not properly paired.\n Skipping these files\n";
		next;
	}
	if($f2 !~ m/_R2_/){
		print STDERR "Warning: Files $f1 and $f2 are not properly paired.\n Skipping these files\n";
		next;
	}
	
	my $infile1  = $Seq_directory.$f1;
	my $infile2  = $Seq_directory.$f2;
	print OUTPUT "java -jar $Trim_path PE -threads $Threads $infile1 $infile2 my_Trimmomatic_OUTPUT/$id1"."_R1_paired.fastq.gz my_Trimmomatic_OUTPUT/$id1"."_R1_unpaired.fastq.gz my_Trimmomatic_OUTPUT/$id2"."_R2_paired.fastq.gz my_Trimmomatic_OUTPUT/$id2"."_R2_unpaired.fastq.gz HEADCROP:$Headcrop LEADING:$Lead TRAILING:$Trail SLIDINGWINDOW:$Sliding_window MINLEN:$Min_length &\n";
	print OUTPUT "PID=\$!\n";
	print OUTPUT "wait \$PID\n";

}

close(OUTPUT);

#go through all of the files in the array, pair them, and then print the appropriate lines

#print STDERR "Printed a total of $Print_count sequences for $Sample_count Samples into the file $outfile\n";
#print STDERR "Printed list of proteins in the order that they were added to alignment into file $outfile2\n ";

#################Begin PERL subroutine definitions#######################################
##########################################################################################3

#subroutine to check if a value is an integer
sub isint($){
  my $val = shift;
  return ($val =~ m/^\d+$/);
}

#subroutine to reverse complement a string of DNA sequence
sub reverse_complement($) {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/NACGTacgt/NTGCAtgca/;
        return $revcomp;
}


# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

# PERL usage function for this program, SNP_density_calc.pl
sub usage
{
	print STDERR << "EOF";

	Name $prg - This Perl script will print a .sh script that can be excuted from the command line
			to process a folder full of gzipped (.gz suffix) paired-end .fastq read files 
			with Trimmomatic. As argument this script takes the directory of gzipped fasta 
			files, the path to the Trimmomatic software and several parameters to use for 
			running Trimmomatic. Naming conventions are expected in the the fastq files 
			in order to process the folder correctly.
			Files must be named as follows:
			Unique-sample-ID_any-text_R(1|2)_any-text.fastq.gz
			The placement of "_" characters is significant. Everything before the first "_" character
			is considered as the unique sample ID. the rest of the filename can follow any convention
			but the string "_R1_" or "_R2_" must be present. This string indentifies members of paired
			read files where R1 = (forward reads) and R2 = (reverse reads). Read files should be placed
			as pairs into a single directory and read files with the same ID, and the R1 and R2 strings
			will be paired in commands generated for the .sh script
			
			NOTE: file names of fastq files must end in .fastq.gz
			
	usage: 	$prg -d directory -p path [-c clip_bases -n threads -l leading_min_score -t trailing_min_score -s sliding-window -m min-length]

	-h	      	:	print this help message
	-d directory   	:	Path to a directory which contains paired end fastq sequence files in .gz format
	-p path		:	Full path to the Trimmomatic software
	-n threads	:	The number of threads to execute under trimmomatic commands (default: 1)
	-c clip_bases	:	The number of bases to trim off the front of the read before proceeding with any other 
				trimming command (default: 10).
	-l min_score	:	The minimum quality score to keep when cutting from the leading edge of 
				a sequence (default: 20)
	-t min_score	:	The minimum quality score to keep when cutting from the trailing edge of 
				a sequence (default: 20)
	-s sliding_window :	The parameters for sliding window trimming (default: 4:20)
	-m min_length	:	The minimum length allowable when keeping a trimmed sequence. (default: 200)

	ex: $prg -d Some_directory/ -p path_to_trimmomatic
		This will pair all files in "Some_directory" into Trimmomatic commands using the default setting above.

EOF
	exit;

}




