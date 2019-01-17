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
my %Dir_hash;

#arrays
my @File_array_paired;
my @File_array_un1;
my @File_array_un2;

#strings


#integers


#declare the options string and associated global variables
#d:directory, -j string to identify joined pe files, -1 string to identify forward unjoined read files, -2 string to identify reverse unjoined read files
my $opt_string='hd:j:1:2:';
my $prg = basename $0;
my %opt;
my $Seq_directory;
my $joined_string;
my $for_string;
my $rev_string;


####process the options and error check input
getopts("$opt_string", \%opt) or usage();
$opt{h} and usage();
usage() unless ($opt{d});
print "Could not read directory: $opt{d}\n" and die(usage()) unless (-e $opt{d});

#set the input files and default parameters
$Seq_directory = $opt{d};
$joined_string= $opt{j} ? $opt{j} : "join";
$for_string= $opt{1} ? $opt{1} : "un1";
$rev_string= $opt{2} ? $opt{2} : "un2";



######inform the user of program settings
print STDERR "Creating backup of the directory: $Seq_directory\n";
print STDERR "All files within the directory $Seq_directory will be renamed with sample id and sorted into 3 directories.\n";
print STDERR "Joined paired-end read files will be identified by the string: $joined_string\n";
print STDERR "Forward unpaired read files will be identified by the string: $for_string\n";
print STDERR "Reverse unpaired read files will be identified by the string: $rev_string\n";


# open up the data directory
opendir(SEQDIR, $Seq_directory) or die "can't open directory: $Seq_directory $!\n";

# if input directory does not end with '/', add it.
#This is for the purposes of printing files later on, no need for the opening of the directory.
if ($Seq_directory !~ /\/$/){
        $Seq_directory .= "/";
}

###loop through all of the directories in the input directory and grab the file names
foreach my $d (readdir (SEQDIR)){
	my $temp_id;
	
	next if ($d =~ /^\./); # skip if this file is default '.' or '..' directory
	next unless (-d $Seq_directory.$d); # skip if this item is not a directory
	if ($d !~ /\/$/){
       		$d .= "/";
	}
	#enter the sample id, as identified from the directory name, into a hash with the directory name as key.
	if($d =~ m/^([A-Za-z0-9]+)-/){
		$Dir_hash{$d}=$1;
	}else{
		print STDERR "Warning! Directory name \'$d\' in directory $Seq_directory does not match naming conventions, skipping this directory.\n";
	}

}

#make a backup of the original directory
system("cp", "-r", $Seq_directory, substr($Seq_directory, 0, -1) . ".bak");


#Initialize the output directories
my $Output_dir1= $Seq_directory . "joined_fastq/";
my $Output_dir2= $Seq_directory . "un1_fastq/";
my $Output_dir3= $Seq_directory . "un2_fastq/";
system("mkdir", $Output_dir1);# or die("Could not create output directory $Output_dir1 : $!\n");
system("mkdir", $Output_dir2);#  or die("Could not create output directory $Output_dir2 : $!\n");
system("mkdir", $Output_dir3);#  or die("Could not create output directory $Output_dir3 : $!\n");



while ((my $key, my $value) = each %Dir_hash){
	#move the files in each folder
	my $move_count = 0;
	print "Could not read directory: $Seq_directory$key\n" and die(usage()) unless (-d $Seq_directory.$key);
	opendir(TEMPDIR, $Seq_directory . $key) or die( "can't open directory $key: $!\n");
	foreach my $f (readdir (TEMPDIR)){
		next if ($f =~ /^\./); # skip if this file is default '.' or '..' directory
		next unless ($f =~ /\.fastq$/); # skip if this file is not in fastq format
		
		#move the joined file
		if($f =~ /$joined_string\.fastq/){
			system("mv", $Seq_directory . $key . $f, $Output_dir1 . $value . "_" . $f);
			$move_count++;
		}
		#move the forward unjoined read file
		if($f =~ /$for_string\.fastq$/){
			system("mv", $Seq_directory . $key . $f, $Output_dir2 . $value . "_" . $f);
			$move_count++;
		}
		#move the joined file
		if($f =~ /$rev_string\.fastq$/){
			system("mv", $Seq_directory . $key . $f, $Output_dir3 . $value . "_" . $f);
			$move_count++;
		}
	}	
	#print a warning if 3 files were not moved
	print STDERR "Warning! Directory with ID \'$value\' does not contain 3 files for transfer" unless ($move_count == 3);
	#remove the old folder, which is now empty
	system("rmdir", $Seq_directory . $key);

}#End for each of %Dir_hash



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

	Name $prg - This Perl script will move and rename to include sample name from all files with the strings 
			"join", "un1", and "un2" from subdirectories of an arguement directory. A backup of the 
			original directory is made as directory.bak. Three directories are made to contain the 
			three file types and are contained within the original arguement directory.
			
	usage: 	$prg -d directory [-j join-string -1 forward-string -2 reverse-string]

	-h	      	:	print this help message
	-d directory   	:	Path to a directory which contains paired end fastq sequence files in .gz format
	-j join-string	:	String used to identify joined read files (default: "join")
	-1 for-string	:	String used to identify unjoined forward read files (default: "un1")
	-2 rev-string	:	String used to identify unjoined forward read files (default: "un2")

	ex: $prg -d Some_directory/ 
		This will move and rename all files in "Some_directory" using the default setting above.

EOF
	exit;

}




