#!/usr/bin/perl -w

#use warnings;
use strict;
use Getopt::Std;
use File::Basename;

sub isint($);
sub reverse_complement($);
sub trim($);
sub usage;


#declare the global variables

#arrays
my @Row_holder;
my @Sample_counts;

#strings
my $Csv_file;
my $Split_char;
#integers
my $Threshold;
my $OTU_count=0;
my $Remove_count=0;
my $Line_count=1;

#declare the options string
my $opt_string='ht:c:s:';
my $prg = basename $0;
my %opt;
my $alphabet;

####process the options and set global variables for options
getopts("$opt_string", \%opt) or usage();
$opt{h} and usage();
usage() unless ($opt{c});
#print "Invalid character separator: $opt{s}\n" and die(usage()) unless (length($opt{s})==1);
print "No such file exists: $opt{c}\n" and die(usage()) unless (-e $opt{c});
if($opt{t}){print "Threshold must be represented as a decimal. Invalid threshold specified: $opt{t}\n" and die(usage()) unless ($opt{t} =~ /^0\.\d+$/);}


######inform the user of program settings
print STDERR "Filtering all reads in .csv spreadsheet: $opt{c} at the threshold of $opt{t} by sample\n";

#set the input files

$Csv_file = $opt{c};
$Threshold = $opt{t} ? $opt{t} : 0.0005;
$Split_char = $opt{s} ? $opt{s} : "\t";

my @basename=split(/\./, basename $Csv_file);
my $outfile = $basename[0] . "_FilteredbySample_" .$Threshold. ".csv";

open(my $CSV_HANDLE, "<:encoding(UTF-8)", $Csv_file) or die("Could not open .csv file $Csv_file: $!\n");
#:encoding(UTF-8)
my @active_header;

#open output file for printing
open(my $OUTFILE, ">:encoding(UTF-8)", $outfile) or die("Could not open output file $outfile: $!\n");


###loop through all of the rows of the .csv and gather information about the samples
###print headers to the new output file
my $header_flag=0;
while(<$CSV_HANDLE>){
	#always assume 2 header lines, skip them for now
	chomp();
	
	if($Line_count<=2){
		$Line_count++;
		$header_flag=1;
		next; 
	}
	#print STDERR "$_\n" and die();
	$Row_holder[$OTU_count]=$_;
	$OTU_count++;
	$Line_count++;
}

print STDERR "A total of $OTU_count OTUs were recognized in the csv file $Csv_file.\n";

#now calculate the totals for each sample

#initialize the $Sample_counts[0] index
$Sample_counts[0]=0;
foreach(@Row_holder){
	my @temp_array = split(/$Split_char/ , $_);
	my $col_num=@temp_array;
	#print STDERR "number of columns is $col_num\n";
	#loop through all the reads in this OTU and add the value to each sample read count
	for(my $col_index=1 ; $col_index < ($col_num-1) ; $col_index++){
		
		if($Sample_counts[$col_index]){
			$Sample_counts[$col_index]+=$temp_array[$col_index];
			#print STDERR "$Sample_counts[$col_index]\n";
			#die();
		}else{
			$Sample_counts[$col_index]=$temp_array[$col_index];
		}
		
	}
}

#now finally loop over the values in the array, and print them to them to the .csv file
#reset the line count and file pointer
seek($CSV_HANDLE,0,0);
$Line_count=1;
while(<$CSV_HANDLE>){
	#always assume 2 header lines, print them
	if($Line_count<=2){
		printf $OUTFILE $_;
		$Line_count++;
		$header_flag=1;
		next; 
	}
	chomp($_);
	my @temp_array = split(/$Split_char/ , $_);
	my $col_num=@temp_array;
	#print STDERR $col_num . "\n" and die();
	#test each read count against the threshold percent of Sample total counts
	#set to zero if reads are removed
	for(my $col_index=1 ; $col_index< ($col_num-1) ; $col_index++){
		#print STDERR $Threshold*$Sample_counts[$col_index] . "\n" and die();
		if($temp_array[$col_index] !=0 && $temp_array[$col_index] < ($Threshold*$Sample_counts[$col_index])){
			$temp_array[$col_index]=0;
			$Remove_count++;
		}
	}

	printf $OUTFILE join("$Split_char", @temp_array) . "\n";
	$Line_count++;
}

print STDERR "A total of ". (scalar @Sample_counts - 1) ." samples were recognized from file $Csv_file\n";
print STDERR "A total of " . ($Line_count-3) . " OTUs were recognized and printed from the $Csv_file\n";
print STDERR "Removed a total of $Remove_count OTU occurences from samples where their total read count was less than $Threshold of total reads in that sample\n";

#################Begin PERL subroutine definitions#######################################
##########################################################################################3


sub isint($){
	my $val= shift;
	return ($val =~ m/^\d+/)
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

# PERL usage function for this program, MES_fasta_insert.pl
sub usage
{
	print STDERR << "EOF";

	Name $prg - Takes BIOM table in .csv format as input and filters out reads from an OTU for a sample
		when the total number of reads in that OTU represent less than a threshold percent value
		specified by the user.


	usage: 	$prg [-t threshold] -c input_csv.csv

	-h	      	:	print this help message
	-t threshold   	:	threshold value (represented as a fraction) at which
				to remove reads from a sample (default: 0.0005)
	-s split char	:	The character that divides fields of the input text file (default: tab (\t))
	-c input_csv	:	the input .csv file (sets inputfile)


	ex: $prg -s \t -c input_csv.csv -t 0.003

EOF
	exit;

}




