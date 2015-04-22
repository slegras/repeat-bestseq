#! /usr/bin/env perl

BEGIN {
	use File::Basename;
	push( @INC, dirname(__FILE__) );
}

use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use lib qw(.);
use Switch;
use findBestSequence;
use POSIX qw(strftime);
use File::Basename;

use config;
use vars qw(%config );
*config=\%config::config;

## Date : 21th April 2015
## Author : Stephanie Le Gras, slegras@igbmc.fr

## Objectives : This script is design to find the ...

my $bed;
my $min;
my $max;
my $workingDir;
my $fasta;
my $help;
my $name;
my $region;
my $genome;
my $clean;

my $num_arg  = scalar @ARGV;

my $result = GetOptions(
	"b|bed=s"      => \$bed,
	"l|min=s" => \$min,
	"u|max=s" => \$max,
	"w|workingDir=s" => \$workingDir,
	"f|fasta=s" => \$fasta,
	"r|region=s" => \$region,
	"n|name=s" => \$name,
	"h|help" => \$help,
	"g|genome=s" => \$genome,
	"c|clean" => \$clean,

);

my $usage = <<END;

Usage: $0 [options] -w DIRNAME -b FILENAME -f FILENAME 

  -w DIRNAME or --workingDir=DIRNAME	 - Working directory (mandatory)
  -b FILENAME or --bed=FILENAME		 - Input bed file (mandatory)
  -f FILENAME or --fasta=FILENAME 	 - Absolute path to fasta file (mandatory)
  -g STRING or --genome=STRING 		 - Name of assembly of genome to use 
  					   (Default: mm9. Only mm9 is available so far)
  -l INT or --min=INT 			 - Minimum size of sequence (Default: 80)
  -u INT or --max=INT 			 - Maximum size of sequence (Default: 120)
  -r STRING or --region=STRING 		 - Region in the fasta sequence to search. 
  					   Can be either "5'", "3'" or entire. (Default: "5'")
  -n STRING or --name=STRING 	 	 - Prefix name for output file (Default: "bestSeq_[date]")
  -c or --clean 			 - Clean older data and rerun a new analysis 
  					   (Default: Unset -> use prior data)
  -h or --help  		 	 - Print this help

END

die "Bad parameter entry\n$usage" if (@ARGV);
die $usage if (! $result);
die $usage if ( $num_arg < 10 );
die $usage if ($help);

### Checking settings
die "No input BED file.\n" unless ($bed);
die "No working directory or not writable.\n" unless ($workingDir and -d $workingDir);
die "No input for required fasta sequence.\n" unless($fasta);

$region = ($region) ? $region : "5'";
$min = ($min) ? $min : 80;
$max = ($max) ? $max : 120;
$genome = ($genome) ? $genome : "mm9";

##################################################################
################################# Creating a tmp dir
my $tmp = "$workingDir/tmp";

if ($clean){
	rmdir $tmp;
}

if(!-d $tmp){
	mkdir $tmp;
}

##################################################################
################################# get prefix for output file
my $prefix;

######## check whether the tool have been run previously.
## If so, previous data are used 
my @previous = glob($tmp.'/'.'bestSeq*');
@previous = map{basename $_} @previous;
@previous = map{ $_ =~ s/(.*)\.\d.*/$1/; $_} @previous;
@previous = &uniq(@previous);

if(@previous > 1){
	die "Previous run of the tool have been detected. Clean them (--clean tag when running the tool)\
	 or set the --name tag.";
}elsif(@previous == 1){
	warn "Previous analysis detected. Will use it to run this analysis.\n";
	warn "Clean tmp directory to run a new analysis.\n\n";
	$name = $previous[0];
}

if($name){
	$prefix = $name;
}else{
	my $date = strftime "%Y-%b-%e-%H%M%S", localtime;
	$prefix = "bestSeq_".$date; 
}

warn "Prefix of output file will be: ".$prefix.".\n"; 

##################################################################
################################# Creating a log file
## Creating the LOG directory if it doesn't exist
if(!-d "$workingDir/LOG"){
	mkdir "$workingDir/LOG";
}

## Getting the log file name
if ( -f "$workingDir/".basename($0).".log"){
	unlink "$workingDir/".basename($0).".log";
}

my $i = 1;
my $logName = "$workingDir/LOG/".basename($0).".$i.log";

while ( -f $logName ){
	$i++;
	$logName = "$workingDir/LOG/".basename($0).".$i.log";
}

symlink $logName, "$workingDir/".basename($0).".log";

open(my $log, ">".$logName) or die "Cannot create log file : $!";


##################################################################
################################# MAIN
print $log "################ Starting Analysis\n";
print $log "##".`date`;
print $log "\n";

######### Reading in the bed files
$bed = abs_path($bed);
my @bed = &readBED($bed, $log); 

######### Generating a .fai file is needed
## We will need the size of the input fasta file
## the size is taken from fai file.
## If the file doesn't exist, generate the file
$fasta = abs_path($fasta);
if(! -s "$fasta.fai"){
	&getFai($fasta, $log);
}

######### Reading in the bed files
for(my $roi=0; $roi <= $#bed; $roi++){

	## getting the size of the fasta sequence
	my $size = `grep "$bed[$roi][0]" $fasta.fai|cut -f 2`;
	chomp $size;

	# End of the region to analyzed cannot be larger than the size of the sequence
	$bed[$roi][2] = $size if($bed[$roi][2] > $size);

	## Determining the coordinate in which to search for the sequence
	my $start = $bed[$roi][1];
	my $end = $bed[$roi][2];
	my $length=$end - $start+1;

	switch ($region){
		case "entire"	{$start=$start; $end=$end }
		case "5'"		{$start=$start; $end=&round($start+($length/2)); if($end > $bed[$roi][2]){$end=$bed[$roi][2] } }
		case "3'"		{$end=$end; $start=&round($start+($length/2)-1) ; if($start < $bed[$roi][1]){$start=$bed[$roi][1] }  }
		else 			{die "Bad entry value for region: $region."}
	}

	## Extracting coordinate of regions to be analyzed
	my @nbBed = glob($tmp.'/'.$prefix.'*.bed');
	if(@nbBed < ($max-$min+1)){
		&extractCoordinates($bed[$roi][0], $min, $max, $start, $end, "$tmp/$prefix", $log);
	}

	## Extracting fasta sequences from coordinate previously computed
	my @nbFa = glob($tmp.'/'.$prefix.'*.fasta');
	if(@nbFa < ($max-$min+1)){	
		&extractSequences($fasta, "$tmp/$prefix", $min, $max, $log);
	}

	## Getting fastq sequences from fasta sequences
	my @nbFq = glob($tmp.'/'.$prefix.'*.fastq');
	if(@nbFq < ($max-$min+1)){	
		&fastaTofastqBatch("$tmp/$prefix", $min, $max, $log);
	}

	## Aligning sequences to the reference genome
	my @nbSam = glob($tmp.'/'.$prefix.'*.sam');
	if(@nbSam < ($max-$min+1)){		
		&mapDataBatch("$tmp/$prefix", $genome, $min, $max, $log)
	}

	## and count the number of match for each sequences.
	my @nbCount = glob($tmp.'/'.$prefix.'*.count');
	if(@nbCount < ($max-$min+1)){	
		&mapToCountsBatch("$tmp/$prefix", $min, $max, $log)
	}

	## compile all data, export a file with all data and a graph using R
	if(! -s "$tmp/$prefix\_final.tsv"){
		&poolData($prefix, $tmp, dirname(__FILE__), $log)
	}
	
}

print $log "################ End of Analysis\n";
print $log "##".`date`;
print $log "\n";









