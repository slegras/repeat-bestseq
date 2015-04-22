package findBestSequence;

use strict;
use warnings;
use Exporter;

our @ISA= qw(Exporter);
our @EXPORT= qw(&poolData &mapToCountsBatch &uniq &readBED &extractCoordinates &getFai &round &extractSequences &fastaTofastqBatch &mapDataBatch);

use config;
use vars qw(%config );
*config=\%config::config;


sub extractCoordinates{
	my ($seq, $min, $max, $start, $end, $prefix, $log) = @_;

	for(my $i=$min; $i<=$max; $i++){
		open(OUT, ">"."$prefix.$i.bed") or die "Cannot create file $prefix.$i.bed: $!";
		for(my $j=$start; $j+$i <= $end; $j++){
			print OUT $seq."\t".$j."\t";
			print OUT $j+$i."\n";
		}
		close OUT;
	}

	return 1;
	
}

sub extractSequence{
	my ($bed, $fasta, $out, $log) = @_;

	## Build up the command line
	my $cmdline = "$config{'bedtools'} getfasta -fi $fasta -bed $bed -fo $out";

	my $title = "Extracting fasta file.";

	&runCmdLine($cmdline, $log, $title);
	return 1;
}

sub extractSequences{
	my ($fasta, $prefix, $min, $max, $log) = @_;

	for(my $i=$min; $i<=$max; $i++){
		extractSequence("$prefix.$i.bed", $fasta, "$prefix.$i.fasta", $log);
	}

	return 1;
}

sub fastaTofastq{

	my ($in, $out, $log)= @_;

	open FILE, "<".$in or die "Cannot open file $in: $!";
	open OUT, ">".$out or die "Cannot create file $out: $!";

	my $header="";
	my ($sequence, $sequence_length, $sequence_quality);
	while(<FILE>) {
	        chomp;
	        if (/^>(.+)/) {
	                if($header ne "") {
	                        print OUT "\@".$header."\n";
	                        print OUT $sequence."\n";
	                        print OUT "+"."\n";
	                        print OUT $sequence_quality."\n";
	                }
	                $header = $1;
			$sequence = "";
			$sequence_length = "";
			$sequence_quality = "";
	        }
		else { 
			$sequence .= $_;
			$sequence_length = length($_); 
			for(my $i=0; $i<$sequence_length; $i++) {$sequence_quality .= "I"} 
		}
	}
	close FILE;
	print OUT "\@".$header."\n";
	print OUT $sequence."\n";
	print OUT "+"."\n";
	print OUT $sequence_quality."\n";

	close OUT;

	return 1;

}

sub fastaTofastqBatch{
	my ($prefix, $min, $max, $log) = @_;

	for(my $i=$min; $i<=$max; $i++){
		fastaTofastq("$prefix.$i.fasta", "$prefix.$i.fastq", $log);
	}

	return 1;

}

sub getFai{
	my ($fasta, $log) = @_;

	## Build up the command line
	my $cmdline = "$config{'samtools'} faidx $fasta";

	my $title = "Generating *.fai file.";

	&runCmdLine($cmdline, $log, $title);
	return 1;

}

sub mapData{
	my ($fastq, $genome, $sam, $log) = @_;

	## Build up the command line
	my $cmdline = "$config{'bowtie'} -a -v 0 -p 4 $genome $fastq -S $sam";

	my $title = "Generating *.fai file.";

	&runCmdLine($cmdline, $log, $title);
	return 1;

}

sub mapDataBatch{
	my ($prefix, $genome, $min, $max, $log) = @_;

	for(my $i=$min; $i<=$max; $i++){
		mapData("$prefix.$i.fastq", $genome, "$prefix.$i.sam", $log);
	}

	return 1;

}

sub mapToCounts{
	my ($sam, $out, $log) = @_;

	my %count;
	open(SAM, "<".$sam) or die "Cannot open the file $sam: $!";
	while(<SAM>){
		next if(/^@/);
		chomp;
		my @tab = split "\t";
		$count{$tab[0]}{'count'}++;
		$count{$tab[0]}{'sequence'} = $tab[9];
		$count{$tab[0]}{'length'} = length($tab[9]);
	}
	close SAM;

	open(OUT, ">".$out) or die "Cannot create output file $out: $!";
	print OUT "Sequence name\tSequence\tCount\n";
	foreach my $key (keys %count){
		print OUT $key."\t".$count{$key}{'sequence'}."\t".$count{$key}{'count'}."\n";
	}
	close OUT;

	return 1;

}

sub mapToCountsBatch{
	my ($prefix, $min, $max, $log) = @_;

	for(my $i=$min; $i<=$max; $i++){
		mapToCounts("$prefix.$i.sam", "$prefix.$i.count", $log);
	}

	return 1;
}

sub poolData{

	my ($prefix, $tmp, $path, $log) = @_;

	## Build up the command line
	my $cmdline = "Rscript $path/findBestSequence.R \"$tmp\" \"$prefix\"";

	my $title = "Compiling data with R.";

	&runCmdLine($cmdline, $log, $title);
	return 1;

}

sub readBED{
	my ($bed, $log) = @_;

	open(BED, "<".$bed) or die "Cannot open BED file $bed: $!";
	my @bed;
	while(<BED>){
		chomp;
		push @{$bed[$#bed+1]}, split "\t"; 
	}
	return @bed;
}

sub round{
 	$_[0] > 0 ? int($_[0] + .5) : -int(-$_[0] + .5)
}


sub runCmdLine {
	my ($cmdline, $log, $title) = @_;

	## Run the command line and output into the log file
	print $log "\n##############\n" ;
	print $log "## Start analysis: ".`date`."\n";
	print $log "## $title. \n\n";
	print $log "$cmdline\n";
	print $log "#<--- Output: --------------------------------------------------\n";

	my $result = `$cmdline 2>&1`;
	print $log $result."\n";
	
	print $log "\n";

	print $log "#--------------------------------------------------------------->\n";
	print $log "## End of analysis: ".`date`."\n";

	return 1;
}

sub uniq { 
	my %seen; 
	grep !$seen{$_}++, @_ 
}

1;
