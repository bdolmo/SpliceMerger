#!/usr/bin/perl

use strict;
use warnings;
use Cwd qw(cwd abs_path);
use Getopt::Long;
use File::Basename;

my $dirname = dirname(__FILE__);
my $vcf;
my $exons = 'ghrc37';
my $genome;
my $bedtools = `which bedtools`;
if($bedtools eq '' || !defined$bedtools) {
	print "Please, introduce a the route of the BEDtools binary\n";
  	$bedtools = <STDIN>;
	chomp $bedtools;
       	print "$bedtools\n";
}

my $keept = 'no';

invocation () if (@ARGV<1 or !GetOptions(
	'genome=s' =>\$genome,
	'vcf=s'    =>\$vcf,
	'exons=s'  =>\$exons,
	'temp=s'   =>\$keept,
         )
);
my $outputDir = $ARGV[0];

my $regions = $exons eq 'ghrc37' ? "$dirname/assembly/exome.ghrc37.strand.bed" : "$dirname/assembly/exome.ghrc38.bed";  

if (!defined($vcf) || $vcf !~/.vcf$/){
	print "#ERROR: VCF not found, please submit a VCF file\n";
	exit;
}
if (!defined($genome)|| $genome !~/.fasta$/){
	print "#ERROR: genome file not found, please submit a genome file in FASTA format\n";
	exit;
}
if (!defined($outputDir)) {
	my @chars = ('a'..'z', 'A'..'Z', '0'..'9');
	$outputDir = join '', map { @chars[rand @chars] } 1 .. 8;
	print "\n · No output directory was defined, $outputDir random directory was created ·\n";
}

sub invocation {
 	print("
Description: SpliceMerger: merging of three splice site predictions for H.sapiens.\n
Usage:  spliceMerger.pl -genome <genome FASTA> -vcf <VCF> -exons ghrc37/ghrc38 OUTPUT_DIR
Params:
	 -genome    STRING	Path to genome file in fasta format.
	 -vcf       STRING	VCF file to be analayzed. Must contain RefSeq or Ensemble transcript IDs.
	 -exons     STRING  	Genome assembly (grch37 or grch38, default = grch37).
	 -temp      STRING      Keep temporary files. yes/no. default = no

2015-2016 Gendiag.exe, SL Barcelona\n\n");
exit;
}

#--------------------------Splicing Predictors---------------------------------#

##### MaxEntScan ######
# 3' splicing sites
my $MaxEntScan3 =
  ( -x "$dirname/bin/MaxEntScan/score3.pl" )
  ? "$dirname/bin/MaxEntScan/score3.pl"
  : die "Can't execute $dirname/bin/MaxEntScan/score3.pl";

# 5' splicing sites
my $MaxEntScan5 =
  ( -x "$dirname/bin/MaxEntScan/score5.pl" )
  ? "$dirname/bin/MaxEntScan/score5.pl"
  : die "Can't execute $dirname/bin/MaxEntScan/score5.pl";

##### FSPLICE ######
my $FSPLICE =
  ( -x "$dirname/bin/FSPLICE/fsplice" )
  ? "$dirname/bin/FSPLICE/fsplice"
  : die "Can't execute $dirname/bin/FSPLICE/fsplice";

##### GeneSplicer ######
my $GeneSplicer =
  ( -x "$dirname/bin/GeneSplicer/bin/linux/genesplicer" )
  ? "$dirname/bin/GeneSplicer/bin/linux/genesplicer"
  : die "Can't execute $dirname/bin/GeneSplicer/bin/linux/genesplicer";

#------------------------------------------------------------------------------#

			#########
			# START #
			#########

my $start_run  = time();

# Displaying parameters
printParameters ($vcf, $outputDir, $genome, $exons);

# We load a Hash of Exon (HoE) coordinates
my $HoE = keepInMemory($regions);
my %HoE = %$HoE;

open (VCF, "<", $vcf) || die "Unable to open $vcf\n";
my $csv = $vcf;
$csv =~s/.vcf/.splicing.csv/;

my $filename = basename($csv);

mkdir $outputDir;

my ($overlap1, $overlap2, $overlap3, $total) = (0, 0, 0, 0);
my ($MES_3variation, $MES_5variation, $MES3out) = (0, 0, 0);
my ($MES5out, $outFSPLICE, $outGeneSplicer) = 0;  
open (CSV, ">", "$outputDir/$filename") || die "Unable to open $outputDir/$filename\n";
open (TEMP,">", "$outputDir/SpliceMerger.tmp")|| die "Unable to open $outputDir/SpliceMerger.tmp\n";
print CSV "CHR\tPOS\tCHANGE\tHGVSC\tGENE\tISOFORM\tEXON\tSTRAND\tDIST to EXON\tTYPE\tMES3'REF score\tMES3' MUT score\tvariation\tMES5'REF score\tMES5' MUT score\tvariation\tFSPLICE REF\tFSPLICE MUT\tGeneSplicer REF\tGeneSplicer MUT\tMaxEntScan Prediction\tFSPLICE Prediction\tGeneSplicer Prediction\n";
print " ----------PROGRESS----------\n";
while (my $line =<VCF>) {
	chomp $line;
	next if $line =~/^#/;
	my @info = split (/\t/, $line);
	my $chr = $info[0];
	my $start = $info[1];
	my $rs  = $info[2];
	my $ref = $info[3];
	my $alt = $info[4];


	print "$ref\t$alt\n";

	# For Gendicall indel annotation we will keep the first indel.
	if ($alt =~/,/) {
		my @tmpAlt = split (/,/, $alt);
		$alt = $tmpAlt[0];
	}
	my @tmp = split (/;/, $info[11]);
	my @geneTmp = grep (/GENE=/ , @tmp);
	my $gene = join ("", @geneTmp);
	$gene =~s/GENE=//;
	my @transTmp = grep (/ENST=/, @tmp);
	my $transcript = join ("", @transTmp);
	$transcript =~s/ENST=//;

	my ($hgvsc) = grep(/HGVSC=/, @tmp);
	if ($hgvsc) {
		$hgvsc=~s/HGVSC=//;
	}
	else {
		$hgvsc = ".";
	}

	my $progress = $total/10;
        my $progressbar = ($progress - int($progress))?'float':'int';
	$progress = $progress * 10;
	print " INFO: $progress variants analyzed\n" if $progressbar eq 'int' && $progress > 0; 

	# Here we retrieve the nearest position to the exon and the type of ss
	my ($distance, $exonSt, $exonEnd, $exon, $nexons, $type, $strand) = getNearestExon ($chr, $start, $gene, $transcript, \%HoE);
	print TEMP "\nDIST:$distance, start of exon: $exonSt, end of exon: $exonEnd, type: $type, strand: $strand\n";

	if(!defined($type)) {
		print " INFO: skipped variant $chr:$start\t ref:$ref/alt:$alt\n";
		next;
	}

	my ($newref, $newalt) = ($ref, $alt); 	

	if ($strand eq '-') {
		$newref = reverse_complement($ref);
		$newalt = reverse_complement($alt);
	}

	my $strandOut = $strand eq '+' ? 'plus' : 'minus';

	if (($strand eq '+' && $type eq '3' && $distance >= -20 && $distance < 3 && $exon > 1)  ||  ($strand eq '+' && $type eq '5' && $distance >= -2 && $distance < 6 && $exon < $nexons) || ($strand eq '-' && $type eq '3' && $distance <= 20 && $distance >= -2 && $exon < $nexons) || ($strand eq '-' && $type eq '5' && $distance <= 2 && $distance >= -6 && $exon < $nexons)) {
  
		# MaxEntScan
		my ($seq3wt , $seq3mut, $seq5wt, $seq5mut, $score3wt, $score3mut, $score5wt, $score5mut) = runMaxEntScan ($chr, $exonSt, $exonEnd, $newref, $newalt, $distance, $type, $strand);

		# FSPLICE
		my ($fspliceWT, $fspliceMUT) = runFSPLICE ($chr, $start, $exonEnd, $newref, $newalt, $distance, $type, $strand);
		my @FsWT  = @$fspliceWT;
		my @FsMUT = @$fspliceMUT;

		# GeneSplicer
		my ($wtGeneSplicer, $mutGeneSplicer) = runGeneSplicer ($chr, $start, $exonEnd, $newref, $newalt, $distance, $type, $strand);
		
		# Scores FSPLICE and GeneSplicer
		$outFSPLICE     = scoreFSPLICE (\@FsWT, \@FsMUT);
		$outGeneSplicer = scoreGeneSplicer ($wtGeneSplicer, $mutGeneSplicer);

		# Here we calculate the % of variation for MaxEnt. Values > 30 % are considered significant
		if ($type == '3') { 
			$MES_3variation = calculateVariation ($score3wt, $score3mut);
			$MES3out = abs($MES_3variation) > 30 ? 1 : 0;
	     print CSV "$chr\t$start\t$ref>$alt\t$hgvsc\t$gene\t$transcript\t$exon\t$strandOut\t$distance\t$type\'ss\t$score3wt\t$score3mut\t$MES_3variation\t\-\t\-\t\-\t@FsWT\t@FsMUT\t$wtGeneSplicer\t$mutGeneSplicer\t$MES3out\t$outFSPLICE\t$outGeneSplicer\n";    
		}
		if ($type == '5') { 
			$MES_5variation = calculateVariation ($score5wt, $score5mut);
			$MES5out = abs($MES_5variation) > 30 ? 1 : 0;
		   print CSV "$chr\t$start\t$ref>$alt\t$hgvsc\t$gene\t$transcript\t$exon\t$strandOut\t$distance\t$type\'ss\t-\t-\t-\t$score5wt\t$score5mut\t$MES_5variation\%\t@FsWT\t@FsMUT\t$wtGeneSplicer\t$mutGeneSplicer\t$MES5out\t$outFSPLICE\t$outGeneSplicer\n";
		}
	}
	else {
		# FSPLICE
		my ($fspliceWT, $fspliceMUT) = runFSPLICE ($chr, $start, $exonEnd, $newref, $newalt, $distance, $type, $strand);
		my @FsWT  = @$fspliceWT;
		my @FsMUT = @$fspliceMUT;

		# GeneSplicer
		my ($wtGeneSplicer, $mutGeneSplicer) = runGeneSplicer ($chr, $start, $exonEnd, $newref, $newalt, $distance, $type, $strand);
	
		$outFSPLICE     = scoreFSPLICE (\@FsWT, \@FsMUT);
		$outGeneSplicer = scoreGeneSplicer ($wtGeneSplicer, $mutGeneSplicer);
		$MES5out = "";
  		$MES3out = "";
		if ($type eq '3') {
			print CSV "$chr\t$start\t$ref>$alt\t$hgvsc\t$gene\t$transcript\t$exon\t$strandOut\t$distance\t$type\'ss\t-\t-\t-\t\-\t\-\t\-\t@FsWT\t@FsMUT\t$wtGeneSplicer\t$mutGeneSplicer\tNA\t$outFSPLICE\t$outGeneSplicer\n";
		}
		if ($type eq '5') {
			print CSV "$chr\t$start\t$ref>$alt\t$hgvsc\t$gene\t$transcript\t$exon\t$strandOut\t$distance\t$type\'ss\t-\t-\t-\t-\t-\t-\t@FsWT\t@FsMUT\t$wtGeneSplicer\t$mutGeneSplicer\tNA\t$outFSPLICE\t$outGeneSplicer\n";
		}
	}
	#exit if $line =~/116283343/;
 my $MESout;
 $MESout = $MES3out if defined($MES3out);
 $MESout = $MES5out if defined($MES5out);
 $total++;
 $overlap3 ++ if ($MESout eq $outFSPLICE  && $MESout eq $outGeneSplicer && $outFSPLICE eq $outGeneSplicer);
 $overlap2 ++ if ($MESout eq $outFSPLICE && $MESout ne $outGeneSplicer || $outFSPLICE eq $outGeneSplicer && $outFSPLICE ne $MESout);	
 $overlap1 ++ if ($MESout ne $outFSPLICE && $MESout  ne $outGeneSplicer && $outFSPLICE ne $outGeneSplicer);
} 
 close (TEMP); 
 if ($keept eq 'no') {
 	unlink (glob ("$outputDir/*.fa"), glob ("$outputDir/*.bed"), glob ("$outputDir/*.out"), glob ("$outputDir/*.mut"), glob ("$outputDir/*.wt"), glob ("$outputDir/*.tmp"));
 }

 my $end_run = time();
 my $run_time = $end_run - $start_run;
 print "\n FINISHED: overall execution time $run_time s\n";
 my $perc_olap1 = sprintf "%2.2f",  100 * $overlap1 / $total;
 my $perc_olap2 = sprintf "%2.2f",  100 * $overlap2 / $total;
 my $perc_olap3 = sprintf "%2.2f",  100 * $overlap3 / $total;
 print "\n ----------SUMMARY----------";
 print "
 No concordant predictors: $overlap1 ($perc_olap1%)
 Two concordant predictors: $overlap2 ($perc_olap2%)
 Three concordant predictors: $overlap3 ($perc_olap3%)\n\n";


			#########
			#  END  #
			#########

##################################################################### 

sub printParameters {
 my ($vcf, $outdir, $genome, $exons) = @_;
 my $datestring = localtime();
 print "\n Program: SpliceMerger (large-scale analysis of splicing sites)\n";
 print " Version: 0.1\n";
 print "  \n   # Date:  $datestring
   # VCF file name: $vcf
   # Output directory: $outdir
   # Genome file: $genome
   # Assembly: $exons\n\n";
}

sub keepInMemory {
 my $regions = shift;
 my %hash = ();
 open (IN, "<", $regions) || die "Unable to open $regions\n";
 while (my $line =<IN>) {
	chomp $line;
	my @info = split (/\t/, $line);
	my $chr = $info[0];
	my $start = $info[1];
	my $end = $info[2];
	my $tmp = $info[3];
	my $strand = $info[5];
	my @tmp  = split (/[;_]/, $tmp);
	my $transcript = $tmp[0];
	my $gene = $tmp[2];
	my $exon = $tmp[3];

	# We load a multidimensional array that contains every gene, transcript and exon
	$hash{$gene}{$transcript}{$exon} = "$chr\t$start\t$end\t$strand";
 }
  return (\%hash);
}

sub getNearestExon {

 my $chr   = shift;
 my $start = shift;
 my $gene  = shift;
 my $transcript = shift;
 my $hash  = shift;
 my %hash  = %$hash;
 my ($distanceL, $distanceR, $type);
 my ($minChr, $minExonSt, $exon, $minExonEnd);
 my $minL = 10e8;
 my $minR = 10e9;
 my @intron = ();
 my @st;
 my @end;
 my $nexons;
 my $j = 0;
 my $i = 0;
 my $strand;
 # Gene
 foreach my $Ge ( keys %hash ) {
	if ($gene eq $Ge) {
		# Transcript
		foreach my $Tr (keys %{$hash{$Ge} } ) {
			if ($transcript eq $Tr) {
				# Exon
				foreach my $Ex ( sort {$a<=>$b} keys %{$hash{$Ge}{$Tr} } ) {
					$nexons = scalar(keys %{$hash{$Ge}{$Tr}});
					my ($chr, $exonSt, $exonEnd) = split ( /\t/, $hash{$Ge}{$Tr}{$Ex});
					$strand = (split /\t/, $hash{$Ge}{$Tr}{$Ex})[3];
					push(@st, $exonSt);
					push(@end, $exonEnd);
					$st[$i] = $exonSt;
					$end[$i] = $exonEnd;
					my $distanceL = $start - $exonSt;
					my $distanceR =	$start - $exonEnd;
					if (abs($distanceL) < abs($minL)) {
						$minL = $distanceL;
						($minChr, $minExonSt, $minExonEnd) = split (/\t/, $hash{$Ge}{$Tr}{$Ex});
						$exon = $Ex;
					}
					if (abs($distanceR) < abs($minR)) {
						$minR = $distanceR;
						($minChr, $minExonSt, $minExonEnd) = split (/\t/, $hash{$Ge}{$Tr}{$Ex});
						$exon = $Ex;
					}
					if (defined ($end[$i-1]) && $st[$i] == $end[$i-1]) {
						last;
					}
				}

				if (abs($minL) < abs($minR)) {

					$type = $strand eq '+' ? '3' : '5';
					return ($minL, $minExonSt, $minExonEnd, $exon, $nexons, $type, $strand);
				}
				else {
					$type = $strand eq '+' ? '5' : '3';
					return ($minR, $minExonSt, $minExonEnd, $exon, $nexons, $type, $strand);

				}
			}
		}
	}
 }
}	

sub scoreFSPLICE {

	my $wildtype = shift;
	my $mutant   = shift;

        my @wildtype = @$wildtype;
	my @mutant   = @$mutant;
 
	# The idea is to compare the WT splicing sites vs MUT splicing sites. If the MUT variants changes the position, type or creates/distroys the WT ss, then it is reported as splicing variant.
	# EXAMPLE
	# WT  -> Acceptor(AG) sites. Treshold      4.175 (90%).        1 P:      51 W:  8.50 Seq: tctgtAGtgtct  Donor(GT) sites. Treshold      6.099 (90%).
	# MUT -> Acceptor(AG) sites. Treshold      4.175 (90%).  Donor(GT) sites. Treshold      6.099 (90%).
	# COMPARISON: MUT lacks the WT splicing site, therefore the MUT variant is classified as a ss variant.

	my $wt  = join ("\t", @wildtype);
	my $mut = join ("\t", @mutant);

	my @WT  = split (/\t/, $wt);
	my @MUT = split (/\t/, $mut);
	my ($donor, $acceptor) = "";

	# grep positions
	my @posWT  = grep(/W/, @WT); 
	my @posMUT = grep(/W/, @MUT);
	
	my $outcome = @posWT eq @posMUT ? 0 : 1;
	
	my @wtN  = grep (/^\d+$/, @WT);
	my @mutN = grep (/^\d+$/, @MUT);

	my %wtN  = map { ($_) => $_ } @wtN;	
	my %mutN = map { ($_) => $_ } @mutN;

	foreach my $keys (keys %wtN ) {
		if (!exists $mutN{$keys}){
			$outcome = 1;
			return $outcome;
		} 
	}
        return $outcome;
}

sub scoreGeneSplicer {

	my $wildtype = shift;
	my $mutant   = shift;
	chomp ($wildtype, $mutant);
        my @wildtype = "";
	my @mutant   = "";
	my $outcome  = 0;

	# Same idea as for FSPLICE. If the MUT variants changes the position, type or creates/distroys the WT ss, then it is reported as splicing variant.
	# EXAMPLE
	# WT  -> 101 102 8.905251 Medium donor
	# MUT -> nothing
	# COMPARISON: MUT lacks the donor site, therefore the MUT variant is classified as a ss variant.
	#print TEMP "GENESPLICER wt:$wildtype\t mut: $mutant\n";

	my @wt  = split (" ", $wildtype);
	my @mut = split (" ", $mutant);

	if (scalar(@wt) != scalar(@mut)) {
		$outcome = 1;
		return $outcome;
	}
	if (defined($wildtype) && defined($mutant) && $wildtype ne '' &&  $mutant ne '') {

		# We grep only the whole numbers (the positions)
		my @WT  = grep (/^\d+$/, @wt);
		my @MUT = grep (/^\d+$/, @mut);

		my %WT  = map { ($_) => $_ } @WT;	
		my %MUT = map { ($_) => $_ } @MUT;

		foreach my $keys (keys %WT) {
			if (!exists $MUT{$keys}){
				$outcome = 1;
				return $outcome;
			} 
		}
	}
	if (!(defined$wildtype) && !(defined$mutant) || $wildtype eq '' && $mutant eq '') {
		$outcome = 0;
		return $outcome;
	}
	if (($wildtype ne '' && $mutant eq '' || defined($wildtype) && !defined($mutant))||($wildtype eq '') && ($mutant ne '') || !defined($wildtype) && defined($mutant)) {
		$outcome = 1;
		return $outcome;
	}
	else {
		$outcome = 0;
		return $outcome;
	}
}

sub calculateVariation {

 my ($scoreWT, $scoreMUT) = @_;
 $scoreWT = '1' if !defined $scoreWT;
 my $sum = $scoreMUT + $scoreWT;
 my $MES_variation;
 if (($scoreMUT < 0) && ($scoreWT > 0)) {
	$sum = $scoreWT + abs$scoreMUT;
	$MES_variation =  sprintf "%2.2f",  100 * $sum / $scoreWT;
 }
 if ($scoreMUT < $scoreWT) {
	$MES_variation = sprintf "%2.2f",  100 * $sum / $scoreWT;
 }
 if (($scoreMUT > 0) && ($scoreWT > 0)) {
	$sum = $scoreWT - $scoreMUT;
	$MES_variation = sprintf "%2.2f",  100 * $sum / $scoreWT;
 }
 if ($scoreWT != 0) {
	 $MES_variation = sprintf "%2.2f",  100 * $sum / $scoreWT;
 }
 else {
	$MES_variation = 0; 
 }
 return ($MES_variation);

}

sub runMaxEntScan {

 my ($chr, $pos, $exonEnd, $ref, $alt, $distance, $ss, $strand)  = @_;

	my ($seq3wt, $seq3mut, $seq5wt, $seq5mut, $score3wt, $score3mut, $score5wt, $score5mut) = (0, 0, 0, 0, 0, 0, 0, 0);

	# for 3' splicing sites, each sequence must be 23 bases long. [20 bases in the intron][3 base in the exon]
	#		Example Fasta File
	#		> dummy1
	#		ttccaaacgaacttttgtAGgga 

	# For Start offset
	my %stOffset  = ( 
	  '3ss' => 20,
	  '5ss' => 3
	);
	# For End offset
	my %endOffset = (
	  '3ss' => 3,
	  '5ss' => 6
	);
	# Type of Splicing site 
	my %type = (
	  '3ss' => 3,
	  '5ss' => 5
	);

	if ($ss == '3') {
		retrieveVariant ($chr, $pos, $exonEnd, $ref, $alt, $stOffset{"3ss"}, $endOffset{"3ss"}, $type{"3ss"}, $distance, $ss, $strand);	

		`perl $MaxEntScan3 $outputDir/3.wt > $outputDir/mes3.WT.out`;
		`perl $MaxEntScan3 $outputDir/3.mut > $outputDir/mes3.MUT.out`;
		open (WT,"<", "$outputDir/mes3.WT.out") || die "Unable to open $outputDir/mes3.WT.out\n";
			while (my $line =<WT>) {
				chomp $line;
				$score3wt = "no score";		
				($seq3wt, $score3wt) = split (/\t/, $line);
			}
		close (WT);

		open (MUT,"<", "$outputDir/mes3.MUT.out") || die "Unable to open $outputDir/mes3.MUT.out\n";
			while (my $line =<MUT>) {
				chomp $line;
				$score3mut = "no score";		
				($seq3mut, $score3mut) = split (/\t/, $line);
			}
		close (MUT);
	}
	# for 5' splicing sites each sequence must be 9 bases long [3 bases in exon][6 bases in intron]. 
	#		Example Fasta File
	#		> dummy1
	#		cagGTAAGT
	
	if ($ss == '5') {
		retrieveVariant ($chr, $pos, $exonEnd, $ref, $alt, $stOffset{"5ss"}, $endOffset{"5ss"}, $type{"5ss"}, $distance, $ss, $strand);

		`perl $MaxEntScan5 $outputDir/5.wt > $outputDir/mes5.WT.out`;
	
		`perl $MaxEntScan5 $outputDir/5.mut > $outputDir/mes5.MUT.out`;
	
		open (WT, "<", "$outputDir/mes5.WT.out") || die "Unable to open $outputDir/mes5.WT.out\n";
			while (my $line =<WT>) {
				chomp $line;
				$score5wt = "no score";		
				($seq5wt, $score5wt) = split (/\t/, $line);
			}
		close (WT);

		open (MUT, "<", "$outputDir/mes5.MUT.out") || die "Unable to open $outputDir/mes5.MUT.out\n";
			while (my $line =<MUT>) {
				chomp $line;	
				$score5mut = "no score";	
				($seq5mut, $score5mut) = split (/\t/, $line);
			}
		close (MUT);
	}
	return ($seq3wt, $seq3mut, $seq5wt, $seq5mut, $score3wt, $score3mut, $score5wt, $score5mut);
}

sub retrieveVariant {

 my ($chr, $start, $exonEnd, $ref, $alt, $stOffset, $endOffset, $type, $distance, $ss, $strand) = @_; 
 my ($deletion, $bed, $fasta, $end, $startN);

 print TEMP "type: $type\t SS: $ss\t distance: $distance\n";

 # 3 Means --> 3' splicing sites
 if ($type eq '3') {
	$bed   = "splice3.bed";
	$fasta = "splice3.fa";
 }
 # 5 Means --> 5' splicing sites
 elsif ($type eq '5')  {
	$bed   = "splice5.bed";
	$fasta = "splice5.fa";
 }
 # For FSPLICE
 elsif ($type eq 'FS') {
	$bed   = 'fsplice.bed';
	$fasta = "wt.fsplice.fa";
 }
 # For GeneSplicer
 elsif ($type eq 'GS') {
	$bed   = 'GeneSplicer.bed';
	$fasta = "wt.GeneSplicer.fa";
 }

 my $SeqLength = $stOffset + $endOffset;
 my $Lalt = length($alt);
 my $Lref = length($ref);
 my $ntd = $alt;	
 $ntd =~s/[-+]//;
 $ntd = length($ntd);    
 my ($MUT, $WT);

	# Here we calculate the coordinates to extract FASTA sequences.
	# Note that when a deletion occurs, we need to extract at end + length of deletion.
	if ($Lref > $Lalt) { 

		$start = $type eq '5' && $strand eq '+' || $type eq '3' && $strand eq '-' ? $exonEnd : $start;

		if (($strand eq '+') && ($type eq '5') && ($distance <= 0)) {
			$start = $exonEnd;
		}
		if (($strand eq '+') && ($type eq '5') && ($distance > 0)) {
			$start = $exonEnd;
		}
		my $Ldeletion  = $Lref - $Lalt;
		$start = $type eq '5' && $strand eq '+' || $type eq '3' && $strand eq '-' ? $exonEnd : $start;
		$startN = $start;

		# Here we define the coordinates to extract the fasta

		$start  = $startN - $stOffset if $type eq '5' && $strand eq '+';
		$start  = $startN - $endOffset -1 if $type eq '5' && $strand eq '-';
		
		print TEMP "stOffset: $stOffset\t endOffset: $endOffset\n";
		print TEMP "Bedtools start1: $start\n";

		$end    = $startN + $endOffset if $type eq '5' && $strand eq '+';
		$end    = $startN + $stOffset -1 if $type eq '5' && $strand eq '-';

		#print TEMP "Bedtools end1: $end\n";

		$start = $startN - $stOffset -1 if $type eq '3' && $strand eq '+';
		$start = $startN - $endOffset if $type eq '3' && $strand eq '-';

		$end    = $startN + $endOffset-1 if $type eq '3' && $strand eq '+';
		$end    = $startN + $stOffset if $type eq '3' && $strand eq '-';

		#$startN   = $start;
		#$start    = $startN - $stOffset;
		#$end      = $startN + $Ldeletion + $endOffset;

		if ($type eq 'FS' || $type eq 'GS') {
			print TEMP "equals to $type\n";
			$start = $startN - $stOffset -1;
			$end = $startN + $endOffset -1;
		}


		print TEMP "Bedtools: $start\t$end\n";
	}
	# For SNPs or insertion
	else { 

		# Here we choose between the start or the end of the exon depending on its strand.
		$start = $type eq '5' && $strand eq '+' || $type eq '3' && $strand eq '-' ? $exonEnd : $start;
		$startN = $start;

		# Here we define the coordinates to extract the fasta

		$start  = $startN - $stOffset if $type eq '5' && $strand eq '+';
		$start  = $startN - $endOffset -1 if $type eq '5' && $strand eq '-';
		
		print TEMP "stOffset: $stOffset\t endOffset: $endOffset\n";
		print TEMP "Bedtools start1: $start\n";

		$end    = $startN + $endOffset if $type eq '5' && $strand eq '+';
		$end    = $startN + $stOffset -1 if $type eq '5' && $strand eq '-';

		#print TEMP "Bedtools end1: $end\n";

		$start = $startN - $stOffset -1 if $type eq '3' && $strand eq '+';
		$start = $startN - $endOffset if $type eq '3' && $strand eq '-';

		$end    = $startN + $endOffset-1 if $type eq '3' && $strand eq '+';
		$end    = $startN + $stOffset if $type eq '3' && $strand eq '-';

		if ($type eq 'FS' || $type eq 'GS') {
			print TEMP "equals to $type\n";
			$start = $startN - $stOffset -1;
			$end = $startN + $endOffset -1;
		}

		print TEMP "Bedtools: $start\t$end\n";
	}
	`echo "$chr\t$start\t$end\t$ref;$alt\t$distance\t$strand\n" >  $outputDir/$bed`;

	my $cmd = "bedtools getfasta -s -fi $genome -bed $outputDir/$bed -fo $outputDir/$fasta";
	system ("$cmd");

	$WT = `grep -v \'>\' $outputDir/$fasta`;
	chomp $WT;
	print TEMP "$ref\t$alt\n";
        my $position;
	if ($type eq '3') {
		if ($strand eq '+') {
			if ($distance <= 0) { # Upstream variant
				$position = length($WT) - abs($distance) - $endOffset;
				print TEMP "$chr, $start, $exonEnd, $ref, $alt, $stOffset, $endOffset, DISTANCE = $distance\tTYPE: $type\tPOSITION:$position\n";
			}
			elsif ($distance > 0 && $distance < 3) { # Exonic variant
				$position = length($WT) + abs($distance) - $endOffset -1;	
				print TEMP "EXONIC 3; strand '+' $chr, $start, $exonEnd, $ref, $alt, $stOffset, $endOffset, DISTANCE = $distance\tTYPE: $type\tPOSITION:$position\n";
			}
		}
		if ($strand eq '-') {
			if ($distance <= 0) { # Exonic variant
				$position = length($WT) - abs($distance) -1;	
				print TEMP "EXONIC 3; strand '-' $chr, $start, $exonEnd, $ref, $alt, $stOffset, $endOffset, DISTANCE = $distance\tTYPE: $type\tPOSITION:$position\n";
			}
			elsif ($distance > 0 ) { # Downstream variant
				$position = length($WT) - abs($distance) - $endOffset;
				print TEMP "$chr, $start, $exonEnd, $ref, $alt, $stOffset, $endOffset, DISTANCE = $distance\tTYPE: $type\tPOSITION:$position\n";	
			}
		}
	}
	elsif ($type eq '5') {
		if ($strand eq '+') {
			if ($distance > 0) { # Downstream variant
				$position = length($WT) + abs($distance) - $endOffset -1;
				print TEMP "$WT, $chr, $start, $exonEnd, $ref, $alt, $stOffset, $endOffset, DISTANCE = $distance\tTYPE: $type\tPOSITION:$position\n";
			}
			else { # Exonic variant
				$position = length($WT) - abs($distance) - $endOffset -1;
				print TEMP "EXONIC 5 $position\n";
				print TEMP "EXONIC 5 $chr, $start, $exonEnd, $ref, $alt, $stOffset, $endOffset, DISTANCE = $distance\tTYPE: $type\tPOSITION:$position\n";
			}
			my $lWT = length($WT);
			if (abs($distance) > 9) {		
			 	$position =  $lWT - $endOffset ;
			}
		}
		if ($strand eq '-') {
			if ($distance > 0) { # Exonic variant
				$position = length($WT) - abs($distance) - $endOffset -1;
				print TEMP "EXONIC 5 $position\n";
				print TEMP "EXONIC 5 $chr, $start, $exonEnd, $ref, $alt, $stOffset, $endOffset, DISTANCE = $distance\tTYPE: $type\tPOSITION:$position\n";
			}
			else { # Upstream variant
				$position = length($WT) + abs($distance) - $endOffset -1;
				print TEMP "EXONIC $position\n";
				print TEMP "$WT, $chr, $start, $exonEnd, $ref, $alt, $stOffset, $endOffset, DISTANCE = $distance\tTYPE: $type\tPOSITION:$position\n";
			}
		}
	}
	else {
		$position = $endOffset if $strand eq '+';
		$position = $endOffset -1 if $strand eq '-';
		print TEMP "$chr, $start, $ref, $exonEnd, $alt, $stOffset, $endOffset, DISTANCE = $distance\tTYPE: $type\tPOSITION:$position\n";
	}
 
	# For SNPs
	if ($Lref == $Lalt) { 
		my @tmp  = split (//, $WT);
		my $ltmp = scalar(@tmp);
		if ($Lref == '1' && $Lalt == '1') {
			print TEMP "REFERENCE: $tmp[$position]\n";
			$tmp[$position] = $alt;
			print TEMP "ALTERNATIVE: $tmp[$position]\n";
			$MUT  = join ("", @tmp);
			print TEMP "WT: $WT\n";
			print TEMP "MUT: $MUT\n";
		}
		else { # NGT - NTC
			$alt =~s/N//;
			$ref =~s/N//;
			my @alt = split ("", $alt);
			for (my $i = 0; $i < @alt; $i++) {
				$tmp[$position+$i] = $alt[$i];
			}
			$MUT  = join ("", @tmp);
			print TEMP "DELINS $chr, $start, $exonEnd, $ref, $alt, $stOffset, $endOffset, DISTANCE = $distance\tTYPE: $type\tPOSITION:$position\n";
			print TEMP "wt: $WT\t mut: $MUT\n";
		}
	}
	# For indels
	else {
		if ($alt !~/N/ || $ref !~/N/) { # Indel does not contain N's
         		print TEMP "NO N's\n";
			if ($alt =~ /\+/ || $Lalt > $Lref) { # insertion

				my $insertion = $alt;
				$insertion =~s/[+-]//;
				#$insertion = "";
				my @alt = split ("", $alt);
				my @ref = split ("", $ref);

				# Here we compare the ref with alt, and we obtain the inserted bases
				print TEMP "reference is $ref\n";
				for (my $i=0; $i<@ref; $i++) {
					if ($insertion eq $ref){
						$insertion = join ("", @alt);
						last;
					}
					if ($alt[$i] eq $ref[$i]) {
						shift @alt;
						$insertion = join ("", @alt);
					}
					else {
						shift @alt;
						$insertion = join ("", @alt);
					}
				}
				my $length = length($insertion);
				my @tmp    = split (//, $WT);
				my $position = length($WT) - abs$distance - $endOffset -1;

				my $first  = substr($WT, 0, $position);
				my $second = substr($WT, $position);
			
				$MUT = $first.$insertion.$second;		

				$MUT  = substr ($MUT,0, $SeqLength);
				$WT   = substr ($WT, 0, $SeqLength);
				print TEMP "INSERTION NO_Ns $chr, $start, $exonEnd, $ref, $alt, $stOffset, $endOffset, DISTANCE = $distance\tTYPE: $type\tPOSITION:$position\n";
				print TEMP "wt: $WT\t mut: $MUT\n";
			}
			elsif ($alt =~ /\-/ || $Lalt < $Lref) { # deletion
				my @ref = split ("", $ref);

				# Here we compare the ref with alt, and we obtain the deleted bases
				for (my $i=0; $i<@ref; $i++) {
					if ($ref[0] eq $alt){
						shift @ref;
						$deletion = join ("", @ref);
						last;
					}
					else {
						$deletion.= $ref[$i];
					}	
				}
				#print TEMP "WT: $WT\n";
				#print TEMP "BEFORE deletion: $deletion\n";
				my $todelete = length($WT) - length($deletion) - abs($distance) - $endOffset;
				#print TEMP "todelete $todelete\n";
				my $first    = substr($WT, 0, $todelete);
				#print TEMP "FIRST: $first\n";
				my $second   = substr($WT, $todelete + length($deletion));
				#print TEMP "SECOND: $second\n";
				my $deleted  = $second;
				$deleted =~s/$ref/$deletion/;


				#if ($strand eq '-' && $type eq '3') {
				#	$todelete = $endOffset + abs($distance) + length($deletion);
				#	$first  = substr($WT, length($deletion), -$todelete);
				#	$second = substr($WT, -$todelete-length($deletion), 0)
				
				#}

				$MUT = $first.$second;
				$WT  = substr($WT, length($deletion));
				print TEMP "DELETION NO_Ns $chr, $start, $exonEnd, $ref, $alt, $stOffset, $endOffset, DISTANCE = $distance\tTYPE: $type\tPOSITION:$position\n";
				print TEMP "wt: $WT\t mut: $MUT\n";
				print TEMP "DELETION: $deleted\n";
			}
		}
		elsif ($alt =~ /^N/ && $ref =~ /^N/) {  # Indel contains N's
         		print TEMP "HAS N\n";
			if ( $Lalt > $Lref) { # Insertion
				print TEMP "Contains N's\n";
				my $insertion = $alt;
				$insertion =~s/[+-N]//;
				$ref =~s/N//;
				if ( $ref eq '') {
					print TEMP "INSERTED:$insertion\n";
					my $length = length($insertion);
					my @tmp    = split (//, $WT);
		         	        my $position = $endOffset -1;
					print TEMP "position to insert: $position\n";

					my $first  = substr($WT, 0, $position);
					my $second = substr($WT, $position);
			
					$MUT = $first.$insertion.$second;		

					$MUT  = substr ($MUT,0, $SeqLength);
					$WT   = substr ($WT, 0, $SeqLength);
					print TEMP "INSERTION Ns $chr, $start, $exonEnd, $ref, $alt, $stOffset, $endOffset, DISTANCE = $distance\tTYPE: $type\tPOSITION:$position\n";
					print TEMP "wt: $WT\t mut: $MUT\n";
				}
				else { # Means delins type insertional
					print TEMP "INSERTED:$insertion\n";
					my $length = length($insertion);
					my @tmp    = split (//, $WT);
					print TEMP "position to insert: $position\n";

					my $first  = substr($WT, 0, $position);
					my $second = substr($WT, $position);
			
					$MUT = $first.$insertion.$second;		

					$MUT = substr ($MUT,0, $SeqLength);
					$WT  = substr ($WT, 0, $SeqLength);
					print TEMP "DELINS-INSERTIONAL Ns $chr, $start, $exonEnd, $ref, $alt, $stOffset, $endOffset, DISTANCE = $distance\tTYPE: $type\tPOSITION:$position\n";
					print TEMP "wt: $WT\t mut: $MUT\n";
				}
			}
			elsif ($Lalt < $Lref) { # deletion
				$ref =~s/N//;
		         	$position --;
				#print TEMP "Contains N's\n";
				if ( $ref eq ' ') {
					print TEMP "WT before modification: $WT\n";
					my $deletion = $ref;
					$deletion =~s/[+-N]//;
					$alt =~s/N//;
					print TEMP "DELETION:$deletion\n";
					my $todelete = length($deletion) + $stOffset;
					my $first    = substr($WT, 0, $stOffset-1);
					my $second   = substr($WT, $todelete);
					my $deleted  = $second;
					$deleted =~s/$ref/$deletion/;
					$MUT = $first.$deleted;
					my $limit = int($endOffset/2);
					$WT  = substr($WT, $limit );
				        $MUT = substr($MUT, $limit );
				        print TEMP "DELETION Ns $chr, $start, $exonEnd, $ref, $alt, $stOffset, $endOffset, DISTANCE = $distance\tTYPE: $type\n";	
					print TEMP "wt: $WT\t mut: $MUT\n";
				}
				else { 
					my $deletion = $ref;
					$deletion =~s/N//;
					print TEMP "WT before modification: $WT\n";
					print TEMP "position before modification: $position\n";
					$alt =~s/N//;
					print TEMP "alt:$alt\n";
					print TEMP "DELETION:$deletion\n";
					if ($type eq '5') {
						my $endoff = $position + length($deletion);
						my $first    = substr($WT, 0, $position+1);
						print TEMP "Position : $position\n";
						my $second   = substr($WT, $endoff);
						print TEMP "first = $first\t second = $second\n";
						print TEMP "alt is: $alt\n";
						$MUT = $first.$alt.$second;
						print TEMP "raw mut : $MUT\n";
				        	print TEMP "DELETION Ns $chr, $start, $exonEnd, $ref, $alt, $stOffset, $endOffset, DISTANCE = $distance\tTYPE: $type\n";	
						$WT = substr($WT, 0, $stOffset + $endOffset);
					}
					elsif ($type eq '3') {
						my $fpos = $position;
						$position  = $position - length($deletion);
						my $endoff = $position + length($deletion);
						print TEMP "endoff = $endoff\n";
						my $first    = substr($WT, 0, $position + 2);
						print TEMP "Position : $position\n";
						my $second   = substr($WT, $endoff + 2);
						print TEMP "first = $first\t second = $second\n";
						$MUT = $first.$alt.$second;
						print TEMP "raw mut : $MUT\n";
						my $limit = int($endOffset/2);
				        	print TEMP "DELETION Ns $chr, $start, $exonEnd, $ref, $alt, $stOffset, $endOffset, DISTANCE = $distance\tTYPE: $type\n";	

						$WT  = substr($WT, $limit + length($deletion) -1);
					}
					else {
						$position  = $position + 1;
						my $endoff = $position + length($deletion);
						my $first    = substr($WT, 0, $position);
						print TEMP "Position : $position\n";
						my $second   = substr($WT, $endoff);
						$MUT = $first.$alt.$second;
						print TEMP "raw mut : $MUT\n";
						my $limit = int($endOffset/2);
						$MUT = substr($MUT, $limit);	
						$WT  = substr($WT, $limit + length($deletion));
					
					}
				        print TEMP "DELINS-DELETIONAL Ns $chr, $start, $exonEnd, $ref, $alt, $stOffset, $endOffset, DISTANCE = $distance\tTYPE: $type\tPOSITION:$position\n";	
					print TEMP "wt: $WT\t mut: $MUT\n";
				}
			}
		}
	}
	if ($type eq '3') {	
		`echo "$WT\n" > $outputDir/3.wt`;
		`echo "$MUT\n" > $outputDir/3.mut`;
	}
	elsif ($type eq '5') {
		`echo "$WT\n" > $outputDir/5.wt`;
		`echo "$MUT\n" > $outputDir/5.mut`;
	}
	elsif ($type eq 'FS') {
		`echo "$WT\n" > $outputDir/wt.fsplice.fa`;
		`echo "$MUT\n" > $outputDir/mut.fsplice.fa`;
	}
	elsif ($type eq 'GS') {
		`echo ">WT\n$WT\n" > $outputDir/wt.GeneSplicer.fa`;
		`echo ">MUT\n$MUT\n" > $outputDir/mut.GeneSplicer.fa`;
	}
}

sub runFSPLICE {

 my ( $chr, $start, $exonEnd, $ref, $alt, $distance, $ss, $strand ) = @_;
	my $dist = 50;

	# For CNVs
	if (length($alt) > 50) {
		$dist = int(length ($alt)) + 50;
	}
	elsif (length($ref) > 50) {
		$dist = int(length ($ref)) +50;
	}
 	my $newStart = $start - $dist;
 	my $newEnd   = $start + $dist-1;
 	my $bed      = "Fsplice.bed";
	my ($stOffset, $endOffset, $type);
	if ($ss eq '3') {
		$stOffset  = $dist ;
		$endOffset = $dist ;
		$type = 'FS';
	}
	else {
		$stOffset  = $dist ;
		$endOffset = $dist ;
		$type = 'FS';
	}
        retrieveVariant ($chr, $start, $exonEnd, $ref, $alt, $stOffset, $endOffset, $type, $distance, $ss, $strand);
 	`$FSPLICE $dirname/bin/FSPLICE/Human.dat $outputDir/wt.fsplice.fa -seq:5 > $outputDir/fsplice.wt.out`;
 	`$FSPLICE $dirname/bin/FSPLICE/Human.dat $outputDir/mut.fsplice.fa -seq:5 > $outputDir/fsplice.mut.out`;
	my $wt =`grep -E 'Acceptor|Donor|P:|W:' $outputDir/fsplice.wt.out`;
	my $mt =`grep -E 'Acceptor|Donor|P:|W:' $outputDir/fsplice.mut.out`;
	chomp ($wt, $mt);
	my @WT  = split (/\n/, $wt);
	my @MUT = split (/\n/, $mt);
 return (\@WT, \@MUT);
}

sub runGeneSplicer {

 my ( $chr, $start, $exonEnd, $ref, $alt, $distance, $ss, $strand ) = @_;
	my $dist = 100;
	
	# For CNVs
	if (length($alt) > 100) {
		$dist = int(length ($alt)) +100;
	}
	elsif (length($ref) > 100) {
		$dist = int(length ($ref)) +100;
	}
 	my $newStart = $start - $dist;
 	my $newEnd   = $start + $dist;
 	my $bed      = "GeneSplicer.bed";
	my $type = 'GS';

        my ($stOffset, $endOffset);

	if ($ss eq '3') {
		$stOffset  = $dist;
		$endOffset = $dist;
	}
	else {
		$stOffset  = $dist;
		$endOffset = $dist;
	}
	#$alt = $strand eq '-' ? &reverse_complement($alt) : $alt;
        retrieveVariant ($chr, $start, $exonEnd, $ref, $alt, $stOffset, $endOffset, $type, $distance, $ss, $strand);

 	my $wt = `$GeneSplicer $outputDir/wt.GeneSplicer.fa $dirname/bin/GeneSplicer/human 2> /dev/null`;
 	my $mt = `$GeneSplicer $outputDir/mut.GeneSplicer.fa $dirname/bin/GeneSplicer/human 2> /dev/null`;

	chomp $wt;
	chomp $mt;
	
	$wt =~ s/\n/ /;
	$mt =~ s/\n/ /;
		
	print TEMP "rungs wt: $wt\t rungs mt: $mt\n"; 

 return ($wt, $mt);
}

sub reverse_complement {
 my $sequence = shift;
 my @seq = split ("", $sequence);
 my $rev = "";
 foreach my $nucleotide (@seq) {
	$rev .= 'N' if $nucleotide eq 'N';
	$rev .= 'A' if $nucleotide eq 'T';
  	$rev .= 'T' if $nucleotide eq 'A';
	$rev .= 'C' if $nucleotide eq 'G';
	$rev .= 'G' if $nucleotide eq 'C';
 }
 return($rev);
}


