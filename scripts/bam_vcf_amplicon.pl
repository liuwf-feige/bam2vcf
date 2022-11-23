#!/usr/bin/perl
#Author: Wanfei Liu <liuwanfei@caas.cn>
#Description: This program can filter bam file and identify variants.
#use strict;
#use warnings;
#use Text::NSP::Measures::2D::Fisher::left;
#use Statistics::Basic qw(:all);

my $version="1.1 version";
my $program=$0;
use Getopt::Long;
my %opts;

my $usage=<<USAGE; #******* Instruction of this program *********#

	Author: Wanfei Liu
	
	Usage: perl $program -b bam_file -e bed_region -s sequence_of_genome -i identity -r read_quality -l read_length -m map_quality -p paried_read -n insert_length -f min_frequency -o out_file.

	   -b: the absolute path of bam file.
	   -e: the absolute path of primers-bioinfo.tsv file.
	   -s: the absolute path of seq file.
	   -i: identity of read (default is 0.95).
	   -r: read quality threshold (default value is 30).
	   -l: read length threshold (default value is 40).
	   -m: map quality threshold (default value is 30).
	   -p: paried read (default values is 1).
	   -n: insert length (default values is 300).
	   -f: minimum frequency of alternative allele (default values is 0.01).
	   -o: the absolute path of output file.

USAGE

#Gather input
GetOptions(\%opts,"b:s","e:s","s:s","i:s","r:s","l:s","m:s","p:s","n:s","f:s","o:s");

#Verify input
if (!defined $opts{b} or !defined $opts{e} or !defined $opts{s} or !defined $opts{o}) {
	die "$usage\n";
}

my $bamfile=$opts{b};
my $bedfile=$opts{e};
my $seqfile=$opts{s};
my $outfile=$opts{o};

#Default parameters
my $identity = (defined $opts{i})?$opts{i}:0.95; #identity_threshold
my $rquality = (defined $opts{r})?$opts{r}:30; #read_quality_threshold
my $read_len = (defined $opts{l})?$opts{l}:40; #read_length_threshold
my $mquality = (defined $opts{m})?$opts{m}:30; #map_quality_threshold
my $readtype = (defined $opts{p})?$opts{p}:1; #paired_read
my $insert_length;#maximum_insert_length_threshold
if (defined $opts{n}) {
	$insert_length=$opts{n};
}elsif ($readtype==1) {
	#evaluate insert length based on the alignment file
	open (STDIN,"samtools view -h $bamfile |")||die("fail to open $bamfile.\n");
	my %insert_start=();
	while (<STDIN>) {
		chomp;
		next if (/^\@/);
		my @list=split /\t/,$_;
		my $read_length=length($list[9]);
		#filter reads by read length
		if ($read_length<$read_len) {
			next;
		}
		#filter reads by map quality
		if ($list[4] eq "255" or $list[4]<$mquality) {
			next;
		}
		#filter read according to the mapping quality marked in flag
		my $bin=unpack("B32",pack("N",$list[1]));
		$bin=sprintf "%012d",$bin;
		my @bin=split //,$bin;
		##each segment properly aligned according to the aligner
		#if ($bin[10] eq "0") {
		#	next;
		#}

		##filter supplementary alignment
		if ($bin[0] eq "1") {
			next;
		}
		##filter quality controls
		if ($bin[2] eq "1") {
			next;
		}
		##filter secondary alignment
		if ($bin[3] eq "1") {
			next;
		}
		#filter reads by insert length
		if ($list[6] eq "=" and $list[8]!=0) {
			if ($list[3]<=$list[7]) {
				$insert_start{"$list[2]:$list[3]\_$list[7]"}=abs($list[8]);
			}else {
				$insert_start{"$list[2]:$list[7]\_$list[3]"}=abs($list[8]);
			}
		}
		last if (scalar(keys %insert_start)>=10000);
	}
	close STDIN;
	my $insert_total=0;
	my @insert_length=();
	foreach my $insert_key (keys %insert_start) {
		$insert_total+=$insert_start{$insert_key};
		push(@insert_length,$insert_start{$insert_key});
		delete $insert_start{$insert_key};
	}
	my $insert_ave=$insert_total/(scalar(@insert_length));
	my $insert_square=0;
	foreach my $insert_key (@insert_length) {
		$insert_square+=(abs($insert_key-$insert_ave))**2;
	}	
	my $insert_sd=sqrt($insert_square/(scalar(@insert_length)-1));
	$insert_length=sprintf "%.0f",$insert_ave+3*($insert_sd);
	%insert_start=();
	@insert_length=();
}
my $minifreq = (defined $opts{f})?$opts{f}:0.01; #minimum_frequency_for_alternative_allele

my $cmdline = "$0";
foreach my $opt (keys %opts) {
	$cmdline.=" -$opt $opts{$opt}";
}

my ($id,%ref,%len,%vcf,%forward,%reverse,%primer,);

my $sample="";
if ($bamfile=~/\//) {
	if ($bamfile=~/^.*\/([^\.]+)\..*$/) {
		$sample=$1;
	}
}else {
	if ($bamfile=~/^([^\.]+)\..*$/) {
		$sample=$1;
	}
}

open (LOG,">$outfile.log")||die("fail to open $outfile.log.\n");

my $time=localtime;
print LOG "Start time: $time.\n";
print LOG "$cmdline\n";
print LOG "Insert length is $insert_length\n" if ($opts{p}==1);

open (IN,"<$bedfile")||die("fail to open $bedfile.\n");
while(<IN>){
	chomp;
	##Chrom	fStart	fEnd	rStart	rEnd	fSeq	Rseq
	#NC_000962	6042	6067	6240	6264	GAACGCGATTCATAGCAGCATCGTG	CCGCAAGCTACTGAAGGACAAGGA	T1P1
	next if (/^\#/);
	my @list=split /\t/,$_;
	$forward{($list[1]+1)}{$list[2]}++;
	$reverse{$list[4]}{$list[3]+1}++;
	$primer{($list[1]+1)}{$list[4]}++;
	$primer{$list[4]}{$list[1]+1}++;
}
close IN;

open (IN,"<$seqfile")||die("fail to open $seqfile.\n");
while(<IN>){
	chomp;
	if (/^>(\S+)/) {
		$id=$1;
		next;
	}
	$ref{$id}.=uc($_);
}
close IN;

foreach my $chr (keys %ref) {
	$len{$chr}=length($ref{$chr});
}

open (OUT,">$outfile")||die("fail to open $outfile.\n");
print OUT "##fileformat=VCFv4.2\n";
print OUT "##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location\">\n";
print OUT "##FILTER=<ID=LowQual,Description=\"Low quality\">\n";
print OUT "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n";
print OUT "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads unsatisfied any filter condition are filtered)\">\n";
print OUT "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print OUT "##INFO=<ID=PDP,Number=1,Type=Integer,Description=\"Primary read depth; reads unsatisfied any read filter condition are filtered\">\n";
print OUT "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">\n";
print OUT "##contig=<ID=NC_000962,length=4411532>\n";
print OUT "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	$sample\n";

#read the alignment file
#my $samheader=&samheader($seqfile);
#print "$samheader";
open (STDIN,"samtools view -h $bamfile |")||die("fail to open $bamfile.\n");
my $chr;
my $depth;
my $genome_len;
my $pos;
my %var=();
my %depth=();
my %readfilter=();
my %read=();
while (<STDIN>) {
	chomp;
	#SRR027895.5259200_1     99      chr1    15443   255     76M     =       15623   256     CTCCAGAGGCCTCAGGTCCAGTCTCTAAAAATATCTCAGGAAGCTGCAGTGGCTGACCATTGCCTTGGACCGCTCT    BCBCCBBABCBBBBBB=BBAB?BBBBBAABBBBBBBBBBA@BBBABBBB?@AAA>?<702>7:;1:+670:9;<==    NM:i:1  NH:i:1
	## Bit    Description                                                Comment                                Value
    ## 0x1    template having multiple segments in sequencing            0: single-end 1: paired end            value: 2^^0 (  1)
    ## 0x2    each segment properly aligned according to the aligner     true only for paired-end alignments    value: 2^^1 (  2)
    ## 0x4    segment unmapped                                           ---                                           ---
    ## 0x8    next segment in the template unmapped                      ---                                           ---
    ## 0x10   SEQ being reverse complemented                                                                    value: 2^^4 ( 16)
    ## 0x20   SEQ of the next segment in the template being reversed                                            value: 2^^5 ( 32)
    ## 0x40   the first segment in the template                          read 1                                 value: 2^^6 ( 64)
    ## 0x80   the last segment in the template                           read 2                                 value: 2^^7 (128)
    ## 0x100  secondary alignment                                        ---                                           ---
    ## 0x200  not passing quality controls                               ---                                           ---
    ## 0x400  PCR or optical duplicate                                   ---                                           ---
	## 0x800  supplementary alignment

	if (/^\@/) {
		print "$_\n";
		next;
	}
	my @list=split /\t/,$_;
	next if (exists $readfilter{$list[0]});

	my ($read_filter,$ref_length,$num,$label,$fend,$rstart);
	if ($readtype==1) {
		($read_filter,$ref_length,$num,$label,$fend,$rstart)=&read_filter($_);
	}else {
		($read_filter,$ref_length,$num,$label,$fend,$rstart)=&read_filter2($_);
	}
	my @num=@{$num};
	my @label=@{$label};
	my @fend=@{$fend};
	my @rstart=@{$rstart};

	if ($read_filter==1) {
		$readfilter{$list[0]}++;
		next;
	}else {
		@{$read{$list[0]}{$_}}=($list[3],$ref_length,\@num,\@label,\@fend,\@rstart);
	}
	if (!defined $chr) {
		$chr=$list[2];
		$depth=0;
		$genome_len=10000;
		$pos=1;
		%var=();
		%depth=();
		%readfilter=();
		%read=();
	}elsif ($chr ne $list[2]) {
		my $start=$genome_len-10000+1;
		my $end=$len{$chr};
		foreach my $seqid (keys %read) {
			my @record=keys %{$read{$seqid}};
			if ($readtype==1) {
				if (scalar(@record)<2) {
					foreach my $record (@record) {
						@{$read{$seqid}{$record}}=();
						delete $read{$seqid}{$record};
					}
					delete $read{$seqid};
				}else {
					#filter unoverlaped reads
					my $map_length=0;my $fragment_length=0;
					foreach my $record (@record) {
						my @info=@{$read{$seqid}{$record}};
						my @read=split /\t/,$record;
						$map_length+=$info[1];
						$fragment_length=abs($read[8]);
					}
					if ($map_length>=$fragment_length) {
						foreach my $record (@record) {
							my @info=@{$read{$seqid}{$record}};
							my @read=split /\t/,$record;
							my @num=@{$info[2]};
							my @label=@{$info[3]};
							my @fend=@{$info[4]};
							my @rstart=@{$info[5]};
							if ($info[0]>=$start and $info[0]<=$end) {
								&var(\@label,\@num,\@fend,\@rstart,$read[9],$read[10],$read[3],$read[2]);
								print "$record\n";
							}
							@{$read{$seqid}{$record}}=();
							delete $read{$seqid}{$record};
						}
						delete $read{$seqid};
					}else {
						foreach my $record (@record) {
							@{$read{$seqid}{$record}}=();
							delete $read{$seqid}{$record};
						}
						delete $read{$seqid};
					}
				}
			}else {
				foreach my $record (@record) {
					my @info=@{$read{$seqid}{$record}};
					my @read=split /\t/,$record;
					my @num=@{$info[2]};
					my @label=@{$info[3]};
					my @fend=@{$info[4]};
					my @rstart=@{$info[5]};
					if ($info[0]>=$start and $info[0]<=$end) {
						&var(\@label,\@num,\@fend,\@rstart,$read[9],$read[10],$read[3],$read[2]);
						print "$record\n";
					}
					@{$read{$seqid}{$record}}=();
					delete $read{$seqid}{$record};
				}
				delete $read{$seqid};
			}
		}
		($depth,$pos)=&var_print($start,$end,$depth,$pos,1);
		$chr=$list[2];
		$depth=0;
		$genome_len=10000;
		$pos=1;
		%var=();
		%depth=();
		%readfilter=();
		%read=();
	}
	if ($list[3]-$insert_length+1>$genome_len) {
		my $start=$genome_len-10000+1;
		my $end=$genome_len;
		foreach my $seqid (keys %read) {
			my @record=keys %{$read{$seqid}};
			if ($readtype==1) {
				my $flag=0;
				if (scalar(@record)<2) {
					foreach my $record (@record) {
						my @info=@{$read{$seqid}{$record}};
						if ($info[0]<=$genome_len) {
							$flag=1;
							@{$read{$seqid}{$record}}=();
							delete $read{$seqid}{$record};
						}
					}
					delete $read{$seqid} if ($flag==1);
				}else {
					#filter unoverlaped reads
					my $map_length=0;my $fragment_length=0;
					foreach my $record (@record) {
						my @info=@{$read{$seqid}{$record}};
						my @read=split /\t/,$record;
						$map_length+=$info[1];
						$fragment_length=abs($read[8]);
					}
					if ($map_length>=$fragment_length) {
						foreach my $record (@record) {
							my @info=@{$read{$seqid}{$record}};
							my @read=split /\t/,$record;
							if ($info[0]<=$genome_len) {
								$flag++;
								my @num=@{$info[2]};
								my @label=@{$info[3]};
								my @fend=@{$info[4]};
								my @rstart=@{$info[5]};
								if ($info[0]>=$start and $info[0]<=$end) {
									&var(\@label,\@num,\@fend,\@rstart,$read[9],$read[10],$read[3],$read[2]);
									print "$record\n";
								}
							}
						}
						if ($flag==2) {
							foreach my $record (@record) {
								@{$read{$seqid}{$record}}=();
								delete $read{$seqid}{$record};
							}
							delete $read{$seqid};
						}
					}else {
						foreach my $record (@record) {
							my @info=@{$read{$seqid}{$record}};
							if ($info[0]<=$genome_len) {
								$flag++;
								@{$read{$seqid}{$record}}=();
								delete $read{$seqid}{$record};
							}
						}
						delete $read{$seqid} if ($flag==2);
					}
				}
			}else {
				foreach my $record (@record) {
					my @info=@{$read{$seqid}{$record}};
					my @read=split /\t/,$record;
					if ($info[0]<=$genome_len) {
						my @num=@{$info[2]};
						my @label=@{$info[3]};
						my @fend=@{$info[4]};
						my @rstart=@{$info[5]};
						if ($info[0]>=$start and $info[0]<=$end) {
							&var(\@label,\@num,\@fend,\@rstart,$read[9],$read[10],$read[3],$read[2]);
							print "$record\n";
						}
						@{$read{$seqid}{$record}}=();
						delete $read{$seqid}{$record};
					}
				}
				delete $read{$seqid};
			}
		}
		($depth,$pos)=&var_print($start,$end,$depth,$pos,0);
		$genome_len+=10000;
	}
}
my $start=$genome_len-10000+1;
my $end=$len{$chr};
foreach my $seqid (keys %read) {
	my @record=keys %{$read{$seqid}};
	if ($readtype==1) {
		if (scalar(@record)<2) {
			foreach my $record (@record) {
				@{$read{$seqid}{$record}}=();
				delete $read{$seqid}{$record};
			}
			delete $read{$seqid};
		}else {
			#filter unoverlaped reads
			my $map_length=0;my $fragment_length=0;
			foreach my $record (@record) {
				my @info=@{$read{$seqid}{$record}};
				my @read=split /\t/,$record;
				$map_length+=$info[1];
				$fragment_length=abs($read[8]);
			}
			if ($map_length>=$fragment_length) {
				foreach my $record (@record) {
					my @info=@{$read{$seqid}{$record}};
					my @read=split /\t/,$record;
					my @num=@{$info[2]};
					my @label=@{$info[3]};
					my @fend=@{$info[4]};
					my @rstart=@{$info[5]};
					if ($info[0]>=$start and $info[0]<=$end) {
						&var(\@label,\@num,\@fend,\@rstart,$read[9],$read[10],$read[3],$read[2]);
						print "$record\n";
					}
					@{$read{$seqid}{$record}}=();
					delete $read{$seqid}{$record};
				}
				delete $read{$seqid};
			}else {
				foreach my $record (@record) {
					@{$read{$seqid}{$record}}=();
					delete $read{$seqid}{$record};
				}
				delete $read{$seqid};
			}
		}
	}else {
		foreach my $record (@record) {
			my @info=@{$read{$seqid}{$record}};
			my @read=split /\t/,$record;
			my @num=@{$info[2]};
			my @label=@{$info[3]};
			my @fend=@{$info[4]};
			my @rstart=@{$info[5]};
			if ($info[0]>=$start and $info[0]<=$end) {
				&var(\@label,\@num,\@fend,\@rstart,$read[9],$read[10],$read[3],$read[2]);
				print "$record\n";
			}
			@{$read{$seqid}{$record}}=();
			delete $read{$seqid}{$record};
		}
		delete $read{$seqid};
	}
}
my ($depth,$pos)=&var_print($start,$end,$depth,$pos,1);
close STDIN;

$time=localtime;
print LOG "End time: $time.\n";
close LOG;

#&read_filter($record)
sub read_filter {
	$record=shift;

	my @num=();
	my @label=();
	my $indel=0;
	my $indel_gap=0;
	my $deletion=0;
	my $del_gap=0;
	my $insertion=0;
	my $ins_gap=0;
	my $alignlen=0;
	my $fstart=0;
	my $fend=0;
	my $rstart=0;
	my $rend=0;
	my @fstart=();
	my @fend=();
	my @rstart=();
	my @rend=();
	my $read_filter=0;
	my $ref_length=0;

	my @list=split /\t/,$record;
	my $read_length=length($list[9]);

	#filter reads by read length
	if ($read_length<$read_len) {
		$read_filter=1;
		return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
	}

	#filter reads by map quality
	if ($list[4] eq "255" or $list[4]<$mquality) {
		$read_filter=1;
		return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
	}
	
	#filter reads by average base quality
	my $qual=0;
	if ($list[10] ne "*") {
		foreach my $basequal (split //,$list[10]) {
			$qual+=ord($basequal)-33;
		}
		$qual=$qual/(length($list[10]));
		if ($qual<$rquality) {
			$read_filter=1;
			return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
		}
	}

	#filter reads by read type
	if ($readtype eq "1") {
		if ($list[6] eq "=") {
			if (abs($list[8])>$insert_length) {
				$read_filter=1;
				return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
			}
		}else {
			$read_filter=1;
			return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
		}
	}

	#filter read according to the mapping quality marked in flag
	my $bin=unpack("B32",pack("N",$list[1]));
	$bin=sprintf "%012d",$bin;
	my @bin=split //,$bin;

	##each segment properly aligned according to the aligner
	#if ($bin[10] eq "0") {
	#	$read_filter=1;
	#	return $read_filter;
	#}

	##filter supplementary alignment
	if ($bin[0] eq "1") {
		$read_filter=1;
		return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
	}

	#filter quality controls
	if ($bin[2] eq "1") {
		$read_filter=1;
		return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
	}

	#filter secondary alignment
	if ($bin[3] eq "1") {
		$read_filter=1;
		return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
	}

	#dissect Optional fields
	my %attr=();
	for (my $i=11;$i<@list;$i++) {
		if ($list[$i]=~/([^\:]+)\:([^\:]+)\:([^\:]+)/) {
			$attr{$1}=$3;
		}
	}

	#dissect cigar
	if ($list[5]=~/S/ or $list[5]=~/H/) {
		$read_filter=1;
		return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
	}
	my @record=($list[5]=~m/(\d+\D+)/g);
	foreach my $record (@record) {
		if ($record=~/(\d+)(\D+)/) {
			push (@num,$1);
			push (@label,$2);
			if ($2 eq "M") {
				$alignlen+=$1;
			}elsif ($2 eq "I") {
				$insertion+=$1;
				$ins_gap++;
				$indel+=$1;
				$indel_gap++;
			}elsif ($2 eq "D") {
				$deletion+=$1;
				$del_gap++;
				$indel+=$1;
				$indel_gap++;
			}
		}
	}
	$ref_length=$alignlen+$deletion;
	
	#filter mismatch
	my $mismatch=$attr{"NM"}-$insertion-$deletion;
	if ($mismatch/$alignlen>(1-$identity)) {
		$read_filter=1;
		return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
	}

	#filter gap
	if (($attr{"NM"})/($alignlen+$indel)>2*(1-$identity)) {
		$read_filter=1;
		return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
	}

	#filter amplicon
	if (!exists $forward{$list[3]} and !exists $reverse{$list[3]-1+$alignlen+$deletion}) {
		$read_filter=1;
		return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
	}

	my $primer_flag=0;
	if (exists $forward{$list[3]}) {
		$fstart=$list[3];
		push(@fstart,$fstart);
		@fend=sort {$a <=> $b} keys %{$forward{$fstart}};
		@rend=sort {$a <=> $b} keys %{$primer{$fstart}};
		foreach my $key (@rend) {
			push (@rstart,sort {$a <=> $b} keys %{$reverse{$key}});
			if ($list[3]-1+abs($list[8])==$key) {
				#filter mismatch and gap for rRNA
				if (exists $amp{$fstart."\t".$key}) {
					if ($mismatch/$alignlen<=(1-$identity)/2 and ($attr{"NM"})/($alignlen+$indel)<=(1-$identity)) {
						$primer_flag=1;
					}
				}else {
					$primer_flag=1;
				}
			}
		}
	}elsif (exists $reverse{$list[3]-1+$alignlen+$deletion}) {
		$rend=$list[3]-1+$alignlen+$deletion;
		push(@rend,$rend);
		@rstart=sort {$a <=> $b} keys %{$reverse{$rend}};
		@fstart=sort {$a <=> $b} keys %{$primer{$rend}};
		foreach my $key (@fstart) {
			push (@fend,sort {$a <=> $b} keys %{$forward{$key}});
			if ($rend-abs($list[8])+1==$key) {
				#filter mismatch and gap for rRNA
				if (exists $amp{$key."\t".$rend}) {
					if ($mismatch/$alignlen<=(1-$identity)/2 and ($attr{"NM"})/($alignlen+$indel)<=(1-$identity)) {
						$primer_flag=1;
					}
				}else {
					$primer_flag=1;
				}
			}
		}
	}
	if ($primer_flag==0) {
		$read_filter=1;
		return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
	}
	return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
}

#&read_filter2($record)
sub read_filter2 {
	$record=shift;

	my @num=();
	my @label=();
	my $indel=0;
	my $indel_gap=0;
	my $deletion=0;
	my $del_gap=0;
	my $insertion=0;
	my $ins_gap=0;
	my $alignlen=0;
	my $fstart=0;
	my $fend=0;
	my $rstart=0;
	my $rend=0;
	my @fstart=();
	my @fend=();
	my @rstart=();
	my @rend=();
	my $read_filter=0;
	my $ref_length=0;

	my @list=split /\t/,$record;
	my $read_length=length($list[9]);

	#filter reads by read length
	if ($read_length<$read_len) {
		$read_filter=1;
		return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
	}

	#filter reads by map quality
	if ($list[4] eq "255" or $list[4]<$mquality) {
		$read_filter=1;
		return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
	}
	
	#filter reads by average base quality
	my $qual=0;
	if ($list[10] ne "*") {
		foreach my $basequal (split //,$list[10]) {
			$qual+=ord($basequal)-33;
		}
		$qual=$qual/(length($list[10]));
		if ($qual<$rquality) {
			$read_filter=1;
			return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
		}
	}

	#filter reads by read type
	if ($readtype eq "1") {
		if ($list[6] eq "=") {
			if (abs($list[8])>$insert_length) {
				$read_filter=1;
				return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
			}
		}else {
			$read_filter=1;
			return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
		}
	}

	#filter read according to the mapping quality marked in flag
	my $bin=unpack("B32",pack("N",$list[1]));
	$bin=sprintf "%012d",$bin;
	my @bin=split //,$bin;

	##each segment properly aligned according to the aligner
	#if ($bin[10] eq "0") {
	#	$read_filter=1;
	#	return $read_filter;
	#}

	##filter supplementary alignment
	if ($bin[0] eq "1") {
		$read_filter=1;
		return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
	}

	#filter quality controls
	if ($bin[2] eq "1") {
		$read_filter=1;
		return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
	}

	#filter secondary alignment
	if ($bin[3] eq "1") {
		$read_filter=1;
		return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
	}

	#dissect Optional fields
	my %attr=();
	for (my $i=11;$i<@list;$i++) {
		if ($list[$i]=~/([^\:]+)\:([^\:]+)\:([^\:]+)/) {
			$attr{$1}=$3;
		}
	}

	#dissect cigar
	if ($list[5]=~/S/ or $list[5]=~/H/) {
		$read_filter=1;
		return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
	}
	my @record=($list[5]=~m/(\d+\D+)/g);
	foreach my $record (@record) {
		if ($record=~/(\d+)(\D+)/) {
			push (@num,$1);
			push (@label,$2);
			if ($2 eq "M") {
				$alignlen+=$1;
			}elsif ($2 eq "I") {
				$insertion+=$1;
				$ins_gap++;
				$indel+=$1;
				$indel_gap++;
			}elsif ($2 eq "D") {
				$deletion+=$1;
				$del_gap++;
				$indel+=$1;
				$indel_gap++;
			}
		}
	}
	$ref_length=$alignlen+$deletion;
	
	#filter mismatch
	my $mismatch=$attr{"NM"}-$insertion-$deletion;
	if ($mismatch/$alignlen>(1-$identity)) {
		$read_filter=1;
		return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
	}

	#filter gap
	if (($attr{"NM"})/($alignlen+$indel)>2*(1-$identity)) {
		$read_filter=1;
		return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
	}

	#filter amplicon
	if (!exists $forward{$list[3]} or !exists $reverse{$list[3]-1+$alignlen+$deletion}) {
		$read_filter=1;
		return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
	}

	my $primer_flag=0;
	if (exists $forward{$list[3]} and exists $reverse{$list[3]-1+$alignlen+$deletion}) {
		$fstart=$list[3];
		push(@fstart,$fstart);
		@fend=sort {$a <=> $b} keys %{$forward{$fstart}};
		$rend=$list[3]-1+$alignlen+$deletion;
		push(@rend,$rend);
		@rstart=sort {$a <=> $b} keys %{$reverse{$rend}};
		#filter mismatch and gap for rRNA
		if (exists $amp{$fstart."\t".$rend}) {
			if ($mismatch/$alignlen<=(1-$identity)/2 and ($attr{"NM"})/($alignlen+$indel)<=(1-$identity)) {
				$primer_flag=1;
			}
		}else {
			$primer_flag=1;
		}
	}
	if ($primer_flag==0) {
		$read_filter=1;
		return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
	}
	return ($read_filter,$ref_length,\@num,\@label,\@fend,\@rstart);
}

#&var(\@label,\@num,\@fend,\@rstart,$seq,$qual,$pos,$chr)
sub var {
	my $label=shift;
	my @label=@{$label};
	my $num=shift;
	my @num=@{$num};
	my $fend=shift;
	my @fend=@{$fend};
	my $rstart=shift;
	my @rstart=@{$rstart};
	my $seq=shift;
	my $qual=shift;
	my $pos=shift;
	my $chr=shift;
	
	my $clip=0;
	my $startr=0;
	my $startv=0;
	for (my $i=0;$i<@label;$i++) {
		if ($label[$i] eq "S") {
			$clip=$num[$i];
		}elsif ($label[$i] eq "M") {
			for (my $j=0;$j<$num[$i];$j++) {
				my $alt=substr($seq,$clip+$startv,1);
				my $var_flag=0;
				for (my $m=0;$m<@fend;$m++) {
					for (my $n=0;$n<@rstart;$n++) {
						if ($pos+$startr>$fend[$m] and $pos+$startr<$rstart[$n]) {
							$var_flag=1;
						}
						last if ($var_flag==1);
					}
					last if ($var_flag==1);
				}
				$var_flag=0 if ((&avequal(substr($qual,$clip+$startv,1)))<$rquality);
				$var{$pos+$startr}{$alt}++ if ($var_flag==1);
				$depth{$pos+$startr}++;
				$startr++;
				$startv++;
			}
		}elsif ($label[$i] eq "N") {
			$startr+=$num[$i];
		}elsif ($label[$i] eq "I") {
			my $alt=substr($seq,$clip+$startv,$num[$i]);
			my $var_flag=0;
			for (my $m=0;$m<@fend;$m++) {
				for (my $n=0;$n<@rstart;$n++) {
					if ($pos+$startr>$fend[$m] and $pos+$startr<$rstart[$n]) {
						$var_flag=1;
					}
					last if ($var_flag==1);
				}
				last if ($var_flag==1);
			}
			$var_flag=0 if ((&avequal(substr($qual,$clip+$startv,$num[$i])))<$rquality);
			$var{$pos+$startr}{"I:$alt"}++ if ($var_flag==1);
			$startv+=$num[$i];
		}elsif ($label[$i] eq "D") {
			my $alt=substr($ref{$chr},$pos-1+$startr,$num[$i]);
			my $var_flag=0;
			for (my $m=0;$m<@fend;$m++) {
				for (my $n=0;$n<@rstart;$n++) {
					if ($pos+$startr>$fend[$m] and $pos+$startr<$rstart[$n]) {
						$var_flag=1;
					}
					last if ($var_flag==1);
				}
				last if ($var_flag==1);
			}
			$var_flag=0 if ((&avequal(substr($qual,$clip+$startv-1,2)))<$rquality);
			$var{$pos+$startr}{"D:$alt"}++ if ($var_flag==1);
			foreach my $posD (($pos+$startr)..($pos+$startr+$num[$i]-1)) {
				$depth{$posD}++;
			}
			$startr+=$num[$i];
		}
	}
}

#&var_print($start,$end,$depth,0|1);#染色体末端:1; 染色体中间:0
sub var_print {
	my $start=shift;
	my $end=shift;
	my $depth=shift;
	my $pos=shift;
	my $region_type=shift;

	my $GT="";
	for (my $loci=$start;$loci<=$end;$loci++) {
		if (!exists $var{$loci}) {
			if (!exists $depth{$loci}) {
				if ($depth!=0) {
					if ($loci-1>=$pos) {
						$GT="0/0";
						print OUT "$chr\t$pos\t\.\t".(uc(substr($ref{$chr},$pos-1,1)))."\t<NON_REF>\t\.\t\.\tEND=".($loci-1).";PDP=$depth\tGT:AD:DP\t$GT:$depth,0:$depth\n";
					}
					$pos=$loci;
					$depth=0;
				}
			}else {
				if ($depth{$loci}!=$depth) {
					if ($loci-1>=$pos) {
						$GT="0/0";
						print OUT "$chr\t$pos\t\.\t".(uc(substr($ref{$chr},$pos-1,1)))."\t<NON_REF>\t\.\t\.\tEND=".($loci-1).";PDP=$depth\tGT:AD:DP\t$GT:$depth,0:$depth\n";
					}
					$pos=$loci;
					$depth=$depth{$loci};
				}
			}
			next;
		}
		my $ref=uc(substr($ref{$chr},$loci-1,1));
		if (!exists $var{$loci}{$ref}) {
			$var{$loci}{$ref}=0;
		}
		my @allele=sort {$var{$loci}{$b}<=>$var{$loci}{$a}} keys %{$var{$loci}};
		my $type=-1;
		my $D="";
		my $total_indel=0;
		my $total_snp=0;
		foreach my $allele (@allele) {
			next if ($allele eq $ref or $var{$loci}{$allele}<3 or $var{$loci}{$allele}/$depth{$loci}<$minifreq);
			if ($allele=~/^D\:(.+)$/) {
				$total_indel+=$var{$loci}{$allele};
				$type=2 if ($type<2);
				if (length($1)>length($D)) {
					$D=$1;
				}
			}elsif ($allele=~/^I\:(.+)$/) {
				$total_indel+=$var{$loci}{$allele};
				$type=1 if ($type<1);
			}else {
				$total_snp+=$var{$loci}{$allele};
				$type=0 if ($type<0);
			}
		}
		$total_snp+=$var{$loci}{$ref} if ($var{$loci}{$ref}>=3 and $var{$loci}{$ref}/$depth{$loci}>=$minifreq);
		my @AD=();
		my @ALT=();
		my $REF="";
		my $GT="";
		my $dp=0;
		if ($type==2) {
			my $type2_flag=0;
			$REF=(uc(substr($ref{$chr},$loci-2,1))).$D;
			@AD=();
			@ALT=();
			$dp=0;
			if ($depth{$loci-1}-$total_indel<3 or ($depth{$loci-1}-$total_indel)/$depth{$loci-1}<$minifreq) {
				push(@AD,0);
			}else {
				push(@AD,$depth{$loci-1}-$total_indel);
				$dp+=$depth{$loci-1}-$total_indel;
			}
			foreach my $allele (@allele) {
				next if (length($allele)==1 or $var{$loci}{$allele}<3 or $var{$loci}{$allele}/$depth{$loci-1}<$minifreq);
				if ($allele=~/^D\:(.+)$/) {
					push(@AD,$var{$loci}{$allele});
					push(@ALT,(uc(substr($ref{$chr},$loci-2,1))).(substr($D,length($1),length($D)-length($1))));
					$dp+=$var{$loci}{$allele};
				}elsif ($allele=~/^I\:(.+)$/) {
					push(@AD,$var{$loci}{$allele});
					push(@ALT,(uc(substr($ref{$chr},$loci-2,1))).$1.$D);
					$dp+=$var{$loci}{$allele};
				}
			}
			if (scalar(@ALT)>0) {
				$type2_flag=1;
				if ($loci-2>=$pos) {
					$GT="0/0";
					print OUT "$chr\t$pos\t\.\t".(uc(substr($ref{$chr},$pos-1,1)))."\t<NON_REF>\t\.\t\.\tEND=".($loci-2).";PDP=$depth\tGT:AD:DP\t$GT:$depth,0:$depth\n";
				}
				if ($depth{$loci-1}-$total_indel>=3 and ($depth{$loci-1}-$total_indel)/$depth{$loci-1}>=$minifreq) {
					$GT="0/1";
				}else {
					$GT="1/1";
				}
				print OUT "$chr\t".($loci-1)."\t\.\t$REF\t".(join ",",@ALT)."\t\.\t\.\tPDP=".($depth{$loci-1}).";DP=$dp\tGT:AD:DP\t$GT:".(join ",",@AD).":$dp\n";
				$pos=$loci;
				$depth=$depth{$loci};
			}
			@AD=();
			@ALT=();
			$dp=0;
			if ($var{$loci}{$ref}<3 or ($var{$loci}{$ref})/$depth{$loci}<$minifreq) {
				push(@AD,0);
			}else {
				push(@AD,$var{$loci}{$ref});
				$dp+=$var{$loci}{$ref};
			}
			foreach my $allele (@allele) {
				next if (length($allele)!=1 or $allele eq $ref or $var{$loci}{$allele}<3 or $var{$loci}{$allele}/$depth{$loci}<$minifreq);
				push(@AD,$var{$loci}{$allele});
				push(@ALT,$allele);
				$dp+=$var{$loci}{$allele};
			}
			if (scalar(@ALT)>0) {
				$type2_flag=1;
				if ($loci-1>=$pos) {
					$GT="0/0";
					print OUT "$chr\t$pos\t\.\t".(uc(substr($ref{$chr},$pos-1,1)))."\t<NON_REF>\t\.\t\.\tEND=".($loci-1).";PDP=$depth\tGT:AD:DP\t$GT:$depth,0:$depth\n";
				}
				if ($var{$loci}{$ref}>=3 and $var{$loci}{$ref}/$depth{$loci}>=$minifreq) {
					$GT="0/1";
				}else {
					$GT="1/1";
				}
				print OUT "$chr\t$loci\t\.\t$ref\t".(join ",",@ALT)."\t\.\t\.\tPDP=$depth{$loci};DP=$dp\tGT:AD:DP\t$GT:".(join ",",@AD).":$dp\n";
				$pos=$loci+1;
				$depth=$depth{$loci+1};
			}
			if ($type2_flag==0) {
				if ($depth{$loci}!=$depth) {
					if ($loci-1>=$pos) {
						$GT="0/0";
						print OUT "$chr\t$pos\t\.\t".(uc(substr($ref{$chr},$pos-1,1)))."\t<NON_REF>\t\.\t\.\tEND=".($loci-1).";PDP=$depth\tGT:AD:DP\t$GT:$depth,0:$depth\n";
					}
					$pos=$loci;
					$depth=$depth{$loci};
				}
			}
		}elsif ($type==1) {
			my $type1_flag=0;
			$REF=uc(substr($ref{$chr},$loci-2,1));
			@AD=();
			@ALT=();
			$dp=0;
			if ($depth{$loci-1}-$total_indel<3 or ($depth{$loci-1}-$total_indel)/$depth{$loci-1}<$minifreq) {
				push(@AD,0);
			}else {
				push(@AD,$depth{$loci-1}-$total_indel);
				$dp+=$depth{$loci-1}-$total_indel;
			}
			foreach my $allele (@allele) {
				next if (length($allele)==1 or $var{$loci}{$allele}<3 or $var{$loci}{$allele}/$depth{$loci-1}<$minifreq);
				if ($allele=~/^I\:(.+)$/) {
					push(@AD,$var{$loci}{$allele});
					push(@ALT,$REF.$1);
					$dp+=$var{$loci}{$allele};
				}
			}
			if (scalar(@ALT)>0) {
				$type1_flag=1;
				if ($loci-2>=$pos) {
					$GT="0/0";
					print OUT "$chr\t$pos\t\.\t".(uc(substr($ref{$chr},$pos-1,1)))."\t<NON_REF>\t\.\t\.\tEND=".($loci-2).";PDP=$depth\tGT:AD:DP\t$GT:$depth,0:$depth\n";
				}
				if ($depth{$loci-1}-$total_indel>=3 and ($depth{$loci-1}-$total_indel)/$depth{$loci-1}>=$minifreq) {
					$GT="0/1";
				}else {
					$GT="1/1";
				}
				print OUT "$chr\t".($loci-1)."\t\.\t$REF\t".(join ",",@ALT)."\t\.\t\.\tPDP=".($depth{$loci-1}).";DP=$dp\tGT:AD:DP\t$GT:".(join ",",@AD).":$dp\n";
				$pos=$loci;
				$depth=$depth{$loci};
			}
			@AD=();
			@ALT=();
			$dp=0;
			if ($var{$loci}{$ref}<3 or ($var{$loci}{$ref})/$depth{$loci}<$minifreq) {
				push(@AD,0);
			}else {
				push(@AD,$var{$loci}{$ref});
				$dp+=$var{$loci}{$ref};
			}
			foreach my $allele (@allele) {
				next if (length($allele)!=1 or $allele eq $ref or $var{$loci}{$allele}<3 or $var{$loci}{$allele}/$depth{$loci}<$minifreq);
				push(@AD,$var{$loci}{$allele});
				push(@ALT,$allele);
				$dp+=$var{$loci}{$allele};
			}
			if (scalar(@ALT)>0) {
				$type1_flag=1;
				if ($loci-1>=$pos) {
					$GT="0/0";
					print OUT "$chr\t$pos\t\.\t".(uc(substr($ref{$chr},$pos-1,1)))."\t<NON_REF>\t\.\t\.\tEND=".($loci-1).";PDP=$depth\tGT:AD:DP\t$GT:$depth,0:$depth\n";
				}
				if ($var{$loci}{$ref}>=3 and $var{$loci}{$ref}/$depth{$loci}>=$minifreq) {
					$GT="0/1";
				}else {
					$GT="1/1";
				}
				print OUT "$chr\t$loci\t\.\t$ref\t".(join ",",@ALT)."\t\.\t\.\tPDP=$depth{$loci};DP=$dp\tGT:AD:DP\t$GT:".(join ",",@AD).":$dp\n";
				$pos=$loci+1;
				$depth=$depth{$loci+1};
			}
			if ($type1_flag==0) {
				if ($depth{$loci}!=$depth) {
					if ($loci-1>=$pos) {
						$GT="0/0";
						print OUT "$chr\t$pos\t\.\t".(uc(substr($ref{$chr},$pos-1,1)))."\t<NON_REF>\t\.\t\.\tEND=".($loci-1).";PDP=$depth\tGT:AD:DP\t$GT:$depth,0:$depth\n";
					}
					$pos=$loci;
					$depth=$depth{$loci};
				}
			}
		}elsif ($type==0) {
			my $type0_flag=0;
			$REF=$ref;
			@AD=();
			@ALT=();
			$dp=0;
			if ($var{$loci}{$ref}<3 or ($var{$loci}{$ref})/$depth{$loci}<$minifreq) {
				push(@AD,0);
			}else {
				push(@AD,$var{$loci}{$ref});
				$dp+=$var{$loci}{$ref};
			}
			foreach my $allele (@allele) {
				next if ($allele eq $ref or $var{$loci}{$allele}<3 or $var{$loci}{$allele}/$depth{$loci}<$minifreq);
				push(@AD,$var{$loci}{$allele});
				push(@ALT,$allele);
				$dp+=$var{$loci}{$allele};
			}
			if (scalar(@ALT)>0) {
				$type0_flag=1;
				if ($loci-1>=$pos) {
					$GT="0/0";
					print OUT "$chr\t$pos\t\.\t".(uc(substr($ref{$chr},$pos-1,1)))."\t<NON_REF>\t\.\t\.\tEND=".($loci-1).";PDP=$depth\tGT:AD:DP\t$GT:$depth,0:$depth\n";
				}
				if ($var{$loci}{$ref}>=3 and $var{$loci}{$ref}/$depth{$loci}>=$minifreq) {
					$GT="0/1";
				}else {
					$GT="1/1";
				}
				print OUT "$chr\t$loci\t\.\t$ref\t".(join ",",@ALT)."\t\.\t\.\tPDP=$depth{$loci};DP=$dp\tGT:AD:DP\t$GT:".(join ",",@AD).":$dp\n";
				$pos=$loci+1;
				$depth=$depth{$loci+1};
			}
			if ($type0_flag==0) {
				if ($depth{$loci}!=$depth) {
					if ($loci-1>=$pos) {
						$GT="0/0";
						print OUT "$chr\t$pos\t\.\t".(uc(substr($ref{$chr},$pos-1,1)))."\t<NON_REF>\t\.\t\.\tEND=".($loci-1).";PDP=$depth\tGT:AD:DP\t$GT:$depth,0:$depth\n";
					}
					$pos=$loci;
					$depth=$depth{$loci};
				}
			}
		}elsif ($type==-1) {
			if ($depth{$loci}!=$depth) {
				if ($loci-1>=$pos) {
					$GT="0/0";
					print OUT "$chr\t$pos\t\.\t".(uc(substr($ref{$chr},$pos-1,1)))."\t<NON_REF>\t\.\t\.\tEND=".($loci-1).";PDP=$depth\tGT:AD:DP\t$GT:$depth,0:$depth\n";
				}
				$pos=$loci;
				$depth=$depth{$loci};
			}
		}
	}
	if ($region_type==1 and $end>=$pos) {
		$GT="0/0";
		print OUT "$chr\t$pos\t\.\t".(uc(substr($ref{$chr},$pos-1,1)))."\t<NON_REF>\t\.\t\.\tEND=".($end).";PDP=$depth\tGT:AD:DP\t$GT:$depth,0:$depth\n";
	}
	for (my $loci=$start;$loci<=$end;$loci++) {
		if (exists $var{$loci}) {
			foreach my $allele (keys %{$var{$loci}}) {
				delete $var{$loci}{$allele};
			}
			delete $var{$loci};
			delete $depth{$loci};
		}
	}
	return ($depth,$pos);
}

sub log10 {
	my $n = shift;
	return log($n)/log(10);
}

sub samheader {
	my $file=shift;

	my $samheader="";
	$samheader.="\@HD	VN:1.0	SO:sorted\n";
	$samheader.="\@PG	ID:bam_vcf.pl	PN:bam_vcf.pl	VN:$version	CL:$cmdline\n";
	open (IN,"<$file")||die("fail to open $file.\n");
	my $id=undef;
	my $str=undef;
	while (<IN>) {
		chomp;
		if (/^>(\S+)/) {
			if (defined $id) {
				$samheader.="\@SQ	SN:$id	LN:".($len{$id})."\n";
			}
			$id=$1;
		}
	}
	close IN;
	$samheader.="\@SQ	SN:$id	LN:".($len{$id})."\n";
	return $samheader;
}

#filter reads by average base quality
sub avequal {
	my $qualstr=shift;
	
	my $qual=0;
	foreach my $basequal (split //,$qualstr) {
		$qual+=ord($basequal)-33;
	}
	$qual=$qual/(length($qualstr));
	return $qual;
}

#perl bam_vcf.pl -b bam_file -s genome_sequence_file -o bam.vcf | samtools view -b - > bam_filter_file &
