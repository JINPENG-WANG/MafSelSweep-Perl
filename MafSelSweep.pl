#!/usr/bin/perl -w
use strict;
use IO::File;
my %parameters=("-infile" => undef,
				"-maf_threshold" => 0.05,
				"-snp_number" => 5,
				"-false_count_ratio" => 0.1,
				"-outfile" => undef );
if(@ARGV%2!=0){
	die "Wrong parameters!\n";
}

while(@ARGV){
	my $p1=shift @ARGV;
	my $p2=shift @ARGV;
	$parameters{$p1}=$p2;
}

my $infile=$parameters{-infile};
my $outfile=$parameters{-outfile};
my $maf_threshold=$parameters{-maf_threshold};
my $snp_number=$parameters{-snp_number};
my $false_count_ratio=$parameters{-false_count_ratio};
if($false_count_ratio=~/(\d+)\/(\d+)/){
	$false_count_ratio=$1/$2;
}
my $fh1=IO::File->new("$infile",'r');
my $fh2=IO::File->new(">$outfile");

my %maf;
my $line_count1=0;
my $header;
while(<$fh1>){
	chomp;
	my $line=$_;
	$line_count1++;
	if($line_count1==1){
		$header=$line;
	}else{
		my @eles=split /\s+/, $line;
		my ($chr,$maf,$pos)=@eles[0,4,6];
		push @{$maf{$chr}{maf}}, $maf;
		push @{$maf{$chr}{pos}}, $pos;
	}
}
$fh2->print("$header\tChromosome\tStart\tEnd\tBlock_size\tSNP_number\n");
my %loc_to_block;
my %block_to_loc;
for my $chr (keys %maf){
	my $string="";
	my $array_count=@{$maf{$chr}{maf}};
	for my $index (0..($array_count-1)){
		my $maf=${$maf{$chr}{maf}}[$index];
		my $tag;
		if($maf eq 'NA'){
			$tag="F";
		}
		elsif($maf<=$maf_threshold){
			$tag="T";
		}else{
			$tag="F";
		}
		$string.=$tag;
	}
	#my %loc_to_block;
	#my %block_to_loc;
	my $block_count=0;
	#my %final_block;
	while($string=~/(T+)/g){
		my $sub_string=$1;
		my $sub_end=pos($string);
		my $sub_length=length($sub_string);
		my $sub_start=$sub_end-$sub_length+1;
		$block_count++;
		for my $do ($sub_start..$sub_end){
			$loc_to_block{$chr}{$do}=$block_count;
		}
		$block_to_loc{$chr}{$block_count}{start}=$sub_start;
		$block_to_loc{$chr}{$block_count}{end}=$sub_end;
		$block_to_loc{$chr}{$block_count}{Fcount}=0;
		$block_to_loc{$chr}{$block_count}{length}=$sub_end-$sub_start+1;
	}
	if($false_count_ratio != 0){
		while($string=~/TFT/g){
			$block_count++;
			my $later_block_start=pos($string);
			my $before_block_end=pos($string)-2;
			my $before_block_count=$loc_to_block{$chr}{$before_block_end};
			my $later_block_count=$loc_to_block{$chr}{$later_block_start};
			$block_to_loc{$chr}{$before_block_count}{Fcount}++;
			my $before_block_start=$block_to_loc{$chr}{$before_block_count}{start};
			my $later_block_end=$block_to_loc{$chr}{$later_block_count}{end};
			$block_to_loc{$chr}{$block_count}{start}=$before_block_start;
			$block_to_loc{$chr}{$block_count}{end}=$later_block_end;
			$block_to_loc{$chr}{$block_count}{length}=$later_block_end-$before_block_start+1;
			$block_to_loc{$chr}{$block_count}{Fcount}=$block_to_loc{$chr}{$before_block_count}{Fcount};
			if(($block_to_loc{$chr}{$block_count}{Fcount}/$block_to_loc{$chr}{$block_count}{length})<=$false_count_ratio){
				for my $do ($before_block_start..$later_block_end){
					$loc_to_block{$chr}{$do}=$block_count;
				}
			}
		}
	}
}

my %lito;
my %wert;
	$fh1=IO::File->new("$infile",'r');
	$line_count1=0;
	my $line_count2=0;
	while(<$fh1>){
		chomp;
		my $line=$_;
		$line_count1++;
		if($line_count1==1){
			$header=$line;
		}else{
			$line_count2++;
			my @eles=split /\s+/, $line;
			my $chr=$eles[0];
			push @{$lito{$chr}}, $chr;
			my $number1=@{$lito{$chr}};

			if($loc_to_block{$chr}{$number1}){
				my $block_count=$loc_to_block{$chr}{$number1};
				my $start=$block_to_loc{$chr}{$block_count}{start};
				my $end=$block_to_loc{$chr}{$block_count}{end};
				my $block_range=$end-$start+1;
				if($block_range>=$snp_number){
					my $start_index=$start-1;
					my $end_index=$end-1;
					my $start_pos=${$maf{$chr}{pos}}[$start_index];
					my $end_pos=${$maf{$chr}{pos}}[$end_index];
					my $block_size=$end_pos-$start_pos+1;
					my $unique_flag="chr$chr-$start_pos";
					unless($wert{$unique_flag}){
						$fh2->print("$line\tchr$chr\t$start_pos\t$end_pos\t$block_size\t$block_range\n");
						$wert{$unique_flag}=1;
					}
					else{
						$fh2->print("$line\tchr$chr\t-\t-\t-\t-\n");
					} 
					
				}
				else{
					$fh2->print("$line\t-\t-\t-\t-\t-\n");
				}
			}else{
				$fh2->print("$line\t-\t-\t-\t-\t-\n");
			}
		}
	}
