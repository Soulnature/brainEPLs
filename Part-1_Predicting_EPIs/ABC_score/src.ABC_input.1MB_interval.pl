use strict;

open INPUT1,"../../results/promoter_region/promoter_region.bed" or die;

open OUTPUT1,">../../results/ABC_input.1MB_interval/ABC_input.1MB_interval.bed" or die;

while(<INPUT1>){
	chomp;
	my @line = split /\t/, $_;
	my $TSS_position = ($line[1]+$line[2])/2;
	my $start = $TSS_position - 1000000;
	if($start<0){
		$start = 0;
	}
	my $end = $TSS_position + 1000000;
	print OUTPUT1 "$line[0]\t$start\t$end\t$line[3]\t$line[4]\t$line[5]\n";
}


`intersectBed -wa -wb -a ../../results/ABC_input.1MB_interval/ABC_input.1MB_interval.bed -b ../../results/enhancer_region/enhancer_region.bed | cut -f 4,10 > ../../results/ABC_input.1MB_interval/ABC_input.1MB_interval.intersect_enhancer.txt`;


open INPUT2,"../../results/ABC_input.1MB_interval/ABC_input.1MB_interval.intersect_enhancer.txt" or die;

my %promoter_id_2_enhancer_id;

while(<INPUT2>){
	chomp;
	my @line = split /\t/, $_;
	if(exists $promoter_id_2_enhancer_id{$line[0]}){
		$promoter_id_2_enhancer_id{$line[0]} = $promoter_id_2_enhancer_id{$line[0]}.";".$line[1];
	}
	else{
		$promoter_id_2_enhancer_id{$line[0]} = $line[1];
	}
}

open OUTPUT2,">../../results/ABC_input.1MB_interval/ABC_input.1MB_interval.promoter_enhancer_pairs.txt" or die;

foreach my $promoter_id (sort keys %promoter_id_2_enhancer_id) {
	my $enhancer_id = $promoter_id_2_enhancer_id{$promoter_id};
	print OUTPUT2 "$promoter_id\t$enhancer_id\n";
}


unlink("../../results/ABC_input.1MB_interval/ABC_input.1MB_interval.bed");
unlink("../../results/ABC_input.1MB_interval/ABC_input.1MB_interval.intersect_enhancer.txt");

