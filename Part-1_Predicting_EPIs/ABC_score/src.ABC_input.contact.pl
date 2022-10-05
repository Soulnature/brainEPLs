use strict;

open INPUT1,"../../results/enhancer_region/enhancer_region.bed" or die;

my %enhancer_id_2_position;

while(<INPUT1>){
	chomp;
	my @line = split /\t/, $_;
	if($line[0] eq $ARGV[0]){
		my $midpoint = int(($line[1]+$line[2])/2);
		$enhancer_id_2_position{$line[3]} = $midpoint;
	}
}

open INPUT2,"../../results/promoter_region/promoter_region.bed" or die;

my %promoter_id_2_chrom;

my %promoter_id_2_position;

while(<INPUT2>){
	chomp;
	my @line = split /\t/, $_;
	if($line[0] eq $ARGV[0]){
		$promoter_id_2_chrom{$line[3]} = $line[0];
		my $midpoint = int(($line[1]+$line[2])/2);
		$promoter_id_2_position{$line[3]} = $midpoint;
	}
}


open INPUT3,"../../results/ABC_input.5MB_interval/ABC_input.5MB_interval.promoter_enhancer_pairs.txt" or die;

open OUTPUT3,">../../results/ABC_input.contact/$ARGV[0].txt" or die;

my $gamma = 0.87;
my $scale = -4.80 + 11.63 * $gamma;

while(<INPUT3>){
	chomp;
	my @line = split /\t/, $_;
	my $chrom = $promoter_id_2_chrom{$line[0]};
	if($chrom eq $ARGV[0]){
		if($line[1] =~ /;/){
			my $promoter_id = $line[0];
			my $promoter_position = $promoter_id_2_position{$promoter_id};
			my @enhancer = split /;/, $line[1];
			foreach my $enhancer_id (@enhancer) {
				my $enhancer_position = $enhancer_id_2_position{$enhancer_id};
				my $dist = abs($promoter_position - $enhancer_position);
				if($dist < 5000){
					$dist = 5000;
				}
				my $log_dist = log($dist + 1);
				my $powerlaw_contact = exp($scale - $gamma * $log_dist);
				print OUTPUT3 "$promoter_id\t$enhancer_id\t$powerlaw_contact\n";
			}
		}
		else{
			my $promoter_id = $line[0];
			my $promoter_position = $promoter_id_2_position{$promoter_id};
			my $enhancer_id = $line[1];
			my $enhancer_position = $enhancer_id_2_position{$enhancer_id};
			my $dist = abs($promoter_position - $enhancer_position);
			if($dist < 5000){
				$dist = 5000;
			}
			my $log_dist = log($dist + 1);
			my $powerlaw_contact = exp($scale - $gamma * $log_dist);
			print OUTPUT3 "$promoter_id\t$enhancer_id\t$powerlaw_contact\n";
		}
	}
}

