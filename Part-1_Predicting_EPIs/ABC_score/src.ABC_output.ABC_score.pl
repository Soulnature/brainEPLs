use strict;

open INPUT1,"/home1/yangyc/project/ImagingGWAS/results/ABC_input.activity/$ARGV[0].txt" or die;

my %enhancer_id_2_activity;

while(<INPUT1>){
	chomp;
	my @line = split /\t/, $_;
	$enhancer_id_2_activity{$line[0]} = $line[1];
}

open INPUT2,"/home1/yangyc/project/ImagingGWAS/results/ABC_input.contact/ABC_input.contact.txt" or die;

my %gene_enhancer_id_2_contact;

while(<INPUT2>){
	chomp;
	my @line = split /\t/, $_;
	$gene_enhancer_id_2_contact{$line[0]."\t".$line[1]} = $line[2];
}

open INPUT3,"/home1/yangyc/project/ImagingGWAS/results/ABC_output.ABC_baseline/$ARGV[0].txt" or die;

my %gene_id_2_ABC_baseline;

while(<INPUT3>){
	chomp;
	my @line = split /\t/, $_;
	$gene_id_2_ABC_baseline{$line[0]} = $line[1];
}


open INPUT4,"/home1/yangyc/project/ImagingGWAS/results/ABC_input.1MB_interval/ABC_input.1MB_interval.promoter_enhancer_pairs.txt" or die;

open OUTPUT4,">/home1/yangyc/project/ImagingGWAS/results/ABC_output.ABC_score/$ARGV[0].txt" or die;

while(<INPUT4>){
	chomp;
	my @line = split /\t/, $_;
	if($line[1] =~ /;/){
		my @enhancer = split /;/, $line[1];
		foreach my $enhancer (@enhancer) {
			my $activity = $enhancer_id_2_activity{$enhancer};
			my $contact = $gene_enhancer_id_2_contact{$line[0]."\t".$enhancer};
			my $ABC_baseline = $gene_id_2_ABC_baseline{$line[0]};
			my $ABC_score;
			if($ABC_baseline == 0) {
				$ABC_score = 0;
			}
			else{
				$ABC_score = ($activity * $contact) / $ABC_baseline;
			}
			print OUTPUT4 "$line[0]\t$enhancer\t$ABC_score\n";
		}
	}
	else{
		my $activity = $enhancer_id_2_activity{$line[1]};
		my $contact = $gene_enhancer_id_2_contact{$line[0]."\t".$line[1]};
		my $ABC_baseline = $gene_id_2_ABC_baseline{$line[0]};
		my $ABC_score;
		if($ABC_baseline == 0) {
			$ABC_score = 0;
		}
		else{
			$ABC_score = ($activity * $contact) / $ABC_baseline;
		}
		print OUTPUT4 "$line[0]\t$line[1]\t$ABC_score\n";
	}
}

