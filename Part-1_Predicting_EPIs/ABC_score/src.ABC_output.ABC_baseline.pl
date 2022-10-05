use strict;

open INPUT1,"/home1/yangyc/project/ImagingGWAS/results/ABC_input.activity/$ARGV[0].txt" or die;

my %enhancer_id_2_activity;

while(<INPUT1>){
	chomp;
	my @line = split /\t/, $_;
	$enhancer_id_2_activity{$line[0]} = $line[1];
}

open INPUT2,"/home1/yangyc/project/ImagingGWAS/results/ABC_input.contact/ABC_input.contact.txt" or die;

my %gene_id_2_ABC_baseline;

while(<INPUT2>){
	chomp;
	my @line = split /\t/, $_;
	my $ABC_score = $enhancer_id_2_activity{$line[1]} * $line[2];
	if(exists $gene_id_2_ABC_baseline{$line[0]}){
		$gene_id_2_ABC_baseline{$line[0]} = $gene_id_2_ABC_baseline{$line[0]} + $ABC_score;
	}
	else{
		$gene_id_2_ABC_baseline{$line[0]} = $ABC_score;
	}
}

open OUTPUT,">/home1/yangyc/project/ImagingGWAS/results/ABC_output.ABC_baseline/$ARGV[0].txt" or die;

foreach my $gene_id (sort keys %gene_id_2_ABC_baseline){
	my $ABC_baseline = $gene_id_2_ABC_baseline{$gene_id};
	print OUTPUT "$gene_id\t$ABC_baseline\n";
}

