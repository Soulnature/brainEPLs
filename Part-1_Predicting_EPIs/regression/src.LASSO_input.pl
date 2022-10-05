use strict;

open INPUT1,"../../results/promoter_signal.matrix/promoter_signal.matrix" or die;

my $header1 = <INPUT1>;

my %promoter_id_2_values;

while(<INPUT1>){
	chomp;
	my @line = split /\t/, $_;
	$promoter_id_2_values{$line[0]} = $_;
}


open INPUT2,"../../results/enhancer_signal.quantile_normalized.matrix/enhancer_signal.quantile_normalized.matrix" or die;

my $header2 = <INPUT2>;

my %enhancer_id_2_values;

while(<INPUT2>){
	chomp;
	my @line = split /\t/, $_;
	$enhancer_id_2_values{$line[0]} = $_;
}


open INPUT3,"../../results/promoter_enhancer_pairs/promoter_enhancer_pairs.txt" or die;

my $idx = 0;

while(<INPUT3>){
	chomp;
	$idx++;
	my @line = split /\t/, $_;
	
	open OUTPUT3_1,">../../results/LASSO_input/LASSO_input.promoter.$idx.txt" or die;
	my $promoter_values = $promoter_id_2_values{$line[0]};
	print OUTPUT3_1 "$header1";
	print OUTPUT3_1 "$promoter_values\n";
	
	open OUTPUT3_2,">../../results/LASSO_input/LASSO_input.enhancer.$idx.txt" or die;
	print OUTPUT3_2 "$header2";
	if($line[1] =~ /;/){
		my @enhancers = split /;/, $line[1];
		foreach my $enhancer (@enhancers) {
			my $enhancer_values = $enhancer_id_2_values{$enhancer};
			print OUTPUT3_2 "$enhancer_values\n";
		}
	}
	else{
		my $enhancer_values = $enhancer_id_2_values{$line[1]};
		print OUTPUT3_2 "$enhancer_values\n";
	}
}

