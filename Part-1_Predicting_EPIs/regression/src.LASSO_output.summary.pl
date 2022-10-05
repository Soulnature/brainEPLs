use strict;

open OUTPUT,">/home1/yangyc/project/imaging_GWAS/results/LASSO_output.summary/LASSO_output.summary.txt" or die;

for(my $i=1;$i<=19874;$i++){
	my $fileExist = -e "/home1/yangyc/project/imaging_GWAS/results/LASSO_output/LASSO_output.$i.txt";
	if( $fileExist ) {
		open INPUT1,"/home1/yangyc/project/imaging_GWAS/results/LASSO_input/LASSO_input.promoter.$i.txt" or die;
		my $header1 = <INPUT1>;
		my $promoter_id;
		while(<INPUT1>){
			chomp;
			my @line = split /\t/, $_;
			$promoter_id = $line[0];
		}
		
		open INPUT2,"/home1/yangyc/project/imaging_GWAS/results/LASSO_output/LASSO_output.$i.txt" or die;
		while(<INPUT2>){
			chomp;
			my @line = split /\t/, $_;
			if($line[1]>0){
				print OUTPUT "$promoter_id\t$line[0]\n";
			}
		}
	}
}

