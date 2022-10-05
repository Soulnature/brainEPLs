use strict;

open INPUT,"../../results/promoter_signal/$ARGV[0].txt" or die;

open OUTPUT,">../../results/promoter_usage/$ARGV[0].txt" or die;

while(<INPUT>){
	chomp;
	my @line = split /\t/, $_;
	if($line[1] > 1){
		print OUTPUT "$line[0]\t1\n";
	}
	else{
		print OUTPUT "$line[0]\t0\n";
	}
}

