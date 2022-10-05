use strict;

open INPUT,"/home1/yangyc/project/ImagingGWAS/results/enhancer_signal.quantile_normalized.matrix/enhancer_signal.quantile_normalized.matrix" or die;

my $header = <INPUT>;
chomp($header);

my @biosample = split /\t/, $header;

my $index;

for(my $i=1;$i<=$#biosample;$i++){
	if($biosample[$i] eq $ARGV[0]){
		$index = $i;
	}
}


my $column_index = $index + 1;

`cut -f 1,$column_index /home1/yangyc/project/ImagingGWAS/results/enhancer_signal.quantile_normalized.matrix/enhancer_signal.quantile_normalized.matrix | grep -v 'enhancer_id' > /home1/yangyc/project/ImagingGWAS/results/ABC_input.activity/$ARGV[0].txt`;

