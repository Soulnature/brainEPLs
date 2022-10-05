use strict;

open INPUT,"../../materials/enhancer_FANTOM/F5.hg38.enhancers.expression.usage.matrix" or die;

my $header = <INPUT>;
chomp($header);
my @header = split /\t/, $header;

my $index;

for(my $i=0;$i<=$#header;$i++){
	if($header[$i] eq $ARGV[1]){
		$index = $i;
	}
}

open OUTPUT,">../../results/enhancer_usage/$ARGV[0].txt" or die;

while(<INPUT>){
	chomp;
	my @line = split /\t/, $_;
	print OUTPUT "$line[0]\t$line[$index]\n";
}
