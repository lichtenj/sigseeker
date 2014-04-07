my $file = shift or die;

open(IN, $file);
while(my $rec = <IN>)
{
	chomp($rec);
	#chr12	70260230	70260576	2	1,2	1	1	Promoter	Rps29	2	70260487	Rps29 Promoter -,Lrr1 Upstream +,
	my @tmp = split(/\t/,$rec);
	print $tmp[0]."\t";
	print $tmp[1]."\t";
	print $tmp[2]."\n";
}
close IN;
