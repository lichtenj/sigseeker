use strict;
use Statistics::Basic qw(:all);
use List::Util qw(max);
use Math::Combinatorics;

my $tempid = pop(@ARGV);
my $output = pop(@ARGV);
my $filter = pop(@ARGV);
my $over = pop(@ARGV);
my @files = @ARGV;

#my $list = join(" ",@files);
my $list;
foreach my $file (@files)
{
	#Reforce Strandedness
	#system('perl add_strandedness.pl TEMP_sorted.bed > TEMP_sorted_stranded.bed');
	if($file =~ /MACS2\_/ && $file !~ /Replicates/)
	{
		open(IN, $file) or die "Cannot open macs2";
		open(OUT, ">".$file.'_macs2.bed');
		while(my $record = <IN>)
		{
			chomp($record);
			my @tmp = split(/\t/,$record);
			print OUT $tmp[0]."\t";
			print OUT $tmp[1]."\t";
			print OUT $tmp[2]."\t";
			print OUT "+\t";
			print OUT $tmp[4]."\n";
		}

		close IN;
		close OUT;
		$list .= $file.'_macs2.bed ';
#		system('perl add_tools.pl '.$file.'_macs2 > '.$file.'_stranded.bed');
	}
	elsif($file =~ /MACS\_/ && $file !~ /Replicates/)
	{
		open(IN, $file) or die "Cannot open macs";
		open(OUT, ">".$file.'_macs.bed');
		while(my $record = <IN>)
		{
			chomp($record);
			my @tmp = split(/\t/,$record);
			print OUT $tmp[0]."\t";
			print OUT $tmp[1]."\t";
			print OUT $tmp[2]."\t";
			print OUT "+\t";
			print OUT $tmp[4]."\n";
		}
		close IN;
		close OUT;
		$list .= $file.'_macs.bed ';
#		system('perl add_tools.pl '.$file.'_macs > '.$file.'_stranded.bed');
	}
	elsif($file =~ /PeaKDEck\_/i && $file !~ /Replicates/)
	{
		open(IN, $file) or die "Cannot open peakdeck";
		open(OUT, ">".$file.'_peakdeck.bed') or die "Cannot open ".$file."_peakdeck for writing";
		while(my $record = <IN>)
		{
			chomp($record);
			my @tmp = split(/\t/,$record);
			print OUT $tmp[0]."\t";
			print OUT $tmp[1]."\t";
			print OUT $tmp[2]."\t";
			print OUT "+\t";
			print OUT $tmp[4]."\n";
		}
		close IN;
		close OUT;
		$list .= $file.'_peakdeck.bed ';
#		system('perl add_tools.pl '.$file.'_peakdeck > '.$file.'_stranded.bed');
	}
	elsif($file =~ /AREMp?\d*\_/i && $file !~ /Replicates/)
	{
		open(IN, $file) or die "Cannot open peakdeck";
		open(OUT, ">".$file.'_arem.bed') or die "Cannot open ".$file."_arem for writing";
		while(my $record = <IN>)
		{
			chomp($record);
			my @tmp = split(/\t/,$record);
			print OUT $tmp[0]."\t";
			print OUT $tmp[1]."\t";
			print OUT $tmp[2]."\t";
			print OUT "+\t";
			print OUT $tmp[4]."\n";
		}
		close IN;
		close OUT;
		$list .= $file.'_arem.bed ';
#		system('perl add_tools.pl '.$file.'_peakdeck > '.$file.'_stranded.bed');
	}
	elsif($file =~ /Replicates\_/)
	{
		open(IN, $file) or die "Cannot open macs";
		open(OUT, ">".$file.'_replicates.bed');
		while(my $record = <IN>)
		{
			chomp($record);
			my @tmp = split(/\t/,$record);
			print OUT $tmp[0]."\t";
			print OUT $tmp[1]."\t";
			print OUT $tmp[2]."\t";
			print OUT "+\t";
			print OUT $tmp[4]."\n";
		}
		close IN;
		close OUT;
		$list .= $file.'_replicates.bed ';
#		system('perl add_tools.pl '.$file.'_macs > '.$file.'_stranded.bed');
	}
	elsif($file =~ /(CISGENOME|ERANGE|SICER|CCAT\-TF|CCAT\-Histone|SWEMBLE|CCATHS\d*|CCATTF\d*|HIDDENDOMAINS)\_/ && $file !~ /Replicates/)
	{
		open(IN, $file) or die "Cannot open cisgenome";
		open(OUT, ">".$file.'_'.$1.'.bed');
		while(my $record = <IN>)
		{
			if($record =~ /^\#/){next;}
			chomp($record);
			my @tmp = split(/\t/,$record);
			print OUT $tmp[0]."\t";
			print OUT $tmp[1]."\t";
			print OUT $tmp[2]."\t";
			print OUT "+\t";
			print OUT $tmp[3]."\n";
		}
		close IN;
		close OUT;
#		system('perl add_tools.pl '.$file.'_cisgenome > '.$file.'_stranded.bed');	
		$list .= $file.'_'.$1.'.bed ';
	}
	else
	{
#		system('perl add_tools.pl '.$file.' > '.$file.'_stranded.bed');
		$list .= $file.' ';
	}
	#	my @tmp = split(/\_/,$file);
	#	$output .= $tmp[-2].'_';
}
chop($list);

my @tmp = split(/\ /,$list);

print "MultiInter\n";
#my $count = 1;
open(OUTMF, ">".$tempid."temp_multi_filter.bed") or die "Cannot open ".$tempid."temp_multi_filter";
#foreach my $list (map { join " ", @$_ } permute(@tmp))
#{
#	system('bedtools multiinter -i '.$list.' > '.$tempid.'temp_multi_'.$count.'.bed');
	system('bedtools multiinter -i '.$list.' > '.$tempid.'temp_multi_temp.bed');
#	my $cmd = 'wc -l '.$tempid.'temp_multi_'.$count.'.bed';
	my $cmd = 'wc -l '.$tempid.'temp_multi_temp.bed';
	my $out = `$cmd`;

	if($out =~ /^0/)
	{
#		print "No multi overlap\n";
		open(OUT,">".$tempid."_FAILED");
		print OUT "temp_multi_filter\n";
		close OUT;
		exit;
	}

#	open(IN, $tempid."temp_multi_".$count.".bed") or die "Cannot open ".$tempid."temp_multi_".$count.".bed";
	open(IN, $tempid."temp_multi_temp.bed") or die "Cannot open ".$tempid."temp_multi_temp.bed";
	while(my $rec = <IN>)
	{
		chomp($rec);
		my @tmp = split(/\t/,$rec);
		if($tmp[3] eq scalar(@files))
		{
			if($tmp[2] - $tmp[1] >= $over)
			{
				print OUTMF $rec."\n";
			}
		}
	}
	close IN;
#
#	system('rm '.$tempid.'temp_multi_'.$count.'.bed');
#
#	$count++;
#}
close OUTMF;

print "Intersect\n";
my $cmdX = 'wc -l '.$tempid.'temp_multi_filter.bed';
my $out = `$cmdX`;
if($out =~ /^0/)
{
	print "No multi filter overlap\n";
	open(OUT,">".$tempid."_FAILED");
	print OUT "temp_multi_filter\n";
	close OUT;
	exit;
}
system('bedtools sort -i '.$tempid.'temp_multi_filter.bed > '.$tempid.'temp_multi_filter_sort.bed');
system('mv '.$tempid.'temp_multi_filter_sort.bed '.$tempid.'temp_multi_filter.bed');

system('bedtools intersect -a '.$tempid.'temp_multi_filter.bed -b '.$list.' -wb > '.$tempid.'temp_inter.bed');
system('bedtools sort -i '.$tempid.'temp_inter.bed > '.$tempid.'temp_inter_sort.bed');
system('mv '.$tempid.'temp_inter_sort.bed '.$tempid.'temp_inter.bed');
open(OUT, ">".$tempid."temp_inter_filter.bed");
open(IN, $tempid."temp_inter.bed");
while(my $rec = <IN>)
{
	chomp($rec);
	my @tmp = split(/\t/,$rec);
	if($tmp[3] eq scalar(@files))
	{
		print OUT $rec."\n";
	}
}
close IN;
close OUT;

system('bedtools sort -i '.$tempid.'temp_inter_filter.bed > '.$tempid.'temp_inter_filter_sort.bed');
system('mv '.$tempid.'temp_inter_filter_sort.bed '.$tempid.'temp_inter_filter.bed');

#print "Merge\n";
my $R;
my $Rcount = 0;
system('bedtools merge -c 5,'.(11 + scalar(@files)).' -o collapse -i '.$tempid.'temp_inter_filter.bed > '.$tempid.'temp_merge.bed');
system('bedtools sort -i '.$tempid.'temp_merge.bed > '.$tempid.'temp_merge_sort.bed');
system('mv '.$tempid.'temp_merge_sort.bed '.$tempid.'temp_merge.bed');

open(OUT, ">".$tempid."temp_merge_filter.bed");
open(IN, $tempid."temp_merge.bed");
while(my $rec = <IN>)
{
	$Rcount++;
	chomp($rec);
	my @tmp = split(/\t/,$rec);
	my @sets = split(/\,/,$tmp[3]);
	my @scores = split(/\,/,$tmp[4]);

	my $hash;
	my $count = 0;
	foreach my $set (@sets)
	{
		if(! $R->{$set}->{$Rcount})
		{
			$R->{$set}->{$Rcount} = $scores[$count];
		}
		$hash->{$set} = 1;
		$count++;
	}
	my $set_data = "";
	foreach my $set (sort {$a <=> $b} keys %$hash)
	{
		$set_data .= @files[$set - 1].',';
	}
	chop($set_data);

	my $max_y = 0;

	my $y_m = mean(@scores);
	my $s = stddev(@scores);

	my $outlier = 0;

	if(scalar(@scores) == scalar(@files) && $outlier == 0)
	{
		print OUT $tmp[0]."\t";
		print OUT $tmp[1]."\t";
		print OUT $tmp[2]."\t";
		print OUT "+";
		print OUT "\t";
		print OUT $y_m;
		print OUT "\t";
		print OUT join("\t",@scores);
		print OUT "\n";
	}
}
close IN;
close OUT;

#system('bedtools sort -i '.$tempid.'temp_merge_filter.bed > '.$tempid.'temp_merge_filter_sort.bed');
#system('mv '.$tempid.'temp_merge_filter_sort.bed '.$tempid.'temp_merge_filter.bed');

my $cmdY = 'wc -l '.$tempid.'temp_merge_filter.bed';
my $out = `$cmdY`;
if($out =~ /^0/)
{
	print "No merge filter overlap\n";
	open(OUT,">".$tempid."_FAILED");
	print OUT "temp_merge_filter\n";
	close OUT;
	exit;
}
system('cp '.$tempid.'temp_merge_filter.bed '.$output.'_'.$over.'_traditionaljoin.bed');

print 'perl SIG_CORR_PIPELINE_3D.pl '.$output.'_'.$over.'_traditionaljoin.bed '.$tempid.' '.$filter."\n";
system('perl SIG_CORR_PIPELINE_3D.pl '.$output.'_'.$over.'_traditionaljoin.bed '.$tempid.' '.$filter);
system('wc -l '.$output.'_'.$over.'_traditionaljoin.bed');
system('wc -l '.$output.'_'.$over.'_quantitative'.$filter.'.bed');

system('rm '.$tempid.'*.bed');


sub normalize
{
	my ($array) = @_;
	my @result;

#	print "array: ";
#	print @$array;
#	print "\n";

	my $max = max @$array;

#	print "\n".$max."\n";

	foreach my $value (@$array)
	{
		push(@result, $value / $max);	
	}

	return @result;
}
