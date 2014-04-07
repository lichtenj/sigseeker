use POSIX;
use strict;

my $set1 = shift or die;
my $set2 = shift or die;
my $output = shift or die;

system('intersectBed -a '.$set1.' -b '.$set2.' -wao > temp_ab.bed');
system('intersectBed -a '.$set2.' -b '.$set1.' -wao > temp_ba.bed');

my $max_peaks1 = 0;
my $max_peaks2 = 0;
my $sets;
my $samples;
my @set_vals;
my @sample_vals;

my $peaks = 0;
open(TMP, "temp_ba.bed");
while(my $record = <TMP>)
{
	chomp($record);

	my @tmp = split(/\t/,$record);

	if($tmp[-1] == 0)
	{
		if($tmp[4] =~ /^\d/)
		{
			push (@set_vals,0);
			push (@sample_vals,$tmp[4]);
		}
		else
		{
			push (@set_vals,0);
			push (@sample_vals,$tmp[3]);
		}
	}
}
close TMP;

open(TMP, "temp_ab.bed");
while(my $record = <TMP>)
{
	chomp($record);

	my @tmp = split(/\t/,$record);

	if($tmp[-1] == 0)
	{
		if($tmp[4] =~ /^\d/)
		{
			push (@set_vals,$tmp[4]);
			push (@sample_vals,0);
		}
		else
		{
			push (@set_vals,$tmp[3]);
			push (@sample_vals,0);
		}
	}

	if($tmp[4] =~ /^\d/)
	{
		if($tmp[9] =~ /^\d/)
		{
			if($tmp[4] != -1 && $tmp[9] != -1 && $tmp[5] ne '.')
			{
				$peaks++;
				if($max_peaks1 < $tmp[4]){$max_peaks1 = $tmp[4];}
				if($max_peaks2 < $tmp[9]){$max_peaks2 = $tmp[9];}

				push (@set_vals,$tmp[4]);
				push (@sample_vals,$tmp[9]);
			}
		}
		else
		{
			if($tmp[4] != -1 && $tmp[8] != -1 && $tmp[5] ne '.')
			{
				$peaks++;
				if($max_peaks1 < $tmp[4]){$max_peaks1 = $tmp[4];}
				if($max_peaks2 < $tmp[8]){$max_peaks2 = $tmp[8];}

				push (@set_vals,$tmp[4]);
				push (@sample_vals,$tmp[8]);
			}
		}
	}
	else
	{
		if($tmp[8] =~ /^\d/)
		{
			if($tmp[3] != -1 && $tmp[8] != -1 && $tmp[4] ne '.')
			{
				$peaks++;
				if($max_peaks1 < $tmp[3]){$max_peaks1 = $tmp[3];}
				if($max_peaks2 < $tmp[8]){$max_peaks2 = $tmp[8];}

				push (@set_vals,$tmp[3]);
				push (@sample_vals,$tmp[8]);
			}
		}
		elsif($tmp[7] =~ /^\d/)
		{
			if($tmp[3] != -1 && $tmp[7] != -1 && $tmp[4] ne '.')
			{
				$peaks++;
				if($max_peaks1 < $tmp[3]){$max_peaks1 = $tmp[3];}
				if($max_peaks2 < $tmp[7]){$max_peaks2 = $tmp[7];}

				push (@set_vals,$tmp[3]);
				push (@sample_vals,$tmp[7]);
			}
		}
#		else
#		{
#
#		}
	}
#	print $tmp[4]."\t".$tmp[9]."\n";
}
close TMP;

my $data = 'NA';
if($max_peaks1 > 0 && $max_peaks2 > 0)
{
	foreach my $val (@set_vals)
	{
	#	if($max_peaks1 == 0){$max_peaks1 = 1;}
	#
		$sets .= ($val / $max_peaks1).',';
	}
	foreach my $val (@sample_vals)
	{
	#	if($max_peaks2 == 0){$max_peaks2 = 1;}
	#
		$samples .= ($val / $max_peaks2).',';
	}
	
	open(R, ">peaks.r");
	print R 'library(ggplot2)'."\n";
	#        print R 'peaks = data.frame('."\n";
    chop($sets);

    print R 'invisible(set <- c('.$sets.'))'."\n";
    chop($samples);
    print R 'invisible(peaks <- c('.$samples.'))'."\n";
#        print R ')'."\n";
	print R 'level <- cor(set,peaks,method="pearson")'."\n";
	print R 'level'."\n";
    print R 'png(file = "'.$output.'.png", width = 1024, height = 640, units = "px", pointsize = 12, bg="transparent")'."\n";

#        print R 'qplot(set, peaks, data = peaks,geom="boxplot",main = "Peak Calling") + aes(ymin=0,ymax='.($max_peaks + 1000).') + scale_y_continuous(breaks=seq(0,'.$max_peaks.',10000)) + geom_jitter(position=position_jitter(w=0.1, h=0.1)) + ylab("#\\ Peaks") + theme(text = element_text(size=20))'."\n";
	$set1 =~ /([\w\-\_]+)\.bed$/;
	my $xlab = $1;
#	print $xlab;
	$set2 =~ /([\w\-\_]+)\.bed$/;
	my $ylab = $1;
	print R 'plot(set~peaks,main=level,ylab="'.$ylab.'",xlab="'.$xlab.'")'."\n";
        print R 'dev.off()'."\n";
        
        close R;
        
#        system('R --vanilla < peaks.r 2> /dev/null > /dev/null');
	$data = `R --vanilla --slave < peaks.r 2> /dev/null`;
}

open(OUT,">".$output.'.tsv');
print OUT $peaks."\t".$max_peaks1."\t".$max_peaks2."\t".$data."\n";
close OUT;