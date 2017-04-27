my $file = shift or die;
my $tempid = shift or die;
my $quant = shift or die;

open(R, ">".$tempid."sigcorr_pipeline.r");

print R 'library(pracma)'."\n";
print R 'library(scatterplot3d)'."\n";

#print $file."\n";
$file =~ /\_([A-Za-z\d\-]+and[A-Za-z\d\-]+)\_peaks/;
#print "\n\n".$1."\n\n";
my @tools = split(/and/,$1);
my $count = -1;
foreach my $tool (reverse @tools)
{
	my $set = $tool.'=c(';

	open(IN, $file) or die "Cannot open $file";
	while(my $record = <IN>)
	{
		chomp($record);
		my @tmp = split(/\t/,$record);

		$tmp[$count] =~ s/\,//g;	
		$set .= $tmp[$count].",";
	}
	close IN;

	chop($set);

	print R $set.')'."\n";

	$count--;
}

#Setting up the dimension assignments
my @subtools = @tools;
my $Y_TOOL = shift(@subtools);
print R 'X <- cbind('.join(",",@subtools).')'."\n";
print R 'Y <- '.$Y_TOOL."\n";
print R 'df = data.frame(X,Y)'."\n";

#Orthogonal Regression
print R 'odr <- odregress(X,Y);odr'."\n";
print R '( cc <- odr$coeff )'."\n";
print R '( res <- odr$resid )'."\n";
print R 'res.qt <- quantile(res, probs = c('.(0.01 * $quant).','.(1 - (0.01 * $quant)).'))'."\n";
print R 'res.qt'."\n";
print R 'res.qt[1]'."\n";
print R 'res.qt[2]'."\n";
print R 'want <- which(res >= res.qt[1] & res <= res.qt[2],arr.ind=TRUE)'."\n";
print R 'unwanted <- which(res < res.qt[1] | res > res.qt[2])'."\n";
#print R 'res'."\n";
#print R 'want'."\n";
#print R 'unwanted'."\n";
print R 'write.csv(df[unwanted,],file="'.$file.'.outliers")'."\n";
##
#print R 'pVAL <- df[want,]'."\n";
##print R 'pVAL[,c(3,1,2)]'."\n";
#print R 'nVAL <- df[unwanted,]'."\n";
##print R 'nVAL[,c(3,1,2)]'."\n";
#
#Visualization
#print R 'png("'.$file.'_CORR_MOD_3D.png",width = 512, height = 512)'."\n";
#print R 's3d <- scatterplot3d('.join(",",@tools).', main="Orthogonal Regression Model")'."\n";
#print R 's3d$points3d(pVAL[,c(3,1,2)], col = "red", pch = 21, bg = "red", cex = 0.8)'."\n";
#print R 's3d$points3d(nVAL[,c(3,1,2)], col = "black", pch = 21, bg = "black", cex = 0.8)'."\n";
#print R 'dev.off()'."\n";

close R;

system('Rscript '.$tempid.'sigcorr_pipeline.r > '.$file.'_Rout 2> '.$file.'_Rerr');
#system('rm '.$tempid."sigcorr_pipeline.r");

##############
my $outliers;

open(IN, $file.'.outliers') or die "Cannot open outliers";
while(my $record = <IN>)
{
	$record =~ s/\"//g;
	chomp($record);
	my @tmp = split(/\,/,$record);

	$outliers->{$tmp[0]} = 1;
}
close IN;

my $count = 1;
open(IN, $file) or die "Cannot open $file";
my $newfile = $file;
$newfile =~ s/traditionaljoin/quantitative$quant/;
open(OUT, ">".$newfile) or die "Cannot open $newfile";
while(my $record = <IN>)
{
	if($outliers->{$count} != 1 || ! $outliers->{$count})
	{
		print OUT $record;
	}
	$count++;
}
close OUT;
close IN;

#if(-e $file.'.outliers')
#{
#	system('rm '.$file.'.outliers');
#}
