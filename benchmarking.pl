use strict;

my $input = shift or die;
my $dir = './';
my $factor = shift or die;

my $positives = $dir.'Positive_classified_peaks_'.$factor.'.txt';
my $negatives = $dir.'Negative_classified_peaks_'.$factor.'.txt';

#print "Positives\n";
my $cmd = 'intersectBed -a '.$positives.' -b '.$input.' -u | wc -l';
my $tp = `$cmd`;
chomp($tp);

#print "Negatives\n";
$cmd = 'intersectBed -a '.$negatives.' -b '.$input.' -u | wc -l';
my $fp = `$cmd`;
chomp($fp);

$cmd = 'cat '.$positives.' | wc -l';
my $mp = `$cmd`;
chomp($mp);

$cmd = 'cat '.$input.' | wc -l';
my $p = `$cmd`;
chomp($p);

$cmd = 'cat '.$negatives.' | wc -l';
my $n = `$cmd`;
chomp($n);

my $tn = $n - $fp;
my $fn = $mp - $tp;

#Total Peaks
print $p."\t";
#print $n."\t";

#print "TP\tFP\tTN\tFN\tSEN\tSPE\n";

#True Positives
print $tp."\t";
#True Negatives
print $tn."\t";
#False Positives
print $fp."\t";
#False Negatives
print $fn."\t";

#if($tp == 0 || $tn == 0)
#{
#	print "0\t0\t0\t0\t0\t0";
#}
#else
#{
	#Sensitivity
	if($tp + $fn != 0)
	{
		print sprintf("%.3f", ($tp / ($tp + $fn)))."\t";
	}
	else
	{
		print "0.000\t";
	}

	#Specificity
	if($tn + $fn != 0)
	{
		print sprintf("%.3f", ($tn / ($fp +$tn)))."\t";
	}
	else
	{
		print "0.000\t";
	}

	#Informedness
	print sprintf("%.3f", (($tn / ($fp +$tn)) + ($tp / ($tp + $fn)) - 1))."\t";

	#Precision
	if($tp + $fp != 0)
	{
		print sprintf("%.3f", ($tp / ($tp + $fp)))."\t";
	}
	else
	{
		print "0.000\t";
	}

#	#Negative Predictive Value
#	if($tn + $fn != 0)
#	{
#		print sprintf("%.3f", ($tn / ($tn + $fn)))."\t";
#	}
#	else
#	{
#		print "0.000\t";
#	}
#
#	#Fall Out
#	if($n != 0)
#	{
#		print sprintf("%.3f", ($fp / $n))."\t";
#	}
#	else
#	{
#		print "0.000\t";
#	}
#
#	#False Discovery Rate
#	if($tp + $fp != 0)
#	{
#		print sprintf("%.3f", ($fp / ($tp + $fp)))."\t";
#	}
#	else
#	{
#		print "0.000\t";
#	}
#
	#Accuracy
	if($p + $n != 0)
	{
		print sprintf("%.3f", (($tp + $tn) / ($p + $n)))."\t";
	}
	else
	{
		print "0.000\t";
	}
#}

	#Jaccard Index Precision
	my $jaccard = jaccard_index_precision($input,$dir.'Positive_classified_peaks_'.$factor.'.txt',$input.'POSITIVE'.$factor);
	print sprintf("%.3f", $jaccard)."\t";

	my $sorenson_coeff = sorenson($input,$dir.'Positive_classified_peaks_'.$factor.'.txt',$input.'POSITIVE'.$factor);
	print sprintf("%.3f", $sorenson_coeff)."\n";

sub jaccard_index_precision
{
	my $a = shift or die;
	my $b = shift or die;
	my $id = shift or die;
	
	#Union
	my $hash;
	if(! -e $id.'_TEMP_INTER')
	{
		system('bedtools intersect -a '.$a.' -b '.$b.' -wao > '.$id.'_TEMP_INTER');
	}
	open(IN, $id.'_TEMP_INTER') or die "Cannot open inter";
	while(my $rec = <IN>)
	{
	        chomp($rec);
	        my @tmp = split(/\t/,$rec);
	        if($tmp[-1] == 0){next;}

		my $count = 0;
		my $offset = 0;
		foreach my $check (@tmp)
		{
			if($check =~ /^chr\w+$/)
			{
				$offset = $count;
			}
			$count++;
		}

	        $hash->{$tmp[0].':'.$tmp[1].'-'.$tmp[2].'_'.$tmp[$offset].':'.$tmp[$offset+1].'-'.$tmp[$offset+2]}->{'intersect'} = $tmp[-1];
	}
	close IN;
#	system('rm '.$id.'_TEMP_INTER');
	
	#Union
	if(! -e $id.'_TEMP_UNI')
	{
		system('bedtools intersect -a '.$a.' -b '.$b.' -wa -wb > '.$id.'_TEMP_UNI');
	}
#print "\n";
	open(IN, $id.'_TEMP_UNI') or die "Cannot open uni";
	while(my $rec = <IN>)
	{
#print $rec;
	        chomp($rec);
	        my @tmp = split(/\t/,$rec);
		my $count = 0;
		my $offset = 0;
		foreach my $check (@tmp)
		{
			if($check =~ /^chr\w+$/)
			{
				$offset = $count;
			}
			$count++;
		}
#print $offset."\n";
	        my $start = 0;
	        my $end = 0;
	        if($tmp[1] < $tmp[$offset + 1])
	        {
	                $start = $tmp[1];
	        }
	        else
	        {
	                $start = $tmp[$offset + 1];
	        }
	        if($tmp[2] < $tmp[$offset + 2])
	        {
	                $end = $tmp[$offset + 2];
	        }
	        else
	        {
	                $end = $tmp[2];
	        }

		#print "\n";
#print $tmp[0].':'.$tmp[1].'-'.$tmp[2].'_'.$tmp[$offset].':'.$tmp[$offset+1].'-'.$tmp[$offset+2]."\n";
	        $hash->{$tmp[0].':'.$tmp[1].'-'.$tmp[2].'_'.$tmp[$offset].':'.$tmp[$offset+1].'-'.$tmp[$offset+2]}->{'union'} = $end - $start;
	}
	close IN;
#	system('rm '.$id.'_TEMP_UNI');
#return 0;
	my $sum = 0;
	foreach my $recid (keys %$hash)
	{
#		print $recid."\n";
#		print $hash->{$recid}->{'intersect'}."\n";
#		print $hash->{$recid}->{'union'}."\n";

	        my $ratio = $hash->{$recid}->{'intersect'} / $hash->{$recid}->{'union'};
	        $sum += $ratio;
	}


	my $jaccard = $sum / scalar(keys %$hash);

#	print "\n";
#	print $sum."\n";
#	print scalar(keys %$hash)."\n";
#	print $jaccard."\n";
	return $jaccard;
}

###########################################

sub sorenson
{
	my $a = shift or die;
	my $b = shift or die;
	my $id = shift or die;
	
	#Intersection
	my $hash;
	if(! -e $id.'_TEMP_SOR')
	{
		system('bedtools intersect -a '.$a.' -b '.$b.' -u > '.$id.'_TEMP_SOR');
	}

	my $inter_cmd = 'wc -l '.$id.'_TEMP_SOR';
	my $inter_cnt = `$inter_cmd`;

	my $a_cmd = 'wc -l '.$a;
	my $a_cnt = `$a_cmd`;

	my $b_cmd = 'wc -l '.$b;
	my $b_cnt = `$b_cmd`;

	$inter_cnt =~ /(\d+)/;
	$inter_cnt = $1;

	$a_cnt =~ /(\d+)/;
	$a_cnt = $1;

	$b_cnt =~ /(\d+)/;
	$b_cnt = $1;

	my $sorenson_coeff = (2 * $inter_cnt) / ( $a_cnt + $b_cnt );

#	print $a_cnt."\n";
#	print $b_cnt."\n";
#	print $inter_cnt."\n";
#	print $sorenson_coeff."\n";
	return $sorenson_coeff;
}

