#!/usr/bin/perl
#BUMP
use strict; 
use warnings;
use Math::Combinatorics;

my $file = './setup.tsv';
my $dir = './';
my $threshold = 1;

#my @alltools = ("CISGENOME","MACS","MACS2","ERANGE","SWEMBLE");
#my @overlaps = (1,2,3,4,5,6,7,8,9,10,100,200,300,400,500,600,700,800,900,1000,1500,2000,2500,3000);
#my @alltools = ("CISGENOME","MACS14","MACS2","ERANGE","SWEMBLE","SICER","HIDDENDOMAINS","CCATHS","CCATTF");

#Transcription Factors
#my @alltools = ("CISGENOME","MACS14","MACS2","ERANGE","SWEMBLE","CCATTF");

#Histones
#my @alltools = ("SICER","HIDDENDOMAINS","CCATHS");

#Methylation
my @alltools = ("CCATHS2","CCATTF2","AREMp2","PeaKDEck","SWEMBLE");

my @overlaps = (1,100);
#for(my $i = 10;$i<100; $i += 10)
for(my $i = 200;$i<3000; $i += 100)
{
	push(@overlaps,$i);
}

#print "Correlations\n";
for(my $quant = 5; $quant <= 15; $quant+=5)
{
	my $outdir = 'Processed_'.$quant.'/';

	foreach my $f (@overlaps)
	{
		for(my $t = 2;$t<=scalar(@alltools);$t++)
		{
			my $combinat = Math::Combinatorics->new(count => $t,data => [@alltools]);
			while(my @tools = $combinat->next_combination)
			{
				@tools = sort @tools;
				open(IN, $file);
				while(my $record = <IN>)
				{
					chomp($record);

					my @entry = split(/\t/,$record);

					my $empty = 0;
					my $cmd2 = 'perl updated_intersectbed.pl ';
					my $cmd3 = 'perl updated_union.pl ';
					foreach my $tool_name (@tools)
					{
						if(! -e $dir.$entry[0].'_'.$entry[1].'_'.$tool_name.'_peaks.bed' || -z $dir.$entry[0].'_'.$entry[1].'_'.$tool_name.'_peaks.bed')
						{
							$empty = 1;
						}

						if(! -e $dir.$entry[0].'_'.$entry[1].'_'.$tool_name.'_peaks_sorted.bed')
						{
							system('bedtools sort -i '.$dir.$entry[0].'_'.$entry[1].'_'.$tool_name.'_peaks.bed > '.$dir.$entry[0].'_'.$entry[1].'_'.$tool_name.'_peaks_sorted.bed');
						}
						$cmd2 .= $dir.$entry[0].'_'.$entry[1].'_'.$tool_name.'_peaks_sorted.bed ';
						$cmd3 .= $dir.$entry[0].'_'.$entry[1].'_'.$tool_name.'_peaks_sorted.bed ';
					}

					my $toolslist = join("and",@tools);

					$cmd2 .= $f.' '.$quant.' '.$outdir.$entry[0].'_'.$entry[1].'_'.$toolslist.'_peaks '.$entry[0].$entry[1].$toolslist.'_'.$f.'_'.$quant;
					$cmd3 .= $f.' '.$quant.' '.$outdir.$entry[0].'_'.$entry[1].'_'.$toolslist.'_peaks '.$entry[0].$entry[1].$toolslist.'_'.$f.'_'.$quant;

					if((! -e $outdir.$entry[0].'_'.$entry[1].'_'.$toolslist.'_peaks_'.$f.'_traditionaljoin.bed' || ! -e $outdir.$entry[0].'_'.$entry[1].'_'.$toolslist.'_peaks_'.$f.'_quantitative'.$quant.'.bed') && $empty != 1)
					{
#						print $cmd2."\n";
					}
				}
				close IN;
			}
		}
	}
}

my @tools = @alltools;
my @combotools = @tools;
print "Benchmarks\n";
for(my $quant = 5; $quant <= 15; $quant+=5)
{
	foreach my $f (@overlaps)
	{
		open(COM,">".$dir."General_benchmark_".$quant."_".$f.".tsv") or die "Cannot print $!";

		for(my $i = 2; $i <= scalar(@combotools); $i++)
		{
			my $combinat = Math::Combinatorics->new(count => $i,data => [@combotools]);
			while(my @tool = $combinat->next_combination)
			{
				@tool = sort @tool;
				my $query = join("and",@tool);
				push(@tools,$query);
			}
		}
		open(OUT, ">".$dir."Specific_benchmark_".$f."_".$quant.".tsv") or die "Cannot print $!";
		open(IN, $file);
		while(my $record = <IN>)
		{
			chomp($record);
			my @entry = split(/\t+/,$record);
			my $factor = $entry[0].'_'.$entry[1];

			foreach my $tool (@tools)
			{
#				print $f."\t".$quant."\t".$sample."\t".$control."\t".$tool."\n";

				#my @tmp = split(/\-/,$entry[0]);
				#my $factor = lc($tmp[1]);
				
				print OUT $entry[0]."_";
				print OUT $entry[1]."\t";
				print OUT $tool."\t";

				my @occs = ($tool =~ /and/g);
				print OUT 1 + scalar(@occs);
				print OUT "\t";	
				if(scalar(@occs) == 0)
				{
					print OUT "Single\t";
					if(! -e $dir.$entry[0].'_'.$entry[1].'_'.$tool.'_peaks.bed' || -z $dir.$entry[0].'_'.$entry[1].'_'.$tool.'_peaks.bed')
					{
						next;
					}
#					print $dir.$entry[0].'_'.$entry[1].'_'.$tool.'_peaks.bed'."\n";
					my $cmd = 'perl benchmarking.pl '.$dir.$entry[0].'_'.$entry[1].'_'.$tool.'_peaks.bed '.$factor;
					my $out;
					$out = `$cmd`;

					print OUT $out;
					print COM $entry[0]."\t";
					print COM $tool."\t";
					print COM 1 + scalar(@occs);
					print COM "\t";
					print COM "NA\t";
					print COM "NA\t";
					print COM "Single\t";
					print COM $out;
				}
				else
				{
					if($quant == 5)
					{
						if(! -e './Processed_'.$quant.'/'.$entry[0].'_'.$entry[1].'_'.$tool.'_peaks_'.$f.'_traditionaljoin.bed' || -z './Processed_'.$quant.'/'.$entry[0].'_'.$entry[1].'_'.$tool.'_peaks_'.$f.'_traditionaljoin.bed')
						{
							next;
						}

						print OUT $entry[0]."\t";
						print OUT $tool."\t";	
						print OUT 1 + scalar(@occs);
						print OUT "\t";
						print OUT "Traditional\t";

						print COM $entry[0]."\t";
						print COM $tool."\t";
						print COM 1 + scalar(@occs);
						print COM "\t";
						print COM $f."\t";
						print COM $quant."\t";
						print COM "Traditional\t";
						my $cmd = 'perl benchmarking.pl ./Processed_'.$quant.'/'.$entry[0].'_'.$entry[1].'_'.$tool.'_peaks_'.$f.'_traditionaljoin.bed '.$factor;
						my $out;
						$out = `$cmd`;
						print OUT $out;	

						print COM $out;
					}

					if(! -e './Processed_'.$quant.'/'.$entry[0].'_'.$entry[1].'_'.$tool.'_peaks_'.$f.'_quantitative'.$quant.'.bed' || -z './Processed_'.$quant.'/'.$entry[0].'_'.$entry[1].'_'.$tool.'_peaks_'.$f.'_quantitative'.$quant.'.bed')
					{
						next;
					}

					print OUT "Novel\t";
					my $cmd = 'perl benchmarking.pl ./Processed_'.$quant.'/'.$entry[0].'_'.$entry[1].'_'.$tool.'_peaks_'.$f.'_quantitative'.$quant.'.bed '.$factor;
					my $out;
					$out = `$cmd`;

					print OUT $out;
					print COM $entry[0]."\t";
					print COM $tool."\t";
					print COM 1 + scalar(@occs);
					print COM "\t";
					print COM $f."\t";
					print COM $quant."\t";
					print COM "Quantitative".$quant."\t";
					print COM $out;
				}
			}
		}
		close IN;
		close OUT;
		print "Done building ".$dir."specific_benchmark_".$f."_".$quant.".tsv\n";
#		system('mail -s "Hello World" lichtenj@gmail.com < '.$dir."specific_benchmark_".$f."_".$quant.".tsv");

		close COM;
	}
}
