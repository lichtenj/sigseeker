#!/usr/bin/perl

use strict;
use Math::Combinatorics;

# [SETUP] This needs to be set to the SigSeeker code directory
use lib '/Users/lichtenbergj/SigSeeker_CODE_05012013';
use SigSeeker;

#Required settings
my $project = shift or die;
my $org = shift or die;
my $peakcalling = "true";
my $encode_threshold = 10000;
my $comparison = shift or die;
my $filtering = shift or die;
my $partitioning = shift or die;
my $visualization = shift or die;
my $expression_correlation = shift or die;
my $conservation = shift or die;
my $conservationfilter = shift or die;
my $peakstatistics = shift or die;
my $email = shift;

#Setting up output directories
if(! -d $SigSeeker::PEAK_REPOSITORY.$project)
{
    system('mkdir '.$SigSeeker::PEAK_REPOSITORY.$project);
}

#Organism Formatting
my $organism;
if($org =~ /mm(\d+)/){$organism = 'mm';}
if($org =~ /hg(\d+)/){$organism = 'hs';}

#Relevant variables
my $setup;
my $setup_opt;
my @partitioning_files;
my @fully_partition_files;
my $encode_filtered;
my $max_peaks = 0;
my $peakcalling_results;

if($peakcalling eq "true")
{
	#Each read directory requires a setup.tsv file of the following format:
	#celltype_name	sample_ids (comma separated)	control_ids (comma separated)
    open(IN, $SigSeeker::READ_UPLOADS.$project.'/setup.tsv') or die "setup.tsv not specified for $project";
    while(my $record = <IN>)
    {
    	chomp($record);
        
    	my @entry = split(/\t+/,$record);
    
    	my @samples = split(/\,/,$entry[1]);
    	my @controls = split(/\,/,$entry[2]);
    
    	my $temp = $entry[0];

        my @breakdown = split(/[\-\_]/,$temp);
        push(@free_order,$breakdown[0]);
        
    	foreach my $sample (@samples)
    	{
    		foreach my $control (@controls)
		{
    			print "Set: ".$sample."\t".$control."\n";

			if($selected_tools->{"BROADPEAK"})
			{
		        	SigSeeker::BROADPEAK($organism,$SigSeeker::READ_UPLOADS.$project.'/'.$project.'_'.$sample.'.BAM', $SigSeeker::PEAK_REPOSITORY.$project.'/'.$entry[0].'_'.$sample.'_BROADPEAK_peaks.bed');
	    	    		system('cp '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$entry[0].'_'.$sample.'_BROADPEAK_peaks.bed '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$entry[0].'_'.$sample.'_'.$control.'_BROADPEAK_peaks.bed');
            			$peakcalling_results->{'BROADPEAK'}->{$sample}->{$control} = $SigSeeker::PEAK_REPOSITORY.$project.'/'.$entry[0].'_'.$sample.'_'.$control.'_BROADPEAK_peaks.bed';
    		    	}

	                if($selected_tools->{"ERANGE"})
        	        {
				$peakcalling_results->{'ERANGE'}->{$sample}->{$control} = SigSeeker::Erange($org,$SigSeeker::READ_UPLOADS.$project.'/'.$project.'_'.$sample.'.BAM', $SigSeeker::READ_UPLOADS.$project.'/'.$project.'_'.$control.'.BAM', $SigSeeker::PEAK_REPOSITORY.$project.'/'.$entry[0].'_'.$sample.'_'.$control.'_ERANGE_peaks.bed');
                	}

                	if($selected_tools->{"SWEMBLE"})
                	{
                    		$peakcalling_results->{'SWEMBLE'}->{$sample}->{$control} = SigSeeker::SWEMBLE($org,$SigSeeker::READ_UPLOADS.$project.'/'.$project.'_'.$sample.'.BAM', $SigSeeker::READ_UPLOADS.$project.'/'.$project.'_'.$control.'.BAM', $SigSeeker::PEAK_REPOSITORY.$project.'/'.$entry[0].'_'.$sample.'_'.$control.'_SWEMBLE_peaks.bed');
                	}

                	if($selected_tools->{"MACS"})
                	{
                    		$peakcalling_results->{'MACS'}->{$sample}->{$control} = SigSeeker::MACS($organism,$SigSeeker::READ_UPLOADS.$project.'/'.$project.'_'.$sample.'.BAM', $SigSeeker::READ_UPLOADS.$project.'/'.$project.'_'.$control.'.BAM', $SigSeeker::PEAK_REPOSITORY.$project.'/'.$entry[0].'_'.$sample.'_'.$control.'_MACS_peaks.bed');
                	}

                	if($selected_tools->{"MACS2"})
                	{
                    		$peakcalling_results->{'MACS2'}->{$sample}->{$control} = SigSeeker::MACS2($organism,$SigSeeker::READ_UPLOADS.$project.'/'.$project.'_'.$sample.'.BAM', $SigSeeker::READ_UPLOADS.$project.'/'.$project.'_'.$control.'.BAM', $SigSeeker::PEAK_REPOSITORY.$project.'/'.$entry[0].'_'.$sample.'_'.$control.'_MACS2_peaks.bed');
                	}

                	if($selected_tools->{"CISGENOME"})
                	{
                    		$peakcalling_results->{'CISGENOME'}->{$sample}->{$control} = SigSeeker::CISGENOME($org,$SigSeeker::READ_UPLOADS.$project.'/'.$project.'_'.$sample.'.BAM', $SigSeeker::READ_UPLOADS.$project.'/'.$project.'_'.$control.'.BAM', $SigSeeker::PEAK_REPOSITORY.$project.'/'.$entry[0].'_'.$sample.'_'.$control.'_CISGENOME_peaks.bed');
                	}

                	if($selected_tools->{"CCAT-TF"})
                	{
                    		$peakcalling_results->{'CCAT-TF'}->{$sample}->{$control} = SigSeeker::CCAT($organism,$SigSeeker::READ_UPLOADS.$project.'/'.$project.'_'.$sample.'.BAM', $SigSeeker::READ_UPLOADS.$project.'/'.$project.'_'.$control.'.BAM', $SigSeeker::PEAK_REPOSITORY.$project.'/'.$entry[0].'_'.$sample.'_'.$control.'_CCAT-TF_peaks.bed','CCAT-TF');
                	}

                	if($selected_tools->{"CCAT-Histone"})
                	{
                    		$peakcalling_results->{'CCAT-Histone'}->{$sample}->{$control} = SigSeeker::CCAT($organism,$SigSeeker::READ_UPLOADS.$project.'/'.$project.'_'.$sample.'.BAM', $SigSeeker::READ_UPLOADS.$project.'/'.$project.'_'.$control.'.BAM', $SigSeeker::PEAK_REPOSITORY.$project.'/'.$entry[0].'_'.$sample.'_'.$control.'_CCAT-Histone_peaks.bed','CCAT-Histone');
                	}

                if($selected_tools->{"SICER"})
                {
                    #push(@partitioning_files,$entry[0].'_'.$sample.'_'.$control.'_SICER_peaks.bed');
                    if(! -e $SigSeeker::PEAK_REPOSITORY.$project.'/'.$entry[0].'_'.$sample.'_'.$control.'_SICER_peaks.bed' || -z $SigSeeker::PEAK_REPOSITORY.$project.'/'.$entry[0].'_'.$sample.'_'.$control.'_SICER_peaks.bed')
        		    {
       				    print "---------------------------------------------------\n";
    	    			print "SICER\n";
    		    		print "---------------------------------------------------\n";

                        if(! -e $SigSeeker::READ_UPLOADS.$project.'/'.$project.'_'.$sample.'.bam.bed' || -z $SigSeeker::READ_UPLOADS.$project.'/'.$project.'_'.$sample.'.bam.bed')
            			{
    						print "Conversion Sample\n";
                            my $cmd = 'bamToBed -i '.$SigSeeker::READ_UPLOADS.$project.'/'.$project.'_'.$sample.'.BAM > '.$SigSeeker::READ_UPLOADS.$project.'/'.$project.'_'.$sample.'.bam.bed';
        					system($cmd);
        				}

                        if(! -e $SigSeeker::READ_UPLOADS.$project.'/'.$project.'_'.$control.'.bam.bed' || -z $SigSeeker::READ_UPLOADS.$project.'/'.$project.'_'.$control.'.bam.bed')
                		{
    						print "Conversion Control\n";
                            my $cmd = 'bamToBed -i '.$SigSeeker::READ_UPLOADS.$project.'/'.$project.'_'.$control.'.BAM > '.$SigSeeker::READ_UPLOADS.$project.'/'.$project.'_'.$control.'.bam.bed';
        					system($cmd);
        				}

        				system('sh '.$SigSeeker::TOOL_REPOSITORY.'SICER/SICER.sh '.$SigSeeker::READ_UPLOADS.$project.' '.$project.'_'.$sample.'.bam.bed '.$project.'_'.$control.'.bam.bed '.$SigSeeker::PEAK_REPOSITORY.$project.' '.$org.' 1 200 150 0.74 600 .01');
                        system('cp '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$project.'_'.$sample.'.bam-W200-G600-FDR.01-island.bed '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$entry[0].'_'.$sample.'_'.$control.'_SICER_peaks.bed');

                        system('rm '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$project.'_'.$sample.'.bam-W200-G600-FDR.01-island.bed');
                        system('rm '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$project.'_'.$sample.'.bam-W200-G600-FDR.01-islandfiltered.bed');
                        system('rm '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$project.'_'.$sample.'.bam-W200-G600-FDR.01-islandfiltered-normalized.bed');
                        system('rm '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$project.'_'.$sample.'.bam-W200-G600-islands-summary-FDR.01');
                        system('rm '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$project.'_'.$sample.'.bam-W200-G600-islands-summary');
                        system('rm '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$project.'_'.$sample.'.bam-W200-G600.scoreisland');
                        system('rm '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$project.'_'.$sample.'.bam-W200-normalized.wig');
                        system('rm '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$project.'_'.$sample.'.bam-W200.graph');
                        system('rm '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$project.'_'.$sample.'.bam-1-removed.bed');
                        system('rm '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$project.'_'.$control.'.bam-1-removed.bed');
    			    }

                    $peakcalling_results->{'SICER'}->{$sample}->{$control} = $SigSeeker::PEAK_REPOSITORY.$project.'/'.$entry[0].'_'.$sample.'_'.$control.'_SICER_peaks.bed';


                    remove_junk($SigSeeker::PEAK_REPOSITORY.$project.'/'.$entry[0].'_'.$sample.'_'.$control.'_SICER_peaks.bed');

                    my $cmd = 'wc -l '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$entry[0].'_'.$sample.'_'.$control.'_SICER_peaks.bed';
                    my $count = `$cmd`;
                    $count =~ /(\d+)/;
                    my $peaks = $1;
                    if($max_peaks < $peaks){$max_peaks = $peaks;}
                    if($peaks >= $encode_threshold)
                    {
                        $setup_opt->{$breakdown[0]}->{$breakdown[1]}->{$sample}->{$control}->{"SICER"} = $peaks;
                    }
                    else
                    {
                        $encode_filtered->{$sample}->{$control} = 1;
                    }
                    $setup->{$breakdown[0]}->{$breakdown[1]}->{$sample}->{$control}->{"SICER"} = $peaks;
                }                                
            }
    	}
    }
    close IN;

	print "Peak Cleaning\n";
    	foreach my $tool (keys %$selected_tools)
    	{
        	foreach my $set (keys %$setup)
        	{
            		foreach my $sub (keys %{$setup->{$set}})
            		{
                		my $subid = $sub;
                		if($sub){$subid = '-'.$sub;}
                
	                	foreach my $sample (keys %{$setup->{$set}->{$sub}})
				{
	                		foreach my $control (keys %{$setup->{$set}->{$sub}->{$sample}})
	                    		{
						my $cmd = 'wc -l '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$set.$subid.'_'.$sample.'_'.$control.'_'.$tool.'_peaks.bed';
						my $count = `$cmd`;
						$count =~ /(\d+)/;
						my $peaks = $1;
						if($max_peaks < $peaks){$max_peaks = $peaks;}
						if($peaks >= $encode_threshold)
						{
							$setup_opt->{$set}->{$sub}->{$sample}->{$control}->{$tool} = $peaks;
						}
						else
						{
							$encode_filtered->{$sample}->{$control} = 1;
						}
						$setup->{$set}->{$sub}->{$sample}->{$control}->{$tool} = $peaks;
					}
				}
			}
		}
	}
	SigSeeker::CleanPeaks($peakcalling_results,$SigSeeker::PEAK_REPOSITORY.$project.'/Called_Peaks.tsv');

    print "Peak Call Visualization\n";

    #Peak Call Visualization
    #my @set_order = ("HSC","CMP","MEP","Mep","GMP","CFUMEG","CFU-Meg","CFUMeg","MEG","Meg","CFUE","CFU-E","ERY","Ery","EB");
    my $order = "";
    #my @set_order = keys %$setup;
    
    my $vis_norm = "true";
    if($vis_norm eq "true")
    {
        my $sets;
        my $subs;
    	my $samples;
    	my $tools;
    	
    	foreach my $set (@free_order)
    	{
            my $bool = 0;
    		foreach my $sub (keys %{$setup->{$set}})
    		{
    			foreach my $sample (keys %{$setup->{$set}->{$sub}})
    			{
                    $bool = 1;
    				my $tool;
    				foreach my $control (keys %{$setup->{$set}->{$sub}->{$sample}})
    				{
    					foreach my $app (keys %$selected_tools)
    					{
                            if($setup->{$set}->{$sub}->{$sample}->{$control}->{$app})
                            {
        						$sets .= '"'.$set.'",';
        						$subs .= '"'.$sub.'",';
        						$tool->{$app} .= $setup->{$set}->{$sub}->{$sample}->{$control}->{$app}.',';			
        						$tools .= '"'.$app.'",';
                            }
    					}
    				}
    				foreach my $app (keys %$tool)
    				{
    					chop($tool->{$app});
    					if(scalar(keys %{$setup->{$set}->{$sub}->{$sample}}) > 1)
    					{
    						$samples .= 'c('.$tool->{$app}.'),';
                        }
    					else
    					{
    						$samples .= $tool->{$app}.',';
                        }
   				}
    			}
    		}
            if($bool == 1)
            {
                $order .= '"'.$set.'",';
            }
    	}
        chop $order;
        
        open(R, ">peaks.r");
        print R 'library(ggplot2)'."\n";
        print R 'peaks = data.frame('."\n";
        chop($sets);
        print R '  set = c('.$sets.'),'."\n";
        chop($subs);
        print R '  sub = c('.$subs.'),'."\n";
        chop($samples);
        print R '  peaks = c('.$samples.'),'."\n";
        chop($tools);
        print R '  tool = c('.$tools.')'."\n";
        print R ')'."\n";
        print R 'png(file = "'.$SigSeeker::PEAK_REPOSITORY.$project.'/Box_Part.png", width = 1024, height = 640, units = "px", pointsize = 12, bg="transparent")'."\n";
        if($order)
        {
            print R 'qplot(set, peaks, data = peaks,geom="boxplot",main = "Peak Calling") + aes(ymin=0,ymax='.($max_peaks + 1000).') + scale_y_continuous(breaks=seq(0,'.$max_peaks.',10000)) + scale_x_discrete(limits=c('.$order.')) + geom_jitter(position=position_jitter(w=0.1, h=0.1)) + ylab("#\\ Peaks") + xlab("Cell Types") + facet_grid(sub~tool) + theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1))'."\n";
            #print R 'qplot(set, peaks, data = peaks,geom="boxplot",main = "Peak Calling") + aes(ymin=0,ymax='.($max_peaks + 1000).') + scale_y_continuous(breaks=seq(0,'.$max_peaks.',10000)) + scale_x_discrete(limits=c('.$order.')) + geom_jitter(position=position_jitter(w=0.1, h=0.1)) + ylab("#\\ Peaks") + xlab("Cell Types") + facet_grid(sub~tool) + theme(text = element_text(size=20))'."\n";
        }
        else
        {
            print R 'qplot(set, peaks, data = peaks,geom="boxplot",main = "Peak Calling") + aes(ymin=0,ymax='.($max_peaks + 1000).') + scale_y_continuous(breaks=seq(0,'.$max_peaks.',10000)) + geom_jitter(position=position_jitter(w=0.1, h=0.1)) + ylab("#\\ Peaks") + xlab("Cell Types") + facet_grid(sub~tool) + theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1))'."\n";
            #print R 'qplot(set, peaks, data = peaks,geom="boxplot",main = "Peak Calling") + aes(ymin=0,ymax='.($max_peaks + 1000).') + scale_y_continuous(breaks=seq(0,'.$max_peaks.',10000)) + geom_jitter(position=position_jitter(w=0.1, h=0.1)) + ylab("#\\ Peaks") + xlab("Cell Types") + facet_grid(sub~tool) + theme(text = element_text(size=20))'."\n";
        }
        print R 'dev.off()'."\n";
        
        close R;
        
        system('R --vanilla < peaks.r 2> /dev/null > /dev/null');
    }

    my $optimized = "true";
    if($optimized eq "true")
    {
        my $order = "";
        my $sets;
        my $subs;
        my $samples;
    	my $tools;
    	
    	foreach my $set (@free_order)
    	{
            my $bool = 0;
    		foreach my $sub (keys %{$setup->{$set}})
    		{
    			foreach my $sample (keys %{$setup->{$set}->{$sub}})
    			{
                    #print scalar(keys %{$encode_filtered->{$sample}})."\t".(scalar(keys %{$setup->{$set}->{$sub}->{$sample}}) * 0.75)."\n";
                    if(scalar(keys %{$encode_filtered->{$sample}}) < scalar(keys %{$setup->{$set}->{$sub}->{$sample}}) * 0.75)
                    {
                        $bool = 1;
                        my $tool;
        				foreach my $control (keys %{$setup->{$set}->{$sub}->{$sample}})
        				{
        					foreach my $app (keys %$selected_tools)
        					{
                                if($setup->{$set}->{$sub}->{$sample}->{$control}->{$app})
                                {
            						$sets .= '"'.$set.'",';
            						$subs .= '"'.$sub.'",';
                            
            						$tool->{$app} .= $setup->{$set}->{$sub}->{$sample}->{$control}->{$app}.',';
            						$tools .= '"'.$app.'",';
                                }
        					}
        				}
        				foreach my $app (keys %$tool)
        				{
        					chop($tool->{$app});
        					if(scalar(keys %{$setup->{$set}->{$sub}->{$sample}}) > 1)
        					{
        						$samples .= 'c('.$tool->{$app}.'),';
                            }
        					else
        					{
        						$samples .= $tool->{$app}.',';
                            }
        				}
                    }
    			}
    		}
            if($bool == 1)
            {
                $order .= '"'.$set.'",';
            }
    	}
        chop $order;
        
        open(R, ">peaks_opt.r");
        print R 'library(ggplot2)'."\n";
        print R 'peaks = data.frame('."\n";
        chop($sets);
        print R '  set = c('.$sets.'),'."\n";
        chop($subs);
        print R '  sub = c('.$subs.'),'."\n";
        chop($samples);
        print R '  peaks = c('.$samples.'),'."\n";
        chop($tools);
        print R '  tool = c('.$tools.')'."\n";
        print R ')'."\n";
        print R 'png(file = "'.$SigSeeker::PEAK_REPOSITORY.$project.'/Box_Part_Optimized.png", width = 1024, height = 640, units = "px", pointsize = 12, bg="transparent")'."\n";
        if($order)
        {
            print R 'qplot(set, peaks, data = peaks,geom="boxplot",main = "Peak Calling") + aes(ymin=0,ymax='.($max_peaks + 1000).') + scale_y_continuous(breaks=seq(0,'.$max_peaks.',10000)) + scale_x_discrete(limits=c('.$order.')) + geom_jitter(position=position_jitter(w=0.1, h=0.1)) + ylab("#\\ Peaks") + xlab("Cell Types") + facet_grid(sub~tool) + theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1))'."\n";
        }
        else
        {
            print R 'qplot(set, peaks, data = peaks,geom="boxplot",main = "Peak Calling") + aes(ymin=0,ymax='.($max_peaks + 1000).') + scale_y_continuous(breaks=seq(0,'.$max_peaks.',10000)) + geom_jitter(position=position_jitter(w=0.1, h=0.1)) + ylab("#\\ Peaks") + xlab("Cell Types") + facet_grid(sub~tool) + theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1))'."\n";
        }
        print R 'dev.off()'."\n";
        
        close R;
        
        system('R --vanilla < peaks_opt.r 2> /dev/null > /dev/null');
    }

    close PEAKS;

    open(OUT, ">".$SigSeeker::PEAK_REPOSITORY.$project.'/Completed_Peaks.txt');
    print OUT 'Completed the comparison of peak predictions.';
    close OUT;
    
    system('mutt -s "Done with Peak Calling" '.$email.' -a '.$SigSeeker::PEAK_REPOSITORY.$project.'/Box_Part.png -a '.$SigSeeker::PEAK_REPOSITORY.$project.'/Box_Part_Optimized.png -a '.$SigSeeker::PEAK_REPOSITORY.$project.'/Called_Peaks.tsv -F /var/.muttrc < '.$SigSeeker::PEAK_REPOSITORY.$project.'/Completed_Peaks.txt');
}

if($comparison eq "true")
{
	my $replicate_files;
    
    	print "Replicate Joining\n";
    	open(OUT,">".$SigSeeker::PEAK_REPOSITORY.$project.'/Replicates.tsv');
    	open(FILES,">".$SigSeeker::PEAK_REPOSITORY.$project."/Files.tsv");
    
    	my $count = 0;
    	foreach my $tool (keys %$selected_tools)
    	{
        	foreach my $set (keys %$setup)
        	{
#            print $set."\n";
        	    	foreach my $sub (keys %{$setup->{$set}})
        	    	{
        	        	my $subid = $sub;
                		if($sub){$subid = '-'.$sub;}
                	
				if(scalar(keys %{$setup->{$set}->{$sub}}) == 1)
				{
					my $peaks;
					#Need to report the original count
	                		foreach my $sample (keys %{$setup->{$set}->{$sub}})
					{
	                			foreach my $control (keys %{$setup->{$set}->{$sub}->{$sample}})
	                		    	{
							$peaks += $setup->{$set}->{$sub}->{$sample}->{$control}->{$tool};
							#Convert the replicate-less file to a pseudo replicate
							system('cp '.$peakcalling_results->{$tool}->{$sample}->{$control}.' '.$SigSeeker::PEAK_REPOSITORY.$project.'/Replicates_'.$set.$subid.'_'.$tool.'_peaks.bed');
						}
					}
	                		print OUT $set."\t".$sub."\t".$tool."\t".$peaks."\n";
	
					#Need to populate replicate_file hash
					$replicate_files->{$count} = $SigSeeker::PEAK_REPOSITORY.$project.'/Replicates_'.$set.$subid.'_'.$tool.'_peaks.bed';
                		    	print FILES $count."\t".$replicate_files->{$count}."\n";
			                push(@partitioning_files, 'Replicates_'.$set.$subid.'_'.$tool.'_peaks.bed');
		                	$count++;
				}
				else
				{
		                	my $cmd = 'multiIntersectBed -cluster -i ';
		                	my $replicates = 0;
		                	foreach my $sample (keys %{$setup->{$set}->{$sub}})
		                	{
		                    		foreach my $control (keys %{$setup->{$set}->{$sub}->{$sample}})
		                    		{
		                            		$replicates++;
		                            		if(-e $SigSeeker::PEAK_REPOSITORY.$project.'/'.$set.$subid.'_'.$sample.'_'.$control.'_'.$tool.'_peaks.bed' && ! -z $SigSeeker::PEAK_REPOSITORY.$project.'/'.$set.$subid.'_'.$sample.'_'.$control.'_'.$tool.'_peaks.bed')
		                            		{
		                                		$cmd .= $SigSeeker::PEAK_REPOSITORY.$project.'/'.$set.$subid.'_'.$sample.'_'.$control.'_'.$tool.'_peaks.bed ';
		                            		}
		                    		}
		                	}
		                	$count++;
		                	$cmd .= ' > '.$SigSeeker::PEAK_REPOSITORY.$project.'/Replicates_'.$set.$subid.'_'.$tool.'_peaks.bed';
		                	if(! -e $SigSeeker::PEAK_REPOSITORY.$project.'/Replicates_'.$set.$subid.'_'.$tool.'_peaks.bed' || -z $SigSeeker::PEAK_REPOSITORY.$project.'/Replicates_'.$set.$subid.'_'.$tool.'_peaks.bed')
		                	{
		                		system($cmd);
		                	}
			               	if(! -z $SigSeeker::PEAK_REPOSITORY.$project.'/Replicates_'.$set.$subid.'_'.$tool.'_peaks.bed')
	       		         	{
	                		    	$replicate_files->{$count} = $SigSeeker::PEAK_REPOSITORY.$project.'/Replicates_'.$set.$subid.'_'.$tool.'_peaks.bed';
	                		    	print FILES $count."\t".$replicate_files->{$count}."\n";
	                		}
	                	
	                		open(REPOUT, ">".$SigSeeker::PEAK_REPOSITORY.$project.'/Replicates_'.$set.$subid.'_'.$tool.'_tmp.bed');
	                		open(REP, $SigSeeker::PEAK_REPOSITORY.$project.'/Replicates_'.$set.$subid.'_'.$tool.'_peaks.bed');
	                		while(my $recrep = <REP>)
	                		{
	                		    	chomp($recrep);
	                		    	my @tmprep = split(/\t/,$recrep);
	                		    	if($tmprep[3] == $replicates)
	                		    	{
	                		        	print REPOUT $recrep."\n";
	                		    	}
	                		}
	                		close REP;
	                		close REPOUT;
	                		system('mv '.$SigSeeker::PEAK_REPOSITORY.$project.'/Replicates_'.$set.$subid.'_'.$tool.'_tmp.bed '.$SigSeeker::PEAK_REPOSITORY.$project.'/Replicates_'.$set.$subid.'_'.$tool.'_peaks.bed');
			
			                push(@partitioning_files, 'Replicates_'.$set.$subid.'_'.$tool.'_peaks.bed');
	        		        
	        	        	$cmd = 'wc -l '.$SigSeeker::PEAK_REPOSITORY.$project.'/Replicates_'.$set.$subid.'_'.$tool.'_peaks.bed';
	        	        	my $count = `$cmd`;
	        	        	$count =~ /(\d+)/;
	        	        	my $peaks = $1;
	        	        	print OUT $set."\t".$sub."\t".$tool."\t".$peaks."\n";
				}        
            		}
        	}
    	}
    	close OUT;

    	print "Replicatation Correlations\n";
    	open(OUT,">".$SigSeeker::PEAK_REPOSITORY.$project.'/ReplicateCorrelations.tsv');
    	foreach my $set (keys %$setup)
    	{
    	    	foreach my $sub (keys %{$setup->{$set}})
    	    	{
    	        	my $subid = $sub;
    		        if($sub){$subid = '-'.$sub;}
            
            		foreach my $sample (keys %{$setup->{$set}->{$sub}})
            		{
                		foreach my $control (keys %{$setup->{$set}->{$sub}->{$sample}})
                		{
                    			#Pairwise Comparison
                    			foreach my $set2 (keys %$setup)
                    			{
                        			foreach my $sub2 (keys %{$setup->{$set2}})
                        			{
                        				my $subid2 = $sub2;
                        	    			if($sub2){$subid2 = '-'.$sub2;}
                        	    
                        	    			foreach my $sample2 (keys %{$setup->{$set2}->{$sub2}})
                        	    			{
                        	        			foreach my $control2 (keys %{$setup->{$set2}->{$sub2}->{$sample2}})
                        	        			{
                                    					foreach my $tool (keys %$selected_tools)
                                   					{
                    								if(! ($set eq $set2 && $sample == $sample2 && $control == $control2))
                    								{
                                            						my $cmd = 'perl '.$SigSeeker::MAIN_DIRECTORY.'correlate.pl '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$set.$subid.'_'.$sample.'_'.$control.'_'.$tool.'_peaks.bed '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$set2.$subid2.'_'.$sample2.'_'.$control2.'_'.$tool.'_peaks.bed '.$SigSeeker::PEAK_REPOSITORY.$project.'/ReplicateCorrelation_'.$set.$subid.'_'.$sample.'_'.$control.'_vs_'.$set2.$subid2.'_'.$sample2.'_'.$control2.'_'.$tool;
                                            						if(! -e $SigSeeker::PEAK_REPOSITORY.$project.'/ReplicateCorrelation_'.$set.$subid.'_'.$sample.'_'.$control.'_vs_'.$set2.$subid2.'_'.$sample2.'_'.$control2.'_'.$tool.'.tsv' || -z $SigSeeker::PEAK_REPOSITORY.$project.'/ReplicateCorrelation_'.$set.$subid.'_'.$sample.'_'.$control.'_vs_'.$set2.$subid2.'_'.$sample2.'_'.$control2.'_'.$tool.'.tsv')
                                            						{
                        									system($cmd);
                                            						}
                                            						open(TMP, $SigSeeker::PEAK_REPOSITORY.$project.'/ReplicateCorrelation_'.$set.$subid.'_'.$sample.'_'.$control.'_vs_'.$set2.$subid2.'_'.$sample2.'_'.$control2.'_'.$tool.'.tsv');
                                            						while(my $record = <TMP>)
                                            						{
                                                						chomp($record);
                        									$record =~ /^(\d+)\t([\d\.]+)\t([\d\.]+)\t\[\d+\]\s*([\-\d\.]+|NA)/;
			                						        print OUT $set.$subid."\t".$set2.$subid2."\t".$sample."\t".$sample2."\t".$control."\t".$control2."\t".$tool."\t".$1."\t".$2."\t".$3."\t".$4."\n";
												last;
                                            						}
                                            						close TMP;
                    								}
                                    					}
                                				}
                            				}
                        			}
                    			}
                		}
            		}
        	}
    	}
    	close OUT;

    	print "Comparisons\n";
    	my $cmd = 'multiIntersectBed -header -names ';
    	foreach my $file (sort {$a <=> $b} keys %$replicate_files)
    	{
    		if(-e $replicate_files->{$file} && ! -z $replicate_files->{$file})
        	{
            		$cmd .= $file.' ';
        	}
    	}
    	$cmd .= '-cluster -i ';
    	foreach my $file (sort {$a <=> $b} keys %$replicate_files)
    	{
        	if(-e $replicate_files->{$file} && ! -z $replicate_files->{$file})
        	{
            		$cmd .= $replicate_files->{$file}.' ';
        	}
    	}

    	$cmd .= '> '.$SigSeeker::PEAK_REPOSITORY.$project.'/Comparisons.bed';

    #print $cmd."\n";
    if(! -e $SigSeeker::PEAK_REPOSITORY.$project.'/Comparisons.bed' || -z $SigSeeker::PEAK_REPOSITORY.$project.'/Comparisons.bed')
    {
        system($cmd);
    }
    
    my $comps;
    open(COMP, $SigSeeker::PEAK_REPOSITORY.$project.'/Comparisons.bed');
    
    my $header = <COMP>;
    chomp($header);
    my @head = split(/\t/,$header);
    
    my $count = 0;
    
    while(my $rec = <COMP>)
    {
        $count++;
        
        chomp($rec);
        my @tmp = split(/\t/,$rec);
        $tmp[4] =~ s/\,/\_/g;
        
        #Remove useless comparisons
        my @files = split(/\_/,$tmp[4]);
        my $bool = 0;
        my $tmp_hash;
        my $tmp_tools;
        foreach my $fid (@files)
        {
            my $file = $replicate_files->{$fid};
            my @tmpdir = split(/\//,$file);
            my @tmpfile = split(/\_/,$tmpdir[-1]);
            $tmp_hash->{$tmpfile[1]}->{$tmpfile[2]} = 1;
            $tmp_tools->{$tmpfile[2]} = 1;
        }
        
#if($count == 10){exit;}

        foreach my $tool (keys %$tmp_tools)
        {
            foreach my $set (keys %$tmp_hash)
            {
#                print ">>>>".$tool."\t".$set."\n";
                if(! $tmp_hash->{$set}->{$tool})
                {
                    $bool = 1;
                }
            }
        }

        #Store usefule entries
        if($bool == 0)
        {
#            print $tmp[4]."\n";
            my @tmpbool = split(/\_/,$tmp[4]);
#            foreach my $fbool (@tmpbool)
#            {
#                print "\t".$replicate_files->{$fbool}."\n";    
#            }
            $comps->{'Comparison_'.$tmp[4]}->{$rec} = 1;
            #print $filenames->{$fn}."\n".$bool."\n";
        }
    }
    close COMP;

    open(COMPOUT,">".$SigSeeker::PEAK_REPOSITORY.$project."/Comparisons.tsv");
    foreach my $fn (keys %$comps)
    {
        print COMPOUT $fn."\t".scalar(keys %{$comps->{$fn}})."\n";
        
        #push(@partitioning_files, $fn.'.bed');
        
        open(COMPBED,">".$SigSeeker::PEAK_REPOSITORY.$project.'/'.$fn.'.bed');
        foreach my $record (keys %{$comps->{$fn}})
        {
            print COMPBED $record."\n";
        }
        close COMPBED;
    }
    
    open(OUT, ">".$SigSeeker::PEAK_REPOSITORY.$project.'/Completed_Comparisons.txt');
    print OUT 'Completed the comparison of peak predictions.';
    close OUT;
    system('mutt -s "Done with Peak Comparison" '.$email.' -a '.$SigSeeker::PEAK_REPOSITORY.$project.'/Comparisons.tsv -a '.$SigSeeker::PEAK_REPOSITORY.$project.'/ReplicateCorrelations.tsv -a '.$SigSeeker::PEAK_REPOSITORY.$project.'/Replicates.tsv -F /var/.muttrc < '.$SigSeeker::PEAK_REPOSITORY.$project.'/Completed_Comparisons.txt');
}

if($filtering eq "true")
{
    print "Filtering\n";
    
    open(OUT,">".$SigSeeker::PEAK_REPOSITORY.$project."/Filters.tsv");
    opendir(DIR, $SIGSEEKER::FILTER_REPOSITORY);
    #print OUT "File\t";
    #while(my $filter = readdir(DIR))
    #{
    #    if($filter =~ /^\./){next;}
    #    print OUT 'Filter_'.$filter."\t";
    #}
    #print OUT "\n";
    
    foreach my $file (@partitioning_files)
    {
        print OUT 'Filter_'.$file."\t";
        my $cmd = 'perl '.$SigSeeker::MAIN_DIRECTORY.'peakfiltering.pl '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$file.' '.$org.' '.$SigSeeker::PEAK_REPOSITORY.$project.'/Filter_'.$file;
        print OUT `$cmd`;
        print OUT "\n";
    }
    close OUT;

    open(OUT, ">".$SigSeeker::PEAK_REPOSITORY.$project.'/Completed_Filtering.txt');
    print OUT 'Completed the filtering of peak predictions.';
    close OUT;
    system('mutt -s "Done with Peak Filtering" '.$email.' -a '.$SigSeeker::PEAK_REPOSITORY.$project.'/Filters.tsv -F /var/.muttrc < '.$SigSeeker::PEAK_REPOSITORY.$project.'/Completed_Filtering.txt');
}

if($partitioning eq "true")
{
    print "Unique Partitioning\n";
    open(OUT,">".$SigSeeker::PEAK_REPOSITORY.$project."/Partitions.tsv");
    open(EXOUT,">".$SigSeeker::PEAK_REPOSITORY.$project."/Expression.tsv");
    
    my $sets;
    my $subs;
    my $peaks;
    my $tools;
    my $parts;
    
    foreach my $file (@partitioning_files)
    {
        print OUT $file."\t";
        my @tmp = split(/\_/,$file);
        my @info = split(/\-/,$tmp[1]);
        
        my $set = $info[0]; #Celltype
        $sets .= '"'.$set.'",';
        my $sub = $info[1]; #Antibody
        $subs .= '"'.$sub.'",';
        my $tool = $tmp[2]; #Tool
        $tools .= '"'.$tool.'",';
        
        my $cmd = "";
        if($file =~ /^Filter/)
        {
            $cmd = 'perl '.$SigSeeker::MAIN_DIRECTORY.'overlap_analysis.pl '.$SigSeeker::PEAK_REPOSITORY.$project.'/Filter_'.$file.' '.$SigSeeker::PEAK_REPOSITORY.$org.'_refFlat.txt '.$SigSeeker::PEAK_REPOSITORY.$project.'/Partition_'.$file;
        }
        else
        {
            $cmd = 'perl '.$SigSeeker::MAIN_DIRECTORY.'overlap_analysis.pl '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$file.' '.$SigSeeker::PEAK_REPOSITORY.$org.'_refFlat.txt '.$SigSeeker::PEAK_REPOSITORY.$project.'/Partition_'.$file;
        }
        my $counts = `$cmd`;
        print OUT $counts;
        chomp($counts);
        $counts =~ s/\t/\,/g;
        $peaks .= $counts.',';
        
        $parts .= '"Upstream","Promoter","RefSeq","Downstream","Intergenic","Generic Upstream","Generic Promoter","Generic RefSeq","Generic Downstream","Generic Intergenic",';
#        if($expression_correlation = "true")
#        {
            my @partitions = ("Upstream","Promoter","RefSeq","Downstream","Intergenic","Generic Upstream","Generic Promoter","Generic RefSeq","Generic Downstream","Generic Intergenic");
            foreach my $partition (@partitions)
            {
                push(@fully_partition_files,'Partition_'.$file.'_'.$partition.'.bed');
#                my $input = $SigSeeker::PEAK_REPOSITORY.$project.'/Partition_'.$file.'_'.$partition.'.bed';
#                my $cmdex = 'perl '.$SigSeeker::MAIN_DIRECTORY.'expression_analysis.pl '.$input.' '.$project.' '.$partition;
#                print EXOUT `$cmdex`;
            }
#        }
    }
    close OUT;
    close EXOUT;
    
    #my @set_order = ("HSC","CMP","MEP","GMP","CFUMEG","MEG","CFUE","ERY");
    my $order = "";
    
    open(R, ">part_peaks.r");
    print R 'library(ggplot2)'."\n";
    print R 'peaks = data.frame('."\n";
    chop($sets);
    print R '  set = c('.$sets.'),'."\n";
    chop($subs);
    print R '  sub = c('.$subs.'),'."\n";
    chop($peaks);
    print R '  samples = c('.$peaks.'),'."\n";
    chop($tools);
    print R '  tool = c('.$tools.'),'."\n";
    chop($parts);
    print R '  part = c('.$parts.')'."\n";
    print R ')'."\n";
    print R 'png(file = "'.$SigSeeker::PEAK_REPOSITORY.$project.'/Partitions.png", width = 1024, height = 1024, units = "px", pointsize = 12, bg="transparent")'."\n";
    if($order)
    {
        print R 'qplot(set, samples, data = sample,geom="bar",main = "Peak Calling",fill=samples) + aes(ymin=0) + scale_x_discrete(limits=c('.$order.')) + ylab("#\\ Peaks") + xlab("Cell Types") + facet_grid(sub~tool~part) + theme(text = element_text(size=20))'."\n";
    }
    else
    {
        print R 'qplot(set, samples, data = peaks,geom="bar",main = "Peak Calling",fill=as.factor(samples)) + aes(ymin=0) + ylab("#\\ Peaks") + xlab("Cell Types") + facet_grid(sub~tool~part) + theme(text = element_text(size=20))'."\n";
    }
    print R 'dev.off()'."\n";
    
    close R;
        
    system('R --vanilla < part_peaks.r 2> /dev/null > /dev/null');

    open(OUT, ">".$SigSeeker::PEAK_REPOSITORY.$project.'/Completed_Partitioning.txt');
    print OUT 'Completed the partitioning of peak predictions.';
    close OUT;
    system('mutt -s "Done with Peak Partitioning" '.$email.' -a '.$SigSeeker::PEAK_REPOSITORY.$project.'/Partitions.png -a '.$SigSeeker::PEAK_REPOSITORY.$project.'/Partitions.tsv -F /var/.muttrc < '.$SigSeeker::PEAK_REPOSITORY.$project.'/Completed_Partitioning.txt');
}

my @conservation_files;
if($conservation eq "true")
{
    print "LiftOver Conservation\n";
    open(OUT,">".$SigSeeker::PEAK_REPOSITORY.$project."/Conservation_LiftOver.tsv");
    foreach my $file (@fully_partition_files)
    {
        print OUT $file."\t";
        my @tmp = split(/\_/,$file);
        my @info = split(/\-/,$tmp[1]);
        
        my $cmd = "";
        if($file =~ /^Filter/)
        {
            if(! -e $SigSeeker::PEAK_REPOSITORY.$project.'/LiftOver_Human_Filter_'.$file || -z $SigSeeker::PEAK_REPOSITORY.$project.'/LiftOver_Human_Filter_'.$file)
            {
                system('perl '.$SigSeeker::MAIN_DIRECTORY.'convertmulti2bed.pl '.$SigSeeker::PEAK_REPOSITORY.$project.'/Filter_'.$file.' > temp.bed');
                $cmd = 'liftOver temp.bed /usr/local/share/mm9.hg19.all.chain '.$SigSeeker::PEAK_REPOSITORY.$project.'/LiftOver_Human_Filter_'.$file.' '.$SigSeeker::PEAK_REPOSITORY.$project.'/LiftOver_Unmapped_Filter_'.$file.' > temp 2> temp2';
                system($cmd);
            }
            push(@conservation_files,$SigSeeker::PEAK_REPOSITORY.$project.'/LiftOver_Human_Filter_'.$file);
            
            my $mapped_count = 'wc -l '.$SigSeeker::PEAK_REPOSITORY.$project.'/LiftOver_Human_Filter_'.$file;
            $mapped_count =~ /^(\d+)\t/;
            $mapped_count = $1;
            my $unmapped_count = 'wc -l '.$SigSeeker::PEAK_REPOSITORY.$project.'/LiftOver_Unmapped_Filter_'.$file;
            $unmapped_count =~ /^(\d+)\t/;
            $unmapped_count = $1;#        }
            
            print OUT $mapped_count."\t";
            print OUT $unmapped_count."\n";
        }
        else
        {
            if(! -e $SigSeeker::PEAK_REPOSITORY.$project.'/LiftOver_Human_'.$file || -z $SigSeeker::PEAK_REPOSITORY.$project.'/LiftOver_Human_'.$file)
            {
                system('perl '.$SigSeeker::MAIN_DIRECTORY.'convertmulti2bed.pl '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$file.' > temp.bed');
                $cmd = 'liftOver temp.bed /usr/local/share/mm9.hg19.all.chain '.$SigSeeker::PEAK_REPOSITORY.$project.'/LiftOver_Human_'.$file.' '.$SigSeeker::PEAK_REPOSITORY.$project.'/LiftOver_Unmapped_'.$file.' > temp 2> temp2';
                system($cmd);
            }
            
            push(@conservation_files,$SigSeeker::PEAK_REPOSITORY.$project.'/LiftOver_Human_'.$file);
            
            my $mapped_count = 'wc -l '.$SigSeeker::PEAK_REPOSITORY.$project.'/LiftOver_Human_'.$file;
            $mapped_count =~ /^(\d+)\t/;
            $mapped_count = $1;
            my $unmapped_count = 'wc -l '.$SigSeeker::PEAK_REPOSITORY.$project.'/LiftOver_Unmapped_'.$file;
            $unmapped_count =~ /^(\d+)\t/;
            $unmapped_count = $1;
            
            print OUT $mapped_count."\t";
            print OUT $unmapped_count."\n";
        }
    }
    close OUT;

    print "Genic Conservation\n";
    open(OUT,">".$SigSeeker::PEAK_REPOSITORY.$project."/Conservation_Genic.tsv");
    foreach my $file (@fully_partition_files)
    {
        print OUT $file."\t";
        my $cmd = "";

        if(! -e $SigSeeker::PEAK_REPOSITORY.$project.'/Genic_Human_'.$file || -z $SigSeeker::PEAK_REPOSITORY.$project.'/Genic_Human_'.$file)
        {
            $cmd = 'perl '.$SigSeeker::MAIN_DIRECTORY.'orthology_lookup.pl '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$file.' '.$SigSeeker::MAIN_DIRECTORY.'Human_Mouse_Orthologs.tsv /usr/local/share/hg19_refFlat.txt '.$SigSeeker::PEAK_REPOSITORY.$project.'/Genic_Human_'.$file;
            #print $cmd."\n";
            system($cmd);
        }
        
        push(@conservation_files,$SigSeeker::PEAK_REPOSITORY.$project.'/Genic_Human_'.$file);
        
        my $cmdcount = 'perl '.$SigSeeker::MAIN_DIRECTORY.'convert_conspartition2genelist.pl '.$SigSeeker::PEAK_REPOSITORY.$project.'/Genic_Human_'.$file.' | wc -l';
        my $mapped_count = `$cmdcount`;
        $mapped_count =~ /^(\d+)/;
    
        my $main_cmdcount = 'perl '.$SigSeeker::MAIN_DIRECTORY.'convert_partition2genelist.pl '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$file.' | wc -l';
        my $main_mapped_count = `$cmdcount`;
        $main_mapped_count =~ /^(\d+)/;
        $main_mapped_count = $1;
        
        print OUT $mapped_count."\t";
        print OUT $main_mapped_count."\n";
    }
    close OUT;
    
    open(OUT, ">".$SigSeeker::PEAK_REPOSITORY.$project.'/Completed_Conservation.txt');
    print OUT 'Completed the conservation analysis of peak predictions.';
    close OUT;
    system('mutt -s "Done with Peak Conservation" '.$email.' -a '.$SigSeeker::PEAK_REPOSITORY.$project.'/Conservation_LiftOver.tsv -a '.$SigSeeker::PEAK_REPOSITORY.$project.'/Conservation_Genic.tsv -F /var/.muttrc < '.$SigSeeker::PEAK_REPOSITORY.$project.'/Completed_Conservation.txt');
}

if($conservationfilter eq "true")
{
    print "Conservation Filter\n";
    open(OUT,">".$SigSeeker::PEAK_REPOSITORY.$project."/Conservation_Filtering.tsv");
    foreach my $file (@conservation_files)
    {
        my $cmd = "";
        if($org eq 'mm9')
        {
            $cmd = 'perl '.$SigSeeker::MAIN_DIRECTORY.'orthology_filter.pl '.$file.' Human';
        }
        else
        {
            $cmd = 'perl '.$SigSeeker::MAIN_DIRECTORY.'orthology_filter.pl '.$file.' Mouse';
        }
        print OUT `$cmd`;
    }
    close OUT;
    open(OUT, ">".$SigSeeker::PEAK_REPOSITORY.$project.'/Completed_Conservation_Filtering.txt');
    print OUT 'Completed the filtering analysis of conserved regions in the peak predictions.';
    close OUT;
    system('mutt -s "Done with Peak Conservation Filtering" '.$email.' -a '.$SigSeeker::PEAK_REPOSITORY.$project.'/Conservation_Filtering.tsv -F /var/.muttrc < '.$SigSeeker::PEAK_REPOSITORY.$project.'/Completed_Conservation_Filtering.txt');
}

if($peakstatistics eq "true")
{
    my $cmd = 'perl '.$SigSeeker::MAIN_DIRECTORY.'peakstat_pipeline.pl '.$SigSeeker::PEAK_REPOSITORY.$project.' /usr/local/share/'.$org.'_refFlat.txt '.$org.' '.$SigSeeker::PEAK_REPOSITORY.$project.'/Peak_Statistics.tsv';

	print $cmd."\n";    
    print `$cmd`;
    
    open(OUT, ">".$SigSeeker::PEAK_REPOSITORY.$project.'/Completed_PeakStatistics.txt');
    print OUT 'Completed the demographic characterization for the peak predictions.';
    close OUT;
    system('mutt -s "Done with Peak Statistics" '.$email.' -a '.$SigSeeker::PEAK_REPOSITORY.$project.'/Peak_Statistics.tsv -F /var/.muttrc < '.$SigSeeker::PEAK_REPOSITORY.$project.'/Completed_PeakStatistics.txt');
}
print "Done\n";

open(OUT, ">".$SigSeeker::PEAK_REPOSITORY.$project.'/Completed_Analysis.txt') or die "Cannot open";
print OUT 'Completed the peak prediction process. You can view the results ';
if($ENV{SERVER_NAME} eq '165.112.60.208')
{
    print OUT '<a href="http://'.$ENV{SERVER_NAME}.'/sigseeker-cgi/SigSeeker_Project_Overview_New.pl?project='.$project_id.'">here</a>';
}
else
{
    print OUT '<a href="http://sigseeker.org/SigSeeker_Project_Overview_New.pl?project='.$project_id.'">here</a>';
}
close OUT;
system('mutt -s "Done with Peak Analysis Pipeline" '.$email.' -F /var/.muttrc < '.$SigSeeker::PEAK_REPOSITORY.$project.'/Completed_Analysis.txt');

sub remove_junk
{
    my $inputfile = shift or die;
    open(RIN, $inputfile);
    open(ROUT,">CLEANTMP") or die "Cannot open temporary cleaning file";
    while(my $record = <RIN>)
    {
        if($record =~ /^chr[\dXY]+\t/)
        {
            print ROUT $record;
        }
    }
    close ROUT;
    close RIN;
    system ('mv CLEANTMP '.$inputfile);
}
