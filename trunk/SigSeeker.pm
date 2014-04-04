package SigSeeker;

use SigSeeker_Help;

$MYSQL_USER = 'root';
$MYSQL_PASS = 'Jen$21513';
$GOOGLE_USER = 'lichtenberg@msseeker.org';
$GOOGLE_PASS = 'tb5013li';
$READ_UPLOADS = '';
$PEAK_REPOSITORY = '';
$RNA_REPOSITORY = '';
$TOOL_REPOSITORY = '';

if($ENV{SERVER_NAME} ne '165.112.60.208')
{
    $READ_UPLOADS = '/media/darklichti/Tython/BAM_Repository/';
    $PEAK_REPOSITORY = '/media/darklichti/Alexandria/Peak_Repository/';
    $RNA_REPOSITORY = '/media/darklichti/Alexandria/RNA_Repository/';
    $TOOL_REPOSITORY = '/var/www/follow/Tool_Peak_Calling_Ensemble/Tools/';
    $FILTER_REPOSITORY = '/media/darklichti/Alexandria/Filters/';
    $MAIN_DIRECTORY = '/home/darklichti/Dropbox/SigSeeker_CODE_05012013/';
    $READ_REPOSITORY = '/media/darklichti/Tython/Read_Repository/'; 
    $SBR_REPOSITORY = '/home/darklichti/Dropbox/SBR/';
}
else
{
    $READ_UPLOADS = '/Volumes/Mac_HD_RAID/BAM_Repository/';    
    $PEAK_REPOSITORY = '/Volumes/Mac_HD_RAID/PEAK_Repository/';
    $RNA_REPOSITORY = '/Volumes/Mac_HD_RAID/RNA_Repository/';
    $TOOL_REPOSITORY = '/Volumes/Mac_HD_RAID/TOOL_Repository/';
    $FILTER_REPOSITORY = '/Volumes/Mac_HD_RAID/FILTER_Repository/';
    $MAIN_DIRECTORY = '/Users/lichtenbergj/Dropbox/SigSeeker_CODE_05012013/';
    $READ_REPOSITORY = '/Volumes/Mac_HD_RAID/Read_Repository/';
    $SBR_REPOSITORY = '/Users/lichtenbergj/Dropbox/SBR/';
}

#Peak Calling Tools

sub PEAKQUENCER
{
    	print "---------------------------------------------------\n";
    	print "Peakquencer\n";
	print "---------------------------------------------------\n";

	my $org = shift or die;
    	my $sample = shift or die;
    	my $control = shift or die;
    	my $peakoutput = shift or die;

    	$sample =~ /\_(\d+)\.BAM$/;
    	my $sample_id = $1;
    	$control =~ /\_(\d+)\.BAM$/;
    	my $control_id = $1;
    
    	my $sample_bedfile = $sample;
    	$sample_bedfile =~ s/\.BAM$/\.bam\.bed/;
    	my $control_bedfile = $control;
    	$control_bedfile =~ s/\.BAM$/\.bam\.bed/;

	print "Peak Calling - ";
    	if(! -e $peakoutput || -z $peakoutput)
    	{
		system('perl /home/darklichti/Dropbox/Peakquencer/simple.pl '.$sample_bedfile.' '.$control_bedfile.' 10 10 '.$peakoutput);
		print "done\n";
		exit;
	}
	else
	{
		print "skipped\n";
	}

	return $peakoutput;
}

sub BROADPEAK
{
    	print "---------------------------------------------------\n";
    	print "BroadPeak\n";
	print "---------------------------------------------------\n";

	my $org = shift or die;
    	my $sample = shift or die;
    	my $peakoutput = shift or die;

    	$sample =~ /\_(\d+)\.BAM$/;
    	my $sample_id = $1;
    
    	my $sample_bedfile = $sample;
    	$sample_bedfile =~ s/\.BAM$/\.bam\.bed/;

    	if(! -e $peakoutput || -z $peakoutput)
    	{
	    	system('sort -k1,1 -k2,2g -o '.$sample_bedfile.'.sorted.bed '.$sample_bedfile);
	    	system('mv '.$sample_bedfile.'.sorted.bed '.$sample_bedfile);
	    	system('genomeCoverageBed -i '.$sample_bedfile.' -g /media/darklichti/Tython/BAM_Repository/'.$org.'.chrom.sizes -bg > '.$sample_bedfile.'graph');
	    	system($SigSeeker::TOOL_REPOSITORY.'BroadPeak/BroadPeak -i '.$sample_bedfile.'graph -m '.$sample_id.'_BROADPEAK -t unsupervised');
	    	system('mv '.$sample_id.'_BROADPEAK/'.$sample_id.'_BROADPEAK_broad_peak_unsupervised/'.$sample_id.'_BROADPEAK_broad_peak_unsupervised.bed '.$peakoutput);
    	}

	return $peakoutput;
}

sub CISGENOME
{
    	print "---------------------------------------------------\n";
	print "CISGENOME\n";
	print "---------------------------------------------------\n";

	my $org = shift or die;
    	my $sample = shift or die;
    	my $control = shift or die;
    	my $peakoutput = shift or die;

    	$sample =~ /\_(\d+)\.BAM$/;
    	my $sample_id = $1;
    	$control =~ /\_(\d+)\.BAM$/;
    	my $control_id = $1;
    
    	my $sample_bedfile = $sample;
    	$sample_bedfile =~ s/\.BAM$/\.bam\.bed/;
    	my $control_bedfile = $control;
    	$control_bedfile =~ s/\.BAM$/\.bam\.bed/;

    	my $sample_alnfile = $sample;
    	$sample_alnfile =~ s/\.BAM$/\.bam\.aln/;
    	my $control_alnfile = $control;
    	$control_alnfile =~ s/\.BAM$/\.bam\.aln/;

	my $fname = $sample_bedfile.'_MACS';
	$fname =~ s/\.bam\.bed/\_$control_id/;

    	if(! -e $sample_alnfile || -z $sample_alnfile)
	{
        	system($SigSeeker::TOOL_REPOSITORY.'cisgenome_project/bin/file_bed2aln -i '.$sample_bedfile.' -o '.$sample_alnfile);
	}

    	if(! -e $control_alnfile || -z $control_alnfile)
	{
        	system($SigSeeker::TOOL_REPOSITORY.'cisgenome_project/bin/file_bed2aln -i '.$control_bedfile.' -o '.$control_alnfile);
	}

    	if(! -e $peakoutput || -z $peakoutput)
    	{
        	open(CGOUT,">"."CG_".$sample_id.'_'.$control_id);
        	print CGOUT $sample_alnfile."\t1\n";
        	print CGOUT $control_alnfile."\t0\n";
        	close CGOUT;

		system($SigSeeker::TOOL_REPOSITORY.'cisgenome_project/bin/seqpeak -i CG_'.$sample_id.'_'.$control_id.' -d TMP -o CG_TMP_'.$sample_id.'_'.$control_id);
        	system($SigSeeker::TOOL_REPOSITORY.'cisgenome_project/bin/file_cod2bed -i TMP/CG_TMP_'.$sample_id.'_'.$control_id.'_peak.cod -o '.$sample_id.'_'.$control_id.'_CISGENOMEunsorted_peaks.bed');
        	system('sort -k1,1 -k2,2g -o '.$peakoutput.' '.$sample_id.'_'.$control_id.'_CISGENOMEunsorted_peaks.bed');
#        	system('rm '.$sample_id.'_'.$control_id.'_CISGENOMEunsorted_peaks.bed');
#        	system('rm CG_'.$sample_id.'_'.$control_id);

		exit;
        	print "done\n";
    	}
    	else
    	{
        	print "skipped\n";    
    	}
    
    	return $peakoutput;
}

sub MACS
{
    print "---------------------------------------------------\n";
    print "MACS 1.4\n";
    print "---------------------------------------------------\n";

	my $org = shift or die;
    my $sample = shift or die;
    my $control = shift or die;
    my $peakoutput = shift or die;

    $sample =~ /\_(\d+)\.BAM$/;
    my $sample_id = $1;
    $control =~ /\_(\d+)\.BAM$/;
    my $control_id = $1;
    
    my $sample_bedfile = $sample;
    $sample_bedfile =~ s/\.BAM$/\.bam\.bed/;
    my $control_bedfile = $control;
    $control_bedfile =~ s/\.BAM$/\.bam\.bed/;

	my $fname = $sample_bedfile.'_MACS';
	$fname =~ s/\.bam\.bed/\_$control_id/;


    print "Peak Calling - ";
    if(! -e $peakoutput || -z $peakoutput)
    {
    	system('macs14 -f BED -t '.$sample_bedfile.' -c '.$control_bedfile.' -p 1e-5 -g '.$org.' -n '.$fname);
    	system('cp '.$fname.'_peaks.bed '.$peakoutput);
    	print "done\n";
    }
    else
    {
        print "skipped\n";    
    }
    return $peakoutput;
}

sub MACS2
{
    print "---------------------------------------------------\n";
    print "MACS 2\n";
    print "---------------------------------------------------\n";

	my $org = shift or die;
    my $sample = shift or die;
    my $control = shift or die;
    my $peakoutput = shift or die;

    $sample =~ /\_(\d+)\.BAM$/;
    my $sample_id = $1;
    $control =~ /\_(\d+)\.BAM$/;
    my $control_id = $1;
    
    my $sample_bedfile = $sample;
    $sample_bedfile =~ s/\.BAM$/\.bam\.bed/;
    my $control_bedfile = $control;
    $control_bedfile =~ s/\.BAM$/\.bam\.bed/;

	my $fname = $sample_bedfile.'_MACS2';
	$fname =~ s/\.bam\.bed/\_$control_id/;


    print "Peak Calling - ";
    if(! -e $peakoutput || -z $peakoutput)
    {
	system('macs2 callpeak -f BED -t '.$sample_bedfile.' -c '.$control_bedfile.' -p 1e-5 -g '.$org.' -n '.$fname);
	system('cp '.$fname.'_peaks.bed '.$peakoutput);
	print "done\n";
    }
    else
    {
        print "skipped\n";    
    }
    return $peakoutput;
}

sub SWEMBLE
{
    print "---------------------------------------------------\n";
    print "SWEMBLE\n";
    print "---------------------------------------------------\n";
    
	my $org = shift or die;
    my $sample = shift or die;
    my $control = shift or die;
    my $peakoutput = shift or die;

    $sample =~ /\_(\d+)\.BAM$/;
    my $sample_id = $1;
    $control =~ /\_(\d+)\.BAM$/;
    my $control_id = $1;
    
    my $sample_bedfile = $sample;
    $sample_bedfile =~ s/\.BAM$/\.bam\.bed/;
    my $control_bedfile = $control;
    $control_bedfile =~ s/\.BAM$/\.bam\.bed/;
    
    SigSeeker::Convert_BAM2BED($sample,$sample_bedfile);
    SigSeeker::Convert_BAM2BED($control,$control_bedfile);

    print "Peak Calling - ";
    if(! -e $peakoutput || -z $peakoutput)
    {
        my $fname = $sample_bedfile.'_SWEMBLE';
        $fname =~ s/\.bam\.bed/\_$control_id/;
#	print $fname."\n";
#	exit;
	#BENCHMARK_6_13_SWEMBLE
        system($TOOL_REPOSITORY.'SWEMBL -i '.$sample_bedfile.' -B -o '.$fname.' -r '.$control_bedfile);
        print "done\n";
        
        print "Output Conversion\n";

        open(OUT,">".$peakoutput) or die "Cannot open $peakoutput";
        open(IN, $fname);
        while(my $record = <IN>)
        {
            chomp($record);
            if($record =~ /^chr/)
            {
                my @field = split(/\t/,$record);
                print OUT $field[0]."\t".$field[1]."\t".$field[2]."\t".$field[6]."\n";
            }
        }
        close OUT;
        close IN;
        #system('rm '.$sample_bedfile.'_TEMP');
    }
    else
    {
        print "skipped\n";    
    }
    return $peakoutput;
}

sub Erange
{
    print "---------------------------------------------------\n";
    print "Erange\n";
    print "---------------------------------------------------\n";
    
	my $org = shift or die;
    my $sample = shift or die;
    my $control = shift or die;
    my $peakoutput = shift or die;
    
    $sample =~ /\_(\d+)\.BAM$/;
    my $sample_id = $1;
    $control =~ /\_(\d+)\.BAM$/;
    my $control_id = $1;

    my $sample_bedfile = $sample;
    $sample_bedfile =~ s/\.BAM$/\.bam\.bed/;
    my $control_bedfile = $control;
    $control_bedfile =~ s/\.BAM$/\.bam\.bed/;
    
    SigSeeker::Convert_BAM2BED($sample,$sample_bedfile);
    SigSeeker::Convert_BAM2BED($control,$control_bedfile);
    
    my $sample_rds = $sample_bedfile;
    $sample_rds =~ s/\.bam\.bed$/\.bam\.rds/;
    my $control_rds = $control_bedfile;
    $control_rds =~ s/\.bam\.bed$/\.bam\.rds/;
    
    my $project = time();
    
    print "Conversion of BED to RDS - ";
    if(! -e $sample_rds || -z $sample_rds)
    {
        system('python '.$TOOL_REPOSITORY.'ERANGE3.3/makerdsfrombed.py '.$project.' '.$sample_bedfile.' '.$sample_rds);
        print "done\n";
    }
    else
    {
        print "skipped\n";    
    }
    
    print "Conversion of BED to RDS - ";
    if(! -e $control_rds || -z $control_rds)
    {
        system('python '.$TOOL_REPOSITORY.'ERANGE3.3/makerdsfrombed.py '.$project.' '.$control_bedfile.' '.$control_rds.' -index');
        print "done\n";
    }
    else
    {
        print "skipped\n";    
    }
    
    print "Peak Calling - ";
    if(! -e $peakoutput || -z $peakoutput)
    {
        my $fname = $sample_rds.'.regions';
        $fname =~ s/\.bam\.rds/\_$control_id\.bam\.rds/;
        #print('python '.$TOOL_REPOSITORY.'ERANGE3.3/findall.py '.$project.' '.$sample_rds.' '.$fname.' -control '.$control_rds.' -listPeak')."\n";
        system('python '.$TOOL_REPOSITORY.'ERANGE3.3/findall.py '.$project.' '.$sample_rds.' '.$fname.' -control '.$control_rds.' -listPeak');
        
        print "Output Conversion from $fname to $peakoutput\n";
        open(OUT,">".$peakoutput);# or die "Cannot open $peakoutput";
        open(IN, $fname);# or die "Cannot open $fname";
        while(my $record = <IN>)
        {
            chomp($record);
            if($record !~ /^\#/)
            {
                my @field = split(/\t/,$record);
                #print $field[1]."\t".$field[2]."\t".$field[3]."\t".$field[11]."\n";
                print OUT $field[1]."\t".$field[2]."\t".$field[3]."\t".$field[11]."\n";
            }
        }
        close OUT;
        close IN;
        #system('rm '.$sample_rds.'.regions');
        system('rm '.$sample_rds);
        system('rm '.$control_rds);

        print "done\n";
    }
    else
    {
        print "skipped\n";    
    }
    return $peakoutput;
}

sub SoleSearch
{
    print "---------------------------------------------------\n";
    print "Sole-Search\n";
    print "---------------------------------------------------\n";

	my $org = shift or die;
    my $sample = shift or die;
    my $control = shift or die;
    my $peakout = shift or die;

    my $sample_bedfile = $sample;
    $sample_bedfile =~ s/\.BAM$/\.bam\.bed/;
    my $control_bedfile = $control;
    $control_bedfile =~ s/\.BAM$/\.bam\.bed/;

    my $project = time();
    my $sample_export = $project.'_'.$sample;
    my $control_export = $project.'_'.$control;
    
    chdir($SigSeeker::TOOL_REPOSITORY.'sole-search');
    
    SigSeeker::Convert_BAM2BED($sample, $sample_bedfile);
    #system('sh parse_bed.sh '.$sample_bedfile.' '.$sample_export);
    
    SigSeeker::Convert_BAM2BED($control, $control_bedfile);
    #system('sh parse_bed.sh '.$control_bedfile.' '.$control_export);
    
    #print "Normalizing Control\n";
    #system('perl normalize_sd.pl -a 0.001 '.$control_export.'_full-length.txt');

    print "Peak Calling - ";
    if(! -e $peakoutput || -z $peakoutput)
    {
        #system('perl Sole-searchV2.pl -t '.$control_export.'.export.txt_tags.txt -c '.$control_export.'.export.txt_full-length.txt_corrected_0.001.sgr -p '.$control_export.'.export.txt_full-length.txt_duplications.gff '.$sample_export.'.export.txt_full-length.txt 2> /dev/null > /dev/null');
        #system('perl convert_gff2bed.pl '.$sample_export.'.export.txt_full-length.txt_signifpeaks.gff > '.$sample_export.'_unsorted_peaks.bed');
        chdir($SigSeeker::READ_UPLOADS);
        
        #print "Output Conversion\n";
        #system('sort -k1,1 -k2,2g -o '.$peakout.' '.$sample_export.'_unsorted_peaks.bed');
        #system('rm '.$sample_export.'_unsorted_peaks.bed');
        
        print "done\n";
    }
    else
    {
        chdir($SigSeeker::READ_UPLOADS);
        print "skipped\n";    
    }
    return $peakoutput;
}

sub CCAT
{
	my $organism = shift or die;    
    	my $sample = shift or die;
    	my $control = shift or die;
    	my $peakoutput = shift or die;
	my $tool = shift or die;
    
    	print "---------------------------------------------------\n";
    	print $tool."\n";
    	print "---------------------------------------------------\n";

    	my $sample_bedfile = $sample;
    	$sample_bedfile =~ s/\.BAM$/\.bam\.bed/;
    	my $control_bedfile = $control;
    	$control_bedfile =~ s/\.BAM$/\.bam\.bed/;

    	$sample =~ /\_(\d+)\.BAM$/;
    	my $sample_id = $1;
    	$control =~ /\_(\d+)\.BAM$/;
    	my $control_id = $1;
    
    	SigSeeker::Convert_BAM2BED($sample,$sample_bedfile);
    	SigSeeker::Convert_BAM2BED($control,$control_bedfile);
    
    	print "Peak Calling - ";
    	if(! -e $peakoutput || -z $peakoutput)
    	{
#  		system('CCAT '.$sample_bedfile.' '.$control_bedfile.' '.$SigSeeker::READ_UPLOADS.$organism.'.chrom.sizes /var/www/follow/Tool_Peak_Calling_Ensemble/Tools/CCAT3.0/example/config_'.$tool.'.txt '.$tool.'_'.$sample_id.'_'.$control_id);

        	open(CCAT, $tool.'_'.$sample_id.'_'.$control_id.'.significant.peak') or die "Cannot open ".$tool."_".$sample_id."_".$control_id.".significant.peak";
        	open(CCAT_OUT, ">".$tool.'_'.$sample_id.'_'.$control_id.'_peaksunsorted.bed') or die "Cannot open".$tool.'_'.$sample_id.'_'.$control_id.'_peaksunsorted.bed' ;
        	while(my $rec_ccat = <CCAT>)
        	{
            		chomp($rec_ccat);
            		my @ccat = split(/\t/,$rec_ccat);
            		print CCAT_OUT $ccat[0]."\t";   #Chromosome
            		print CCAT_OUT $ccat[2]."\t";   #Region Start
            		print CCAT_OUT $ccat[3]."\t";   #Region End
			print CCAT_OUT $ccat[7];        #Fold change score
	        	print CCAT_OUT "\n";
        	}
        	close CCAT;
        	close CCAT_OUT;
        
        	system('sort -k1,1 -k2,2g -o '.$peakoutput.' '.$tool.'_'.$sample_id.'_'.$control_id.'_peaksunsorted.bed');
		system('rm '.$tool.'_'.$sample_id.'_'.$control_id.'_peaksunsorted.bed');
		system('rm '.$tool.'_'.$sample_id.'_'.$control_id.'_significant.peak');

		print "done\n";
    	}
	else
	{
		print "skipped\n";
    	}

	return $peakoutput;
}

#Cleaning
sub CleanPeaks
{
    my $peakcalling_results = shift or die;
    my $peak_report = shift or die;
    
    print "Cleaning Peaks\n";
    
    open(PEAKS, ">".$peak_report) or die "Cannot open $peak_report\n";
    foreach my $file (keys %$peakcalling_results)
    {
        foreach my $sample (keys %{$peakcalling_results->{$file}})
        {
            foreach my $control (keys %{$peakcalling_results->{$file}->{$sample}})
            {
                #print $file."\t";
                #print $sample."\t";
                #print $control."\t";
                my @path = split(/\//, $peakcalling_results->{$file}->{$sample}->{$control});
                my $fname = $path[-1];
                my @field = split(/\_/,$fname);
                my @breakdown = split(/\-/,$field[-5]);
                
                open(RIN, $peakcalling_results->{$file}->{$sample}->{$control});
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
                system ('mv CLEANTMP '.$peakcalling_results->{$file}->{$sample}->{$control});
                
                my $cmd = 'wc -l '.$peakcalling_results->{$file}->{$sample}->{$control};
                my $count = `$cmd`;
                $count =~ /(\d+)/;
                my $peaks = $1;
                #print $peaks."\n";
                if($max_peaks < $peaks){$max_peaks = $peaks;}
                #if($peaks >= $encode_threshold)
                #{
                    #$setup_opt->{$breakdown[0]}->{$breakdown[1]}->{$sample}->{$control}->{"MACS2"} = $peaks;
                #}
                #else
                #{
                    #$encode_filtered->{$sample}->{$control} = 1;
                #}
                #$setup->{$breakdown[0]}->{$breakdown[1]}->{$sample}->{$control}->{"MACS2"} = $peaks;
                
                print PEAKS $breakdown[0]."\t".$breakdown[1]."\t".$sample."\t".$control."\t".$field[-2]."\t".$peaks."\n";
            }
        }
    }
    close PEAKS;
}

#Conversions

sub Convert_BAM2BED
{
    my $bamfile = shift or die;
    my $bedfile = shift or die;
    
    print "Conversion of BAM to BED - ";  
    
    if(! -e $bedfile || -z $bedfile)
	{
		print "done\n";
        my $cmd = 'bamToBed -i '.$bamfile.' > '.$bedfile;
		system($cmd);
	}
    else
    {
    	print "skipped\n";    
    }
}

sub Convert_CCAT2BED
{
    my $input = shift or die;
    my $tool = shift or die;
    
    open(CCAT, $input);
    open(CCAT_OUT, ">".$SigSeeker::PEAK_REPOSITORY.$project.'/'.$entry[0].'_'.$sample.'_'.$control.'_'.$tool.'_peaksunsorted.bed');
    while(my $rec_ccat = <CCAT>)
    {
        chomp($rec_ccat);
        my @ccat = split(/\t/,$rec_ccat);
        print CCAT_OUT $ccat[0]."\t";   #Chromosome
        print CCAT_OUT $ccat[2]."\t";   #Region Start
        print CCAT_OUT $ccat[3]."\t";   #Region End
        print CCAT_OUT $ccat[7];        #Fold change score
        print CCAT_OUT "\n";
    }
    close CCAT;
    close CCAT_OUT;
    
    system('sort -k1,1 -k2,2g -o '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$entry[0].'_'.$sample.'_'.$control.'_'.$tool.'_peaks.bed '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$entry[0].'_'.$sample.'_'.$control.'_'.$tool.'_peaksunsorted.bed');
    system('rm '.$SigSeeker::PEAK_REPOSITORY.$project.'/'.$entry[0].'_'.$sample.'_'.$control.'_'.$tool.'_peaksunsorted.bed');
}

#Interface Controls

sub PeakPlus
{
	my $project = shift or die;
    my $current = shift or die;
    
	if($current eq 'Quality Control'){print '<a href="'.$ENV{SERVER}.'/cgi-bin/SigSeeker2/SigSeeker_Project_Overview_QualityControl.pl?project='.$project.'">'};
    print 'Quality Control';
    if($current eq 'Quality Control'){print '</a>';}
	print '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;';
    
	if($current eq 'Overview'){print '<a href="'.$ENV{SERVER}.'/cgi-bin/SigSeeker2/SigSeeker_Project_Overview_Peaks.pl?project='.$project.'">';}
    print 'Overview';
    if($current eq 'Overview'){print '</a>';}
	print '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;';
    
	if($current eq 'Partitions'){print '<a href="'.$ENV{SERVER}.'/cgi-bin/SigSeeker2/SigSeeker_Partitioning.pl?project='.$project.'">';}
    print 'Partitions';
    if($current eq 'Partitions'){print '</a>';}
	print '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;';
    
	if($current eq 'Partition Visualization'){print '<a href="'.$ENV{SERVER}.'/cgi-bin/SigSeeker2/SigSeeker_Partitioning_Visualization.pl?project='.$project.'">';}
    print 'Partition Visualization';
    if($current eq 'Partition Visualization'){print '</a>';}
	print '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;';
    
	if($current eq 'Comparisons'){print '<a href="'.$ENV{SERVER}.'/cgi-bin/SigSeeker2/SigSeeker_Project_Overview_PeakComparisons.pl?project='.$project.'">';}
    print 'Comparisons';
    if($current eq 'Comparisons'){print '</a>';}
	print '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;';
    
	if($current eq 'Background Models'){print '<a href="'.$ENV{SERVER}.'/cgi-bin/SigSeeker2/SigSeeker_Project_Overview_PeakModels.pl?project='.$project.'">';}
    print 'Background Models';
    if($current eq 'Background Models'){print '</a>';}
	print '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;';
    
	if($current eq 'Filtering'){print '<a href="'.$ENV{SERVER}.'/cgi-bin/SigSeeker2/SigSeeker_Project_Overview_PeakConservation.pl?project='.$project.'">';}
    print 'Filtering';
    if($current eq 'Filtering'){print '</a>';}
	print '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;';
    
	print '<hr>';
}

sub FilterPlus
{
	my $project = shift or die;

	print '<a href="'.$ENV{SERVER}.'/cgi-bin/SigSeeker2/SigSeeker_Project_Overview_PeakConservation.pl?project='.$project.'&filter=conserved">Conservation</a>';
	print '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;';
	print '<a href="'.$ENV{SERVER}.'/cgi-bin/SigSeeker2/SigSeeker_Project_Overview_PeakConservation.pl?project='.$project.'&filter=cpgisland">CpG Islands</a>';
	print '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;';
	print '<hr>';
}

sub subheader
{
    my $project = shift or die;
    my $current = shift or die;
    
#	if($current ne 'Quality Control'){print '<a href="'.$ENV{SERVER}.'/sigseeker-cgi/SigSeeker_Project_Overview_QualityControl.pl?project='.$project.'">'};
#   print 'Quality Control';
#   if($current ne 'Quality Control'){print '</a>';}
#	print '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;';
    
	if($current ne 'Overview'){print '<a href="'.$ENV{SERVER}.'/sigseeker-cgi/SigSeeker_Project_Overview_New.pl?project='.$project.'">';}
    print 'Overview';
    if($current ne 'Overview'){print '</a>';}
	print '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;';
    
	if($current ne 'Comparisons'){print '<a href="'.$ENV{SERVER}.'/sigseeker-cgi/SigSeeker_Project_Overview_NewComparisons.pl?project='.$project.'">';}
    print 'Comparisons';
    if($current ne 'Comparisons'){print '</a>';}
	print '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;';

    if($current ne 'Correlations'){print '<a href="'.$ENV{SERVER}.'/sigseeker-cgi/SigSeeker_Project_Overview_NewCorrelations.pl?project='.$project.'">';}
    print 'Correlations';
    if($current ne 'Correlations'){print '</a>';}
	print '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;';
    
    if($current ne 'Partitions'){print '<a href="'.$ENV{SERVER}.'/sigseeker-cgi/SigSeeker_Project_Overview_NewPartitions.pl?project='.$project.'">';}
    print 'Partitions';
    if($current ne 'Partitions'){print '</a>';}
	print '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;';
    
    
    if($current ne 'Metrics'){print '<a href="'.$ENV{SERVER}.'/sigseeker-cgi/SigSeeker_Project_Overview_NewMetrics.pl?project='.$project.'">';}
    print 'Metrics';
    if($current ne 'Metrics'){print '</a>';}
    print '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;';
    
#	if($current ne 'Background Models'){print '<a href="'.$ENV{SERVER}.'/sigseeker-cgi/SigSeeker_Project_Overview_PeakModels.pl?project='.$project.'">';}
#    print 'Background Models';
#    if($current ne 'Background Models'){print '</a>';}
#	print '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;';
  
#    if($current ne 'Filtering'){print '<a href="'.$ENV{SERVER}.'/sigseeker-cgi/SigSeeker_Project_Overview_NewFiltering.pl?project='.$project.'">';}
#    print 'Filtering';
#    if($current ne 'Filtering'){print '</a>';}
#    print '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;';

    print '<br/><br/>';
    
    if($current eq 'Overview')
    {
        print '<a href="'.$ENV{SERVER}.'/sigseeker-cgi/SigSeeker_Project_Overview_New.pl?project='.$project.'">Data</a>';
    }
    else
    {
        print '<a href="'.$ENV{SERVER}.'/sigseeker-cgi/SigSeeker_Project_Overview_New'.$current.'.pl?project='.$project.'">Data</a>';
    }
    print '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;';
    if($current ne "Correlations" && $current ne "Metrics")
    {
        if($current eq 'Overview')
        {
            print '<a href="'.$ENV{SERVER}.'/sigseeker-cgi/SigSeeker_Project_Overview_New.pl?project='.$project.'&visualization=1">Visualization</a>';
        }
        else
        {
            print '<a href="'.$ENV{SERVER}.'/sigseeker-cgi/SigSeeker_Project_Overview_New'.$current.'.pl?project='.$project.'&visualization=1">Visualization</a>';
        }
        print '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;';
    }
    print '<hr>';
}

sub header
{
	my $title = shift;

	print "Content-type: text/html\n\n";
	print '<html>'."\n";
	print '<head>'."\n";
	print '<title>SigSeeker+</title>'."\n";
	print '<script src="'.$ENV->{SERVER}.'/sorttable.js"></script>'."\n";
    print '<script src="'.$ENV->{SERVER}.'/table.js"></script>'."\n";
	if($title eq "Status")
	{
		print '<meta http-equiv="refresh" content="60" >';
	}
	print '<link rel="stylesheet" type="text/css" href="'.$ENV->{SERVER}.'/SigSeeker.css">';
	print '</head>'."\n";
	print "<body>"."\n";
	print "<SCRIPT>"."\n";
	print 'function showStuff(show_id,hide_id) {'."\n";
	print 'document.getElementById(hide_id).style.display = \'block\';'."\n";
	print 'document.getElementById(show_id).style.display = \'none\';'."\n";
	print '}'."\n";

	print 'function hideStuff(show_id,hide_id)'."\n";
    print '{'."\n";
	print '     document.getElementById(hide_id).style.display = \'none\';'."\n";
	print '     document.getElementById(show_id).style.display = \'block\';'."\n";
	print '}'."\n";

	print 'function displayResult(target_type,select_type)'."\n";
	print '{'."\n";
	print '		var x=document.getElementById(target_type);'."\n";
	print '		var selectedArray = new Array();'."\n";
	print '		var selObj = document.getElementById(select_type);'."\n";
	print '		var i;'."\n";
	print '		var count = 0;'."\n";
	print '		for (i=0; i < selObj.options.length; i++)'."\n";
	print '		{'."\n";
	print '			if (selObj.options[i].selected)'."\n";
	print '			{'."\n";
	print '				selectedArray[count] = selObj.options[i].value;'."\n";
	print '				count++;'."\n";
	print '				var option=document.createElement("option");'."\n";
	print '				option.text=selObj.options[i].text;'."\n";
	print '				option.value=selObj.options[i].value;'."\n";
	print '				option.selected=true;'."\n";
	print '				try'."\n";
	print '				{'."\n";
	print '					x.add(option,x.options[null]);'."\n";
	print '				}'."\n";
	print '				catch (e)'."\n";
	print '				{'."\n";
	print '					x.add(option,null);'."\n";
	print '				}'."\n";
	print '			}'."\n";
	print '		}'."\n";
	print '}'."\n";
    
    print 'function formcheck()'."\n";
    print '{'."\n";
    print '     var elements = document.forms["Main"].elements;'."\n";
    print '     for (i=0; i < elements.length; i++)'."\n";
    print '     {'."\n";
    print '         var x=document.forms["Main"].elements[i].value;'."\n";
    print '         if (x==null || x=="")'."\n";
    print '         {'."\n";
    print '             if(document.forms["Main"].elements[i].className == "TBR")'."\n";
    print '             {'."\n";
    print '                 var name = document.forms["Main"].elements[i].name;'."\n";
    print '                 var colClass = document.forms["Main"].elements[i].className;'."\n";
    print '                 alert(name + " " + colClass + " must be filled out");'."\n";
    print '                 return false;'."\n";
    print '             }'."\n";
    print '         }'."\n";
    print '     }'."\n";
    
#    print '     var x=document.forms["Main"]["projectname"].value;'."\n";
#    print '     if (x==null || x=="")'."\n";
#    print '     {'."\n";
#    print '         alert("First name must be filled out");'."\n";
#    print '         return false;'."\n";
#    print '     }'."\n";
    print '}'."\n";

	print '</SCRIPT>'."\n";

	print "<FORM name=\"jump1\">";
	print '<table style="width:100%">';
	print "<tr>";
	print '<td style="width:75px">';
	print '<a href="./SigSeeker.pl"><img src="http://'.$ENV{'SERVER_NAME'}.'/Images/SigSeeker.png" width=75 alt="Cannot find image"/></a>';
	print "</td><td align=\"left\" valign=\"middle\">&nbsp;Navigation:<BR>";
	print "<select name=\"myjumpbox\" OnChange=\"location.href=jump1.myjumpbox.options[selectedIndex].value\">";
	print "<option selected>Please Select...";
	print "<option value=\"./SigSeeker.pl\">Home";
    
#    print '<optgroup label="Quality Control">';
#    print "<option value=\"./SigSeeker_QualityControl.pl\">Quality Control";
#    print '</optgroup>';

    print '<optgroup label="Basic">';
	print "<option value=\"./SigSeeker_QualityControl.pl\">Quality Control";
#	print "<option value=\"./SigSeeker_ReadMapping.pl\">Read Mapping";
    print '</optgroup>';

    print '<optgroup label="Analysis">';
	print "<option value=\"./SigSeeker_RNAAnalysis.pl\">RNA Analysis";
	print "<option value=\"./SigSeeker_PeakCalling.pl\">Peak Calling";
    print '</optgroup>';

    print '<optgroup label="Overview">';
	print "<option value=\"./SigSeeker_Project_Overview.pl\">Projects";
    print '</optgroup>';
    
    print '<optgroup label="User Management">';
    print "<option value=\"./SigSeeker_User.pl\">User Generation";
    print '</optgroup>';

    print "</select>";
	print "</td>";
	print "</tr>";
	print '</table>';
	print "</FORM>";

	print '<table style="border-spacing:0;width:100%">';
	print '<tr align="center">';
	print '<td style="width:25px;">';
	print '</td>';
	print '<td>';
	print '<p style="font-weight:bold"><a href="./SigSeeker_Examples.pl">Examples</a></p>';
	print '</td>';
	print '<td>';
	print '<p style="font-weight:bold"><a href="./SigSeeker_HELP.pl">Documentation</a></p>';
	print '</td>';
	print '<td>';
	print '<p style="font-weight:bold"><a href="http://code.google.com/p/sigseeker">Downloads</a></p>';
	print '</td>';
	print '</tr>';
	print "</table>";
	print "<hr>";

	#Courtesy of SimplytheBest.net - http://simplythebest.net/scripts/
	print '<div id="overDiv" style="position:absolute; visibility:hide; z-index:1;">';
	print '</div>';
	print '<script LANGUAGE="JavaScript" SRC="../../overlib.js"></script>';
	
	if($title)
	{
		print "<h1>".$title."</h1>";
		print '<p style="width: 800px;">'.$SigSeeker_HELP::helphash->{$title}->{'Documentation'}.'<p>';
	}
}

sub footer
{
	print "<hr>";
	print "Copyright <a href=\"http://msseeker.org\">Jens Lichtenberg</a>";
	print "</body>";

	return;
}

1;
