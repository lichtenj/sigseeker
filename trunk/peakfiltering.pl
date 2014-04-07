use strict;
use lib '/home/darklichti/Dropbox/SigSeeker_CODE_05012013';
use lib '/Users/lichtenbergj/SigSeeker_CODE_05012013';
use SigSeeker;

my $file = shift or die;
my $organism = shift or die;
my $id = shift or die;

#open(IN, $file);
#open(OUT, ">".$id.'_modified.bed') or die "Cannot write modified BED file".$id.'_modified.bed';
#while(my $record = <IN>)
#{
#	my @field = split(/\t/,$record);
#	print OUT $field[0]."\t";
#	print OUT $field[1]."\t";
#	print OUT $field[2]."\n";
#}
#close OUT;
#close IN;

opendir(DIR, $SIGSEEKER::FILTER_REPOSITORY);
while(my $filter = readdir(DIR))
{
    if($filter =~ /^\./){next;}
#	print "\t".$filter."\n";
#    my $filter_id = $filter;
#    $filter_id =~ /\/([\w\_\-\.]+)$/;
#    $filter_id = $1;
	my $output = $id.'_'.$filter.'.bed';
	if(! -e $output || -z $output)
	{
#		my $cmd = 'intersectBed -a '.$id.'_modified.bed -b '.$SIGSEEKER::FILTER_REPOSITORY.$filter.' > '.$output.' 2> /dev/null';
    	my $cmd = 'intersectBed -a '.$file.' -b '.$SIGSEEKER::FILTER_REPOSITORY.$filter.' > '.$output.' 2> /dev/null';
		system($cmd);
#        print $cmd."\n";
	}
    my $cmd = 'wc -l '.$output;
    my $count = `$cmd`;
    $count =~ /^(\d+)/;
    print $1.' ';
}
