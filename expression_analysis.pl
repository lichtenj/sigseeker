#!/usr/bin/perl

use strict;
use Math::Combinatorics;
use lib '/home/darklichti/Dropbox/SigSeeker_CODE_05012013';
use lib '/Users/lichtenbergj/SigSeeker_CODE_05012013';use SigSeeker;

my $partitioned_peaklist = shift or die;
my $projectname = shift or die;
my $partition = shift or die;

my $genes="";

open(IN, $partitioned_peaklist);
while(my $record = <IN>)
{
    chomp($record);
    my @fl = split(/\t/,$record);

    $genes .= $fl[-4].',';
}
close IN;

chop($genes);
#print $genes."\n";
my $bodyname = $partitioned_peaklist;
$bodyname =~ /([\.\-\_\w]+)$/;
$bodyname = $1;
open(OUT, '>'.$SigSeeker::SBR_REPOSITORY.'SigSeeker_Correlations/'.$projectname.'_expression_'.$bodyname.'_'.$partition.'.txt') or die "Cannot open $!";

print OUT '$examples->{Title} = "'.$projectname.': Expression Correlation ('.$partition.')";'."\n";
print OUT '$examples->{Description} = "This breakdown illustrates which of the genes, with a peak in their '.$partition.' are expressed for a selected set of cell types.";'."\n";
print OUT '$examples->{Type} = "GeneCentric";'."\n";
#print OUT '$examples->{celltypes}->{2} = 1;'."\n";
#print OUT '$examples->{celltypes}->{3} = 1;'."\n";
#print OUT '$examples->{celltypes}->{4} = 1;'."\n";
#print OUT '$examples->{celltypes}->{5} = 1;'."\n";
#print OUT '$examples->{celltypes}->{7} = 1;'."\n";
print OUT '$examples->{celltypes}->{8} = 1;'."\n";
#print OUT '$examples->{celltypes}->{10} = 1;'."\n";
print OUT '$examples->{celltypes}->{14} = 1;'."\n";
print OUT '$examples->{input_genes} = "'.$genes.'";'."\n";
print OUT '$examples->{analysis_threshold} = 0.05;'."\n";

close OUT;

print 'http://72.83.3.85/sbr-cgi/GENE_ANALYSIS.pl?example_file=./SigSeeker_Correlations/'.$projectname.'_expression_'.$bodyname.'_'.$partition.'.txt';
print "\n";