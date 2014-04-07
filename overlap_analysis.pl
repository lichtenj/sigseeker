use strict;
use POSIX qw(ceil);

my $file = shift or die;
my $refFlat = shift or die;
my $output = shift or die;

#print "In\n";

#$file =~ s/\|/\\\|/g;
#print $file."\n\n";

#$output =~ s/\|/\\\|/g;
#print $output."\n\n";

my $coordinates;
my $genes;

my $count = 0;
open(IN, $refFlat) or die "Cannot open refFlat at $refFlat";
while(my $record = <IN>)
{
	$count++;
	chomp($record);
	my @field = split(/\t/,$record);
	$coordinates->{$count}->{SYMBOL} = $field[0];
	$coordinates->{$count}->{ID} = $field[1];
	$coordinates->{$count}->{CHR} = $field[2];
	$coordinates->{$count}->{STR} = $field[3];
	$coordinates->{$count}->{TSS} = $field[4];
	$coordinates->{$count}->{TES} = $field[5];

	if($field[4] > $field[5])
	{
		print $record."\n";
	}
}
close IN;
if(! -e $output.'_RefSeq.bed' || -z $output.'_RefSeq.bed' || ! -e $output.'_GenericRefSeq.bed' || -z $output.'_GenericRefSeq.bed')
{
    open(OUTG, ">".$output."_GenericIntergenic.bed") or die "Cannot open ".$output."_GenericIntergenic.bed file $!";
    open(OUTGU, ">".$output."_GenericUpstream.bed") or die "Cannot open ".$output."_GenericUpstream.bed file $!";
    open(OUTGP, ">".$output."_GenericPromoter.bed") or die "Cannot open ".$output."_GenericPromoter.bed file $!";
    open(OUTGR, ">".$output."_GenericRefSeq.bed") or die "Cannot open ".$output."_GenericRefSeq.bed file $!";
    open(OUTGD, ">".$output."_GenericDownstream.bed") or die "Cannot open ".$output."_GenericDownstream.bed file $!";
    
    open(OUT, ">".$output."_Intergenic.bed") or die "Cannot open ".$output."_Intergenic.bed file $!";
    open(OUTU, ">".$output."_Upstream.bed") or die "Cannot open ".$output."_Upstream.bed file $!";
    open(OUTP, ">".$output."_Promoter.bed") or die "Cannot open ".$output."_Promoter.bed file $!";
    open(OUTR, ">".$output."_RefSeq.bed") or die "Cannot open ".$output."_RefSeq.bed file $!";
    open(OUTD, ">".$output."_Downstream.bed") or die "Cannot open ".$output."_Downstream.bed file $!";
    
    open(IN, $file) or die "Cannot open input file $file";
    while(my $record = <IN>)
    {
    	chomp($record);
    	my @field = split(/\t/,$record);
    #	if($field[6] >= $min)
    #	{
            my $bool = 0;
    		my $max_left;
    		my $gene_left;
    		my $partition_left;
    
    		my $max_right;
    		my $gene_right;
    		my $partition_right;
    
    		my @starts;
    		my @ends;
    
    		my @partitions;
    			
    		foreach my $id (keys %$coordinates)
    		{
    			if($field[0] ne $coordinates->{$id}->{CHR})
    			{
    				next;
    			}
    
    			#Downstream
    			if($coordinates->{$id}->{STR} eq "+" && ($field[1] >= $coordinates->{$id}->{TES} && $field[1] <= $coordinates->{$id}->{TES} + 10000))
    			{
                    print OUTGD $record."\n";
                    $bool = 1;
                    
    				push @starts,$coordinates->{$id}->{TES};
    				push @ends,$coordinates->{$id}->{TES} + 10000;
    
    				push @partitions, $coordinates->{$id}->{SYMBOL}." Downstream +";
    #				print "Downstream +\t";
    				if(!$max_left || $max_left > $coordinates->{$id}->{TES})
    				{
    					$max_left = $coordinates->{$id}->{TES};
    					$gene_left = $coordinates->{$id}->{SYMBOL};
    					$partition_left = "Downstream";
    				}
    				if(!$max_right || $max_right < $coordinates->{$id}->{TES} + 10000)
    				{
    					$max_right = $coordinates->{$id}->{TES} + 10000;
    					$gene_right = $coordinates->{$id}->{SYMBOL};
    					$partition_right = "Downstream";
    				}				
    			}
    			if($coordinates->{$id}->{STR} eq "-" && ($field[1] <= $coordinates->{$id}->{TSS} && $field[1] >= $coordinates->{$id}->{TSS} - 10000))
    			{
                    print OUTGD $record."\n";
                    $bool = 1;
                    
    				push @starts,$coordinates->{$id}->{TSS} - 10000;
    				push @ends,$coordinates->{$id}->{TSS};
    
    				push @partitions, $coordinates->{$id}->{SYMBOL}." Downstream -";
    #				$partitions++;
    #				print "Downstream -\t";
    				if(!$max_left || $max_left > $coordinates->{$id}->{TSS} - 10000)
    				{
    					$max_left = $coordinates->{$id}->{TSS} - 10000;
    					$gene_left = $coordinates->{$id}->{SYMBOL};
    					$partition_left = "Downstream";
    				}
    				if(!$max_right || $max_right < $coordinates->{$id}->{TSS})
    				{
    					$max_right = $coordinates->{$id}->{TSS};
    					$gene_right = $coordinates->{$id}->{SYMBOL};
    					$partition_right = "Downstream";
    				}				
    			}
    
    			#Upstream
    			if($coordinates->{$id}->{STR} eq "+" && ($field[1] <= $coordinates->{$id}->{TSS} - 1000 && $field[1] >= $coordinates->{$id}->{TSS} - 10000))
    			{
                    print OUTGU $record."\n";
                    $bool = 1;
                    
    				push @starts,$coordinates->{$id}->{TSS} - 10000;
    				push @ends,$coordinates->{$id}->{TSS} - 1000;
    
    				push @partitions, $coordinates->{$id}->{SYMBOL}." Upstream +";
    #				$partitions++;
    #				print "Upstream +\t";
    				if(!$max_left || $max_left > $coordinates->{$id}->{TSS} - 1000)
    				{
    					$max_left = $coordinates->{$id}->{TSS} - 1000;
    					$gene_left = $coordinates->{$id}->{SYMBOL};
    					$partition_left = "Upstream";
    				}
    				if(!$max_right || $max_right < $coordinates->{$id}->{TSS} - 10000)
    				{
    					$max_right = $coordinates->{$id}->{TSS} - 10000;
    					$gene_right = $coordinates->{$id}->{SYMBOL};
    					$partition_right = "Upstream";
    				}				
    			}
    			if($coordinates->{$id}->{STR} eq "-" && ($field[1] >= $coordinates->{$id}->{TES} + 1000 && $field[1] <= $coordinates->{$id}->{TES} + 10000))
    			{
                    print OUTGU $record."\n";
                    $bool = 1;
                    
    				push @starts,$coordinates->{$id}->{TES} + 1000;
    				push @ends,$coordinates->{$id}->{TES} + 10000;
    
    				push @partitions, $coordinates->{$id}->{SYMBOL}." Upstream -";
    #				$partitions++;
    #				print "Upstream -\t";
    				if(!$max_left || $max_left > $coordinates->{$id}->{TES} + 1000)
    				{
    					$max_left = $coordinates->{$id}->{TES} + 1000;
    					$gene_left = $coordinates->{$id}->{SYMBOL};
    					$partition_left = "Upstream";
    				}
    				if(!$max_right || $max_right < $coordinates->{$id}->{TES} + 10000)
    				{
    					$max_right = $coordinates->{$id}->{TES} + 10000;
    					$gene_right = $coordinates->{$id}->{SYMBOL};
    					$partition_right = "Upstream";
    				}				
    			}
    
    			#Promoter
    			if($coordinates->{$id}->{STR} eq "+" && ($field[1] <= $coordinates->{$id}->{TSS} + 50 && $field[1] >= $coordinates->{$id}->{TSS} - 1000))
    			{
                    print OUTGP $record."\n";
                    $bool = 1;
                    
    				push @starts,$coordinates->{$id}->{TSS} - 1000;
    				push @ends,$coordinates->{$id}->{TSS} + 50;
    
    				push @partitions, $coordinates->{$id}->{SYMBOL}." Promoter +";
    #				$partitions++;
    #				print "Promoter +\t";
    				if(!$max_left || $max_left > $coordinates->{$id}->{TSS} + 50)
    				{
    					$max_left = $coordinates->{$id}->{TSS} + 50;
    					$gene_left = $coordinates->{$id}->{SYMBOL};
    					$partition_left = "Promoter";
    				}
    				if(!$max_right || $max_right < $coordinates->{$id}->{TSS} - 1000)
    				{
    					$max_right = $coordinates->{$id}->{TSS} - 1000;
    					$gene_right = $coordinates->{$id}->{SYMBOL};
    					$partition_right = "Promoter";
    				}				
    			}
    			if($coordinates->{$id}->{STR} eq "-" && ($field[1] >= $coordinates->{$id}->{TES} - 50 && $field[1] <= $coordinates->{$id}->{TES} + 1000))
    			{
                    print OUTGP $record."\n";
                    $bool = 1;
                    
    				push @starts,$coordinates->{$id}->{TES} - 50;
    				push @ends,$coordinates->{$id}->{TES} + 1000;
    
    				push @partitions, $coordinates->{$id}->{SYMBOL}." Promoter -";
    #				$partitions++;
    #				print "Promoter -\t";
    				if(!$max_left || $max_left > $coordinates->{$id}->{TES} - 50)
    				{
    					$max_left = $coordinates->{$id}->{TES} - 50;
    					$gene_left = $coordinates->{$id}->{SYMBOL};
    					$partition_left = "Promoter";
    				}
    				if(!$max_right || $max_right < $coordinates->{$id}->{TES} + 1000)
    				{
    					$max_right = $coordinates->{$id}->{TES} + 1000;
    					$gene_right = $coordinates->{$id}->{SYMBOL};
    					$partition_right = "Promoter";
    				}				
    			}
    
    			#RefSeq
    			if($coordinates->{$id}->{STR} eq "+" && ($field[1] >= $coordinates->{$id}->{TSS} + 51 && $field[1] <= $coordinates->{$id}->{TES}))
    			{
                    print OUTGR $record."\n";
                    $bool = 1;
                    
    				push @starts,$coordinates->{$id}->{TSS} + 51;
    				push @ends,$coordinates->{$id}->{TES};
    
    				push @partitions, $coordinates->{$id}->{SYMBOL}." RefSeq +";
    #				$partitions++;
    #				print "RefSeq +\t";
    				if(!$max_left || $max_left > $coordinates->{$id}->{TSS} + 51)
    				{
    					$max_left = $coordinates->{$id}->{TSS} + 51;
    					$gene_left = $coordinates->{$id}->{SYMBOL};
    					$partition_left = "RefSeq";
    				}
    				if(!$max_right || $max_right < $coordinates->{$id}->{TES})
    				{
    					$max_right = $coordinates->{$id}->{TES};
    					$gene_right = $coordinates->{$id}->{SYMBOL};
    					$partition_right = "RefSeq";
    				}				
    			}
    			if($coordinates->{$id}->{STR} eq "-" && ($field[1] <= $coordinates->{$id}->{TES} - 51 && $field[1] >= $coordinates->{$id}->{TSS}))
    			{
                    print OUTGR $record."\n";
                    $bool = 1;
                    
    				push @starts,$coordinates->{$id}->{TSS};
    				push @ends,$coordinates->{$id}->{TES} - 51;
    
    				push @partitions, $coordinates->{$id}->{SYMBOL}." RefSeq -";
    #				$partitions++;
    #				print "RefSeq -\t";
    				if(!$max_left || $max_left > $coordinates->{$id}->{TSS})
    				{
    					$max_left = $coordinates->{$id}->{TSS};
    					$gene_left = $coordinates->{$id}->{SYMBOL};
    					$partition_left = "RefSeq";
    				}
    				if(!$max_right || $max_right < $coordinates->{$id}->{TES} - 51)
    				{
    					$max_right = $coordinates->{$id}->{TES} - 51;
    					$gene_right = $coordinates->{$id}->{SYMBOL};
    					$partition_right = "RefSeq";
    				}				
    			}
    		}

            if($bool == 0)
            {
                print OUTG $record."\n";
            }

    		my @ts = sort @starts;
    		my @te = sort @ends;
    
    		my $middle = ceil( abs($ts[0] + $te[0]) / 2);
    #		my $middle = ceil( abs($field[4] + $field[5]) / 2);
    		#print $middle."\t";
    		if($field[1] < $middle && scalar(@partitions) > 0)
    		{
                if($partition_left eq "Upstream")
                {
        			print OUTU $record."\t";
        			print OUTU $partition_left."\t".$gene_left."\t";
        			print OUTU scalar(@partitions)."\t";
        			print OUTU $middle."\t";
        			foreach my $partition (@partitions)
        			{
					$partition =~ /(\w+)/;
					$gl->{UPSTREAM}->{$1} = 1;
        				print OUTU $partition.",";
        			}
        			print OUTU "\n";
                }
                elsif($partition_left eq "Promoter")
                {
            		print OUTP $record."\t";
        			print OUTP $partition_left."\t".$gene_left."\t";
        			print OUTP scalar(@partitions)."\t";
        			print OUTP $middle."\t";
        			foreach my $partition (@partitions)
        			{
					$partition =~ /(\w+)/;
					$gl->{PROMOTER}->{$1} = 1;
        				print OUTP $partition.",";
        			}
        			print OUTP "\n";
                }
                elsif($partition_left eq "RefSeq")
                {
            		print OUTR $record."\t";
        			print OUTR $partition_left."\t".$gene_left."\t";
        			print OUTR scalar(@partitions)."\t";
        			print OUTR $middle."\t";
        			foreach my $partition (@partitions)
        			{
					$partition =~ /(\w+)/;
					$gl->{REFSEQ}->{$1} = 1;
        				print OUTR $partition.",";
        			}
        			print OUTR "\n";
                }
                elsif($partition_left eq "Downstream")
                {
            		print OUTD $record."\t";
        			print OUTD $partition_left."\t".$gene_left."\t";
        			print OUTD scalar(@partitions)."\t";
        			print OUTD $middle."\t";
        			foreach my $partition (@partitions)
        			{
					$partition =~ /(\w+)/;
					$gl->{DOWNSTREAM}->{$1} = 1;
        				print OUTD $partition.",";
        			}
        			print OUTD "\n";
                }
                else
        	    {
    			    print OUT $record."\n";
    		    }
            }
    		elsif($field[1] >= $middle && scalar(@partitions) > 0)
    		{
    			if($partition_right eq "Upstream")
                {
            		print OUTU $record."\t";
        			print OUTU $partition_right."\t".$gene_right."\t";
        			print OUTU scalar(@partitions)."\t";
        			print OUTU $middle."\t";
        			foreach my $partition (@partitions)
        			{
					$partition =~ /(\w+)/;
					$gl->{UPSTREAM}->{$1} = 1;
        				print OUTU $partition.",";
        			}
        			print OUTU "\n";
                }
                elsif($partition_right eq "Promoter")
                {
            		print OUTP $record."\t";
        			print OUTP $partition_right."\t".$gene_right."\t";
        			print OUTP scalar(@partitions)."\t";
        			print OUTP $middle."\t";
        			foreach my $partition (@partitions)
        			{
					$partition =~ /(\w+)/;
					$gl->{PROMOTER}->{$1} = 1;
        				print OUTP $partition.",";
        			}
        			print OUTP "\n";
                }
                elsif($partition_right eq "RefSeq")
                {
            		print OUTR $record."\t";
        			print OUTR $partition_right."\t".$gene_right."\t";
        			print OUTR scalar(@partitions)."\t";
        			print OUTR $middle."\t";
        			foreach my $partition (@partitions)
        			{
					$partition =~ /(\w+)/;
					$gl->{REFSEQ}->{$1} = 1;
        				print OUTR $partition.",";
        			}
        			print OUTR "\n";
                }
                elsif($partition_right eq "Downstream")
                {
            		print OUTD $record."\t";
        			print OUTD $partition_right."\t".$gene_right."\t";
        			print OUTD scalar(@partitions)."\t";
        			print OUTD $middle."\t";
        			foreach my $partition (@partitions)
        			{
					$partition =~ /(\w+)/;
					$gl->{DOWNSTREAM}->{$1} = 1;
        				print OUTD $partition.",";
        			}
        			print OUTD "\n";
                }
                else
        	    {
    			    print OUT $record."\n";
    		    }
    		}
    		else
    		{
    			print OUT $record."\n";
    		}
    #	}
    }
    close IN;
    close OUT;
    close OUTU;
    close OUTP;
    close OUTR;
    close OUTD;
    
    close OUTG;
    close OUTGU;
    close OUTGP;
    close OUTGR;
    close OUTGD;
}

open(OGLU, ">".$output."_Upstream_GL.txt") or die "Cannot open ".$output."_Upstream_GL.bed file $!";
open(OGLP, ">".$output."_Promoter_GL.txt") or die "Cannot open ".$output."_Promoter_GL.bed file $!";
open(OGLR, ">".$output."_RefSeq_GL.txt") or die "Cannot open ".$output."_RefSeq_GL.bed file $!";
open(OGLD, ">".$output."_Downstream_GL.txt") or die "Cannot open ".$output."_Downstream_GL.bed file $!";
foreach my $id (keys %$gl->{UPSTREAM})
{
	print OGLU $id."\n";
}
foreach my $id (keys %$gl->{PROMOTER})
{
	print OGLP $id."\n";
}
foreach my $id (keys %$gl->{REFSEQ})
{
	print OGLR $id."\n";
}
foreach my $id (keys %$gl->{DOWNSTREAM})
{
	print OGLD $id."\n";
}
close OGLU;
close OGLP;
close OGLR;
close OGLD;

#$output =~ s/\|/\\\|/g;
my $cmd = 'wc -l '.$output.'_Intergenic.bed';
#$cmd =~ s/\|/\\\|/g;
#print $cmd."\n";
my $count = `$cmd`;
$count =~ /(\d+)/;
print $1."\t";

$cmd = 'wc -l '.$output.'_Upstream.bed';
#$cmd =~ s/\|/\\\|/g;
$count = `$cmd`;
$count =~ /(\d+)/;
print $1."\t";

$cmd = 'wc -l '.$output.'_Promoter.bed';
#$cmd =~ s/\|/\\\|/g;
$count = `$cmd`;
$count =~ /(\d+)/;
print $1."\t";

$cmd = 'wc -l '.$output.'_RefSeq.bed';
#$cmd =~ s/\|/\\\|/g;
$count = `$cmd`;
$count =~ /(\d+)/;
print $1."\t";

$cmd = 'wc -l '.$output.'_Downstream.bed';
#$cmd =~ s/\|/\\\|/g;
$count = `$cmd`;
$count =~ /(\d+)/;
print $1."\t";

#$output =~ s/\|/\\\|/g;
my $cmd = 'wc -l '.$output.'_GenericIntergenic.bed';
#$cmd =~ s/\|/\\\|/g;
#print $cmd."\n";
my $count = `$cmd`;
$count =~ /(\d+)/;
print $1."\t";

$cmd = 'wc -l '.$output.'_GenericUpstream.bed';
#$cmd =~ s/\|/\\\|/g;
$count = `$cmd`;
$count =~ /(\d+)/;
print $1."\t";

$cmd = 'wc -l '.$output.'_GenericPromoter.bed';
#$cmd =~ s/\|/\\\|/g;
$count = `$cmd`;
$count =~ /(\d+)/;
print $1."\t";

$cmd = 'wc -l '.$output.'_GenericRefSeq.bed';
#$cmd =~ s/\|/\\\|/g;
$count = `$cmd`;
$count =~ /(\d+)/;
print $1."\t";

$cmd = 'wc -l '.$output.'_GenericDownstream.bed';
#$cmd =~ s/\|/\\\|/g;
$count = `$cmd`;
$count =~ /(\d+)/;
print $1."\n";

sub closest 
{ 
        my $find = shift; 
        my $closest = shift; 
        for my $num (@_) 
        { 
                if (abs($num - $find) < abs($closest - $find)) 
                { 
                        $closest = $num; 
                } 
        } 
        return $closest; 
}
