#!/usr/bin/perl

my $usage= "
Prints out list of commands for launcher_creator.py
in order to recover .fastq files from .bam files
The arguments are:
1: glob to fastq files
2: awk code for searching tabs and printing reads'
";

if (!$ARGV[0]) { die $usage;}
my $glob=$ARGV[0];
if (!$ARGV[1]) { die $usage;}
my $ref=$ARGV[1];


opendir THIS, ".";
my @fqs=grep /$glob/,readdir THIS;
my $outname="";

foreach $fqf (@fqs) {
        if ($ARGV[2]) {
                my @parts=split('_',$fqf);
                $outname=$parts[$ARGV[1]-1].".fastq";
        }
	else { $outname=$fqf.".fastq";}
        print "samtools view $fqf | awk $ref > $outname\n";
}
