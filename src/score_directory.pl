#!/usr/bin/perl

use strict;
use Cwd qw(realpath);
use File::Basename;

# Globals

my $SCORE_PREDICTION = dirname(realpath($0))."/score_prediction";
my @dir_names = ();

############################################################
# Version()
#
# Display program version.
############################################################
 
sub Version 
{
    print STDERR "ScoreDirectory version 1.00 - Compare two directories of RNA secondary structure predictions\n";
    print STDERR "\n";
    print STDERR "Written by Chuong B. Do\n";
    exit(0);
}

############################################################
# Usage()
#
# Print program usage.
############################################################

sub Usage
{
    print STDERR "\n";
    print STDERR "Usage: ".basename(realpath($0))." [protein|rna] [OPTION]... REF_DIR TEST_DIR\n";
    print STDERR "\n";
    print STDERR "       where [OPTION]...     is a list of zero or more optional arguments\n";
    print STDERR "             REF_DIR         is the name of the reference directory\n";
    print STDERR "             TEST_DIR        is the name of the test directory\n";
    print STDERR "\n";
    print STDERR "Miscellaneous arguments:\n";
    print STDERR "  --version                  display program version information\n";
    print STDERR "\n";
    exit(1);
}

############################################################
# ParseParameters()
#
# Parse program parameters.
############################################################

sub ParseParameters 
{
    if (@ARGV < 3)
    { 
	Usage(); 
    }
    @dir_names = ();

    if ($ARGV[0] ne "protein" && $ARGV[0] ne "rna")
    {
	print STDERR "ERROR: First argument should be \"protein\" or \"rna\".\n";
	exit(1);
    }
    else
    {
	$SCORE_PREDICTION .= " ".$ARGV[0];
    }
    
    for (my $argno = 1; $argno < @ARGV; $argno++)
    {
	if ($ARGV[$argno] eq "--version")
	{
	    Version();
	} 
	elsif ($ARGV[$argno] eq "--core")
	{
	    $SCORE_PREDICTION .= " --core";
	}
	else 
	{
	    push(@dir_names, $ARGV[$argno]);
	}
    }

    if (@dir_names != 2)
    {
	print STDERR "ERROR: Incorrect number of directory names.\n";
	exit(1);
    }
}

############################################################
# main()
#
# Main program.
############################################################

ParseParameters();

my $ref_dir_name = $dir_names[0];
my $test_dir_name = $dir_names[1];

my @ref_filenames = `ls $ref_dir_name`;
my @test_filenames = `ls $test_dir_name`;

my $count = 0;
my $Q_total = 0;
my $fM_total = 0;
my $sens_total = 0;
my $ppv_total = 0;

foreach my $test_file (@test_filenames)
{
    chomp $test_file;

    # check for reference file

    if (-e "$ref_dir_name/$test_file")
    {
	my @ret = `$SCORE_PREDICTION $ref_dir_name/$test_file $test_dir_name/$test_file`;
	
	# now, parse results from score_prediction
	
	$ret[0] =~ /Q=([e\.0-9+-]+); fM=([e\.0-9+-]+); sens=([e\.0-9+-]+); ppv=([e\.0-9+-]+)/;
	
	$count++;
	$Q_total += $1;
	$fM_total += $2;
	$sens_total += $3;
	$ppv_total += $4;
    }
    else
    {
	print STDERR "ERROR: Unable to find test file \"$test_file\" in reference directory.\n";
    }
}

$Q_total /= $count;
$fM_total /= $count;
$sens_total /= $count;
$ppv_total /= $count;

# print results

print "ref=$ref_dir_name; test=$test_dir_name; n=$count; Q=$Q_total; fM=$fM_total; sens=$sens_total; ppv=$ppv_total;\n";
