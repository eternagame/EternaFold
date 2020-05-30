#!/usr/bin/perl

use strict;
use Cwd qw(realpath);
use File::Basename;

# Globals

my $SCORE_DIRECTORY = dirname(realpath($0))."/score_directory.pl";
my @dir_names = ();

############################################################
# Version()
#
# Display program version.
############################################################
 
sub Version 
{
    print STDERR "ROCarea version 1.00 - Compute ROC area for a set of RNA secondary structure predictions\n";
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
# CompareDirectoryNames()
#
# Compare two directory names for sorting.
############################################################

sub CompareDirectoryNames
{
    $a =~ /gamma=([.0-9]*)/;
    my $gamma_a = $1;
    $b =~ /gamma=([.0-9]*)/;
    my $gamma_b = $1;

    if ($gamma_a < $gamma_b) {
	return -1;
    } elsif ($gamma_a > $gamma_b) {
	return 1;
    } else {
	return 0;
    }
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
	$SCORE_DIRECTORY .= " ".$ARGV[0];
    }
    
    for (my $argno = 1; $argno < @ARGV; $argno++)
    {
	if ($ARGV[$argno] eq "--version")
	{
	    Version();
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
chomp(my @test_subdirs = `ls $test_dir_name`);

my @sens = ();
my @ppv = ();

# compute score for each subdirectory

for (my $i = 0; $i < @test_subdirs; $i++)
{
    print ".";
    my @ret = `$SCORE_DIRECTORY $ref_dir_name $test_dir_name/$test_subdirs[$i]`;
    $ret[0] =~ /sens=([e\.0-9+-]+); ppv=([e\.0-9+-]+)/;
    push(@sens, $1);
    push(@ppv, $2);
}

print "\n";

# sort scores by increasing PPV (and decreasing sensitivity)

for (my $i = 0; $i < @ppv; $i++)
{
    for (my $j = $i+1; $j < @ppv; $j++)
    {
	if ($ppv[$j] < $ppv[$i] || ($ppv[$j] == $ppv[$i] && $sens[$j] > $sens[$i]))
	{
	    ($sens[$i], $sens[$j]) = ($sens[$j], $sens[$i]);
	    ($ppv[$i], $ppv[$j]) = ($ppv[$j], $ppv[$i]);
	    ($test_subdirs[$i], $test_subdirs[$j]) = ($test_subdirs[$j], $test_subdirs[$i]);
	}
    }
}

# compute ROC area

my $area = 0;
my $prev_sens = $sens[0];
my $prev_ppv = 0;

for (my $i = 0; $i < @ppv; $i++)
{
    $area += ($ppv[$i] - $prev_ppv) * ($prev_sens + $sens[$i]) / 2.0;
    $prev_sens = $sens[$i];
    $prev_ppv = $ppv[$i];

    chomp(my $count = `ls $test_dir_name/$test_subdirs[$i] | wc | awk '{print \$1;}'`);
    print "ref=$ref_dir_name; test=$test_dir_name/$test_subdirs[$i]; n=$count; sens=$sens[$i]; ppv=$ppv[$i];\n";
}

print "\nroc=$area;\n";
