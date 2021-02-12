#!/usr/bin/perl

use strict;

sub WriteDefaultValues 
{
    my ($filename, $function_name) = @_;
    open (INFILE, "<$filename");

    my @names = ();
    my @values = ();
    
    while (my $line = <INFILE>)
    {
        $line =~ /(\S+)\s+(\S+)/;
        push (@names, $1);
        push (@values, $2);
    }
    my $length = @names;

    close INFILE;

    print OUTFILE "/////////////////////////////////////////////////////////////////\n";
    print OUTFILE "// $function_name()\n";
    print OUTFILE "//\n";
    print OUTFILE "// Retrieve default parameter values.\n";
    print OUTFILE "/////////////////////////////////////////////////////////////////\n";
    print OUTFILE "\n";
    print OUTFILE "std::vector<RealT> $function_name()\n";
    print OUTFILE "{\n";
    print OUTFILE "    RealT values[] =\n";
    print OUTFILE "    {\n";
    for (my $i = 0; $i < @values; $i++)
    {
	my $value = sprintf("%-20s", sprintf("%.10lf%s", $values[$i], ($i == $#values ? " " : ",")));
	print OUTFILE "        $value // $names[$i]\n";
    }
    print OUTFILE "    };\n\n";
    print OUTFILE "    return std::vector<RealT>(values, values + ".scalar(@values).");\n";
    print OUTFILE "}\n\n";
}

sub main 
{
    if (@ARGV != 3)
    {
        print STDERR "Usage: perl MakeDefaults.pl PARAMFILE_COMPLEMENTARY PARAMFILE_NONCOMPLEMENTARY PARAMFILE_PROFILE\n";
        exit(1);
    }

    open(OUTFILE, ">Defaults.hpp");
    print OUTFILE "#pragma once \n\n";
    WriteDefaultValues($ARGV[0], "GetDefaultComplementaryValues");
    print OUTFILE "\n";
    WriteDefaultValues($ARGV[1], "GetDefaultNoncomplementaryValues");
    print OUTFILE "\n";
    WriteDefaultValues($ARGV[2], "GetDefaultProfileValues");
    close OUTFILE;
}

main();
