#!/usr/bin/perl
use strict;

# input file on command line: perl smatparse.pl ifile.in
my $iFile = $ARGV[0];
#variables
my @line = "";
my @line2 = "";
my $test = 0;
my $oFile = "stot.txt";

# open input file
open(FILE, '<', $iFile) or die "Could not open $iFile:  $!";

# open output file
open(TEMP, '>', $oFile) or die "Could not open $oFile:  $!";

# read input file line by line
while (<FILE>) {
   #get rid of new line
   chomp;
   # split by white spaces
   @line = split(/ +/);

   # if this is the real or imaginary part of the S matrix
   # then check if the first number is an integer
   # construct rows of the S matrix  
   if ($test == 1 || $test == 2){
   if ($line[1] =~ /^-?\d+$/){
      #print "$line[1]\n";
      print TEMP "\n";
      $test = 1;
      }
      if ($line[2] eq 'matrix,' || $line[2] eq 'matrix,imag:'){
         # do nothing, don't print 'S matrix,imag:'
	 }
      else{
         print TEMP "$_";
         }
      }

   # get the energy and momentum for the given phase shift
   # tells the code above what part of the S matrix you have
   if ($line[2] eq 'matrix,real:'){
      print TEMP "$line2[4] $line2[7]"; 
      print "$line2[4] $line2[7] \n";
      $test = 1;
      }
   elsif ($line[2] eq 'matrix,imag:'){
      #print TEMP "\n";
      $test = 2;
      }
   elsif($line[2] eq 'matrix,'){
      print TEMP "\n";
      $test = 0;
      }
   @line2 = @line;
}