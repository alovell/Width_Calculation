#!/usr/bin/perl
use strict;
use warnings;

# takes in the potentials and separates out the channels into different files

# define a trim command
sub ltrim { my $s = shift; $s =~ s/^\s+//;       return $s };

# reads input file from command line
my $iFile = $ARGV[0];
my @oFile;

# variables

# open input file
open (FILE, '<', $iFile) or die "Could not opne $iFile:  $!";

# read input line by line and match patterns
while (<FILE>) {
   chomp;
   # split line based on white spaces
   my @line = split(/ +/);
   # check for 3body header
   if ($_ =~ /THREE BODY POTENTIAL/) {
      #print "found!\n";
      my $string1 = <FILE>; #blank line
      my $string2 = <FILE>; #header line V(i,j)
      my @line = split(/V/,$string2);
      #print "$line[1] $line[2] $line[3]\n";
      # $line[1] - Ro, $line[2] - fm $line[3] & up V(i,j)
      # split V(i,j) to get i_j for file name
      for (my $j=1; $j <= 9; $j++) {
         my $v1 = ltrim(substr($line[$j],1,2));
	 my $v2 = ltrim(substr($line[$j],4,2));
	 #print "$v1 $v2\n";
	 #print "j= $j line= $line[$j] v1= $v1 v2= $v2\n";
	 #print "substr($line[$j],0) substr($line[$j],1) substr($line[$j],2) substr($line[$j],3) substr($line[$j],4)\n";
         $oFile[$j] = $v1 . "_" . $v2 . ".pot";
      }
      for (my $i=0; $i <= 399; $i++) {
         my $string = <FILE>;
	 chomp($string);
         #print "$string";
	 my @line = split(/ +/,$string);
	 #print "$line[0] $line[1] $line[2] $line[3] $line[4]\n";
	 for (my $j=1; $j <=9; $j++) {
	    open(TEMP, '>>', $oFile[$j]) or die "Could not open file :$!";
	    print TEMP "$line[1] $line[$j+1]\n";
	    close(TEMP);
	 }
      }
   }
   
}