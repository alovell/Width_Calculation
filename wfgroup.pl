#!/usr/bin/perl
use strict;
use warnings;

# takes in fort.11 and groups the wave functions based on quantum numbers

#read input file from command line
#probably fort.11
my $iFile = $ARGV[0];

#variables
my @oFile; #output file
my $count = 0; #counter
my $poleE = 0; #pole energy
my @ival; #channel number
my @norm; #norm
my @nval; #channel label
my @Kval; #K
my $ncount; #number of channels with norm > 0.01
my @nstates; #holds channel numbers for above
my $tempnstates; #same as above, changes for each iteration

#open input file
open (FILE, '<', $iFile) or die "Could not open $iFile:  $!";

#read input file and match patterns
while (<FILE>) {
   chomp;
   #split line based on spaces
   my @line = split(/ +/);
   if ($_ =~ /# Bound /) {
      $poleE = $line[5];
      #print "$poleE\n";
   }
   elsif ($_ =~ /#    i: /) { 
      # rest the counter, nstates
      $count = -1;
      $ncount = -1;
   }
   elsif ($_ =~/#  /){
      $count++;
      $ival[$count] = $line[1];
      #remove the semicolon
      $ival[$count] =~ tr/://d;
      $norm[$count] = $line[2];
      $nval[$count] = $line[3];
      $Kval[$count] = $line[4];
      if ($norm[$count] > 0.01){
         $ncount++;
	 $nstates[$ncount] = $count;
         print "$ival[$count] $norm[$count] $nval[$count] $Kval[$count] \n";
         # create a file for all of the wave function
         $oFile[$count] = $ival[$count] . "K" . ".wf";
         print "$oFile[$count]\n";
	 # put relivant numbers in output file - useful for width calc
         open (TEMP, '>>', $oFile[$count]) or die "Could not open wave function file :$!";
         print TEMP "$poleE $norm[$count] $nval[$count] $Kval[$count] \n";
	 close (TEMP);
      }
   }
   elsif ($_ =~/&/){
      # reset the counter if there is an & 
      $ncount = -1;
   }
   else {
      for (my $i=0; $i <= $ncount; $i++){
         # print the wave functions to correct files (one for each set of quantum numbers
	 #print "$oFile[$nstates[$i]] $i";
	 open (TEMP, '>>', $oFile[$nstates[$i]]) or die "Could not open wave function file :$!";
	 $tempnstates = $nstates[$i];
	 if ($tempnstates > 9){
	    $tempnstates = $tempnstates % 10;
	 }
	 print TEMP "$line[1] $line[$tempnstates+2]\n";
	 print "$_ $nstates[$i] $tempnstates\n";
	 close(TEMP);
	 #print "\n";
      }
   }
}   