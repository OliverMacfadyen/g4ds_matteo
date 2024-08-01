#!/usr/bin/perl -w
# Author: Stefano Davini <stefano.davini@gssi.infn.it>
#
# This perl script is used to generate placement files for PMT positions in the LSV;
# it is intended to be used in the design phase of DS-20k
#

use strict;
use warnings;
use Math::Trig;
use Math::Trig ':pi';

# global default variables
my $filename = 'DS20kPMTNeutronVeto.dat';
my $pmt_type = 2;    # 1-> 8 inches, 2-> 20 inches
my @theta_deg;       # azimuthal angle in degree
my @row_nums;  
# define basic placements depending on the PMT type
# numbers given from Yury Suvorov (private conversation)
# 8 inches
my @theta_deg1 = (20, 35, 50, 65, 80, 100, 115, 130, 145, 160); 
my @row_nums1  = (15, 45, 50, 60, 80,  80,  60,  50,  45,  15);
# 20 inches
my @theta_deg2 = (20, 40, 60, 80, 100, 120, 140, 160);
my @row_nums2  = ( 5, 10, 15, 20,  20,  15,  10,   5);

# user input parsing
foreach my $arg (@ARGV){
  if (($arg =~ /^help$/i) || ($arg =~/^usage$/i)){
    #usage();
    exit;
}
elsif ($arg =~ /^pmt_type\=(\d+)/){
  $pmt_type = $1;
  if ($pmt_type==1){
    @theta_deg = @theta_deg1;  
    @row_nums  = @row_nums1; 
  }
  elsif ($pmt_type==2){
    @theta_deg = @theta_deg2;
    @row_nums  = @row_nums2;
  }
  else {die " ERROR: pmt_type = $pmt_type not defined \n";}
}
elsif ($arg =~ /^filename\=(.*)/){
  $filename = $1;
}
# future development: add possibility to change theta_deg and row_nums with user input
else {
  die " ERROR: unknown parameter $arg \n";
  }
}


# user input processing




# create the file with PMT placement info.
# the format is the following:
# 1st number is the total number of PMTs
# then 2 colums, with angles theta and phi in radians
my $FILE;
open  $FILE, ">", $filename or die " ERROR: can't open $filename \n";

my $tot_num_pmts = 0;
foreach my $num (@row_nums){
  $tot_num_pmts += $num;
}
print $FILE "$tot_num_pmts"."\n";

my $row = 0;
foreach my $theta (@theta_deg){
  my $theta_rad  = deg2rad($theta);
  my $phi_rad    = 0.;     # starting angle, will be increased in next loop
  my $row_num    = $row_nums[$row];
  my $phi_step   = ($row_num!=0) ? (2.*pi/$row_num) : 0.; # (avoids division by zero)
  for (my $i = 0; $i<$row_num; $i++){
    print $FILE "$theta_rad \t $phi_rad"."\n";
    $phi_rad += $phi_step;
  }
  $row++;
} 

close $FILE;
