#!/usr/bin/perl

$trainingfile = $ARGV[0];
$nevent = $ARGV[1];

print "Processing $nevent events from the $ARGV[0] training file...\n"; 


$path = $ENV{PWD};

if($path =~ /(.*)g4ds10.*/) {
  $newpath = "${1}/g4ds10";

}

print "$newpath\n";
system("g++ ${newpath}/tools/g4rootered_prereco.C -o ${newpath}/tools/g4rooter_prereco `root-config --cflags --glibs`");
print "g4rooter_prereco compiled\n";
system("${newpath}/tools/g4rooter_prereco $trainingfile nevents=$nevent ${newpath}/Linux-g++/covariance.root");
system("mv ${newpath}/Linux-g++/fitxMDF.cxx ${newpath}/tools/reco/");
system("mv ${newpath}/Linux-g++/fityMDF.cxx ${newpath}/tools/reco/");
system("mv ${newpath}/Linux-g++/pca.C ${newpath}/tools/reco");
system("g++ ${newpath}/tools/g4rootered_reco.C -o ${newpath}/tools/g4rooter_reco `root-config --cflags --glibs`");

print "...covariance matrix and multi-fitter scripts produced and copied in ../tools\n";
print "g4rooter_reco compiled. \n";
print "\nHave fun!!\n";

