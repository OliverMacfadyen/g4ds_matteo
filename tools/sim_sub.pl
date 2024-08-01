#!/usr/bin/perl -w
## Author: Stefano Davini <stefano.davini@ge.infn.it>
##
## This perl script is used to generate a mac file 
## for g4ds and submit a similation job
##

use Cwd;
use strict;
use warnings;

my $workdir 	 = getcwd;
my $tooldir 	 = $workdir.'/../tools/';
my $out_base; 				# name of output directory
my $file_base; 				# base name for mac, fil, and root files
my $verbosity	 = 1;			# verbosity
my $wrph	 = 0;			# write photons?
my $wrdau	 = 0;			# write daughters?
my $wrdep	 = 0;			# write deposits?
my $checkoverlap = 0;			# check overlap?
my $seed	 = 12345678;		# random seed
my $nev     	 = 100;			# number of events to simulate
my $detector	 = 'ds50';		# detector to simulate
my $config	 = 1;			# detector configuration
my $tmb	    	 = 0.05;		# TMB volume fraction
my $gen;				# primary particle generator
my $pos;				# position generator
my $energy  	 = 100;			# energy in keV (if using g4gun)
my $energyfile;				# energy file (used to generate neutrons from volumes)
my $x		 = 0.;			# position in x (m) (if using sphere_radius)
my $y		 = 0.;			# position in y (m)
my $z 		 = 0.; 			# position in z (m)
my $vetoyf 	 = 1.; 			# veto yield factor
my $interactive  = 0;			# run in interactive?
my $batch        = 0;			# submit job on batch?
my $cluster	 = 'cnaf';		# cluster farm
my $run_g4ds 	 = 1;			# run g4ds?
my $run_g4rooter = 0;			# run g4rooter?
my $run_g4rooter_full = 0;		# run g4rooter_full?
my $executable   = '';			# an additional command to execute after all (usually analysis)
my $cleardir     = 0;			# clean output directory?
my $no_rewrite   = 0;

unless (-x $workdir."/g4ds"){
	die " ERROR: g4ds is not installed in $workdir \n";	
}


foreach my $arg (@ARGV){
	if (($arg =~ /^help$/i) || ($arg =~/^usage$/i)){
		usage();
		exit;
	}
	#global options
	elsif ($arg =~ /^out\=(.*)/i){
		$out_base = $1;
		print "INFO: output directory = $out_base\n";
	}
	elsif ($arg =~ /^file\=(.*)/i){
		$file_base = $1;
		print "INFO: file name = $file_base\n";
	}
	elsif ($arg =~ /^verbosity\=(\d+)/i){
		$verbosity = $1;
		print "INFO: verbosity = $verbosity\n";
	}
	elsif ($arg =~ /^wrph\=(\d+)/i){
		$wrph = $1;
		if ($wrph > 1) { $wrph=1}
		print "INFO: write photons = $wrph\n";
	}
	elsif ($arg =~ /^wrdep\=(\d+)/i){
		$wrdep = $1;
		if ($wrdep > 1) { $wrdep=1}
		print "INFO: write deposits = $wrdep\n";
	}
	elsif ($arg =~ /^wrdau\=(\d+)/i){
		$wrdau = $1;
		if ($wrdau > 1) { $wrdau=1}
		print "INFO: write daughters = $wrdau\n";
	}
	elsif ($arg =~ /^gen\=(.*)/i){
		$gen = lc $1;
		if (($gen eq 'electron')    || ($gen eq 'e-')) { $gen='e-'}
		elsif (($gen eq 'positron') || ($gen eq 'e+')) { $gen='e+'}
		elsif (($gen eq 'gamma')    || ($gen eq 'g'))  { $gen='gamma'}
		elsif (($gen eq 'proton')   || ($gen eq 'p'))  { $gen='proton'}
		elsif (($gen eq 'neutron')  || ($gen eq 'n'))  { $gen='neutron'}
		elsif (($gen eq 'alpha')    || ($gen eq 'a'))  { $gen='alpha'}
		elsif (($gen eq 'ambe')     || ($gen eq lc 'AmBeSource')) {$gen='AmBeSource'}
		elsif (($gen eq 'n_u238'))   { $gen='n_u238'; $energyfile='1';}
		elsif (($gen eq 'n_th232'))  { $gen='n_th232'; $energyfile='1';}
		else {die " ERROR: unknown generator $gen \n";}
		print "INFO: generator = $gen\n";
	}
	elsif ($arg =~ /^pos\=(.*)/i){
		$pos = lc $1;
		if    ( ($pos eq 'sphere') || ($pos eq 'sphere_radius')) {$pos = 'sphere_radius'}
		elsif ( ($pos eq 'liquidargon')  || ($pos eq 'lar')) {$pos = 'liquidargon'}
		elsif ( ($pos eq 'liquidscintillator') || ($pos eq 'ls')) {$pos = 'liquidscintillator'}
		elsif ( ($pos eq 'bgd_teflon')   || ($pos eq 'teflon')) {$pos = 'bgd_teflon'}
		elsif ( ($pos eq 'bgd_cryostats') || ($pos eq 'cryostats')) {$pos = 'bgd_cryostats'}
		elsif ( ($pos eq 'bgd_sipm') || ($pos eq 'sipm')) {$pos = 'bgd_sipm'}
		elsif ( ($pos eq 'bgd_fused_silica') || ($pos eq 'fused_silica')) {$pos = 'bgd_fused_silica'}
		elsif ( ($pos eq 'bgd_pmt_photocathode') || ($pos eq 'pmt_photocathode')) {$pos = 'bgd_pmt_photocathode'}
		elsif ( ($pos eq 'bgd_pmt_stem') || ($pos eq 'pmt_stem')) {$pos = 'bgd_pmt_stem'}
		elsif ( ($pos eq 'bgd_rings') || ($pos eq 'rings')) {$pos = 'bgd_rings'}
		elsif ( ($pos eq 'bgd_grid') || ($pos eq 'grid')) {$pos = 'bgd_grid'}
		else {die " ERROR: unknown spatial generator $pos \n";}
		print "INFO: Spatial generator = $pos \n";
        } 
	elsif ($arg =~ /^energy\=(.*)/i){
		$energy = $1;
		if (lc $energy eq 'file') {
			$energyfile = 1;
			print "INFO: Energy from file (pos must be set to a meaningful value) \n";
		}
		else {
			if ($energy<=0) {die " ERROR: kinetic energy must be positive \n";}
			print "INFO: Energy = $energy keV \n";
		}
        }
	elsif ($arg =~ /^x\=(.*)/i){
		$x = $1;
		print "INFO: x-coordinate of primary event = $x cm\n";
	}
	elsif ($arg =~ /^y\=(.*)/i){
		$y = $1;
		print "INFO: y-coordinate of primary event = $y cm\n";
	}
	elsif ($arg =~ /^z\=(.*)/i){
		$z = $1;
		print "INFO: z-coordinate of primary event = $z cm\n";
	}
	elsif ($arg =~ /^nev\=(\d+)/i){
		$nev = $1;
		if ($nev<=0) {$nev=1;}
		print "INFO: number of events to simulate = $nev\n";
	}
	elsif ($arg =~ /^seed\=(\d+)/i){
		$seed = $1;
		if ($seed<=0) {$seed=12345678;}
		print "INFO: hep random seed = $seed\n";
	}
	elsif ($arg =~ /^detector\=(.*)/i){
		$detector = $1;
		if ( (lc $detector eq 'ds50') ) {$config = 1;}
		elsif ( (lc $detector eq 'ds20k') ) {$config = 10;}
		else {die " ERROR: unknown detector $detector \n"; }
		print "INFO: detector = $detector, configuration = $config \n";
        }
	elsif ($arg =~ /^tmb\=(.*)/i){
		$tmb = $1;
		if (($tmb<0) || ($tmb>1)) {die " ERROR: TMB fraction must be between 0 and 1, not $tmb \n";}
		print "INFO: TMB fraction = $tmb\n";
        }
	elsif ($arg =~ /^vetoyieldfactor\=(.*)/i){
		$vetoyf = $1;
		if (($vetoyf<0)) {die " ERROR: veto yield factor fraction must be >=0, not $vetoyf \n";}
		print "INFO: veto yield factor = $vetoyf\n";
        }
	elsif ($arg =~ /^run_g4ds\=(\d+)/i){
		$run_g4ds = $1;
        }
	elsif ($arg =~ /^run_g4rooter\=(\d+)/i){
		$run_g4rooter = $1;
	}
	elsif ($arg =~ /^run_g4rooter_full\=(\d+)/i){
		$run_g4rooter_full = $1;
	}
	elsif ($arg =~ /^exec\=(.*)/i){
		my $exec = $1;
		if ((lc $exec eq 'interactive')) {$interactive=1; $batch=0; }
		elsif ((lc $exec eq 'batch')) {$interactive=0; $batch=1; }
		else {die " ERROR: Unknown execution mode $exec \n"; }
	}
	elsif ($arg =~ /^cleardir\=(\d+)/i){
		$cleardir = $1;
	}
	elsif ($arg =~ /^executable\=(.*)/i){
		$executable = $1;
		print " INFO: additional executable: $executable \n";
	}
	else {
		die " ERROR: unknown parameter $arg \n";
	}
}

unless ($out_base) {
	die " ERROR: you need to specify an output folder (out=)\n";
}

if (($out_base eq '.') || ($out_base eq $workdir)){
	die " ERROR: the output directory can not be the working directory\n";
}

if ($run_g4rooter){
	unless (-x $tooldir."/g4rooter"){
		die " ERROR: g4rooter is not installed\n";	
	}
}

if ($run_g4rooter_full){
	unless (-x $tooldir."/g4rooter_full"){
		die " ERROR: g4rooter_full is not installed\n";	
	}
}

if ($cleardir) {
	print "INFO: removing all files in $out_base \n";
	system("rm -rf $out_base") ;
}

# check: existence of output folder; if not exist, create it
unless (-d $out_base) {
	mkdir $out_base or die " ERROR: i can't create OUTPUT DIRECTORY $out_base \n";
	print "INFO: created directory $out_base \n";
}

unless (defined($gen)){
	die " ERROR: you need to specify a primary particle generator (gen=)\n";
}

unless (defined($pos)){
	die " ERROR: you need to specify a spatial generator (pos=)\n";
}

if ($energyfile){
	if ($pos eq 'bgd_teflon') {
		if ($gen eq 'n_u238') { $energyfile = "../data/physics/U238Teflon.dat"; $gen='neutron';}
		elsif ($gen eq 'n_th232') { $energyfile = "../data/physics/Th232Teflon.dat"; $gen='neutron';}
		else { $energyfile = "../data/physics/U238Teflon.dat"; }
	}
	elsif ($pos eq 'bgd_cryostats') {
		if ($gen eq 'n_u238') { $energyfile = "../data/physics/U238StainlessSteel.dat"; $gen='neutron';}
                elsif ($gen eq 'n_th232') { $energyfile = "../data/physics/Th232StainlessSteel.dat"; $gen='neutron';}
                else { $energyfile = "../data/physics/U238U235Th232StainlessSteel.dat"; }
	}
	elsif ($pos eq 'bgd_fused_silica') {
		if ($gen eq 'n_u238') { $energyfile = "../data/physics/U238FusedSilica.dat"; $gen='neutron';}
                elsif ($gen eq 'n_th232') { $energyfile = "../data/physics/Th232FusedSilica.dat"; $gen='neutron';}
                else { $energyfile = "../data/physics/U238FusedSilica.dat"; }
	}
	else {die " ERROR: energy file generator not defined for $pos \n";}
}


unless (defined($file_base)){
	$file_base = $gen;
	print "INFO: file name automatically set to $file_base\n";
}



my $macfilename  = $file_base.'.mac';
my $filfilename  = $file_base.'.fil';
my $rootfilename = $file_base.'.root';
my $shfilename   = $file_base.'.sh';
my $outfilename  = $file_base.'.out';
my $errfilename  = $file_base.'.err';

# change $out_base path to absolute path
chdir $out_base;
$out_base = getcwd;

my $basefilepath = $out_base.'/'.$file_base;
my $macfilepath  = $out_base.'/'.$macfilename;
my $filfilepath  = $out_base.'/'.$filfilename;
my $rootfilepath = $out_base.'/'.$rootfilename;
my $shfilepath   = $out_base.'/'.$shfilename;
my $outfilepath  = $out_base.'/'.$outfilename;
my $errfilepath  = $out_base.'/'.$errfilename;

if (-e $macfilepath){
	print " WARNING: a mac file with the same path already exists: $macfilepath; the file will be overwritten; \n";
}
if (-e $filfilepath){
	if ($no_rewrite) { die " INFO: a fil file with the same path already exists: $filfilepath. Exiting; \n"}
	print " WARNING: a fil file with the same path already exists: $filfilepath;\n";
}

my $eventcounter = ($interactive) ? ($nev/10) : ($nev/100);
if ($eventcounter<1) {$eventcounter = 1;}

# create macfile
write_mac_file();

# create sh script file
write_sh_file();
system("chmod u+x $shfilepath");

# run the sh script
if ($interactive) {
	system("$shfilepath");
}
elsif ($batch) {
	my $batchcmd;
	if ($cluster eq 'cnaf'){
		$batchcmd = "bsub -q darkside -o $outfilepath -e $errfilepath $shfilepath ";
	} 
	system("$batchcmd");
}


print "Exiting sim_sub.\n";


# routine to write the sh file
sub write_sh_file{
	# g4rooter command 
	my $g4rc = ($run_g4rooter) ? './../tools/g4rooter '.$filfilepath."\n" : "";
	# run g4rooter_full command
	my $g4rfc = ($run_g4rooter_full) ? './../tools/g4rooter_full '.$filfilepath."\n" : "";

	my $SH; # macfile handle
		open $SH, ">", $shfilepath or die " ERROR: I can't open $shfilepath \n";
	print $SH " \n";
	print $SH "cd $workdir \n";
	if ($run_g4ds) {print $SH "./g4ds $macfilepath \n"}
	print $SH "$g4rc";
	print $SH "$g4rfc";
	print $SH "$executable";
	close $SH;
}

# routine to write the mac file
sub write_mac_file{
	print "INFO: Creating macfile $macfilepath \n";
	my $MAC; # macfile handle
		open $MAC, ">", $macfilepath or die " ERROR: I can't open $macfilepath \n";
	print $MAC "# macfile created by sim_sub.pl \n \n";

	print $MAC "/ds/manager/log routine \n";
	print $MAC "/ds/manager/verbosity $verbosity \n";
	print $MAC "/ds/manager/checkoverlap $checkoverlap \n";
	print $MAC "/ds/manager/eventcounter $eventcounter \n";
	print $MAC "/ds/manager/writephotons $wrph \n";
	print $MAC "/ds/manager/writedeposits $wrdep \n";
	print $MAC "/ds/manager/writedaughters $wrdau \n";
	print $MAC "/ds/manager/TMB_fraction $tmb \n";
	print $MAC "/ds/detector/configuration $config \n";
	if ($config == 10 ){
		print $MAC "/ds/detector/ds20lsv_detector 1 \n";
		print $MAC "/ds/detector/ds20lsv_diameter 8 m \n";
	}
	print $MAC "/run/filename $basefilepath\n";
	print $MAC "/run/heprandomseed $seed \n";
	if (($gen eq 'AmBeSource') || ($gen eq 'neutron')){
		print $MAC "/ds/physics/hadronic_list HP \n";
	}
	else {
		print $MAC "/ds/physics/hadronic_list none \n";
	}
	print $MAC "/ds/physics/em_list livermore \n";
	print $MAC "/ds/physics/optics 3 \n";
	print $MAC "/ds/physics/LSoptics 2 \n";
	print $MAC "/ds/detector/vetoyieldfactor $vetoyf \n";
	print $MAC "/ds/detector/ExtLarScintillating 0 \n";
	print $MAC "/ds/physics/killS1S2 1 \n";
	print $MAC "/ds/physics/killS2 1 \n";
	print $MAC "/run/initialize \n";
	if (($gen eq 'e-') || ($gen eq 'e+') || ($gen eq 'gamma') || ($gen eq 'proton') || ($gen eq 'neutron') || ($gen eq 'alpha')){
		print $MAC "/ds/generator/select G4Gun\n";
		print $MAC "/ds/generator/particle $gen\n";
		unless (defined($energyfile)){
			print $MAC "/ds/generator/energy $energy keV\n";
		}
		else { print $MAC "/ds/generator/energyfile $energyfile \n";}
	}
	elsif ($gen eq 'AmBeSource'){
		print $MAC "/ds/generator/select AmBeSource\n";
		print $MAC "/ds/generator/AmBe/source neutron0G\n";
		print $MAC "/ds/generator/AmBe/disable gamma\n";
	}
	if ($pos eq 'sphere_radius'){
		print $MAC "/ds/generator/sphere_radius 1. mm\n";
		print $MAC "/ds/generator/set_center $x $y $z cm\n";
	}
	elsif ($pos eq 'liquidargon'){
		print $MAC "/ds/generator/liquidargon 1 \n";
	}
	elsif ($pos eq 'liquidscintillator'){
		print $MAC "/ds/generator/liquidscintillator 1 \n";
	}
	elsif ($pos eq 'bgd_teflon'){
		print $MAC "/ds/generator/bgd_teflon 1 \n";
	}
	elsif ($pos eq 'bgd_cryostats'){
		print $MAC "/ds/generator/bgd_cryostats 1 \n";
	}
	print $MAC "/run/beamOn $nev\n";
	close $MAC;

}

sub usage{
  print "\n $0 [options]\n\n";
  print" options:\n";
  print" out=\t\t\t Output folder\n";
  print" gen=[]\t\t\t Particle generator [gamma, e+, e-, p, n, alpha, AmBe, ... ]\n";
  print" pos=[]\t\t\t Position generator [sphere, LAr, LS, teflon, ... ]\n";
  print" file=\t\t\t Optional: Output file base names\n";
  print" nev=\t\t\t Optional: Number of events to simulate  (default $nev)\n";
  print" seed=\t\t\t Optional: Random seed  (default $seed)\n";
  print" detector=\t\t Optional: Detector [ds50, ds20k]  (default $detector)\n";
  print" energy=\t\t Optional: Energy of the particle in keV  (default $energy keV)\n";
  print" x=\t\t\t Optional: X-position in m, for sphere generator  (default $x m)\n";
  print" y=\t\t\t Optional: Y-position in m, for sphere generator  (default $y m)\n";
  print" z=\t\t\t Optional: Z-position in m, for sphere generator  (default $z m)\n";
  print" vetoyieldfactor=\t Optional: Light yield factor for veto  (default $vetoyf)\n";
  print" run_g4ds=[0/1]\t\t Optional: Run g4ds?  (default $run_g4ds)\n";
  print" run_g4rooter=[0/1]\t Optional: Run g4rooter?  (default $run_g4rooter)\n";
  print" run_g4rooter_full=[0/1] Optional: Run g4rooter_full?  (default $run_g4rooter_full)\n";
  print" exec=[] \t\t Optional: How to execute [interactive, batch]\n";
  print" cleardir=[0/1] \t Optional: Remove all files in the output folder  (default $cleardir)\n";
  print" executable= \t\t Optional: Run an executable after?\n";
  #print"\t v=[v]\t\t Optional: V  (default $v)\n";
  print" verbosity=\t\t Optional: Verbosity  (default $verbosity)\n";
  print" wrph=[0/1]\t\t Optional: Write photons?  (default $wrph)\n";
  print" wrdep=[0/1]\t\t Optional: Write deposits?  (default $wrdep)\n";
  print" wrdau=[0/1]\t\t Optional: Write daughter?  (default $wrdau)\n";
  print"\n";
}
