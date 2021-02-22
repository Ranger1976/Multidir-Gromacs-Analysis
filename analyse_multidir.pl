#!/usr/bin/perl -w

use strict;
use Getopt::ArgParse;
use Cwd;
use List::Util qw(max);
use Statistics::Basic qw(:all);
use Statistics::PointEstimation;
# use Data::Dumper::Perltidy;

my $cwd=getcwd;

my $parser = Getopt::ArgParse->new_parser(
    help => 'Analyses trajectories from gromacs',
    description =>  
    "$0 is used to run typical analysis of multiple Gromacs runs, where \n".
    "each run is kept in its own directory and all runs are supposed to represent\n".
    "the same system.\n".
    "\n".
    "Available choices for analysis (select with --what):\n".
    "rmsd (1)\n".
    "rmsf (1)\n".
    "rgyr (1)\n".
    "sasa (1)\n".
    "hbond (protein-protein) (1)\n".
    "cluster (1)\n".
    "hbond_pl (hbond protein-ligand) \n".
    "contact_pl (protein-ligand contacts, within 0.4nm)\n".
    "dipole (of the protein)\n".
    "\n".
    "(1) run by default when running $0 --what all_basic \n".
    "If you run the program with --traj_prepare [group], analysis will not be run.\n".
    "Instead, trajectories will be generated that include only this group. This \n".
    "can be useful to speed up calculations afterwards. For example, you can run\n".
    "$0 --traj_prepare protein --deffnm md\n".
    "In this case, in each directory a trjectory named md_protein.xtc will be \n".
    "generated that includes only the protein, and you can thereafter run analysis\n".
    "with these files using the flag --deffnm md_protein\n",
    epilog      => 'Written by Ran Friedman, Linnaeus University, Kalmar, Sweden',
);

$parser->add_args(
    ['--dir', default => "sim", help => "generic directory name. If this is e.g., sim, will look for sim0, sim1 etc", metavar=>"dir"],
    ['--directories', type => 'Array', split =>',', 
    help => "a list of directories to check. --dir is ignored if --directories is set.", metavar=>"sim0,sim1,sim2"],
    ['--deffnm', default => "md", help => "generic name for files in each dir", metavar=>"name"],
    ['--what', type => 'Array', split =>',', default => "all_basic", 
    help => "what to analyse (choose one or more, comma-separated and no spaces in-between)", 
    metavar=>"all all_basic rmsd rmsf rgyr sasa hbond cluster hbond_pl contact_pl dipole"],
    ['--ligand_group',help => "name of ligand group for analysis", metavar=>"drug"],
    ['--rmsd_ref', default => "em.tpr", help => "reference structure for rmsd", metavar=>"tpr or structure file"],
    ['--rmsf_trj_length', default => "2.0", help => "length of trajectory to run RMSD on", metavar=>"[0-9].[0-9]"],
    ['--traj_prepare', required => 0, 
    help => "Do not run analysis. Instead, prepare trajectory with the given index group", metavar=>"protein"],
    ['-n', default => "index.ndx", help => "index file", metavar=>"file.ndx"],
    ['-v', type => 'Bool', help => "verbose"],
);

my $args = $parser->parse_args();
my $dir = $args->dir;
my $directories_ref = $args->directories;
my @directories=@$directories_ref;
my $deffnm = $args->deffnm;
my $what_ref = $args->what;
my @w_analysis = @$what_ref;
my $ligand_group = $args->ligand_group;
my $traj_group = $args->traj_prepare;
my $rmsd_ref = $args->rmsd_ref;
my $rmsf_trj_length = $args->rmsf_trj_length;
my $ndx_file = $args->n;
my $b_verbose= $args->v;

my %required_files = (
    "index" => $ndx_file,
    "rmsd reference" => $rmsd_ref
);

if ((grep (/^all$/,@w_analysis) || grep (/hbond_pl/,@w_analysis) || grep (/contact_pl/,@w_analysis)) && 
    ! $ligand_group=~/\w+/)
{ 
    die "Protein-ligand analysis requested but ligand group is not given.\n".
    "Set ligand group with --ligand_group\n";
}

if ($b_verbose) {
    print "$0 will be run with the following arguments\n";
    if (defined $directories[0]) {
        print "--directories ";
        print join(" ",@directories),"\n";
    }
    else {
        print "--dir $dir\n";
    }
    print "--deffnm $deffnm\n";
    if (defined $traj_group) {
        $traj_group=lc($traj_group);
        print "--traj_prepare $traj_group\n";
    }
    else {
        print "--what ";
        print join(" ",@w_analysis),"\n";
        # Do we need to use the ligand group?
        my @ar = grep (/^all$/,@w_analysis);
        push (@ar, grep (/hbond_pl/,@w_analysis));
        push (@ar, grep (/contact_pl/,@w_analysis));
        if (scalar @ar > 0 && defined ($ligand_group)) {
            # @ar has values, which means that either "all" or "hbond_pl" and/or "contact_pl" is set
            print "--ligand_group $ligand_group\n";
        }   
    }
    print "--rmsd_ref $rmsd_ref\n";
    print "--rmsf_trj_length $rmsf_trj_length\n";
    print "--ndx_file $ndx_file\n";
}

# check if a generic directory name is used, in which case put all directories
# in @directories
if (! defined($directories[0])) {
    # use the value stored in $dir
    opendir ((my $dh), $cwd) || die "Cannot open dir $cwd: $!\n";
    my @match = grep {/^$dir\d+$/} readdir($dh);
    close($dh);
    push @directories,@match;
}

# check if it is possible to open all directories
die "No directory with name $dir\[0-9\]\* is found\n" unless (defined($directories[0]));

foreach my $dir (@directories) {
    opendir ((my $dh), "$cwd/$dir") || die "Cannot open dir $dir: $!\n";
    close ($dh);
    # check that the directory includes the relevant files
    foreach my $filetype (".xtc",".tpr") {
        die "File $deffnm$filetype not found in directory $dir" unless (-e "$cwd/$dir/$deffnm$filetype");
    }
    unless (-d "$cwd/$dir/Analysis") {
        mkdir ("$cwd/$dir/Analysis",0755) || die "Cannot generate directory $cwd/$dir/Analysis for Analysis:$!";
    }
}

# check if necessary files exist
foreach my $file (keys %required_files) {
    unless (-e $required_files{$file}) {
        die "$file file $required_files{$file} is missing\n";
    } 
}

if (defined $traj_group) {
    my $ndx_grp = get_index_group_number ($traj_group,$ndx_file);
    die "Cannot get index group for $traj_group in file $ndx_file\n" unless (defined $ndx_grp); 
    run("trjconv",$ndx_grp,$traj_group);
    exit(0);
}

if ( (grep /^all/, @w_analysis) || (grep /^rmsd$/, @w_analysis)) {
    print "Will run RMSD\n";
    run("rmsd");
}
if ( (grep /^all/, @w_analysis) || (grep /^rmsf$/, @w_analysis)) {
    print "Will run RMSF\n";
    run("rmsf");
}
if ( (grep /^all/, @w_analysis) || (grep /^rgyr$/, @w_analysis)) {
    print "Will run gmx gyrate\n";
    run("rgyr");
}
if ( (grep /^all/, @w_analysis) || (grep /^sasa$/, @w_analysis)) {
    print "Will run sasa analysis\n";
    run("sasa");
}
if ( (grep /^all/, @w_analysis) || (grep /^hbond$/, @w_analysis)) {
    print "Will run hydrogen bond analysis\n";
    run("hbond");
}
if ( (grep /^all/, @w_analysis) || (grep /^cluster$/, @w_analysis)) {
    print "Will run cluster analysis\n";
    run("cluster");
}
if ( (grep /^all$/, @w_analysis) || (grep /^hbond_pl$/, @w_analysis) || (grep /^contact_pl$/, @w_analysis)) {
    # need to get the index group number for the ligand
    # if the user gave a number, we just use it.
    # otherwise, we need to find the index group based on the name
    unless ($ligand_group =~ /^\d+$/) {
        my $ndx_grp = get_index_group_number($ligand_group,$ndx_file);
        die "Cannot get index group for $ligand_group in file $ndx_file\n" unless (defined $ndx_grp); 
        $ligand_group = $ndx_grp;
    }  
}
if ( (grep /^all$/, @w_analysis) || (grep /^hbond_pl$/, @w_analysis)) {
    print "Will run protein-ligand hydrogen bond analysis\n";
    run("hbond_pl");
}
if ( (grep /^all$/, @w_analysis) || (grep /^contact_pl$/, @w_analysis)) {
    print "Will run protein-ligand contact analysis\n";
    run("contact_pl");
}
if ( (grep /^all$/, @w_analysis) || (grep /^dipole$/, @w_analysis)) {
    print "Will calculate the protein dipole\n";
    run("dipole");
}


# subs

sub get_index_group_number {
    my ($group,$ndxfile)=@_;
    my $cmd ="gmx check -n $ndxfile 2>& 1";
    my @output=`$cmd `;
    foreach $_ (@output) {
        # example line that would match below
        #    0  System                 70152       1   70152
        if (/\s*(\d+)\s+(\w+)\s+\d+\s+\d+\s+\d+/) {
            if ($group eq lc($2)) {
                return $1;
            }
        }
    }
    return undef; 
}

sub run {
    my $what = $_[0];
    my $cmd;
    my $newtraj;
    my $max_rmsd;
    my @file_list;

    if ($what eq "trjconv") {
        die "cannot run trjconv\n" unless (defined $_[1] && defined $_[2]);
        my $group_number = $_[1];
        my $group_name = $_[2];
        $newtraj=$deffnm."_".$group_name;
        $cmd = "echo $group_number | gmx trjconv -f $deffnm -s $deffnm -o  $newtraj -n $cwd/$ndx_file > /dev/null 2>&1 ";
    }
    elsif ($what eq "rmsd") {
        $cmd = 'echo "3 3" | gmx rms '. "-f $deffnm.xtc -s $cwd/$rmsd_ref  -n $cwd/$ndx_file -tu ns -o Analysis/rmsd_ca > /dev/null 2>&1";
    }
    elsif ($what eq "rmsf") {
        print "RMSF per residue will be runs in chunks of $rmsf_trj_length ns\n";
        print "Only non-hydrogen atoms are considered.\n";
        print "Output will be written to file rmsf_all.dat\n";
    }
    elsif ($what eq "rgyr") {
        $cmd = "echo 1 | gmx gyrate -f $deffnm -s $deffnm -p -o Analysis/gyrate > /dev/null 2>&1";
    }
    elsif ($what eq "sasa") {
        $cmd = " echo 1 | gmx sasa -f $deffnm -s $deffnm -tu ns -o Analysis/area > /dev/null 2>&1";
    }
    elsif ($what eq "hbond") {
        $cmd = "echo \"1 1\" | gmx hbond -f $deffnm -s $deffnm -tu ns -num Analysis/hbnum > /dev/null 2>&1";
    }
    elsif ($what eq "cluster") {
        $cmd = "echo \"3 1\" | gmx cluster -f md_protein -s md_protein -cutoff 0.15 -method gromos ".
        " -g Analysis/cluster -sz Analysis/clust-size -ntr Analysis/clust-trans -clid Analysis/clust-id.xvg ".
        " -cl Analysis/clusters.xtc -o Analysis/rmsd-clust -dist Analysis/rmsd-dist > /dev/null 2>&1";
    }
    elsif ($what eq "hbond_pl") {
        $cmd = "echo \"1 $ligand_group\" | gmx hbond -f $deffnm -s $deffnm -tu ns -n $cwd/$ndx_file -num ". 
        "Analysis/hbnum_protein_ligand > /dev/null 2>&1";
    }
    elsif ($what eq "contact_pl") {
        $cmd = "echo \"1 $ligand_group\" | gmx mindist -f $deffnm -s $deffnm -tu ns -n $cwd/$ndx_file ".
        " -od Analysis/mindist_protein_ligand -on Analysis/numcont_protein_ligand -or Analysis/mindistres_protein-ligand ".
        "-d 0.4 > /dev/null 2>&1";
    }
    elsif ($what eq "dipole") {
        $cmd = "echo 1 | gmx dipoles -f $deffnm -s $deffnm -o Analysis/mtot -eps Analysis/epsilon ".
        "-a Analysis/dipole_avg -d Analysis/dipdist > Analysis/gmx_dipoles.log 2>&1"; 
    }
    else {
        print "Cannot run $what\n";
        return;
    }
    foreach my $dir (@directories) {
        chdir "$cwd/$dir";

        if ($what eq "rmsf") {
            # Per residue RMSF is not run on the whole trajectory but in chunks of $rmsf_trj_length ns
            # we should first extract the length of the trajectory
            my ($trj_start,$trj_end)=trajectory_get_start_end($deffnm);
            for (my $t=$trj_start; $t<=$trj_end-$rmsf_trj_length;$t+=$rmsf_trj_length ) {
                # run rmsf from $t to $t+=$rmsf_trj_length
                # frames at start/end are coming twice, but if rmsf_trj_length is long enough this won't matter
                my $begin=$t;
                my $end=$t+$rmsf_trj_length;
                my $rmsf_outfile="Analysis/rmsf_$begin"."_".$end;
                $begin*=1000; # gmx rmsf wants ps, not ns
                $end*=1000; # same here
                $cmd = "echo 2 | gmx rmsf -f $deffnm -s $deffnm -res -b $begin -e $end -o $rmsf_outfile > /dev/null 2>&1";
                unless (system($cmd) == 0) {
                    print "Could not run $cmd\n";
                    return;
                }
                push(@file_list,"$cwd/$dir/$rmsf_outfile.xvg");
            }
        }
        else {
            unless (system($cmd) == 0) {
                print "Could not run $cmd\n";
                return;
            }
        }

        if ($what eq "trjconv") {
            unless (-e "$newtraj.tpr") {
                print "Linking $deffnm.tpr to $newtraj.tpr\n";
                $cmd = "ln -s $deffnm.tpr $newtraj.tpr";
                system($cmd);
            } 
        }

        if ($what eq "rmsd") {
            $max_rmsd=-1E32;
            my @rmsd_vals=xvg2ar("Analysis/rmsd_ca.xvg");
            $max_rmsd = $max_rmsd < max(@rmsd_vals) ? max(@rmsd_vals) : $max_rmsd;
        }

        if ($what eq "rgyr") {
            push (@file_list,"$cwd/$dir/Analysis/gyrate.xvg");
        }

        if ($what eq "sasa") {
            push (@file_list,"$cwd/$dir/Analysis/area.xvg");
        }

        if ($what eq "hbond") {
            push (@file_list,"$cwd/$dir/Analysis/hbnum.xvg");
        }

        if ($what eq "cluster") {
            push (@file_list,"$cwd/$dir/Analysis/clusters.xtc");
        }

        if ($what eq "hbond_pl") {
            push (@file_list,"$cwd/$dir/Analysis/hbnum_protein_ligand.xvg");
        }

        if ($what eq "contact_pl") {
            push (@file_list,"$cwd/$dir/Analysis/numcont_protein_ligand.xvg");
        }

        if ($what eq "dipole") {
            push (@file_list,"$cwd/$dir/Analysis/mtot.xvg");
        }

        chdir $cwd;
    }

    # output / further analysis
    if ($what eq "rmsd") {
        printf "The maximal RMSD was %.2f nm\n",$max_rmsd;
    }

    if ($what eq "rmsf") {
        my $output;
        my @residues = xvg2ar($file_list[0],1);
        my $nres=scalar @residues;
        my @all_values;
        my @rmsf_values;
        foreach my $rmsf_file(@file_list) {
            @rmsf_values=xvg2ar($rmsf_file,2);
            push(@all_values,@rmsf_values);
        }
        for (my $i=0;$i<$nres;$i++) {
            my @ar;
            # push all values that correspond to residue $residues[i] to @ar
            for (my $j=$i;$j<scalar @all_values;$j+=$nres) {
                push(@ar,$all_values[$j]);
            }
            my $avg=mean(@ar);
            my $stddev=stddev(@ar);
            $output.=sprintf("%10d %12.4f %12.4f %12.4f\n",$residues[$i],$avg-$stddev,$avg,$avg+$stddev);
        }
        unless (open (RMSF_ALL,">rmsf_all.dat")) {
            print "Error - cannot open file rmsf_all.dat: $!\n";
            return;
        }
        print RMSF_ALL $output;
        close (RMSF_ALL);
    }

    if ($what eq "rgyr") {
        my (@rgyr,@pc1,@pc2,@pc3);
        foreach my $file(@file_list) {
            push @rgyr,xvg2ar($file,2);
            push @pc1,xvg2ar($file,3);
            push @pc2,xvg2ar($file,4);
            push @pc3,xvg2ar($file,5);
        }
        my $rg=mean(@rgyr);
        my $pc1=mean(@pc1);
        my $pc2=mean(@pc2);
        my $pc3=mean(@pc3);
        printf "Stats for Rgyr below, mean +/- standard errors\n";
        printf "Gyration radius %.2f +/- %.2f PC1 %.2f +/- %.2f PC2 %.2f +/- %.2f PC3 %.2f +/- %.2f\n",$rg,stddev(@rgyr)/sqrt(scalar @rgyr),
        $pc1,stddev(@pc1)/sqrt(scalar @pc1),$pc2,stddev(@pc2)/sqrt(scalar @pc2),$pc3,stddev(@pc3)/sqrt(scalar @pc3);
        my $b=1.5*$pc3*$pc3 - 0.5*$rg*$rg;
        my $c=$pc2*$pc2-$pc1*$pc1;
        my $k2=($b*$b-0.75*$c*$c)/($rg**4);
        printf ("Asphericity %.2f acylindricity %.2f relative shape anisotropy %.2f\n",$b,$c,$k2);
    }

    if ($what eq "sasa") {
        my @all_values;
        foreach my $file(@file_list) {
            push @all_values, xvg2ar($file,2);
        }
        printf "Statistics for SASA\n";
        my $stat = new Statistics::PointEstimation;
        $stat->set_significance(95); #set the significance(confidence) level to 95%
        $stat->add_data(@all_values);
        my ($mean, $std_error, $lower_clm, $upper_clm)=($stat->mean(), $stat->standard_error(), $stat->lower_clm(), $stat->upper_clm());
        printf ("SASA %.2f +/- %.2f CI (0.95) %.2f-%.2f \n", $mean, $std_error, $lower_clm, $upper_clm);
    }

    if ($what =~ /hbond/ || $what eq "contact_pl" || $what eq "dipole") {
        my @all_values;
        my $field;
        if ($what =~ /hbond/ || $what eq "contact_pl") {
            $field=2;
        }
        if ($what eq "dipole") {
            $field=5;
        }
        foreach my $file(@file_list) {
            push @all_values, xvg2ar($file,$field);
        }
        if ($what =~ /hbond/) {
            print "Statistics for H-bonds\n";
        }
        elsif ($what eq "contact_pl") {
            print "Statistics for the number of contacts\n";
        }
        elsif ($what eq "dipole") {
            print "Statistics for total dipole\n";
        }
        my $stat = new Statistics::PointEstimation;
        $stat->set_significance(95); #set the significance(confidence) level to 95%
        $stat->add_data(@all_values);
        my ($mean, $std_error, $lower_clm, $upper_clm)=($stat->mean(), $stat->standard_error(), $stat->lower_clm(), $stat->upper_clm());
        printf ("N %.2f +/- %.2f CI (0.95) %.2f-%.2f \n", $mean, $std_error, $lower_clm, $upper_clm);
    }

    if ($what eq "cluster") {
        #concatantae cluster files
        $cmd = "gmx trjcat -f @file_list -o all_clust > /dev/null 2>&1";
        unless (system($cmd) == 0) {
            print "Could not run $cmd\n";
            return;
        }
        $cmd = "echo \"3 1\" | gmx cluster -f all_clust -s $directories[0]/$deffnm -cutoff 0.2 -method gromos -sz > /dev/null 2>&1";
        unless (system($cmd) == 0) {
            print "Could not run $cmd\n";
            return;
        }
        my @clusters=xvg2ar("clust-size.xvg",2);
        print "Number of clusters: ",scalar(@clusters),"\n";
        my $nc=0;
        foreach my $clust_size (@clusters) {
            if ($clust_size>1) {
                $nc++;
            }
        }
        print "Number of clusters with >1 members ",$nc,"\n";
    }

}

# xvg analysis

sub xvg2ar {
    # Use: xvg2ar(file.xvg,n);
    # returns an array with all values in field n
    my ($xvg,$field)=@_;
    unless (-e $xvg) {
        print "File $xvg is not found in xvg2ar. Expect problems!\n";
        return;
    }
    my @lines;
    my $line;
    my @ar;
    my @return;

    $field--;
    open(XVG,$xvg);
    @lines = <XVG>;
    close(XVG);
    while ($line=shift(@lines)) {
	    next if ($line=~/[\@\#\&]/);
	    @ar=split(/\s+/,$line);
	    if ($ar[0] eq '') {
	        shift (@ar); # get rid of a leading empty space
	    }
	    push(@return,$ar[$field]);
    }
    return @return;
}

# trajectories
sub trajectory_get_start_end {
    # Input: trajectory file
    # Output: start and end, in nanoseconds
    my $trjfile=shift;
    my $tmpfile="tmp".int(rand(10000));
    my $cmd = "gmx check -f $trjfile > $tmpfile 2>&1";
    system($cmd);
    open(TMPFILE,$tmpfile);
    my @lines= <TMPFILE>;
    close(TMPFILE);
    @lines = grep (/Reading/,@lines);
    unlink $tmpfile;
    # @lines should now have two lines, the first is for the first time and the last is for the last time
    # e.g.
    # Reading frame       0 time    0.000
    # Last frame       1000 time 10000.000
    $_=$lines[0];
    my @ar=split;
    my $firsttime=$ar[-1]/1000; # ps → ns
    $_=$lines[-1];
    @ar=split;
    my $lasttime=$ar[-1]/1000;
    return($firsttime,$lasttime);
}