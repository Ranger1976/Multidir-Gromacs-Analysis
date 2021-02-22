# Multidir-Gromacs-Analysis

This tool is used for analysis of multiple protein MD runs of a similar system that are stored in multiple directories. For example, directories sim0 sim1 and sim2 each contain an MD simulation of the same system, which you want to analyse simoltaneously.

# How to run 

In the following, it is assumed that you have:
* directories named sim0, sim1, sim2 etc. each containing a file called md.xtc which is a trjectory file, and a .tpr file called md.tpr
* an index file located in index.ndx, located in the directory from which you run the script
* a reference structure file for RMSD calculations em.tpr, located in the directory from which you run the script

In a typical run you first want to convert the trajectories so that you only have the protein. 
./analyse_multidir.pl --traj_prepare protein --deffnm md --dir sim -n index.ndx 

Then, you want to run complete analysis on the protein simulations:
./analyse_multidir.pl --deffnm md --dir sim -n index.ndx --rmsd_ref em.tpr --what all_basic

This will, if all works well, produce some output and generate .xvg files under sim*/Analysis

./analyse_multidir.pl -h 
or 
perl ./analyse_multidir.pl -h 

to see a list of options and get some info.

# Requirements to run this tool
* Multiple directories, each containing a .xvg and .tpr file of the same name (e.g., md.xtc and md.tpr), same length and corresponding to the same system
* An index file on the top level directory (need to have the index group relevant to run trjconv - usally protein, CA as group 3 and a ligand group if relevant)
* A reference file for RMSD calculations, on the top directory, such as em.tpr

# Requirement for the tool to execute in the first place
* Perl installed (almost always available if you use Linux)
and the following Perl packages:
* Getopt::ArgParse
* Cwd
* List::Util
* Statistics::Basic
* Statistics::PointEstimation

The easiest way to get these is from CPAN. Usually, you'll have CPAN installed as well, but check your distro. For installation of tools see https://www.linuxcloudvps.com/blog/how-to-install-perl-modules-using-cpan-on-linux/ or just google on how to install packages with CPAN.

