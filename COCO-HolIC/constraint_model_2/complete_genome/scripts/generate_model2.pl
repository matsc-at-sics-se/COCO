#!/usr/bin/perl
## will generate the CLP eclispe program for constraint model 2 based on a given 
## interaction matrix
##
## argument 1: the interaction matrix
## argument 2: the value to scale the interaction frequency by (in order to convert it to
## an integer
##
## Kimberly MacKay June 8, 2016
## updated: July 28, 2017; removal of residual and redundant code
##
## license: This work is licensed under the Creative Commons Attribution-NonCommercial-
## ShareAlike 3.0 Unported License. To view a copy of this license, visit 
## http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 
## PO Box 1866, Mountain View, CA 94042, USA.

use strict;
use warnings;
use List::MoreUtils 'uniq';

## check to ensure two arguments was passed in
die "ERROR: must pass in two argumnets." if @ARGV != 2;

my $HiC_file = $ARGV[0];
my $scale = $ARGV[1]; # for s.pombe work this was set to 1000

###########################################################################################
##	parse the interaction matrix
###########################################################################################

## open the interaction matrix file
open HIC, "$HiC_file" or die "ERROR: $HiC_file could not be opened.";
chomp(my @matrix_file = <HIC>);
close HIC;

## $#matrix_file will be the size of the whole-genome contact map since the file also 
## contains a header line
my $num_variables = $#matrix_file;

my @frequencies;

## for each line after the header line
for(my $row = 1; $row <= $#matrix_file; $row++)
{
	## split the line
	my @matrix_line = split /\t/, $matrix_file[$row];
	
	## loop through the entire file to extract the frequencies
	## note: we only have to extract one half of the matrix since it is symmetric
	## along the diagonal
	for(my $col = $row; $col <= $num_variables; $col++)
	{
		## adjusts NA's to 0's
		if($matrix_line[$col] =~ "NA")
		{
			$frequencies[$row][$col] = 0;
		}
		else
		{
			# convert the frequency to an integer to improve speed
			$frequencies[$row][$col] = int($matrix_line[$col]*$scale);
		}
	}
}

###########################################################################################
##	printing the program
###########################################################################################

print "% Load the relevant libraries \n";
print ":- lib(gfd). \n\n";

print "% alldifferent_except(Vars) is true if  each term in Vars is \n";
print "% pairwise different from every other term, or has a value \n";
print "% of zero \n";
print "alldifferent_except(Vars) :- \n";
	print "\t% the list Vars has length N \n";
	print "\tlength(Vars, N), \n\n";

	print "\t% for each pair of terms in Vars, check if they are \n";
	print "\t% different or zero \n";
	print "\t( for(I,1,N), param(Vars, N) do \n";
		print "\t\t( for(J, I, N), param(Vars, I) do \n\n";
			
			print "\t\t\t% The variable I is the list position of the element X \n";
			print "\t\t\t% in Vars \n";
			print "\t\t\telement(I, Vars, X), \n\n";
			
			print "\t\t\t% The variable J is the list position of the element Y \n";
			print "\t\t\t% in Vars \n"; 
			print "\t\t\telement(J, Vars, Y), \n\n";

			print "\t\t\t% If I and J do not correspond to the same position in \n";
			print "\t\t\t% Vars \n";
			print "\t\t\t(I #\\= J) -> \n";
				print "\t\t\t\t% constrain the elements X and Y to be pairwise \n";
				print "\t\t\t\t% different or have one (or both) equal zero \n";
				print "\t\t\t\tX #\\= Y or (X #= 0 or Y #= 0) \n";
			print "\t\t\t; \n";
				print "\t\t\t\ttrue \n";
			print "\t\t\t) \n";
		print "\t\t). \n\n";

print "% maximize(Rows, Freqs)  is true if the elements in Row are all \n";
print "% different or zero, the elements in Freqs are bound to zero \n";
print "% or the associated rounded and scaled integer value (based on \n";
print "% the interaction frequency from the whole-genome contact \n";
print "% map), and the elements in Freqs represent the subset of \n";
print "% rounded and scaled interaction frequencies that have maximum \n";
print "% sum \n";
print "maximize(Rows, Freqs) :- \n\n";

	print "\t% Variable Declarations: \n";
	print "\t% The list Rows has one variable for each row of the \n";
	print "\t% whole-genome contact map \n";
	print "\tRows = [";
	for(my $i = 1; $i <= $num_variables; $i++)
	{
		# account for the last variable, it won't have a comma
		if($i == $num_variables)
		{
			print "Row".$i."],\n\n";
		}
		else
		{
			print "Row".$i.", ";
		}
	}
	
	print "\t% The list Freqs has one variable for each row of the \n";
	print "\t% whole-genome contact map	\n";
	print "\tFreqs = [";
	for(my $i = 1; $i <= $num_variables; $i++)
	{
		# account for the last variable, it won't have a comma
		if($i == $num_variables)
		{
			print "Freq".$i."],\n\n";
		}
		else
		{
			print "Freq".$i.", ";
		}
	}

	print "\t% Representation of the Genome: \n";
	print "\t% Each Row term can assume a value based on interacting bin \n";
	print "\t% indices where `0' represents an interaction not being \n";
	print "\t% selected and a non-zero value (ranging from 1 to N) \n";
	print "\t% represents which genomic bin is involved in the selected \n";
	print "\t% interaction \n";
	for(my $i = 1; $i <= $num_variables; $i++)
	{
		# account for the last number, it won't have a ..
		if($i == $num_variables)
		{
			print "\tRow".$i." :: [0, ".$num_variables."], \n\n";
		}
		else
		{
			print "\tRow".$i." :: [0, ".$i."..".$num_variables."], \n";
		}
	}
	
	print "\t% Each frequency term can assume either the rounded \n";
	print "\t% and scaled integer value (based on the corresponding \n";
	print "\t% interaction frequency from the whole-genome contact \n";
	print "\t% map) or a value of `0' where `0'  represents an \n";
	print "\t% interaction not being selected \n";
	for(my $i = 1; $i <= $num_variables; $i++)
	{
	
		print "\tFreq".$i." :: [0, ";
	
		# note we only have to loop through half of the matrix since it is symmetric
		for(my $j = $i; $j <= $num_variables; $j++)
		{
			# account for the last number, it won't have a comma
			if($j == $num_variables)
			{
				print $frequencies[$i][$j]."], \n";
			}
			else
			{
				print $frequencies[$i][$j].", ";
			}
		}
	}

	print "\n\t% Constraints: \n";
	print "\t% Each pair of corresponding (Row<i>, Freq<i>) variables \n";
	print "\t% must assume dependent values based on data from the \n";
	print "\t% whole-genome contact map; A (Row, Freq) pair ground to \n";
	print "\t% (0,0) encodes that nothing is chosen \n";
	for(my $i = 1; $i <= $num_variables; $i++)
	{
		print "\n\t% Row ".$i." Constraints \n";
		# print the "no" case
		print "\t((Row".$i." #= 0) and (Freq".$i." #= 0)) or \n";
		
		# note we only have to loop through half of the matrix since it is symmetric
		for(my $j = $i; $j <= $num_variables; $j++)
		{
			# account for the last number, it won't have an or
			if($j == $num_variables)
			{
				print "\t((Row".$i." #= ".$j.") and (Freq".$i." #=  ".$frequencies[$i][$j].")), \n"; 
			}
			else
			{
				print "\t((Row".$i." #= ".$j.") and (Freq".$i." #=  ".$frequencies[$i][$j].")) or\n";
			}
		}
	}

	print "\n\t% All of the values assumed by the Row<i> variables must be \n";
	print "\t% all different or zero to ensure each genomic bin is \n";
	print "\t% involved in only 1 interaction; multiple zeros are allowed \n";
	print "\talldifferent_except(Rows), \n\n";

	print "\t% Optimize: the sum of the selected interaction \n";
	print "\t% frequencies is maximal (it is the minimum of \n";
	print "\t% the additive inverse of the sum) \n";
	print "\tCost #= -sum(Freqs), \n";
	print "\tsearch(Freqs, 0, input_order, indomain_max, bb_min(Cost), []). \n\n";
