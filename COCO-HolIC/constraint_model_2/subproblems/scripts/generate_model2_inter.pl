#!/usr/bin/perl
## will generate the CLP eclispe sub-program for inter-interactions based on constraint 
## model 2 and a given interaction matrix when using the divide-and-conquer approach
##
## argument 1: the interaction matrix
## argument 2: the value to scale the interaction frequency by (in order to convert it to
## an integer
## argument 3: the first chromosome fo interest (refered to as chr1 throughout the rest of the code)
## argument 4: chr1 start index
## argument 5: chr1 stop index
## argument 6: the second chromosome fo interest (refered to as chr2 throughout the rest of the code)
## argument 7: chr2 start index
## argument 8: chr2 stop index
## argument 9: size of chr2 (number of bins)
##
## Kimberly MacKay June 13, 2016
## license: This work is licensed under the Creative Commons Attribution-NonCommercial-
## ShareAlike 3.0 Unported License. To view a copy of this license, visit 
## http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 
## PO Box 1866, Mountain View, CA 94042, USA.

use strict;
use warnings;
use List::MoreUtils 'uniq';

## check to ensure nine arguments was passed in
die "ERROR: must pass in nine argumnets." if @ARGV != 9;

my $HiC_file = $ARGV[0];
my $scale = $ARGV[1]; # for s.pombe work this was set to 1000

## info pertaining to the first chromsome of interest
my $chr1 = $ARGV[2];
my $chr1_start = $ARGV[3];
my $chr1_stop = $ARGV[4];

## infor pertaining to the second chromosome of interest
my $chr2 = $ARGV[5];
my $chr2_start = $ARGV[6];
my $chr2_stop = $ARGV[7];
my $chr2_size = $ARGV[8];


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
	for(my $col = 1; $col <= $row; $col++)
	{
		## adjusts NA's to 0's
		if($matrix_line[$col] =~ "NA")
		{
			$frequencies[$row][$col] = 0;
		}
		else
		{
			## trying a smaller integer size to see if it improves speed
			$frequencies[$row][$col] = int($matrix_line[$col]*$scale);
		}	
	}
}

###########################################################################################
##	extract the non-zero rows
###########################################################################################

## loop through the submatrix to determine the non-zero rows
my %non_zero;
my @non_zero_row_index;
my %only_zero_rows;

for(my $row = $chr2_start; $row <= $chr2_stop; $row++)
{
	my $only_zero = 1;
	
	for(my $col = $chr1_start; $col <= $chr1_stop; $col++)
	{
		##check to see if it is non-zero
		if($frequencies[$row][$col] != 0)
		{
			$only_zero = 0;
			$non_zero{$row}{$col} = $col;
		}
	}
	
	if($only_zero)
	{
		$only_zero_rows{$row} = $row;
	}
}

@non_zero_row_index = (sort {$a <=> $b} keys %non_zero);

## print out the non-zero rows
my $out_filename = "../chr".$chr1."_chr".$chr2."_non_zero_bins.txt";
open OUT, '>', "$out_filename" or die "ERROR: $out_filename could not be opened.";

for(my $i = 0; $i <= $#non_zero_row_index; $i++)
{
	print OUT "$non_zero_row_index[$i]\n";
}

close OUT;


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
		
print "% maximize(RowFile, FreqFile, Non_Zero_Rows)  is true if the elements in Row are all \n";
print "% different or zero, the elements in Freqs are bound to zero \n";
print "% or the associated rounded and scaled integer value (based on \n";
print "% the interaction frequency from the whole-genome contact \n";
print "% map), and the elements in Freqs represent the subset of \n";
print "% rounded and scaled interaction frequencies that have maximum \n";
print "% sum. RowFile will be the name of the selected columns for each \n";
print "% row in the solution set. FreqFile will be the corresponding interaction \n";
print "% frequencies for the row,column pairs in the solution set. Non_Zero_Rows is an \n";
print "% output variable \n";
print "maximize(RowFile, FreqFile, Non_Zero_Rows) :-\n\n";

	print "\t% Variable Declarations: \n";
	print "\t% The list Non_Zero_Rows has one variable for each non-zero row of the \n";
	print "\t% whole-genome contact map \n";
	print "\tNon_Zero_Rows = [";
	for(my $i = 0; $i <= $#non_zero_row_index; $i++)
	{
		## account for the last variable, it won't have a comma
		if($i == $#non_zero_row_index)
		{
			print "Chr".$chr1."_Chr".$chr2."_Row".$non_zero_row_index[$i]."],\n\n";
		}
		else
		{
			print "Chr".$chr1."_Chr".$chr2."_Row".$non_zero_row_index[$i].", ";
		}
	}
	
	print "\t% The list Freqs has one variable for each non-zero row of the \n";
	print "\t% whole-genome contact map	\n";
	print "\tFreqs = [";
	for(my $i = 0; $i <= $#non_zero_row_index; $i++)
	{
		## account for the last variable, it won't have a comma
		if($i == $#non_zero_row_index)
		{
			print "Chr".$chr1."_Chr".$chr2."_Freq".$non_zero_row_index[$i]."],\n\n";
		}
		else
		{
			print "Chr".$chr1."_Chr".$chr2."_Freq".$non_zero_row_index[$i].", ";
		}
	}
	
	print "\t% Representation of the Genome: \n";
	print "\t% Each Row term can assume a value based on interacting bin \n";
	print "\t% indices where `0' represents an interaction not being \n";
	print "\t% selected and a non-zero value (ranging from 1 to N) \n";
	print "\t% represents which genomic bin is involved in the selected \n";
	print "\t% interaction \n";
	for(my $i = 0; $i <= $#non_zero_row_index; $i++)
	{
		my @domain;
		my @unique_values;
		
		for(my $col = $chr1_start; $col <= $chr1_stop; $col++)
		{
			##check to see if it is non-zero
			if($frequencies[$non_zero_row_index[$i]][$col] != 0)
			{
				push @domain, $col;
			}	
		}

		@unique_values = uniq @domain;

		print "\tChr".$chr1."_Chr".$chr2."_Row".$non_zero_row_index[$i]." :: [0";
		for(my $u = 0; $u <= $#unique_values; $u++)
		{
			print ", ".$unique_values[$u];	
		}	
		print "],\n";
	}
	
	print "\n\t% Each frequency term can assume either the rounded \n";
	print "\t% and scaled integer value (based on the corresponding \n";
	print "\t% interaction frequency from the whole-genome contact \n";
	print "\t% map) or a value of `0' where `0'  represents an \n";
	print "\t% interaction not being selected \n";
	for(my $i = 0; $i <= $#non_zero_row_index; $i++)
	{
		my @domain;
		my @unique_values;

		for(my $col = $chr1_start; $col <= $chr1_stop; $col++)
		{
			##check to see if it is non-zero
			if($frequencies[$non_zero_row_index[$i]][$col] != 0)
			{
				push @domain, $frequencies[$non_zero_row_index[$i]][$col];
			}	
		}
		
		@unique_values = uniq @domain;
	
		print "\tChr".$chr1."_Chr".$chr2."_Freq".$non_zero_row_index[$i]." :: [0";
		for(my $u = 0; $u <= $#unique_values; $u++)
		{
			print ", ".$unique_values[$u];	
		}	
		print "],\n";	
	}

	print "\n\t% Constraints: \n";
	print "\t% Each pair of corresponding (Row<i>, Freq<i>) variables \n";
	print "\t% must assume dependent values based on data from the \n";
	print "\t% whole-genome contact map; A (Row, Freq) pair ground to \n";
	print "\t% (0,0) encodes that nothing is chosen \n";
	for(my $row = $chr2_start; $row <= $chr2_stop; $row++)
	{		
		## if non-zero values exist
		if(defined $non_zero{$row})
		{
			## first print the non-zero options
			for my $nz_col (sort keys $non_zero{$row})
			{	
				print "\t((Chr".$chr1."_Chr".$chr2."_Row".$row." #= ".$non_zero{$row}{$nz_col}.") and (Chr".$chr1."_Chr".$chr2."_Freq".$row." #= ".$frequencies[$row][$non_zero{$row}{$nz_col}].")) or \n";			
			}			
			# print the "0" case
			print "\t((Chr".$chr1."_Chr".$chr2."_Freq".$row." #= 0) and (Chr".$chr1."_Chr".$chr2."_Row".$row." #= 0)),\n\n";
		}
	}		

	print "\n\t% All of the values assumed by the Row<i> variables must be \n";
	print "\t% all different or zero to ensure each genomic bin is \n";
	print "\t% involved in only 1 interaction; multiple zeros are allowed \n";
	print "\talldifferent_except(Non_Zero_Rows), \n\n";
	print "\tatmost(".($#non_zero_row_index-1).", Non_Zero_Rows, 0),\n";

	print "\t% Optimize: the sum of the selected interaction \n";
	print "\t% frequencies is maximal (it is the minimum of \n";
	print "\t% the additive inverse of the sum) \n";
	print "\tCost #= -sum(Freqs), \n";
	print "\tsearch(Freqs, 0, input_order, indomain_max, bb_min(Cost), []), \n\n";

	print "\t%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
	print "\t%% 	Output the results\n";
	print "\t%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";

	print "\topen(FreqFile, 'write', FREQ_OUT),\n";
	print "\t%%list the frequencies\n";
	print "\t(foreach(X,Freqs),\n";
	print "\t\tparam(FREQ_OUT) do\n";
	print "\t\t\tget_domain_as_list(X, DomList),\n";
	print "\t\t\t\t(foreach(Y,DomList),\n";
	print "\t\t\t\t\tparam(FREQ_OUT) do\n";
	print "\t\t\t\t\t\twrite(FREQ_OUT, Y),\n";
	print "\t\t\t\t\t\twrite(FREQ_OUT, ' ')\n";
	print "\t\t\t\t),\n";
	print "\t\t\twrite(FREQ_OUT, \"\\n\")\n";
	print "\t),\n";
	print "\tclose(FREQ_OUT),\n\n";
	
	print "\t%% list the potential rows\n";
	print "\topen(RowFile, 'write', ROW_OUT),\n";
	print "\t(foreach(X,Non_Zero_Rows),\n";
	print "\t\tparam(ROW_OUT) do\n";
	print "\t\t\tget_domain_as_list(X, DomList),\n";
	print "\t\t\t\t(foreach(Y,DomList),\n";
	print "\t\t\t\t\tparam(ROW_OUT) do\n";
	print "\t\t\t\t\t\twrite(ROW_OUT, Y),\n";
	print "\t\t\t\t\t\twrite(ROW_OUT, ' ')\n";
	print "\t\t\t\t),\n";
	print "\t\t\t\twrite(ROW_OUT, \"\\n\")\n";
	print "\t),\n";
	print "\tclose(ROW_OUT).\n";
	
