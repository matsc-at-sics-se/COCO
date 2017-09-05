% license: This work is licensed under the Creative Commons Attribution-NonCommercial-
% ShareAlike 3.0 Unported License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 
% PO Box 1866, Mountain View, CA 94042, USA.

% Load the relevant libraries 
:- lib(gfd).

% alldifferent_except(Vars) is true if  each term in Vars is 
% pairwise different from every other term, or has a value 
% of zero 
alldifferent_except(Vars) :-
	% the list Vars has length N 
	length(Vars, N),

	% for each pair of terms in Vars, check if they are 
	% different or zero 
	( for(I,1,N), param(Vars, N) do
		( for(J, I, N), param(Vars, I) do
			
			% The variable I is the list position of the element X 
			% in Vars 
			element(I, Vars, X),
			
			% The variable J is the list position of the element Y 
			% in Vars 
			element(J, Vars, Y),

			% If I and J do not correspond to the same position in 
			% Vars 
			(I #\= J) ->
				% constrain the elements X and Y to be pairwise 
				% different or have one (or both) equal zero 
				X #\= Y or (X #= 0 or Y #= 0)
			;
				true
			)
		).

% maximize(RowFile, FreqFile, Non_Zero_Rows)  is true if the elements in Row are all 
% different or zero, the elements in Freqs are bound to zero 
% or the associated rounded and scaled integer value (based on 
% the interaction frequency from the whole-genome contact 
% map), and the elements in Freqs represent the subset of 
% rounded and scaled interaction frequencies that have maximum 
% sum. RowFile will be the name of the selected columns for each 
% row in the solution set. FreqFile will be the corresponding interaction 
% frequencies for the row,column pairs in the solution set. Non_Zero_Rows is an 
% output variable 
maximize(RowFile, FreqFile, Non_Zero_Rows) :-

	% Variable Declarations: 
	% The list Non_Zero_Rows has one variable for each non-zero row of the 
	% whole-genome contact map 
	Non_Zero_Rows = [Chr1_Chr2_Row562, Chr1_Chr2_Row563, Chr1_Chr2_Row564, Chr1_Chr2_Row565, Chr1_Chr2_Row566, Chr1_Chr2_Row567, Chr1_Chr2_Row568, Chr1_Chr2_Row569, Chr1_Chr2_Row570, Chr1_Chr2_Row571, Chr1_Chr2_Row572, Chr1_Chr2_Row573, Chr1_Chr2_Row574, Chr1_Chr2_Row575, Chr1_Chr2_Row576, Chr1_Chr2_Row579, Chr1_Chr2_Row711, Chr1_Chr2_Row712, Chr1_Chr2_Row713, Chr1_Chr2_Row714, Chr1_Chr2_Row715, Chr1_Chr2_Row716, Chr1_Chr2_Row717, Chr1_Chr2_Row724, Chr1_Chr2_Row725, Chr1_Chr2_Row726, Chr1_Chr2_Row727, Chr1_Chr2_Row728, Chr1_Chr2_Row729, Chr1_Chr2_Row730, Chr1_Chr2_Row731, Chr1_Chr2_Row732, Chr1_Chr2_Row733, Chr1_Chr2_Row734, Chr1_Chr2_Row735, Chr1_Chr2_Row736, Chr1_Chr2_Row772, Chr1_Chr2_Row995, Chr1_Chr2_Row996, Chr1_Chr2_Row997, Chr1_Chr2_Row998, Chr1_Chr2_Row999, Chr1_Chr2_Row1000, Chr1_Chr2_Row1001, Chr1_Chr2_Row1002, Chr1_Chr2_Row1003, Chr1_Chr2_Row1004, Chr1_Chr2_Row1005, Chr1_Chr2_Row1006, Chr1_Chr2_Row1007, Chr1_Chr2_Row1008],

	% The list Freqs has one variable for each non-zero row of the 
	% whole-genome contact map	
	Freqs = [Chr1_Chr2_Freq562, Chr1_Chr2_Freq563, Chr1_Chr2_Freq564, Chr1_Chr2_Freq565, Chr1_Chr2_Freq566, Chr1_Chr2_Freq567, Chr1_Chr2_Freq568, Chr1_Chr2_Freq569, Chr1_Chr2_Freq570, Chr1_Chr2_Freq571, Chr1_Chr2_Freq572, Chr1_Chr2_Freq573, Chr1_Chr2_Freq574, Chr1_Chr2_Freq575, Chr1_Chr2_Freq576, Chr1_Chr2_Freq579, Chr1_Chr2_Freq711, Chr1_Chr2_Freq712, Chr1_Chr2_Freq713, Chr1_Chr2_Freq714, Chr1_Chr2_Freq715, Chr1_Chr2_Freq716, Chr1_Chr2_Freq717, Chr1_Chr2_Freq724, Chr1_Chr2_Freq725, Chr1_Chr2_Freq726, Chr1_Chr2_Freq727, Chr1_Chr2_Freq728, Chr1_Chr2_Freq729, Chr1_Chr2_Freq730, Chr1_Chr2_Freq731, Chr1_Chr2_Freq732, Chr1_Chr2_Freq733, Chr1_Chr2_Freq734, Chr1_Chr2_Freq735, Chr1_Chr2_Freq736, Chr1_Chr2_Freq772, Chr1_Chr2_Freq995, Chr1_Chr2_Freq996, Chr1_Chr2_Freq997, Chr1_Chr2_Freq998, Chr1_Chr2_Freq999, Chr1_Chr2_Freq1000, Chr1_Chr2_Freq1001, Chr1_Chr2_Freq1002, Chr1_Chr2_Freq1003, Chr1_Chr2_Freq1004, Chr1_Chr2_Freq1005, Chr1_Chr2_Freq1006, Chr1_Chr2_Freq1007, Chr1_Chr2_Freq1008],
	
	% Representation of the Genome: 
	% Each Row term can assume a value based on interacting bin 
	% indices where `0' represents an interaction not being 
	% selected and a non-zero value (ranging from 1 to N) 
	% represents which genomic bin is involved in the selected 
	% interaction 
	Chr1_Chr2_Row562 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 462, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557],
	Chr1_Chr2_Row563 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557],
	Chr1_Chr2_Row564 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557],
	Chr1_Chr2_Row565 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557],
	Chr1_Chr2_Row566 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557],
	Chr1_Chr2_Row567 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557],
	Chr1_Chr2_Row568 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557],
	Chr1_Chr2_Row569 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 543, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557],
	Chr1_Chr2_Row570 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557],
	Chr1_Chr2_Row571 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 557],
	Chr1_Chr2_Row572 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 557],
	Chr1_Chr2_Row573 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 546, 547, 548, 549, 550, 551, 552, 553, 555, 556, 557],
	Chr1_Chr2_Row574 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 549, 550, 551, 552, 553],
	Chr1_Chr2_Row575 :: [0, 6, 7, 8, 9, 10, 11, 12, 548, 550],
	Chr1_Chr2_Row576 :: [0, 11],
	Chr1_Chr2_Row579 :: [0, 6],
	Chr1_Chr2_Row711 :: [0, 370, 374],
	Chr1_Chr2_Row712 :: [0, 373, 374, 380, 381],
	Chr1_Chr2_Row713 :: [0, 373, 374, 381, 382, 383, 384],
	Chr1_Chr2_Row714 :: [0, 372, 373, 374, 380, 381, 384, 385, 386],
	Chr1_Chr2_Row715 :: [0, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 385, 386],
	Chr1_Chr2_Row716 :: [0, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387],
	Chr1_Chr2_Row717 :: [0, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386],
	Chr1_Chr2_Row724 :: [0, 364, 365, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390],
	Chr1_Chr2_Row725 :: [0, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390],
	Chr1_Chr2_Row726 :: [0, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390],
	Chr1_Chr2_Row727 :: [0, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391],
	Chr1_Chr2_Row728 :: [0, 367, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387],
	Chr1_Chr2_Row729 :: [0, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387],
	Chr1_Chr2_Row730 :: [0, 371, 372, 373, 374, 380, 381, 382, 383, 385, 386, 387],
	Chr1_Chr2_Row731 :: [0, 372, 373, 374, 380, 381, 382, 383, 385, 386],
	Chr1_Chr2_Row732 :: [0, 372, 373, 374, 381, 382],
	Chr1_Chr2_Row733 :: [0, 373, 381, 382, 385],
	Chr1_Chr2_Row734 :: [0, 381, 387, 388],
	Chr1_Chr2_Row735 :: [0, 373, 381],
	Chr1_Chr2_Row736 :: [0, 373],
	Chr1_Chr2_Row772 :: [0, 555, 556, 557],
	Chr1_Chr2_Row995 :: [0, 6, 7, 11, 557],
	Chr1_Chr2_Row996 :: [0, 6, 7, 8, 9, 10, 11],
	Chr1_Chr2_Row997 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14],
	Chr1_Chr2_Row998 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 549, 550, 556, 557],
	Chr1_Chr2_Row999 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557],
	Chr1_Chr2_Row1000 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557],
	Chr1_Chr2_Row1001 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 545, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557],
	Chr1_Chr2_Row1002 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557],
	Chr1_Chr2_Row1003 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 20, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557],
	Chr1_Chr2_Row1004 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 543, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557],
	Chr1_Chr2_Row1005 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557],
	Chr1_Chr2_Row1006 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557],
	Chr1_Chr2_Row1007 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557],
	Chr1_Chr2_Row1008 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557],

	% Each frequency term can assume either the rounded 
	% and scaled integer value (based on the corresponding 
	% interaction frequency from the whole-genome contact 
	% map) or a value of `0' where `0'  represents an 
	% interaction not being selected 
	Chr1_Chr2_Freq562 :: [0, 3, 4, 6, 5, 1, 2],
	Chr1_Chr2_Freq563 :: [0, 4, 5, 6, 3, 1, 2],
	Chr1_Chr2_Freq564 :: [0, 4, 5, 2, 1, 3],
	Chr1_Chr2_Freq565 :: [0, 4, 6, 5, 1, 2, 3],
	Chr1_Chr2_Freq566 :: [0, 3, 4, 5, 1, 2],
	Chr1_Chr2_Freq567 :: [0, 3, 4, 1, 2],
	Chr1_Chr2_Freq568 :: [0, 2, 1],
	Chr1_Chr2_Freq569 :: [0, 2, 1],
	Chr1_Chr2_Freq570 :: [0, 2, 3, 1],
	Chr1_Chr2_Freq571 :: [0, 2, 3, 1],
	Chr1_Chr2_Freq572 :: [0, 1, 2],
	Chr1_Chr2_Freq573 :: [0, 1, 2],
	Chr1_Chr2_Freq574 :: [0, 1],
	Chr1_Chr2_Freq575 :: [0, 1],
	Chr1_Chr2_Freq576 :: [0, 1],
	Chr1_Chr2_Freq579 :: [0, 1],
	Chr1_Chr2_Freq711 :: [0, 1],
	Chr1_Chr2_Freq712 :: [0, 1],
	Chr1_Chr2_Freq713 :: [0, 1],
	Chr1_Chr2_Freq714 :: [0, 1],
	Chr1_Chr2_Freq715 :: [0, 1],
	Chr1_Chr2_Freq716 :: [0, 1, 2],
	Chr1_Chr2_Freq717 :: [0, 1, 2],
	Chr1_Chr2_Freq724 :: [0, 1, 2, 3],
	Chr1_Chr2_Freq725 :: [0, 1, 2, 3],
	Chr1_Chr2_Freq726 :: [0, 1, 2, 3],
	Chr1_Chr2_Freq727 :: [0, 1, 2],
	Chr1_Chr2_Freq728 :: [0, 1],
	Chr1_Chr2_Freq729 :: [0, 1, 2],
	Chr1_Chr2_Freq730 :: [0, 1],
	Chr1_Chr2_Freq731 :: [0, 1],
	Chr1_Chr2_Freq732 :: [0, 1],
	Chr1_Chr2_Freq733 :: [0, 1],
	Chr1_Chr2_Freq734 :: [0, 1],
	Chr1_Chr2_Freq735 :: [0, 1],
	Chr1_Chr2_Freq736 :: [0, 1],
	Chr1_Chr2_Freq772 :: [0, 1],
	Chr1_Chr2_Freq995 :: [0, 1],
	Chr1_Chr2_Freq996 :: [0, 1],
	Chr1_Chr2_Freq997 :: [0, 1],
	Chr1_Chr2_Freq998 :: [0, 1],
	Chr1_Chr2_Freq999 :: [0, 2, 1],
	Chr1_Chr2_Freq1000 :: [0, 3, 2, 1],
	Chr1_Chr2_Freq1001 :: [0, 2, 1],
	Chr1_Chr2_Freq1002 :: [0, 3, 2, 1],
	Chr1_Chr2_Freq1003 :: [0, 3, 2, 1],
	Chr1_Chr2_Freq1004 :: [0, 4, 5, 6, 3, 2, 1],
	Chr1_Chr2_Freq1005 :: [0, 5, 6, 4, 2, 1, 3],
	Chr1_Chr2_Freq1006 :: [0, 5, 6, 7, 2, 1, 3, 4],
	Chr1_Chr2_Freq1007 :: [0, 5, 6, 7, 2, 1, 3, 4],
	Chr1_Chr2_Freq1008 :: [0, 5, 6, 8, 4, 2, 1, 3],

	% Constraints: 
	% Each pair of corresponding (Row<i>, Freq<i>) variables 
	% must assume dependent values based on data from the 
	% whole-genome contact map; A (Row, Freq) pair ground to 
	% (0,0) encodes that nothing is chosen 
	((Chr1_Chr2_Row562 #= 10) and (Chr1_Chr2_Freq562 #= 3)) or 
	((Chr1_Chr2_Row562 #= 11) and (Chr1_Chr2_Freq562 #= 1)) or 
	((Chr1_Chr2_Row562 #= 12) and (Chr1_Chr2_Freq562 #= 1)) or 
	((Chr1_Chr2_Row562 #= 13) and (Chr1_Chr2_Freq562 #= 1)) or 
	((Chr1_Chr2_Row562 #= 14) and (Chr1_Chr2_Freq562 #= 1)) or 
	((Chr1_Chr2_Row562 #= 15) and (Chr1_Chr2_Freq562 #= 1)) or 
	((Chr1_Chr2_Row562 #= 462) and (Chr1_Chr2_Freq562 #= 3)) or 
	((Chr1_Chr2_Row562 #= 545) and (Chr1_Chr2_Freq562 #= 1)) or 
	((Chr1_Chr2_Row562 #= 546) and (Chr1_Chr2_Freq562 #= 1)) or 
	((Chr1_Chr2_Row562 #= 547) and (Chr1_Chr2_Freq562 #= 1)) or 
	((Chr1_Chr2_Row562 #= 548) and (Chr1_Chr2_Freq562 #= 2)) or 
	((Chr1_Chr2_Row562 #= 549) and (Chr1_Chr2_Freq562 #= 2)) or 
	((Chr1_Chr2_Row562 #= 550) and (Chr1_Chr2_Freq562 #= 3)) or 
	((Chr1_Chr2_Row562 #= 551) and (Chr1_Chr2_Freq562 #= 4)) or 
	((Chr1_Chr2_Row562 #= 552) and (Chr1_Chr2_Freq562 #= 5)) or 
	((Chr1_Chr2_Row562 #= 553) and (Chr1_Chr2_Freq562 #= 6)) or 
	((Chr1_Chr2_Row562 #= 554) and (Chr1_Chr2_Freq562 #= 5)) or 
	((Chr1_Chr2_Row562 #= 555) and (Chr1_Chr2_Freq562 #= 3)) or 
	((Chr1_Chr2_Row562 #= 556) and (Chr1_Chr2_Freq562 #= 4)) or 
	((Chr1_Chr2_Row562 #= 557) and (Chr1_Chr2_Freq562 #= 4)) or 
	((Chr1_Chr2_Row562 #= 6) and (Chr1_Chr2_Freq562 #= 3)) or 
	((Chr1_Chr2_Row562 #= 7) and (Chr1_Chr2_Freq562 #= 4)) or 
	((Chr1_Chr2_Row562 #= 8) and (Chr1_Chr2_Freq562 #= 6)) or 
	((Chr1_Chr2_Row562 #= 9) and (Chr1_Chr2_Freq562 #= 5)) or 
	((Chr1_Chr2_Freq562 #= 0) and (Chr1_Chr2_Row562 #= 0)),

 	((Chr1_Chr2_Row563 #= 10) and (Chr1_Chr2_Freq563 #= 3)) or 
	((Chr1_Chr2_Row563 #= 11) and (Chr1_Chr2_Freq563 #= 1)) or 
	((Chr1_Chr2_Row563 #= 12) and (Chr1_Chr2_Freq563 #= 1)) or 
	((Chr1_Chr2_Row563 #= 13) and (Chr1_Chr2_Freq563 #= 1)) or 
	((Chr1_Chr2_Row563 #= 14) and (Chr1_Chr2_Freq563 #= 1)) or 
	((Chr1_Chr2_Row563 #= 15) and (Chr1_Chr2_Freq563 #= 1)) or 
	((Chr1_Chr2_Row563 #= 543) and (Chr1_Chr2_Freq563 #= 1)) or 
	((Chr1_Chr2_Row563 #= 544) and (Chr1_Chr2_Freq563 #= 1)) or 
	((Chr1_Chr2_Row563 #= 545) and (Chr1_Chr2_Freq563 #= 1)) or 
	((Chr1_Chr2_Row563 #= 546) and (Chr1_Chr2_Freq563 #= 1)) or 
	((Chr1_Chr2_Row563 #= 547) and (Chr1_Chr2_Freq563 #= 1)) or 
	((Chr1_Chr2_Row563 #= 548) and (Chr1_Chr2_Freq563 #= 2)) or 
	((Chr1_Chr2_Row563 #= 549) and (Chr1_Chr2_Freq563 #= 2)) or 
	((Chr1_Chr2_Row563 #= 550) and (Chr1_Chr2_Freq563 #= 3)) or 
	((Chr1_Chr2_Row563 #= 551) and (Chr1_Chr2_Freq563 #= 5)) or 
	((Chr1_Chr2_Row563 #= 552) and (Chr1_Chr2_Freq563 #= 5)) or 
	((Chr1_Chr2_Row563 #= 553) and (Chr1_Chr2_Freq563 #= 6)) or 
	((Chr1_Chr2_Row563 #= 554) and (Chr1_Chr2_Freq563 #= 5)) or 
	((Chr1_Chr2_Row563 #= 555) and (Chr1_Chr2_Freq563 #= 4)) or 
	((Chr1_Chr2_Row563 #= 556) and (Chr1_Chr2_Freq563 #= 4)) or 
	((Chr1_Chr2_Row563 #= 557) and (Chr1_Chr2_Freq563 #= 3)) or 
	((Chr1_Chr2_Row563 #= 6) and (Chr1_Chr2_Freq563 #= 4)) or 
	((Chr1_Chr2_Row563 #= 7) and (Chr1_Chr2_Freq563 #= 5)) or 
	((Chr1_Chr2_Row563 #= 8) and (Chr1_Chr2_Freq563 #= 6)) or 
	((Chr1_Chr2_Row563 #= 9) and (Chr1_Chr2_Freq563 #= 6)) or 
	((Chr1_Chr2_Freq563 #= 0) and (Chr1_Chr2_Row563 #= 0)),

 	((Chr1_Chr2_Row564 #= 10) and (Chr1_Chr2_Freq564 #= 4)) or 
	((Chr1_Chr2_Row564 #= 11) and (Chr1_Chr2_Freq564 #= 2)) or 
	((Chr1_Chr2_Row564 #= 12) and (Chr1_Chr2_Freq564 #= 2)) or 
	((Chr1_Chr2_Row564 #= 13) and (Chr1_Chr2_Freq564 #= 1)) or 
	((Chr1_Chr2_Row564 #= 14) and (Chr1_Chr2_Freq564 #= 1)) or 
	((Chr1_Chr2_Row564 #= 15) and (Chr1_Chr2_Freq564 #= 1)) or 
	((Chr1_Chr2_Row564 #= 16) and (Chr1_Chr2_Freq564 #= 1)) or 
	((Chr1_Chr2_Row564 #= 544) and (Chr1_Chr2_Freq564 #= 1)) or 
	((Chr1_Chr2_Row564 #= 545) and (Chr1_Chr2_Freq564 #= 1)) or 
	((Chr1_Chr2_Row564 #= 546) and (Chr1_Chr2_Freq564 #= 1)) or 
	((Chr1_Chr2_Row564 #= 547) and (Chr1_Chr2_Freq564 #= 2)) or 
	((Chr1_Chr2_Row564 #= 548) and (Chr1_Chr2_Freq564 #= 2)) or 
	((Chr1_Chr2_Row564 #= 549) and (Chr1_Chr2_Freq564 #= 2)) or 
	((Chr1_Chr2_Row564 #= 550) and (Chr1_Chr2_Freq564 #= 3)) or 
	((Chr1_Chr2_Row564 #= 551) and (Chr1_Chr2_Freq564 #= 4)) or 
	((Chr1_Chr2_Row564 #= 552) and (Chr1_Chr2_Freq564 #= 4)) or 
	((Chr1_Chr2_Row564 #= 553) and (Chr1_Chr2_Freq564 #= 4)) or 
	((Chr1_Chr2_Row564 #= 554) and (Chr1_Chr2_Freq564 #= 3)) or 
	((Chr1_Chr2_Row564 #= 555) and (Chr1_Chr2_Freq564 #= 3)) or 
	((Chr1_Chr2_Row564 #= 556) and (Chr1_Chr2_Freq564 #= 2)) or 
	((Chr1_Chr2_Row564 #= 557) and (Chr1_Chr2_Freq564 #= 3)) or 
	((Chr1_Chr2_Row564 #= 6) and (Chr1_Chr2_Freq564 #= 4)) or 
	((Chr1_Chr2_Row564 #= 7) and (Chr1_Chr2_Freq564 #= 4)) or 
	((Chr1_Chr2_Row564 #= 8) and (Chr1_Chr2_Freq564 #= 5)) or 
	((Chr1_Chr2_Row564 #= 9) and (Chr1_Chr2_Freq564 #= 4)) or 
	((Chr1_Chr2_Freq564 #= 0) and (Chr1_Chr2_Row564 #= 0)),

 	((Chr1_Chr2_Row565 #= 10) and (Chr1_Chr2_Freq565 #= 4)) or 
	((Chr1_Chr2_Row565 #= 11) and (Chr1_Chr2_Freq565 #= 1)) or 
	((Chr1_Chr2_Row565 #= 12) and (Chr1_Chr2_Freq565 #= 1)) or 
	((Chr1_Chr2_Row565 #= 13) and (Chr1_Chr2_Freq565 #= 1)) or 
	((Chr1_Chr2_Row565 #= 14) and (Chr1_Chr2_Freq565 #= 1)) or 
	((Chr1_Chr2_Row565 #= 15) and (Chr1_Chr2_Freq565 #= 1)) or 
	((Chr1_Chr2_Row565 #= 16) and (Chr1_Chr2_Freq565 #= 1)) or 
	((Chr1_Chr2_Row565 #= 544) and (Chr1_Chr2_Freq565 #= 1)) or 
	((Chr1_Chr2_Row565 #= 545) and (Chr1_Chr2_Freq565 #= 1)) or 
	((Chr1_Chr2_Row565 #= 546) and (Chr1_Chr2_Freq565 #= 1)) or 
	((Chr1_Chr2_Row565 #= 547) and (Chr1_Chr2_Freq565 #= 2)) or 
	((Chr1_Chr2_Row565 #= 548) and (Chr1_Chr2_Freq565 #= 2)) or 
	((Chr1_Chr2_Row565 #= 549) and (Chr1_Chr2_Freq565 #= 2)) or 
	((Chr1_Chr2_Row565 #= 550) and (Chr1_Chr2_Freq565 #= 3)) or 
	((Chr1_Chr2_Row565 #= 551) and (Chr1_Chr2_Freq565 #= 5)) or 
	((Chr1_Chr2_Row565 #= 552) and (Chr1_Chr2_Freq565 #= 5)) or 
	((Chr1_Chr2_Row565 #= 553) and (Chr1_Chr2_Freq565 #= 4)) or 
	((Chr1_Chr2_Row565 #= 554) and (Chr1_Chr2_Freq565 #= 4)) or 
	((Chr1_Chr2_Row565 #= 555) and (Chr1_Chr2_Freq565 #= 3)) or 
	((Chr1_Chr2_Row565 #= 556) and (Chr1_Chr2_Freq565 #= 3)) or 
	((Chr1_Chr2_Row565 #= 557) and (Chr1_Chr2_Freq565 #= 2)) or 
	((Chr1_Chr2_Row565 #= 6) and (Chr1_Chr2_Freq565 #= 4)) or 
	((Chr1_Chr2_Row565 #= 7) and (Chr1_Chr2_Freq565 #= 4)) or 
	((Chr1_Chr2_Row565 #= 8) and (Chr1_Chr2_Freq565 #= 6)) or 
	((Chr1_Chr2_Row565 #= 9) and (Chr1_Chr2_Freq565 #= 5)) or 
	((Chr1_Chr2_Freq565 #= 0) and (Chr1_Chr2_Row565 #= 0)),

 	((Chr1_Chr2_Row566 #= 10) and (Chr1_Chr2_Freq566 #= 3)) or 
	((Chr1_Chr2_Row566 #= 11) and (Chr1_Chr2_Freq566 #= 1)) or 
	((Chr1_Chr2_Row566 #= 12) and (Chr1_Chr2_Freq566 #= 1)) or 
	((Chr1_Chr2_Row566 #= 13) and (Chr1_Chr2_Freq566 #= 1)) or 
	((Chr1_Chr2_Row566 #= 14) and (Chr1_Chr2_Freq566 #= 1)) or 
	((Chr1_Chr2_Row566 #= 15) and (Chr1_Chr2_Freq566 #= 1)) or 
	((Chr1_Chr2_Row566 #= 16) and (Chr1_Chr2_Freq566 #= 1)) or 
	((Chr1_Chr2_Row566 #= 545) and (Chr1_Chr2_Freq566 #= 1)) or 
	((Chr1_Chr2_Row566 #= 546) and (Chr1_Chr2_Freq566 #= 1)) or 
	((Chr1_Chr2_Row566 #= 547) and (Chr1_Chr2_Freq566 #= 1)) or 
	((Chr1_Chr2_Row566 #= 548) and (Chr1_Chr2_Freq566 #= 2)) or 
	((Chr1_Chr2_Row566 #= 549) and (Chr1_Chr2_Freq566 #= 2)) or 
	((Chr1_Chr2_Row566 #= 550) and (Chr1_Chr2_Freq566 #= 3)) or 
	((Chr1_Chr2_Row566 #= 551) and (Chr1_Chr2_Freq566 #= 4)) or 
	((Chr1_Chr2_Row566 #= 552) and (Chr1_Chr2_Freq566 #= 4)) or 
	((Chr1_Chr2_Row566 #= 553) and (Chr1_Chr2_Freq566 #= 4)) or 
	((Chr1_Chr2_Row566 #= 554) and (Chr1_Chr2_Freq566 #= 3)) or 
	((Chr1_Chr2_Row566 #= 555) and (Chr1_Chr2_Freq566 #= 2)) or 
	((Chr1_Chr2_Row566 #= 556) and (Chr1_Chr2_Freq566 #= 2)) or 
	((Chr1_Chr2_Row566 #= 557) and (Chr1_Chr2_Freq566 #= 2)) or 
	((Chr1_Chr2_Row566 #= 6) and (Chr1_Chr2_Freq566 #= 3)) or 
	((Chr1_Chr2_Row566 #= 7) and (Chr1_Chr2_Freq566 #= 4)) or 
	((Chr1_Chr2_Row566 #= 8) and (Chr1_Chr2_Freq566 #= 5)) or 
	((Chr1_Chr2_Row566 #= 9) and (Chr1_Chr2_Freq566 #= 5)) or 
	((Chr1_Chr2_Freq566 #= 0) and (Chr1_Chr2_Row566 #= 0)),

 	((Chr1_Chr2_Row567 #= 10) and (Chr1_Chr2_Freq567 #= 3)) or 
	((Chr1_Chr2_Row567 #= 11) and (Chr1_Chr2_Freq567 #= 1)) or 
	((Chr1_Chr2_Row567 #= 12) and (Chr1_Chr2_Freq567 #= 1)) or 
	((Chr1_Chr2_Row567 #= 13) and (Chr1_Chr2_Freq567 #= 1)) or 
	((Chr1_Chr2_Row567 #= 14) and (Chr1_Chr2_Freq567 #= 1)) or 
	((Chr1_Chr2_Row567 #= 16) and (Chr1_Chr2_Freq567 #= 1)) or 
	((Chr1_Chr2_Row567 #= 545) and (Chr1_Chr2_Freq567 #= 1)) or 
	((Chr1_Chr2_Row567 #= 546) and (Chr1_Chr2_Freq567 #= 1)) or 
	((Chr1_Chr2_Row567 #= 547) and (Chr1_Chr2_Freq567 #= 1)) or 
	((Chr1_Chr2_Row567 #= 548) and (Chr1_Chr2_Freq567 #= 1)) or 
	((Chr1_Chr2_Row567 #= 549) and (Chr1_Chr2_Freq567 #= 2)) or 
	((Chr1_Chr2_Row567 #= 550) and (Chr1_Chr2_Freq567 #= 2)) or 
	((Chr1_Chr2_Row567 #= 551) and (Chr1_Chr2_Freq567 #= 3)) or 
	((Chr1_Chr2_Row567 #= 552) and (Chr1_Chr2_Freq567 #= 3)) or 
	((Chr1_Chr2_Row567 #= 553) and (Chr1_Chr2_Freq567 #= 3)) or 
	((Chr1_Chr2_Row567 #= 554) and (Chr1_Chr2_Freq567 #= 2)) or 
	((Chr1_Chr2_Row567 #= 555) and (Chr1_Chr2_Freq567 #= 2)) or 
	((Chr1_Chr2_Row567 #= 556) and (Chr1_Chr2_Freq567 #= 1)) or 
	((Chr1_Chr2_Row567 #= 557) and (Chr1_Chr2_Freq567 #= 1)) or 
	((Chr1_Chr2_Row567 #= 6) and (Chr1_Chr2_Freq567 #= 3)) or 
	((Chr1_Chr2_Row567 #= 7) and (Chr1_Chr2_Freq567 #= 3)) or 
	((Chr1_Chr2_Row567 #= 8) and (Chr1_Chr2_Freq567 #= 4)) or 
	((Chr1_Chr2_Row567 #= 9) and (Chr1_Chr2_Freq567 #= 4)) or 
	((Chr1_Chr2_Freq567 #= 0) and (Chr1_Chr2_Row567 #= 0)),

 	((Chr1_Chr2_Row568 #= 10) and (Chr1_Chr2_Freq568 #= 2)) or 
	((Chr1_Chr2_Row568 #= 11) and (Chr1_Chr2_Freq568 #= 2)) or 
	((Chr1_Chr2_Row568 #= 12) and (Chr1_Chr2_Freq568 #= 1)) or 
	((Chr1_Chr2_Row568 #= 13) and (Chr1_Chr2_Freq568 #= 1)) or 
	((Chr1_Chr2_Row568 #= 14) and (Chr1_Chr2_Freq568 #= 1)) or 
	((Chr1_Chr2_Row568 #= 15) and (Chr1_Chr2_Freq568 #= 1)) or 
	((Chr1_Chr2_Row568 #= 16) and (Chr1_Chr2_Freq568 #= 1)) or 
	((Chr1_Chr2_Row568 #= 18) and (Chr1_Chr2_Freq568 #= 1)) or 
	((Chr1_Chr2_Row568 #= 544) and (Chr1_Chr2_Freq568 #= 1)) or 
	((Chr1_Chr2_Row568 #= 545) and (Chr1_Chr2_Freq568 #= 1)) or 
	((Chr1_Chr2_Row568 #= 546) and (Chr1_Chr2_Freq568 #= 1)) or 
	((Chr1_Chr2_Row568 #= 547) and (Chr1_Chr2_Freq568 #= 1)) or 
	((Chr1_Chr2_Row568 #= 548) and (Chr1_Chr2_Freq568 #= 1)) or 
	((Chr1_Chr2_Row568 #= 549) and (Chr1_Chr2_Freq568 #= 1)) or 
	((Chr1_Chr2_Row568 #= 550) and (Chr1_Chr2_Freq568 #= 2)) or 
	((Chr1_Chr2_Row568 #= 551) and (Chr1_Chr2_Freq568 #= 2)) or 
	((Chr1_Chr2_Row568 #= 552) and (Chr1_Chr2_Freq568 #= 2)) or 
	((Chr1_Chr2_Row568 #= 553) and (Chr1_Chr2_Freq568 #= 1)) or 
	((Chr1_Chr2_Row568 #= 554) and (Chr1_Chr2_Freq568 #= 1)) or 
	((Chr1_Chr2_Row568 #= 555) and (Chr1_Chr2_Freq568 #= 1)) or 
	((Chr1_Chr2_Row568 #= 556) and (Chr1_Chr2_Freq568 #= 1)) or 
	((Chr1_Chr2_Row568 #= 557) and (Chr1_Chr2_Freq568 #= 1)) or 
	((Chr1_Chr2_Row568 #= 6) and (Chr1_Chr2_Freq568 #= 2)) or 
	((Chr1_Chr2_Row568 #= 7) and (Chr1_Chr2_Freq568 #= 2)) or 
	((Chr1_Chr2_Row568 #= 8) and (Chr1_Chr2_Freq568 #= 2)) or 
	((Chr1_Chr2_Row568 #= 9) and (Chr1_Chr2_Freq568 #= 2)) or 
	((Chr1_Chr2_Freq568 #= 0) and (Chr1_Chr2_Row568 #= 0)),

 	((Chr1_Chr2_Row569 #= 10) and (Chr1_Chr2_Freq569 #= 2)) or 
	((Chr1_Chr2_Row569 #= 11) and (Chr1_Chr2_Freq569 #= 1)) or 
	((Chr1_Chr2_Row569 #= 12) and (Chr1_Chr2_Freq569 #= 1)) or 
	((Chr1_Chr2_Row569 #= 13) and (Chr1_Chr2_Freq569 #= 1)) or 
	((Chr1_Chr2_Row569 #= 14) and (Chr1_Chr2_Freq569 #= 1)) or 
	((Chr1_Chr2_Row569 #= 16) and (Chr1_Chr2_Freq569 #= 1)) or 
	((Chr1_Chr2_Row569 #= 543) and (Chr1_Chr2_Freq569 #= 1)) or 
	((Chr1_Chr2_Row569 #= 545) and (Chr1_Chr2_Freq569 #= 1)) or 
	((Chr1_Chr2_Row569 #= 546) and (Chr1_Chr2_Freq569 #= 1)) or 
	((Chr1_Chr2_Row569 #= 547) and (Chr1_Chr2_Freq569 #= 1)) or 
	((Chr1_Chr2_Row569 #= 548) and (Chr1_Chr2_Freq569 #= 1)) or 
	((Chr1_Chr2_Row569 #= 549) and (Chr1_Chr2_Freq569 #= 1)) or 
	((Chr1_Chr2_Row569 #= 550) and (Chr1_Chr2_Freq569 #= 2)) or 
	((Chr1_Chr2_Row569 #= 551) and (Chr1_Chr2_Freq569 #= 1)) or 
	((Chr1_Chr2_Row569 #= 552) and (Chr1_Chr2_Freq569 #= 2)) or 
	((Chr1_Chr2_Row569 #= 553) and (Chr1_Chr2_Freq569 #= 1)) or 
	((Chr1_Chr2_Row569 #= 554) and (Chr1_Chr2_Freq569 #= 1)) or 
	((Chr1_Chr2_Row569 #= 555) and (Chr1_Chr2_Freq569 #= 1)) or 
	((Chr1_Chr2_Row569 #= 556) and (Chr1_Chr2_Freq569 #= 1)) or 
	((Chr1_Chr2_Row569 #= 557) and (Chr1_Chr2_Freq569 #= 1)) or 
	((Chr1_Chr2_Row569 #= 6) and (Chr1_Chr2_Freq569 #= 2)) or 
	((Chr1_Chr2_Row569 #= 7) and (Chr1_Chr2_Freq569 #= 2)) or 
	((Chr1_Chr2_Row569 #= 8) and (Chr1_Chr2_Freq569 #= 2)) or 
	((Chr1_Chr2_Row569 #= 9) and (Chr1_Chr2_Freq569 #= 1)) or 
	((Chr1_Chr2_Freq569 #= 0) and (Chr1_Chr2_Row569 #= 0)),

 	((Chr1_Chr2_Row570 #= 10) and (Chr1_Chr2_Freq570 #= 2)) or 
	((Chr1_Chr2_Row570 #= 11) and (Chr1_Chr2_Freq570 #= 2)) or 
	((Chr1_Chr2_Row570 #= 12) and (Chr1_Chr2_Freq570 #= 1)) or 
	((Chr1_Chr2_Row570 #= 13) and (Chr1_Chr2_Freq570 #= 1)) or 
	((Chr1_Chr2_Row570 #= 14) and (Chr1_Chr2_Freq570 #= 1)) or 
	((Chr1_Chr2_Row570 #= 15) and (Chr1_Chr2_Freq570 #= 1)) or 
	((Chr1_Chr2_Row570 #= 544) and (Chr1_Chr2_Freq570 #= 1)) or 
	((Chr1_Chr2_Row570 #= 545) and (Chr1_Chr2_Freq570 #= 1)) or 
	((Chr1_Chr2_Row570 #= 546) and (Chr1_Chr2_Freq570 #= 1)) or 
	((Chr1_Chr2_Row570 #= 547) and (Chr1_Chr2_Freq570 #= 1)) or 
	((Chr1_Chr2_Row570 #= 548) and (Chr1_Chr2_Freq570 #= 1)) or 
	((Chr1_Chr2_Row570 #= 549) and (Chr1_Chr2_Freq570 #= 1)) or 
	((Chr1_Chr2_Row570 #= 550) and (Chr1_Chr2_Freq570 #= 1)) or 
	((Chr1_Chr2_Row570 #= 551) and (Chr1_Chr2_Freq570 #= 2)) or 
	((Chr1_Chr2_Row570 #= 552) and (Chr1_Chr2_Freq570 #= 2)) or 
	((Chr1_Chr2_Row570 #= 553) and (Chr1_Chr2_Freq570 #= 2)) or 
	((Chr1_Chr2_Row570 #= 554) and (Chr1_Chr2_Freq570 #= 1)) or 
	((Chr1_Chr2_Row570 #= 555) and (Chr1_Chr2_Freq570 #= 1)) or 
	((Chr1_Chr2_Row570 #= 556) and (Chr1_Chr2_Freq570 #= 1)) or 
	((Chr1_Chr2_Row570 #= 557) and (Chr1_Chr2_Freq570 #= 1)) or 
	((Chr1_Chr2_Row570 #= 6) and (Chr1_Chr2_Freq570 #= 2)) or 
	((Chr1_Chr2_Row570 #= 7) and (Chr1_Chr2_Freq570 #= 2)) or 
	((Chr1_Chr2_Row570 #= 8) and (Chr1_Chr2_Freq570 #= 3)) or 
	((Chr1_Chr2_Row570 #= 9) and (Chr1_Chr2_Freq570 #= 2)) or 
	((Chr1_Chr2_Freq570 #= 0) and (Chr1_Chr2_Row570 #= 0)),

 	((Chr1_Chr2_Row571 #= 10) and (Chr1_Chr2_Freq571 #= 2)) or 
	((Chr1_Chr2_Row571 #= 11) and (Chr1_Chr2_Freq571 #= 1)) or 
	((Chr1_Chr2_Row571 #= 12) and (Chr1_Chr2_Freq571 #= 1)) or 
	((Chr1_Chr2_Row571 #= 13) and (Chr1_Chr2_Freq571 #= 1)) or 
	((Chr1_Chr2_Row571 #= 14) and (Chr1_Chr2_Freq571 #= 1)) or 
	((Chr1_Chr2_Row571 #= 15) and (Chr1_Chr2_Freq571 #= 1)) or 
	((Chr1_Chr2_Row571 #= 16) and (Chr1_Chr2_Freq571 #= 1)) or 
	((Chr1_Chr2_Row571 #= 544) and (Chr1_Chr2_Freq571 #= 1)) or 
	((Chr1_Chr2_Row571 #= 545) and (Chr1_Chr2_Freq571 #= 1)) or 
	((Chr1_Chr2_Row571 #= 546) and (Chr1_Chr2_Freq571 #= 1)) or 
	((Chr1_Chr2_Row571 #= 547) and (Chr1_Chr2_Freq571 #= 1)) or 
	((Chr1_Chr2_Row571 #= 548) and (Chr1_Chr2_Freq571 #= 1)) or 
	((Chr1_Chr2_Row571 #= 549) and (Chr1_Chr2_Freq571 #= 1)) or 
	((Chr1_Chr2_Row571 #= 550) and (Chr1_Chr2_Freq571 #= 1)) or 
	((Chr1_Chr2_Row571 #= 551) and (Chr1_Chr2_Freq571 #= 2)) or 
	((Chr1_Chr2_Row571 #= 552) and (Chr1_Chr2_Freq571 #= 2)) or 
	((Chr1_Chr2_Row571 #= 553) and (Chr1_Chr2_Freq571 #= 1)) or 
	((Chr1_Chr2_Row571 #= 554) and (Chr1_Chr2_Freq571 #= 1)) or 
	((Chr1_Chr2_Row571 #= 555) and (Chr1_Chr2_Freq571 #= 1)) or 
	((Chr1_Chr2_Row571 #= 557) and (Chr1_Chr2_Freq571 #= 1)) or 
	((Chr1_Chr2_Row571 #= 6) and (Chr1_Chr2_Freq571 #= 2)) or 
	((Chr1_Chr2_Row571 #= 7) and (Chr1_Chr2_Freq571 #= 2)) or 
	((Chr1_Chr2_Row571 #= 8) and (Chr1_Chr2_Freq571 #= 2)) or 
	((Chr1_Chr2_Row571 #= 9) and (Chr1_Chr2_Freq571 #= 3)) or 
	((Chr1_Chr2_Freq571 #= 0) and (Chr1_Chr2_Row571 #= 0)),

 	((Chr1_Chr2_Row572 #= 10) and (Chr1_Chr2_Freq572 #= 1)) or 
	((Chr1_Chr2_Row572 #= 11) and (Chr1_Chr2_Freq572 #= 1)) or 
	((Chr1_Chr2_Row572 #= 12) and (Chr1_Chr2_Freq572 #= 1)) or 
	((Chr1_Chr2_Row572 #= 13) and (Chr1_Chr2_Freq572 #= 1)) or 
	((Chr1_Chr2_Row572 #= 14) and (Chr1_Chr2_Freq572 #= 1)) or 
	((Chr1_Chr2_Row572 #= 15) and (Chr1_Chr2_Freq572 #= 1)) or 
	((Chr1_Chr2_Row572 #= 545) and (Chr1_Chr2_Freq572 #= 1)) or 
	((Chr1_Chr2_Row572 #= 546) and (Chr1_Chr2_Freq572 #= 1)) or 
	((Chr1_Chr2_Row572 #= 547) and (Chr1_Chr2_Freq572 #= 1)) or 
	((Chr1_Chr2_Row572 #= 548) and (Chr1_Chr2_Freq572 #= 1)) or 
	((Chr1_Chr2_Row572 #= 549) and (Chr1_Chr2_Freq572 #= 1)) or 
	((Chr1_Chr2_Row572 #= 550) and (Chr1_Chr2_Freq572 #= 1)) or 
	((Chr1_Chr2_Row572 #= 551) and (Chr1_Chr2_Freq572 #= 1)) or 
	((Chr1_Chr2_Row572 #= 552) and (Chr1_Chr2_Freq572 #= 1)) or 
	((Chr1_Chr2_Row572 #= 553) and (Chr1_Chr2_Freq572 #= 1)) or 
	((Chr1_Chr2_Row572 #= 554) and (Chr1_Chr2_Freq572 #= 1)) or 
	((Chr1_Chr2_Row572 #= 557) and (Chr1_Chr2_Freq572 #= 1)) or 
	((Chr1_Chr2_Row572 #= 6) and (Chr1_Chr2_Freq572 #= 1)) or 
	((Chr1_Chr2_Row572 #= 7) and (Chr1_Chr2_Freq572 #= 2)) or 
	((Chr1_Chr2_Row572 #= 8) and (Chr1_Chr2_Freq572 #= 1)) or 
	((Chr1_Chr2_Row572 #= 9) and (Chr1_Chr2_Freq572 #= 2)) or 
	((Chr1_Chr2_Freq572 #= 0) and (Chr1_Chr2_Row572 #= 0)),

 	((Chr1_Chr2_Row573 #= 10) and (Chr1_Chr2_Freq573 #= 1)) or 
	((Chr1_Chr2_Row573 #= 11) and (Chr1_Chr2_Freq573 #= 1)) or 
	((Chr1_Chr2_Row573 #= 12) and (Chr1_Chr2_Freq573 #= 1)) or 
	((Chr1_Chr2_Row573 #= 13) and (Chr1_Chr2_Freq573 #= 1)) or 
	((Chr1_Chr2_Row573 #= 14) and (Chr1_Chr2_Freq573 #= 1)) or 
	((Chr1_Chr2_Row573 #= 15) and (Chr1_Chr2_Freq573 #= 1)) or 
	((Chr1_Chr2_Row573 #= 546) and (Chr1_Chr2_Freq573 #= 1)) or 
	((Chr1_Chr2_Row573 #= 547) and (Chr1_Chr2_Freq573 #= 1)) or 
	((Chr1_Chr2_Row573 #= 548) and (Chr1_Chr2_Freq573 #= 1)) or 
	((Chr1_Chr2_Row573 #= 549) and (Chr1_Chr2_Freq573 #= 1)) or 
	((Chr1_Chr2_Row573 #= 550) and (Chr1_Chr2_Freq573 #= 1)) or 
	((Chr1_Chr2_Row573 #= 551) and (Chr1_Chr2_Freq573 #= 1)) or 
	((Chr1_Chr2_Row573 #= 552) and (Chr1_Chr2_Freq573 #= 1)) or 
	((Chr1_Chr2_Row573 #= 553) and (Chr1_Chr2_Freq573 #= 1)) or 
	((Chr1_Chr2_Row573 #= 555) and (Chr1_Chr2_Freq573 #= 1)) or 
	((Chr1_Chr2_Row573 #= 556) and (Chr1_Chr2_Freq573 #= 1)) or 
	((Chr1_Chr2_Row573 #= 557) and (Chr1_Chr2_Freq573 #= 1)) or 
	((Chr1_Chr2_Row573 #= 6) and (Chr1_Chr2_Freq573 #= 1)) or 
	((Chr1_Chr2_Row573 #= 7) and (Chr1_Chr2_Freq573 #= 1)) or 
	((Chr1_Chr2_Row573 #= 8) and (Chr1_Chr2_Freq573 #= 1)) or 
	((Chr1_Chr2_Row573 #= 9) and (Chr1_Chr2_Freq573 #= 2)) or 
	((Chr1_Chr2_Freq573 #= 0) and (Chr1_Chr2_Row573 #= 0)),

 	((Chr1_Chr2_Row574 #= 10) and (Chr1_Chr2_Freq574 #= 1)) or 
	((Chr1_Chr2_Row574 #= 11) and (Chr1_Chr2_Freq574 #= 1)) or 
	((Chr1_Chr2_Row574 #= 12) and (Chr1_Chr2_Freq574 #= 1)) or 
	((Chr1_Chr2_Row574 #= 13) and (Chr1_Chr2_Freq574 #= 1)) or 
	((Chr1_Chr2_Row574 #= 14) and (Chr1_Chr2_Freq574 #= 1)) or 
	((Chr1_Chr2_Row574 #= 549) and (Chr1_Chr2_Freq574 #= 1)) or 
	((Chr1_Chr2_Row574 #= 550) and (Chr1_Chr2_Freq574 #= 1)) or 
	((Chr1_Chr2_Row574 #= 551) and (Chr1_Chr2_Freq574 #= 1)) or 
	((Chr1_Chr2_Row574 #= 552) and (Chr1_Chr2_Freq574 #= 1)) or 
	((Chr1_Chr2_Row574 #= 553) and (Chr1_Chr2_Freq574 #= 1)) or 
	((Chr1_Chr2_Row574 #= 6) and (Chr1_Chr2_Freq574 #= 1)) or 
	((Chr1_Chr2_Row574 #= 7) and (Chr1_Chr2_Freq574 #= 1)) or 
	((Chr1_Chr2_Row574 #= 8) and (Chr1_Chr2_Freq574 #= 1)) or 
	((Chr1_Chr2_Row574 #= 9) and (Chr1_Chr2_Freq574 #= 1)) or 
	((Chr1_Chr2_Freq574 #= 0) and (Chr1_Chr2_Row574 #= 0)),

 	((Chr1_Chr2_Row575 #= 10) and (Chr1_Chr2_Freq575 #= 1)) or 
	((Chr1_Chr2_Row575 #= 11) and (Chr1_Chr2_Freq575 #= 1)) or 
	((Chr1_Chr2_Row575 #= 12) and (Chr1_Chr2_Freq575 #= 1)) or 
	((Chr1_Chr2_Row575 #= 548) and (Chr1_Chr2_Freq575 #= 1)) or 
	((Chr1_Chr2_Row575 #= 550) and (Chr1_Chr2_Freq575 #= 1)) or 
	((Chr1_Chr2_Row575 #= 6) and (Chr1_Chr2_Freq575 #= 1)) or 
	((Chr1_Chr2_Row575 #= 7) and (Chr1_Chr2_Freq575 #= 1)) or 
	((Chr1_Chr2_Row575 #= 8) and (Chr1_Chr2_Freq575 #= 1)) or 
	((Chr1_Chr2_Row575 #= 9) and (Chr1_Chr2_Freq575 #= 1)) or 
	((Chr1_Chr2_Freq575 #= 0) and (Chr1_Chr2_Row575 #= 0)),

 	((Chr1_Chr2_Row576 #= 11) and (Chr1_Chr2_Freq576 #= 1)) or 
	((Chr1_Chr2_Freq576 #= 0) and (Chr1_Chr2_Row576 #= 0)),

 	((Chr1_Chr2_Row579 #= 6) and (Chr1_Chr2_Freq579 #= 1)) or 
	((Chr1_Chr2_Freq579 #= 0) and (Chr1_Chr2_Row579 #= 0)),

 	((Chr1_Chr2_Row711 #= 370) and (Chr1_Chr2_Freq711 #= 1)) or 
	((Chr1_Chr2_Row711 #= 374) and (Chr1_Chr2_Freq711 #= 1)) or 
	((Chr1_Chr2_Freq711 #= 0) and (Chr1_Chr2_Row711 #= 0)),

 	((Chr1_Chr2_Row712 #= 373) and (Chr1_Chr2_Freq712 #= 1)) or 
	((Chr1_Chr2_Row712 #= 374) and (Chr1_Chr2_Freq712 #= 1)) or 
	((Chr1_Chr2_Row712 #= 380) and (Chr1_Chr2_Freq712 #= 1)) or 
	((Chr1_Chr2_Row712 #= 381) and (Chr1_Chr2_Freq712 #= 1)) or 
	((Chr1_Chr2_Freq712 #= 0) and (Chr1_Chr2_Row712 #= 0)),

 	((Chr1_Chr2_Row713 #= 373) and (Chr1_Chr2_Freq713 #= 1)) or 
	((Chr1_Chr2_Row713 #= 374) and (Chr1_Chr2_Freq713 #= 1)) or 
	((Chr1_Chr2_Row713 #= 381) and (Chr1_Chr2_Freq713 #= 1)) or 
	((Chr1_Chr2_Row713 #= 382) and (Chr1_Chr2_Freq713 #= 1)) or 
	((Chr1_Chr2_Row713 #= 383) and (Chr1_Chr2_Freq713 #= 1)) or 
	((Chr1_Chr2_Row713 #= 384) and (Chr1_Chr2_Freq713 #= 1)) or 
	((Chr1_Chr2_Freq713 #= 0) and (Chr1_Chr2_Row713 #= 0)),

 	((Chr1_Chr2_Row714 #= 372) and (Chr1_Chr2_Freq714 #= 1)) or 
	((Chr1_Chr2_Row714 #= 373) and (Chr1_Chr2_Freq714 #= 1)) or 
	((Chr1_Chr2_Row714 #= 374) and (Chr1_Chr2_Freq714 #= 1)) or 
	((Chr1_Chr2_Row714 #= 380) and (Chr1_Chr2_Freq714 #= 1)) or 
	((Chr1_Chr2_Row714 #= 381) and (Chr1_Chr2_Freq714 #= 1)) or 
	((Chr1_Chr2_Row714 #= 384) and (Chr1_Chr2_Freq714 #= 1)) or 
	((Chr1_Chr2_Row714 #= 385) and (Chr1_Chr2_Freq714 #= 1)) or 
	((Chr1_Chr2_Row714 #= 386) and (Chr1_Chr2_Freq714 #= 1)) or 
	((Chr1_Chr2_Freq714 #= 0) and (Chr1_Chr2_Row714 #= 0)),

 	((Chr1_Chr2_Row715 #= 368) and (Chr1_Chr2_Freq715 #= 1)) or 
	((Chr1_Chr2_Row715 #= 369) and (Chr1_Chr2_Freq715 #= 1)) or 
	((Chr1_Chr2_Row715 #= 370) and (Chr1_Chr2_Freq715 #= 1)) or 
	((Chr1_Chr2_Row715 #= 371) and (Chr1_Chr2_Freq715 #= 1)) or 
	((Chr1_Chr2_Row715 #= 372) and (Chr1_Chr2_Freq715 #= 1)) or 
	((Chr1_Chr2_Row715 #= 373) and (Chr1_Chr2_Freq715 #= 1)) or 
	((Chr1_Chr2_Row715 #= 374) and (Chr1_Chr2_Freq715 #= 1)) or 
	((Chr1_Chr2_Row715 #= 380) and (Chr1_Chr2_Freq715 #= 1)) or 
	((Chr1_Chr2_Row715 #= 381) and (Chr1_Chr2_Freq715 #= 1)) or 
	((Chr1_Chr2_Row715 #= 382) and (Chr1_Chr2_Freq715 #= 1)) or 
	((Chr1_Chr2_Row715 #= 383) and (Chr1_Chr2_Freq715 #= 1)) or 
	((Chr1_Chr2_Row715 #= 385) and (Chr1_Chr2_Freq715 #= 1)) or 
	((Chr1_Chr2_Row715 #= 386) and (Chr1_Chr2_Freq715 #= 1)) or 
	((Chr1_Chr2_Freq715 #= 0) and (Chr1_Chr2_Row715 #= 0)),

 	((Chr1_Chr2_Row716 #= 367) and (Chr1_Chr2_Freq716 #= 1)) or 
	((Chr1_Chr2_Row716 #= 368) and (Chr1_Chr2_Freq716 #= 1)) or 
	((Chr1_Chr2_Row716 #= 369) and (Chr1_Chr2_Freq716 #= 1)) or 
	((Chr1_Chr2_Row716 #= 370) and (Chr1_Chr2_Freq716 #= 1)) or 
	((Chr1_Chr2_Row716 #= 371) and (Chr1_Chr2_Freq716 #= 1)) or 
	((Chr1_Chr2_Row716 #= 372) and (Chr1_Chr2_Freq716 #= 1)) or 
	((Chr1_Chr2_Row716 #= 373) and (Chr1_Chr2_Freq716 #= 1)) or 
	((Chr1_Chr2_Row716 #= 374) and (Chr1_Chr2_Freq716 #= 1)) or 
	((Chr1_Chr2_Row716 #= 380) and (Chr1_Chr2_Freq716 #= 1)) or 
	((Chr1_Chr2_Row716 #= 381) and (Chr1_Chr2_Freq716 #= 2)) or 
	((Chr1_Chr2_Row716 #= 382) and (Chr1_Chr2_Freq716 #= 1)) or 
	((Chr1_Chr2_Row716 #= 383) and (Chr1_Chr2_Freq716 #= 1)) or 
	((Chr1_Chr2_Row716 #= 384) and (Chr1_Chr2_Freq716 #= 1)) or 
	((Chr1_Chr2_Row716 #= 385) and (Chr1_Chr2_Freq716 #= 1)) or 
	((Chr1_Chr2_Row716 #= 386) and (Chr1_Chr2_Freq716 #= 1)) or 
	((Chr1_Chr2_Row716 #= 387) and (Chr1_Chr2_Freq716 #= 1)) or 
	((Chr1_Chr2_Freq716 #= 0) and (Chr1_Chr2_Row716 #= 0)),

 	((Chr1_Chr2_Row717 #= 369) and (Chr1_Chr2_Freq717 #= 1)) or 
	((Chr1_Chr2_Row717 #= 370) and (Chr1_Chr2_Freq717 #= 1)) or 
	((Chr1_Chr2_Row717 #= 371) and (Chr1_Chr2_Freq717 #= 1)) or 
	((Chr1_Chr2_Row717 #= 372) and (Chr1_Chr2_Freq717 #= 1)) or 
	((Chr1_Chr2_Row717 #= 373) and (Chr1_Chr2_Freq717 #= 1)) or 
	((Chr1_Chr2_Row717 #= 374) and (Chr1_Chr2_Freq717 #= 2)) or 
	((Chr1_Chr2_Row717 #= 380) and (Chr1_Chr2_Freq717 #= 2)) or 
	((Chr1_Chr2_Row717 #= 381) and (Chr1_Chr2_Freq717 #= 2)) or 
	((Chr1_Chr2_Row717 #= 382) and (Chr1_Chr2_Freq717 #= 2)) or 
	((Chr1_Chr2_Row717 #= 383) and (Chr1_Chr2_Freq717 #= 1)) or 
	((Chr1_Chr2_Row717 #= 384) and (Chr1_Chr2_Freq717 #= 1)) or 
	((Chr1_Chr2_Row717 #= 385) and (Chr1_Chr2_Freq717 #= 1)) or 
	((Chr1_Chr2_Row717 #= 386) and (Chr1_Chr2_Freq717 #= 1)) or 
	((Chr1_Chr2_Freq717 #= 0) and (Chr1_Chr2_Row717 #= 0)),

 	((Chr1_Chr2_Row724 #= 364) and (Chr1_Chr2_Freq724 #= 1)) or 
	((Chr1_Chr2_Row724 #= 365) and (Chr1_Chr2_Freq724 #= 1)) or 
	((Chr1_Chr2_Row724 #= 367) and (Chr1_Chr2_Freq724 #= 1)) or 
	((Chr1_Chr2_Row724 #= 368) and (Chr1_Chr2_Freq724 #= 1)) or 
	((Chr1_Chr2_Row724 #= 369) and (Chr1_Chr2_Freq724 #= 1)) or 
	((Chr1_Chr2_Row724 #= 370) and (Chr1_Chr2_Freq724 #= 1)) or 
	((Chr1_Chr2_Row724 #= 371) and (Chr1_Chr2_Freq724 #= 2)) or 
	((Chr1_Chr2_Row724 #= 372) and (Chr1_Chr2_Freq724 #= 2)) or 
	((Chr1_Chr2_Row724 #= 373) and (Chr1_Chr2_Freq724 #= 2)) or 
	((Chr1_Chr2_Row724 #= 374) and (Chr1_Chr2_Freq724 #= 3)) or 
	((Chr1_Chr2_Row724 #= 380) and (Chr1_Chr2_Freq724 #= 3)) or 
	((Chr1_Chr2_Row724 #= 381) and (Chr1_Chr2_Freq724 #= 3)) or 
	((Chr1_Chr2_Row724 #= 382) and (Chr1_Chr2_Freq724 #= 2)) or 
	((Chr1_Chr2_Row724 #= 383) and (Chr1_Chr2_Freq724 #= 2)) or 
	((Chr1_Chr2_Row724 #= 384) and (Chr1_Chr2_Freq724 #= 2)) or 
	((Chr1_Chr2_Row724 #= 385) and (Chr1_Chr2_Freq724 #= 1)) or 
	((Chr1_Chr2_Row724 #= 386) and (Chr1_Chr2_Freq724 #= 1)) or 
	((Chr1_Chr2_Row724 #= 387) and (Chr1_Chr2_Freq724 #= 1)) or 
	((Chr1_Chr2_Row724 #= 388) and (Chr1_Chr2_Freq724 #= 1)) or 
	((Chr1_Chr2_Row724 #= 389) and (Chr1_Chr2_Freq724 #= 1)) or 
	((Chr1_Chr2_Row724 #= 390) and (Chr1_Chr2_Freq724 #= 1)) or 
	((Chr1_Chr2_Freq724 #= 0) and (Chr1_Chr2_Row724 #= 0)),

 	((Chr1_Chr2_Row725 #= 365) and (Chr1_Chr2_Freq725 #= 1)) or 
	((Chr1_Chr2_Row725 #= 366) and (Chr1_Chr2_Freq725 #= 1)) or 
	((Chr1_Chr2_Row725 #= 367) and (Chr1_Chr2_Freq725 #= 1)) or 
	((Chr1_Chr2_Row725 #= 368) and (Chr1_Chr2_Freq725 #= 1)) or 
	((Chr1_Chr2_Row725 #= 369) and (Chr1_Chr2_Freq725 #= 1)) or 
	((Chr1_Chr2_Row725 #= 370) and (Chr1_Chr2_Freq725 #= 1)) or 
	((Chr1_Chr2_Row725 #= 371) and (Chr1_Chr2_Freq725 #= 2)) or 
	((Chr1_Chr2_Row725 #= 372) and (Chr1_Chr2_Freq725 #= 2)) or 
	((Chr1_Chr2_Row725 #= 373) and (Chr1_Chr2_Freq725 #= 2)) or 
	((Chr1_Chr2_Row725 #= 374) and (Chr1_Chr2_Freq725 #= 3)) or 
	((Chr1_Chr2_Row725 #= 380) and (Chr1_Chr2_Freq725 #= 3)) or 
	((Chr1_Chr2_Row725 #= 381) and (Chr1_Chr2_Freq725 #= 3)) or 
	((Chr1_Chr2_Row725 #= 382) and (Chr1_Chr2_Freq725 #= 2)) or 
	((Chr1_Chr2_Row725 #= 383) and (Chr1_Chr2_Freq725 #= 2)) or 
	((Chr1_Chr2_Row725 #= 384) and (Chr1_Chr2_Freq725 #= 2)) or 
	((Chr1_Chr2_Row725 #= 385) and (Chr1_Chr2_Freq725 #= 2)) or 
	((Chr1_Chr2_Row725 #= 386) and (Chr1_Chr2_Freq725 #= 1)) or 
	((Chr1_Chr2_Row725 #= 387) and (Chr1_Chr2_Freq725 #= 1)) or 
	((Chr1_Chr2_Row725 #= 388) and (Chr1_Chr2_Freq725 #= 1)) or 
	((Chr1_Chr2_Row725 #= 389) and (Chr1_Chr2_Freq725 #= 1)) or 
	((Chr1_Chr2_Row725 #= 390) and (Chr1_Chr2_Freq725 #= 1)) or 
	((Chr1_Chr2_Freq725 #= 0) and (Chr1_Chr2_Row725 #= 0)),

 	((Chr1_Chr2_Row726 #= 366) and (Chr1_Chr2_Freq726 #= 1)) or 
	((Chr1_Chr2_Row726 #= 367) and (Chr1_Chr2_Freq726 #= 1)) or 
	((Chr1_Chr2_Row726 #= 368) and (Chr1_Chr2_Freq726 #= 1)) or 
	((Chr1_Chr2_Row726 #= 369) and (Chr1_Chr2_Freq726 #= 1)) or 
	((Chr1_Chr2_Row726 #= 370) and (Chr1_Chr2_Freq726 #= 1)) or 
	((Chr1_Chr2_Row726 #= 371) and (Chr1_Chr2_Freq726 #= 1)) or 
	((Chr1_Chr2_Row726 #= 372) and (Chr1_Chr2_Freq726 #= 2)) or 
	((Chr1_Chr2_Row726 #= 373) and (Chr1_Chr2_Freq726 #= 2)) or 
	((Chr1_Chr2_Row726 #= 374) and (Chr1_Chr2_Freq726 #= 2)) or 
	((Chr1_Chr2_Row726 #= 380) and (Chr1_Chr2_Freq726 #= 3)) or 
	((Chr1_Chr2_Row726 #= 381) and (Chr1_Chr2_Freq726 #= 2)) or 
	((Chr1_Chr2_Row726 #= 382) and (Chr1_Chr2_Freq726 #= 2)) or 
	((Chr1_Chr2_Row726 #= 383) and (Chr1_Chr2_Freq726 #= 2)) or 
	((Chr1_Chr2_Row726 #= 384) and (Chr1_Chr2_Freq726 #= 1)) or 
	((Chr1_Chr2_Row726 #= 385) and (Chr1_Chr2_Freq726 #= 2)) or 
	((Chr1_Chr2_Row726 #= 386) and (Chr1_Chr2_Freq726 #= 1)) or 
	((Chr1_Chr2_Row726 #= 387) and (Chr1_Chr2_Freq726 #= 1)) or 
	((Chr1_Chr2_Row726 #= 388) and (Chr1_Chr2_Freq726 #= 1)) or 
	((Chr1_Chr2_Row726 #= 389) and (Chr1_Chr2_Freq726 #= 1)) or 
	((Chr1_Chr2_Row726 #= 390) and (Chr1_Chr2_Freq726 #= 1)) or 
	((Chr1_Chr2_Freq726 #= 0) and (Chr1_Chr2_Row726 #= 0)),

 	((Chr1_Chr2_Row727 #= 365) and (Chr1_Chr2_Freq727 #= 1)) or 
	((Chr1_Chr2_Row727 #= 366) and (Chr1_Chr2_Freq727 #= 1)) or 
	((Chr1_Chr2_Row727 #= 367) and (Chr1_Chr2_Freq727 #= 1)) or 
	((Chr1_Chr2_Row727 #= 368) and (Chr1_Chr2_Freq727 #= 1)) or 
	((Chr1_Chr2_Row727 #= 369) and (Chr1_Chr2_Freq727 #= 1)) or 
	((Chr1_Chr2_Row727 #= 370) and (Chr1_Chr2_Freq727 #= 1)) or 
	((Chr1_Chr2_Row727 #= 371) and (Chr1_Chr2_Freq727 #= 1)) or 
	((Chr1_Chr2_Row727 #= 372) and (Chr1_Chr2_Freq727 #= 1)) or 
	((Chr1_Chr2_Row727 #= 373) and (Chr1_Chr2_Freq727 #= 2)) or 
	((Chr1_Chr2_Row727 #= 374) and (Chr1_Chr2_Freq727 #= 2)) or 
	((Chr1_Chr2_Row727 #= 380) and (Chr1_Chr2_Freq727 #= 2)) or 
	((Chr1_Chr2_Row727 #= 381) and (Chr1_Chr2_Freq727 #= 2)) or 
	((Chr1_Chr2_Row727 #= 382) and (Chr1_Chr2_Freq727 #= 2)) or 
	((Chr1_Chr2_Row727 #= 383) and (Chr1_Chr2_Freq727 #= 1)) or 
	((Chr1_Chr2_Row727 #= 384) and (Chr1_Chr2_Freq727 #= 1)) or 
	((Chr1_Chr2_Row727 #= 385) and (Chr1_Chr2_Freq727 #= 1)) or 
	((Chr1_Chr2_Row727 #= 386) and (Chr1_Chr2_Freq727 #= 1)) or 
	((Chr1_Chr2_Row727 #= 387) and (Chr1_Chr2_Freq727 #= 1)) or 
	((Chr1_Chr2_Row727 #= 388) and (Chr1_Chr2_Freq727 #= 1)) or 
	((Chr1_Chr2_Row727 #= 389) and (Chr1_Chr2_Freq727 #= 1)) or 
	((Chr1_Chr2_Row727 #= 390) and (Chr1_Chr2_Freq727 #= 1)) or 
	((Chr1_Chr2_Row727 #= 391) and (Chr1_Chr2_Freq727 #= 1)) or 
	((Chr1_Chr2_Freq727 #= 0) and (Chr1_Chr2_Row727 #= 0)),

 	((Chr1_Chr2_Row728 #= 367) and (Chr1_Chr2_Freq728 #= 1)) or 
	((Chr1_Chr2_Row728 #= 369) and (Chr1_Chr2_Freq728 #= 1)) or 
	((Chr1_Chr2_Row728 #= 370) and (Chr1_Chr2_Freq728 #= 1)) or 
	((Chr1_Chr2_Row728 #= 371) and (Chr1_Chr2_Freq728 #= 1)) or 
	((Chr1_Chr2_Row728 #= 372) and (Chr1_Chr2_Freq728 #= 1)) or 
	((Chr1_Chr2_Row728 #= 373) and (Chr1_Chr2_Freq728 #= 1)) or 
	((Chr1_Chr2_Row728 #= 374) and (Chr1_Chr2_Freq728 #= 1)) or 
	((Chr1_Chr2_Row728 #= 380) and (Chr1_Chr2_Freq728 #= 1)) or 
	((Chr1_Chr2_Row728 #= 381) and (Chr1_Chr2_Freq728 #= 1)) or 
	((Chr1_Chr2_Row728 #= 382) and (Chr1_Chr2_Freq728 #= 1)) or 
	((Chr1_Chr2_Row728 #= 383) and (Chr1_Chr2_Freq728 #= 1)) or 
	((Chr1_Chr2_Row728 #= 384) and (Chr1_Chr2_Freq728 #= 1)) or 
	((Chr1_Chr2_Row728 #= 385) and (Chr1_Chr2_Freq728 #= 1)) or 
	((Chr1_Chr2_Row728 #= 386) and (Chr1_Chr2_Freq728 #= 1)) or 
	((Chr1_Chr2_Row728 #= 387) and (Chr1_Chr2_Freq728 #= 1)) or 
	((Chr1_Chr2_Freq728 #= 0) and (Chr1_Chr2_Row728 #= 0)),

 	((Chr1_Chr2_Row729 #= 366) and (Chr1_Chr2_Freq729 #= 1)) or 
	((Chr1_Chr2_Row729 #= 367) and (Chr1_Chr2_Freq729 #= 1)) or 
	((Chr1_Chr2_Row729 #= 368) and (Chr1_Chr2_Freq729 #= 1)) or 
	((Chr1_Chr2_Row729 #= 369) and (Chr1_Chr2_Freq729 #= 1)) or 
	((Chr1_Chr2_Row729 #= 370) and (Chr1_Chr2_Freq729 #= 1)) or 
	((Chr1_Chr2_Row729 #= 371) and (Chr1_Chr2_Freq729 #= 1)) or 
	((Chr1_Chr2_Row729 #= 372) and (Chr1_Chr2_Freq729 #= 1)) or 
	((Chr1_Chr2_Row729 #= 373) and (Chr1_Chr2_Freq729 #= 1)) or 
	((Chr1_Chr2_Row729 #= 374) and (Chr1_Chr2_Freq729 #= 1)) or 
	((Chr1_Chr2_Row729 #= 380) and (Chr1_Chr2_Freq729 #= 1)) or 
	((Chr1_Chr2_Row729 #= 381) and (Chr1_Chr2_Freq729 #= 2)) or 
	((Chr1_Chr2_Row729 #= 382) and (Chr1_Chr2_Freq729 #= 1)) or 
	((Chr1_Chr2_Row729 #= 383) and (Chr1_Chr2_Freq729 #= 1)) or 
	((Chr1_Chr2_Row729 #= 384) and (Chr1_Chr2_Freq729 #= 1)) or 
	((Chr1_Chr2_Row729 #= 385) and (Chr1_Chr2_Freq729 #= 1)) or 
	((Chr1_Chr2_Row729 #= 386) and (Chr1_Chr2_Freq729 #= 1)) or 
	((Chr1_Chr2_Row729 #= 387) and (Chr1_Chr2_Freq729 #= 1)) or 
	((Chr1_Chr2_Freq729 #= 0) and (Chr1_Chr2_Row729 #= 0)),

 	((Chr1_Chr2_Row730 #= 371) and (Chr1_Chr2_Freq730 #= 1)) or 
	((Chr1_Chr2_Row730 #= 372) and (Chr1_Chr2_Freq730 #= 1)) or 
	((Chr1_Chr2_Row730 #= 373) and (Chr1_Chr2_Freq730 #= 1)) or 
	((Chr1_Chr2_Row730 #= 374) and (Chr1_Chr2_Freq730 #= 1)) or 
	((Chr1_Chr2_Row730 #= 380) and (Chr1_Chr2_Freq730 #= 1)) or 
	((Chr1_Chr2_Row730 #= 381) and (Chr1_Chr2_Freq730 #= 1)) or 
	((Chr1_Chr2_Row730 #= 382) and (Chr1_Chr2_Freq730 #= 1)) or 
	((Chr1_Chr2_Row730 #= 383) and (Chr1_Chr2_Freq730 #= 1)) or 
	((Chr1_Chr2_Row730 #= 385) and (Chr1_Chr2_Freq730 #= 1)) or 
	((Chr1_Chr2_Row730 #= 386) and (Chr1_Chr2_Freq730 #= 1)) or 
	((Chr1_Chr2_Row730 #= 387) and (Chr1_Chr2_Freq730 #= 1)) or 
	((Chr1_Chr2_Freq730 #= 0) and (Chr1_Chr2_Row730 #= 0)),

 	((Chr1_Chr2_Row731 #= 372) and (Chr1_Chr2_Freq731 #= 1)) or 
	((Chr1_Chr2_Row731 #= 373) and (Chr1_Chr2_Freq731 #= 1)) or 
	((Chr1_Chr2_Row731 #= 374) and (Chr1_Chr2_Freq731 #= 1)) or 
	((Chr1_Chr2_Row731 #= 380) and (Chr1_Chr2_Freq731 #= 1)) or 
	((Chr1_Chr2_Row731 #= 381) and (Chr1_Chr2_Freq731 #= 1)) or 
	((Chr1_Chr2_Row731 #= 382) and (Chr1_Chr2_Freq731 #= 1)) or 
	((Chr1_Chr2_Row731 #= 383) and (Chr1_Chr2_Freq731 #= 1)) or 
	((Chr1_Chr2_Row731 #= 385) and (Chr1_Chr2_Freq731 #= 1)) or 
	((Chr1_Chr2_Row731 #= 386) and (Chr1_Chr2_Freq731 #= 1)) or 
	((Chr1_Chr2_Freq731 #= 0) and (Chr1_Chr2_Row731 #= 0)),

 	((Chr1_Chr2_Row732 #= 372) and (Chr1_Chr2_Freq732 #= 1)) or 
	((Chr1_Chr2_Row732 #= 373) and (Chr1_Chr2_Freq732 #= 1)) or 
	((Chr1_Chr2_Row732 #= 374) and (Chr1_Chr2_Freq732 #= 1)) or 
	((Chr1_Chr2_Row732 #= 381) and (Chr1_Chr2_Freq732 #= 1)) or 
	((Chr1_Chr2_Row732 #= 382) and (Chr1_Chr2_Freq732 #= 1)) or 
	((Chr1_Chr2_Freq732 #= 0) and (Chr1_Chr2_Row732 #= 0)),

 	((Chr1_Chr2_Row733 #= 373) and (Chr1_Chr2_Freq733 #= 1)) or 
	((Chr1_Chr2_Row733 #= 381) and (Chr1_Chr2_Freq733 #= 1)) or 
	((Chr1_Chr2_Row733 #= 382) and (Chr1_Chr2_Freq733 #= 1)) or 
	((Chr1_Chr2_Row733 #= 385) and (Chr1_Chr2_Freq733 #= 1)) or 
	((Chr1_Chr2_Freq733 #= 0) and (Chr1_Chr2_Row733 #= 0)),

 	((Chr1_Chr2_Row734 #= 381) and (Chr1_Chr2_Freq734 #= 1)) or 
	((Chr1_Chr2_Row734 #= 387) and (Chr1_Chr2_Freq734 #= 1)) or 
	((Chr1_Chr2_Row734 #= 388) and (Chr1_Chr2_Freq734 #= 1)) or 
	((Chr1_Chr2_Freq734 #= 0) and (Chr1_Chr2_Row734 #= 0)),

 	((Chr1_Chr2_Row735 #= 373) and (Chr1_Chr2_Freq735 #= 1)) or 
	((Chr1_Chr2_Row735 #= 381) and (Chr1_Chr2_Freq735 #= 1)) or 
	((Chr1_Chr2_Freq735 #= 0) and (Chr1_Chr2_Row735 #= 0)),

 	((Chr1_Chr2_Row736 #= 373) and (Chr1_Chr2_Freq736 #= 1)) or 
	((Chr1_Chr2_Freq736 #= 0) and (Chr1_Chr2_Row736 #= 0)),

 	((Chr1_Chr2_Row772 #= 555) and (Chr1_Chr2_Freq772 #= 1)) or 
	((Chr1_Chr2_Row772 #= 556) and (Chr1_Chr2_Freq772 #= 1)) or 
	((Chr1_Chr2_Row772 #= 557) and (Chr1_Chr2_Freq772 #= 1)) or 
	((Chr1_Chr2_Freq772 #= 0) and (Chr1_Chr2_Row772 #= 0)),

 	((Chr1_Chr2_Row995 #= 11) and (Chr1_Chr2_Freq995 #= 1)) or 
	((Chr1_Chr2_Row995 #= 557) and (Chr1_Chr2_Freq995 #= 1)) or 
	((Chr1_Chr2_Row995 #= 6) and (Chr1_Chr2_Freq995 #= 1)) or 
	((Chr1_Chr2_Row995 #= 7) and (Chr1_Chr2_Freq995 #= 1)) or 
	((Chr1_Chr2_Freq995 #= 0) and (Chr1_Chr2_Row995 #= 0)),

 	((Chr1_Chr2_Row996 #= 10) and (Chr1_Chr2_Freq996 #= 1)) or 
	((Chr1_Chr2_Row996 #= 11) and (Chr1_Chr2_Freq996 #= 1)) or 
	((Chr1_Chr2_Row996 #= 6) and (Chr1_Chr2_Freq996 #= 1)) or 
	((Chr1_Chr2_Row996 #= 7) and (Chr1_Chr2_Freq996 #= 1)) or 
	((Chr1_Chr2_Row996 #= 8) and (Chr1_Chr2_Freq996 #= 1)) or 
	((Chr1_Chr2_Row996 #= 9) and (Chr1_Chr2_Freq996 #= 1)) or 
	((Chr1_Chr2_Freq996 #= 0) and (Chr1_Chr2_Row996 #= 0)),

 	((Chr1_Chr2_Row997 #= 10) and (Chr1_Chr2_Freq997 #= 1)) or 
	((Chr1_Chr2_Row997 #= 11) and (Chr1_Chr2_Freq997 #= 1)) or 
	((Chr1_Chr2_Row997 #= 12) and (Chr1_Chr2_Freq997 #= 1)) or 
	((Chr1_Chr2_Row997 #= 13) and (Chr1_Chr2_Freq997 #= 1)) or 
	((Chr1_Chr2_Row997 #= 14) and (Chr1_Chr2_Freq997 #= 1)) or 
	((Chr1_Chr2_Row997 #= 6) and (Chr1_Chr2_Freq997 #= 1)) or 
	((Chr1_Chr2_Row997 #= 7) and (Chr1_Chr2_Freq997 #= 1)) or 
	((Chr1_Chr2_Row997 #= 8) and (Chr1_Chr2_Freq997 #= 1)) or 
	((Chr1_Chr2_Row997 #= 9) and (Chr1_Chr2_Freq997 #= 1)) or 
	((Chr1_Chr2_Freq997 #= 0) and (Chr1_Chr2_Row997 #= 0)),

 	((Chr1_Chr2_Row998 #= 10) and (Chr1_Chr2_Freq998 #= 1)) or 
	((Chr1_Chr2_Row998 #= 11) and (Chr1_Chr2_Freq998 #= 1)) or 
	((Chr1_Chr2_Row998 #= 12) and (Chr1_Chr2_Freq998 #= 1)) or 
	((Chr1_Chr2_Row998 #= 13) and (Chr1_Chr2_Freq998 #= 1)) or 
	((Chr1_Chr2_Row998 #= 14) and (Chr1_Chr2_Freq998 #= 1)) or 
	((Chr1_Chr2_Row998 #= 15) and (Chr1_Chr2_Freq998 #= 1)) or 
	((Chr1_Chr2_Row998 #= 549) and (Chr1_Chr2_Freq998 #= 1)) or 
	((Chr1_Chr2_Row998 #= 550) and (Chr1_Chr2_Freq998 #= 1)) or 
	((Chr1_Chr2_Row998 #= 556) and (Chr1_Chr2_Freq998 #= 1)) or 
	((Chr1_Chr2_Row998 #= 557) and (Chr1_Chr2_Freq998 #= 1)) or 
	((Chr1_Chr2_Row998 #= 6) and (Chr1_Chr2_Freq998 #= 1)) or 
	((Chr1_Chr2_Row998 #= 7) and (Chr1_Chr2_Freq998 #= 1)) or 
	((Chr1_Chr2_Row998 #= 8) and (Chr1_Chr2_Freq998 #= 1)) or 
	((Chr1_Chr2_Row998 #= 9) and (Chr1_Chr2_Freq998 #= 1)) or 
	((Chr1_Chr2_Freq998 #= 0) and (Chr1_Chr2_Row998 #= 0)),

 	((Chr1_Chr2_Row999 #= 10) and (Chr1_Chr2_Freq999 #= 2)) or 
	((Chr1_Chr2_Row999 #= 11) and (Chr1_Chr2_Freq999 #= 1)) or 
	((Chr1_Chr2_Row999 #= 12) and (Chr1_Chr2_Freq999 #= 1)) or 
	((Chr1_Chr2_Row999 #= 13) and (Chr1_Chr2_Freq999 #= 1)) or 
	((Chr1_Chr2_Row999 #= 14) and (Chr1_Chr2_Freq999 #= 1)) or 
	((Chr1_Chr2_Row999 #= 15) and (Chr1_Chr2_Freq999 #= 1)) or 
	((Chr1_Chr2_Row999 #= 546) and (Chr1_Chr2_Freq999 #= 1)) or 
	((Chr1_Chr2_Row999 #= 547) and (Chr1_Chr2_Freq999 #= 1)) or 
	((Chr1_Chr2_Row999 #= 548) and (Chr1_Chr2_Freq999 #= 1)) or 
	((Chr1_Chr2_Row999 #= 549) and (Chr1_Chr2_Freq999 #= 1)) or 
	((Chr1_Chr2_Row999 #= 550) and (Chr1_Chr2_Freq999 #= 1)) or 
	((Chr1_Chr2_Row999 #= 551) and (Chr1_Chr2_Freq999 #= 1)) or 
	((Chr1_Chr2_Row999 #= 552) and (Chr1_Chr2_Freq999 #= 1)) or 
	((Chr1_Chr2_Row999 #= 553) and (Chr1_Chr2_Freq999 #= 1)) or 
	((Chr1_Chr2_Row999 #= 554) and (Chr1_Chr2_Freq999 #= 1)) or 
	((Chr1_Chr2_Row999 #= 555) and (Chr1_Chr2_Freq999 #= 1)) or 
	((Chr1_Chr2_Row999 #= 556) and (Chr1_Chr2_Freq999 #= 1)) or 
	((Chr1_Chr2_Row999 #= 557) and (Chr1_Chr2_Freq999 #= 1)) or 
	((Chr1_Chr2_Row999 #= 6) and (Chr1_Chr2_Freq999 #= 2)) or 
	((Chr1_Chr2_Row999 #= 7) and (Chr1_Chr2_Freq999 #= 2)) or 
	((Chr1_Chr2_Row999 #= 8) and (Chr1_Chr2_Freq999 #= 2)) or 
	((Chr1_Chr2_Row999 #= 9) and (Chr1_Chr2_Freq999 #= 1)) or 
	((Chr1_Chr2_Freq999 #= 0) and (Chr1_Chr2_Row999 #= 0)),

 	((Chr1_Chr2_Row1000 #= 10) and (Chr1_Chr2_Freq1000 #= 2)) or 
	((Chr1_Chr2_Row1000 #= 11) and (Chr1_Chr2_Freq1000 #= 2)) or 
	((Chr1_Chr2_Row1000 #= 12) and (Chr1_Chr2_Freq1000 #= 1)) or 
	((Chr1_Chr2_Row1000 #= 13) and (Chr1_Chr2_Freq1000 #= 1)) or 
	((Chr1_Chr2_Row1000 #= 14) and (Chr1_Chr2_Freq1000 #= 1)) or 
	((Chr1_Chr2_Row1000 #= 15) and (Chr1_Chr2_Freq1000 #= 1)) or 
	((Chr1_Chr2_Row1000 #= 544) and (Chr1_Chr2_Freq1000 #= 1)) or 
	((Chr1_Chr2_Row1000 #= 545) and (Chr1_Chr2_Freq1000 #= 1)) or 
	((Chr1_Chr2_Row1000 #= 546) and (Chr1_Chr2_Freq1000 #= 1)) or 
	((Chr1_Chr2_Row1000 #= 547) and (Chr1_Chr2_Freq1000 #= 1)) or 
	((Chr1_Chr2_Row1000 #= 548) and (Chr1_Chr2_Freq1000 #= 1)) or 
	((Chr1_Chr2_Row1000 #= 549) and (Chr1_Chr2_Freq1000 #= 1)) or 
	((Chr1_Chr2_Row1000 #= 550) and (Chr1_Chr2_Freq1000 #= 2)) or 
	((Chr1_Chr2_Row1000 #= 551) and (Chr1_Chr2_Freq1000 #= 1)) or 
	((Chr1_Chr2_Row1000 #= 552) and (Chr1_Chr2_Freq1000 #= 1)) or 
	((Chr1_Chr2_Row1000 #= 553) and (Chr1_Chr2_Freq1000 #= 1)) or 
	((Chr1_Chr2_Row1000 #= 554) and (Chr1_Chr2_Freq1000 #= 1)) or 
	((Chr1_Chr2_Row1000 #= 555) and (Chr1_Chr2_Freq1000 #= 2)) or 
	((Chr1_Chr2_Row1000 #= 556) and (Chr1_Chr2_Freq1000 #= 1)) or 
	((Chr1_Chr2_Row1000 #= 557) and (Chr1_Chr2_Freq1000 #= 1)) or 
	((Chr1_Chr2_Row1000 #= 6) and (Chr1_Chr2_Freq1000 #= 3)) or 
	((Chr1_Chr2_Row1000 #= 7) and (Chr1_Chr2_Freq1000 #= 2)) or 
	((Chr1_Chr2_Row1000 #= 8) and (Chr1_Chr2_Freq1000 #= 2)) or 
	((Chr1_Chr2_Row1000 #= 9) and (Chr1_Chr2_Freq1000 #= 2)) or 
	((Chr1_Chr2_Freq1000 #= 0) and (Chr1_Chr2_Row1000 #= 0)),

 	((Chr1_Chr2_Row1001 #= 10) and (Chr1_Chr2_Freq1001 #= 2)) or 
	((Chr1_Chr2_Row1001 #= 11) and (Chr1_Chr2_Freq1001 #= 1)) or 
	((Chr1_Chr2_Row1001 #= 12) and (Chr1_Chr2_Freq1001 #= 1)) or 
	((Chr1_Chr2_Row1001 #= 13) and (Chr1_Chr2_Freq1001 #= 1)) or 
	((Chr1_Chr2_Row1001 #= 14) and (Chr1_Chr2_Freq1001 #= 1)) or 
	((Chr1_Chr2_Row1001 #= 15) and (Chr1_Chr2_Freq1001 #= 1)) or 
	((Chr1_Chr2_Row1001 #= 16) and (Chr1_Chr2_Freq1001 #= 1)) or 
	((Chr1_Chr2_Row1001 #= 545) and (Chr1_Chr2_Freq1001 #= 1)) or 
	((Chr1_Chr2_Row1001 #= 547) and (Chr1_Chr2_Freq1001 #= 1)) or 
	((Chr1_Chr2_Row1001 #= 548) and (Chr1_Chr2_Freq1001 #= 1)) or 
	((Chr1_Chr2_Row1001 #= 549) and (Chr1_Chr2_Freq1001 #= 1)) or 
	((Chr1_Chr2_Row1001 #= 550) and (Chr1_Chr2_Freq1001 #= 2)) or 
	((Chr1_Chr2_Row1001 #= 551) and (Chr1_Chr2_Freq1001 #= 2)) or 
	((Chr1_Chr2_Row1001 #= 552) and (Chr1_Chr2_Freq1001 #= 2)) or 
	((Chr1_Chr2_Row1001 #= 553) and (Chr1_Chr2_Freq1001 #= 2)) or 
	((Chr1_Chr2_Row1001 #= 554) and (Chr1_Chr2_Freq1001 #= 2)) or 
	((Chr1_Chr2_Row1001 #= 555) and (Chr1_Chr2_Freq1001 #= 1)) or 
	((Chr1_Chr2_Row1001 #= 556) and (Chr1_Chr2_Freq1001 #= 1)) or 
	((Chr1_Chr2_Row1001 #= 557) and (Chr1_Chr2_Freq1001 #= 1)) or 
	((Chr1_Chr2_Row1001 #= 6) and (Chr1_Chr2_Freq1001 #= 2)) or 
	((Chr1_Chr2_Row1001 #= 7) and (Chr1_Chr2_Freq1001 #= 2)) or 
	((Chr1_Chr2_Row1001 #= 8) and (Chr1_Chr2_Freq1001 #= 2)) or 
	((Chr1_Chr2_Row1001 #= 9) and (Chr1_Chr2_Freq1001 #= 2)) or 
	((Chr1_Chr2_Freq1001 #= 0) and (Chr1_Chr2_Row1001 #= 0)),

 	((Chr1_Chr2_Row1002 #= 10) and (Chr1_Chr2_Freq1002 #= 2)) or 
	((Chr1_Chr2_Row1002 #= 11) and (Chr1_Chr2_Freq1002 #= 2)) or 
	((Chr1_Chr2_Row1002 #= 12) and (Chr1_Chr2_Freq1002 #= 1)) or 
	((Chr1_Chr2_Row1002 #= 13) and (Chr1_Chr2_Freq1002 #= 1)) or 
	((Chr1_Chr2_Row1002 #= 14) and (Chr1_Chr2_Freq1002 #= 1)) or 
	((Chr1_Chr2_Row1002 #= 15) and (Chr1_Chr2_Freq1002 #= 1)) or 
	((Chr1_Chr2_Row1002 #= 16) and (Chr1_Chr2_Freq1002 #= 1)) or 
	((Chr1_Chr2_Row1002 #= 546) and (Chr1_Chr2_Freq1002 #= 1)) or 
	((Chr1_Chr2_Row1002 #= 547) and (Chr1_Chr2_Freq1002 #= 1)) or 
	((Chr1_Chr2_Row1002 #= 548) and (Chr1_Chr2_Freq1002 #= 1)) or 
	((Chr1_Chr2_Row1002 #= 549) and (Chr1_Chr2_Freq1002 #= 1)) or 
	((Chr1_Chr2_Row1002 #= 550) and (Chr1_Chr2_Freq1002 #= 1)) or 
	((Chr1_Chr2_Row1002 #= 551) and (Chr1_Chr2_Freq1002 #= 1)) or 
	((Chr1_Chr2_Row1002 #= 552) and (Chr1_Chr2_Freq1002 #= 2)) or 
	((Chr1_Chr2_Row1002 #= 553) and (Chr1_Chr2_Freq1002 #= 1)) or 
	((Chr1_Chr2_Row1002 #= 554) and (Chr1_Chr2_Freq1002 #= 1)) or 
	((Chr1_Chr2_Row1002 #= 555) and (Chr1_Chr2_Freq1002 #= 1)) or 
	((Chr1_Chr2_Row1002 #= 556) and (Chr1_Chr2_Freq1002 #= 2)) or 
	((Chr1_Chr2_Row1002 #= 557) and (Chr1_Chr2_Freq1002 #= 2)) or 
	((Chr1_Chr2_Row1002 #= 6) and (Chr1_Chr2_Freq1002 #= 3)) or 
	((Chr1_Chr2_Row1002 #= 7) and (Chr1_Chr2_Freq1002 #= 2)) or 
	((Chr1_Chr2_Row1002 #= 8) and (Chr1_Chr2_Freq1002 #= 2)) or 
	((Chr1_Chr2_Row1002 #= 9) and (Chr1_Chr2_Freq1002 #= 2)) or 
	((Chr1_Chr2_Freq1002 #= 0) and (Chr1_Chr2_Row1002 #= 0)),

 	((Chr1_Chr2_Row1003 #= 10) and (Chr1_Chr2_Freq1003 #= 2)) or 
	((Chr1_Chr2_Row1003 #= 11) and (Chr1_Chr2_Freq1003 #= 2)) or 
	((Chr1_Chr2_Row1003 #= 12) and (Chr1_Chr2_Freq1003 #= 1)) or 
	((Chr1_Chr2_Row1003 #= 13) and (Chr1_Chr2_Freq1003 #= 1)) or 
	((Chr1_Chr2_Row1003 #= 14) and (Chr1_Chr2_Freq1003 #= 1)) or 
	((Chr1_Chr2_Row1003 #= 15) and (Chr1_Chr2_Freq1003 #= 1)) or 
	((Chr1_Chr2_Row1003 #= 16) and (Chr1_Chr2_Freq1003 #= 1)) or 
	((Chr1_Chr2_Row1003 #= 17) and (Chr1_Chr2_Freq1003 #= 1)) or 
	((Chr1_Chr2_Row1003 #= 20) and (Chr1_Chr2_Freq1003 #= 1)) or 
	((Chr1_Chr2_Row1003 #= 545) and (Chr1_Chr2_Freq1003 #= 1)) or 
	((Chr1_Chr2_Row1003 #= 546) and (Chr1_Chr2_Freq1003 #= 1)) or 
	((Chr1_Chr2_Row1003 #= 547) and (Chr1_Chr2_Freq1003 #= 1)) or 
	((Chr1_Chr2_Row1003 #= 548) and (Chr1_Chr2_Freq1003 #= 1)) or 
	((Chr1_Chr2_Row1003 #= 549) and (Chr1_Chr2_Freq1003 #= 1)) or 
	((Chr1_Chr2_Row1003 #= 550) and (Chr1_Chr2_Freq1003 #= 1)) or 
	((Chr1_Chr2_Row1003 #= 551) and (Chr1_Chr2_Freq1003 #= 1)) or 
	((Chr1_Chr2_Row1003 #= 552) and (Chr1_Chr2_Freq1003 #= 2)) or 
	((Chr1_Chr2_Row1003 #= 553) and (Chr1_Chr2_Freq1003 #= 2)) or 
	((Chr1_Chr2_Row1003 #= 554) and (Chr1_Chr2_Freq1003 #= 1)) or 
	((Chr1_Chr2_Row1003 #= 555) and (Chr1_Chr2_Freq1003 #= 2)) or 
	((Chr1_Chr2_Row1003 #= 556) and (Chr1_Chr2_Freq1003 #= 2)) or 
	((Chr1_Chr2_Row1003 #= 557) and (Chr1_Chr2_Freq1003 #= 2)) or 
	((Chr1_Chr2_Row1003 #= 6) and (Chr1_Chr2_Freq1003 #= 3)) or 
	((Chr1_Chr2_Row1003 #= 7) and (Chr1_Chr2_Freq1003 #= 2)) or 
	((Chr1_Chr2_Row1003 #= 8) and (Chr1_Chr2_Freq1003 #= 2)) or 
	((Chr1_Chr2_Row1003 #= 9) and (Chr1_Chr2_Freq1003 #= 2)) or 
	((Chr1_Chr2_Freq1003 #= 0) and (Chr1_Chr2_Row1003 #= 0)),

 	((Chr1_Chr2_Row1004 #= 10) and (Chr1_Chr2_Freq1004 #= 3)) or 
	((Chr1_Chr2_Row1004 #= 11) and (Chr1_Chr2_Freq1004 #= 2)) or 
	((Chr1_Chr2_Row1004 #= 12) and (Chr1_Chr2_Freq1004 #= 2)) or 
	((Chr1_Chr2_Row1004 #= 13) and (Chr1_Chr2_Freq1004 #= 1)) or 
	((Chr1_Chr2_Row1004 #= 14) and (Chr1_Chr2_Freq1004 #= 1)) or 
	((Chr1_Chr2_Row1004 #= 15) and (Chr1_Chr2_Freq1004 #= 1)) or 
	((Chr1_Chr2_Row1004 #= 16) and (Chr1_Chr2_Freq1004 #= 1)) or 
	((Chr1_Chr2_Row1004 #= 543) and (Chr1_Chr2_Freq1004 #= 1)) or 
	((Chr1_Chr2_Row1004 #= 545) and (Chr1_Chr2_Freq1004 #= 1)) or 
	((Chr1_Chr2_Row1004 #= 546) and (Chr1_Chr2_Freq1004 #= 1)) or 
	((Chr1_Chr2_Row1004 #= 547) and (Chr1_Chr2_Freq1004 #= 1)) or 
	((Chr1_Chr2_Row1004 #= 548) and (Chr1_Chr2_Freq1004 #= 1)) or 
	((Chr1_Chr2_Row1004 #= 549) and (Chr1_Chr2_Freq1004 #= 2)) or 
	((Chr1_Chr2_Row1004 #= 550) and (Chr1_Chr2_Freq1004 #= 2)) or 
	((Chr1_Chr2_Row1004 #= 551) and (Chr1_Chr2_Freq1004 #= 4)) or 
	((Chr1_Chr2_Row1004 #= 552) and (Chr1_Chr2_Freq1004 #= 4)) or 
	((Chr1_Chr2_Row1004 #= 553) and (Chr1_Chr2_Freq1004 #= 4)) or 
	((Chr1_Chr2_Row1004 #= 554) and (Chr1_Chr2_Freq1004 #= 3)) or 
	((Chr1_Chr2_Row1004 #= 555) and (Chr1_Chr2_Freq1004 #= 3)) or 
	((Chr1_Chr2_Row1004 #= 556) and (Chr1_Chr2_Freq1004 #= 3)) or 
	((Chr1_Chr2_Row1004 #= 557) and (Chr1_Chr2_Freq1004 #= 2)) or 
	((Chr1_Chr2_Row1004 #= 6) and (Chr1_Chr2_Freq1004 #= 4)) or 
	((Chr1_Chr2_Row1004 #= 7) and (Chr1_Chr2_Freq1004 #= 5)) or 
	((Chr1_Chr2_Row1004 #= 8) and (Chr1_Chr2_Freq1004 #= 6)) or 
	((Chr1_Chr2_Row1004 #= 9) and (Chr1_Chr2_Freq1004 #= 5)) or 
	((Chr1_Chr2_Freq1004 #= 0) and (Chr1_Chr2_Row1004 #= 0)),

 	((Chr1_Chr2_Row1005 #= 10) and (Chr1_Chr2_Freq1005 #= 4)) or 
	((Chr1_Chr2_Row1005 #= 11) and (Chr1_Chr2_Freq1005 #= 2)) or 
	((Chr1_Chr2_Row1005 #= 12) and (Chr1_Chr2_Freq1005 #= 2)) or 
	((Chr1_Chr2_Row1005 #= 13) and (Chr1_Chr2_Freq1005 #= 2)) or 
	((Chr1_Chr2_Row1005 #= 14) and (Chr1_Chr2_Freq1005 #= 2)) or 
	((Chr1_Chr2_Row1005 #= 15) and (Chr1_Chr2_Freq1005 #= 1)) or 
	((Chr1_Chr2_Row1005 #= 16) and (Chr1_Chr2_Freq1005 #= 1)) or 
	((Chr1_Chr2_Row1005 #= 17) and (Chr1_Chr2_Freq1005 #= 1)) or 
	((Chr1_Chr2_Row1005 #= 544) and (Chr1_Chr2_Freq1005 #= 1)) or 
	((Chr1_Chr2_Row1005 #= 545) and (Chr1_Chr2_Freq1005 #= 1)) or 
	((Chr1_Chr2_Row1005 #= 546) and (Chr1_Chr2_Freq1005 #= 1)) or 
	((Chr1_Chr2_Row1005 #= 547) and (Chr1_Chr2_Freq1005 #= 1)) or 
	((Chr1_Chr2_Row1005 #= 548) and (Chr1_Chr2_Freq1005 #= 2)) or 
	((Chr1_Chr2_Row1005 #= 549) and (Chr1_Chr2_Freq1005 #= 2)) or 
	((Chr1_Chr2_Row1005 #= 550) and (Chr1_Chr2_Freq1005 #= 2)) or 
	((Chr1_Chr2_Row1005 #= 551) and (Chr1_Chr2_Freq1005 #= 4)) or 
	((Chr1_Chr2_Row1005 #= 552) and (Chr1_Chr2_Freq1005 #= 4)) or 
	((Chr1_Chr2_Row1005 #= 553) and (Chr1_Chr2_Freq1005 #= 4)) or 
	((Chr1_Chr2_Row1005 #= 554) and (Chr1_Chr2_Freq1005 #= 4)) or 
	((Chr1_Chr2_Row1005 #= 555) and (Chr1_Chr2_Freq1005 #= 3)) or 
	((Chr1_Chr2_Row1005 #= 556) and (Chr1_Chr2_Freq1005 #= 3)) or 
	((Chr1_Chr2_Row1005 #= 557) and (Chr1_Chr2_Freq1005 #= 3)) or 
	((Chr1_Chr2_Row1005 #= 6) and (Chr1_Chr2_Freq1005 #= 5)) or 
	((Chr1_Chr2_Row1005 #= 7) and (Chr1_Chr2_Freq1005 #= 5)) or 
	((Chr1_Chr2_Row1005 #= 8) and (Chr1_Chr2_Freq1005 #= 6)) or 
	((Chr1_Chr2_Row1005 #= 9) and (Chr1_Chr2_Freq1005 #= 6)) or 
	((Chr1_Chr2_Freq1005 #= 0) and (Chr1_Chr2_Row1005 #= 0)),

 	((Chr1_Chr2_Row1006 #= 10) and (Chr1_Chr2_Freq1006 #= 5)) or 
	((Chr1_Chr2_Row1006 #= 11) and (Chr1_Chr2_Freq1006 #= 2)) or 
	((Chr1_Chr2_Row1006 #= 12) and (Chr1_Chr2_Freq1006 #= 2)) or 
	((Chr1_Chr2_Row1006 #= 13) and (Chr1_Chr2_Freq1006 #= 2)) or 
	((Chr1_Chr2_Row1006 #= 14) and (Chr1_Chr2_Freq1006 #= 2)) or 
	((Chr1_Chr2_Row1006 #= 15) and (Chr1_Chr2_Freq1006 #= 1)) or 
	((Chr1_Chr2_Row1006 #= 16) and (Chr1_Chr2_Freq1006 #= 1)) or 
	((Chr1_Chr2_Row1006 #= 17) and (Chr1_Chr2_Freq1006 #= 1)) or 
	((Chr1_Chr2_Row1006 #= 545) and (Chr1_Chr2_Freq1006 #= 1)) or 
	((Chr1_Chr2_Row1006 #= 546) and (Chr1_Chr2_Freq1006 #= 1)) or 
	((Chr1_Chr2_Row1006 #= 547) and (Chr1_Chr2_Freq1006 #= 1)) or 
	((Chr1_Chr2_Row1006 #= 548) and (Chr1_Chr2_Freq1006 #= 2)) or 
	((Chr1_Chr2_Row1006 #= 549) and (Chr1_Chr2_Freq1006 #= 1)) or 
	((Chr1_Chr2_Row1006 #= 550) and (Chr1_Chr2_Freq1006 #= 3)) or 
	((Chr1_Chr2_Row1006 #= 551) and (Chr1_Chr2_Freq1006 #= 4)) or 
	((Chr1_Chr2_Row1006 #= 552) and (Chr1_Chr2_Freq1006 #= 4)) or 
	((Chr1_Chr2_Row1006 #= 553) and (Chr1_Chr2_Freq1006 #= 5)) or 
	((Chr1_Chr2_Row1006 #= 554) and (Chr1_Chr2_Freq1006 #= 4)) or 
	((Chr1_Chr2_Row1006 #= 555) and (Chr1_Chr2_Freq1006 #= 4)) or 
	((Chr1_Chr2_Row1006 #= 556) and (Chr1_Chr2_Freq1006 #= 3)) or 
	((Chr1_Chr2_Row1006 #= 557) and (Chr1_Chr2_Freq1006 #= 3)) or 
	((Chr1_Chr2_Row1006 #= 6) and (Chr1_Chr2_Freq1006 #= 5)) or 
	((Chr1_Chr2_Row1006 #= 7) and (Chr1_Chr2_Freq1006 #= 5)) or 
	((Chr1_Chr2_Row1006 #= 8) and (Chr1_Chr2_Freq1006 #= 6)) or 
	((Chr1_Chr2_Row1006 #= 9) and (Chr1_Chr2_Freq1006 #= 7)) or 
	((Chr1_Chr2_Freq1006 #= 0) and (Chr1_Chr2_Row1006 #= 0)),

 	((Chr1_Chr2_Row1007 #= 10) and (Chr1_Chr2_Freq1007 #= 5)) or 
	((Chr1_Chr2_Row1007 #= 11) and (Chr1_Chr2_Freq1007 #= 2)) or 
	((Chr1_Chr2_Row1007 #= 12) and (Chr1_Chr2_Freq1007 #= 2)) or 
	((Chr1_Chr2_Row1007 #= 13) and (Chr1_Chr2_Freq1007 #= 2)) or 
	((Chr1_Chr2_Row1007 #= 14) and (Chr1_Chr2_Freq1007 #= 2)) or 
	((Chr1_Chr2_Row1007 #= 15) and (Chr1_Chr2_Freq1007 #= 1)) or 
	((Chr1_Chr2_Row1007 #= 16) and (Chr1_Chr2_Freq1007 #= 1)) or 
	((Chr1_Chr2_Row1007 #= 17) and (Chr1_Chr2_Freq1007 #= 1)) or 
	((Chr1_Chr2_Row1007 #= 544) and (Chr1_Chr2_Freq1007 #= 1)) or 
	((Chr1_Chr2_Row1007 #= 545) and (Chr1_Chr2_Freq1007 #= 1)) or 
	((Chr1_Chr2_Row1007 #= 546) and (Chr1_Chr2_Freq1007 #= 1)) or 
	((Chr1_Chr2_Row1007 #= 547) and (Chr1_Chr2_Freq1007 #= 2)) or 
	((Chr1_Chr2_Row1007 #= 548) and (Chr1_Chr2_Freq1007 #= 2)) or 
	((Chr1_Chr2_Row1007 #= 549) and (Chr1_Chr2_Freq1007 #= 2)) or 
	((Chr1_Chr2_Row1007 #= 550) and (Chr1_Chr2_Freq1007 #= 3)) or 
	((Chr1_Chr2_Row1007 #= 551) and (Chr1_Chr2_Freq1007 #= 4)) or 
	((Chr1_Chr2_Row1007 #= 552) and (Chr1_Chr2_Freq1007 #= 4)) or 
	((Chr1_Chr2_Row1007 #= 553) and (Chr1_Chr2_Freq1007 #= 5)) or 
	((Chr1_Chr2_Row1007 #= 554) and (Chr1_Chr2_Freq1007 #= 4)) or 
	((Chr1_Chr2_Row1007 #= 555) and (Chr1_Chr2_Freq1007 #= 4)) or 
	((Chr1_Chr2_Row1007 #= 556) and (Chr1_Chr2_Freq1007 #= 4)) or 
	((Chr1_Chr2_Row1007 #= 557) and (Chr1_Chr2_Freq1007 #= 4)) or 
	((Chr1_Chr2_Row1007 #= 6) and (Chr1_Chr2_Freq1007 #= 5)) or 
	((Chr1_Chr2_Row1007 #= 7) and (Chr1_Chr2_Freq1007 #= 6)) or 
	((Chr1_Chr2_Row1007 #= 8) and (Chr1_Chr2_Freq1007 #= 7)) or 
	((Chr1_Chr2_Row1007 #= 9) and (Chr1_Chr2_Freq1007 #= 6)) or 
	((Chr1_Chr2_Freq1007 #= 0) and (Chr1_Chr2_Row1007 #= 0)),

 	((Chr1_Chr2_Row1008 #= 10) and (Chr1_Chr2_Freq1008 #= 4)) or 
	((Chr1_Chr2_Row1008 #= 11) and (Chr1_Chr2_Freq1008 #= 2)) or 
	((Chr1_Chr2_Row1008 #= 12) and (Chr1_Chr2_Freq1008 #= 2)) or 
	((Chr1_Chr2_Row1008 #= 13) and (Chr1_Chr2_Freq1008 #= 2)) or 
	((Chr1_Chr2_Row1008 #= 14) and (Chr1_Chr2_Freq1008 #= 2)) or 
	((Chr1_Chr2_Row1008 #= 15) and (Chr1_Chr2_Freq1008 #= 1)) or 
	((Chr1_Chr2_Row1008 #= 16) and (Chr1_Chr2_Freq1008 #= 1)) or 
	((Chr1_Chr2_Row1008 #= 17) and (Chr1_Chr2_Freq1008 #= 1)) or 
	((Chr1_Chr2_Row1008 #= 545) and (Chr1_Chr2_Freq1008 #= 1)) or 
	((Chr1_Chr2_Row1008 #= 546) and (Chr1_Chr2_Freq1008 #= 1)) or 
	((Chr1_Chr2_Row1008 #= 547) and (Chr1_Chr2_Freq1008 #= 2)) or 
	((Chr1_Chr2_Row1008 #= 548) and (Chr1_Chr2_Freq1008 #= 2)) or 
	((Chr1_Chr2_Row1008 #= 549) and (Chr1_Chr2_Freq1008 #= 2)) or 
	((Chr1_Chr2_Row1008 #= 550) and (Chr1_Chr2_Freq1008 #= 3)) or 
	((Chr1_Chr2_Row1008 #= 551) and (Chr1_Chr2_Freq1008 #= 5)) or 
	((Chr1_Chr2_Row1008 #= 552) and (Chr1_Chr2_Freq1008 #= 5)) or 
	((Chr1_Chr2_Row1008 #= 553) and (Chr1_Chr2_Freq1008 #= 5)) or 
	((Chr1_Chr2_Row1008 #= 554) and (Chr1_Chr2_Freq1008 #= 5)) or 
	((Chr1_Chr2_Row1008 #= 555) and (Chr1_Chr2_Freq1008 #= 4)) or 
	((Chr1_Chr2_Row1008 #= 556) and (Chr1_Chr2_Freq1008 #= 4)) or 
	((Chr1_Chr2_Row1008 #= 557) and (Chr1_Chr2_Freq1008 #= 4)) or 
	((Chr1_Chr2_Row1008 #= 6) and (Chr1_Chr2_Freq1008 #= 5)) or 
	((Chr1_Chr2_Row1008 #= 7) and (Chr1_Chr2_Freq1008 #= 6)) or 
	((Chr1_Chr2_Row1008 #= 8) and (Chr1_Chr2_Freq1008 #= 8)) or 
	((Chr1_Chr2_Row1008 #= 9) and (Chr1_Chr2_Freq1008 #= 6)) or 
	((Chr1_Chr2_Freq1008 #= 0) and (Chr1_Chr2_Row1008 #= 0)),

	% All of the values assumed by the Row<i> variables must be 
	% all different or zero to ensure each genomic bin is 
	% involved in only 1 interaction; multiple zeros are allowed
	alldifferent_except(Non_Zero_Rows),
	atmost(49, Non_Zero_Rows, 0),

	% Optimize: the sum of the selected interaction 
	% frequencies is maximal (it is the minimum of 
	% the additive inverse of the sum
	Cost #= -sum(Freqs),
	search(Freqs, 0, input_order, indomain_max, bb_min(Cost), []),
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% 	Output the results
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	open(FreqFile, 'write', FREQ_OUT),
	%%list the frequencies
	(foreach(X,Freqs),
		param(FREQ_OUT) do
			get_domain_as_list(X, DomList),
				(foreach(Y,DomList),
					param(FREQ_OUT) do
						write(FREQ_OUT, Y),
						write(FREQ_OUT, ' ')
				),
			write(FREQ_OUT, "\n")
	),
	close(FREQ_OUT),

	%% list the potential rows
	open(RowFile, 'write', ROW_OUT),
	(foreach(X,Non_Zero_Rows),
		param(ROW_OUT) do
			get_domain_as_list(X, DomList),
				(foreach(Y,DomList),
					param(ROW_OUT) do
						write(ROW_OUT, Y),
						write(ROW_OUT, ' ')
				),
				write(ROW_OUT, "\n")
	),
	close(ROW_OUT).
