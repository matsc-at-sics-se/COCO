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
	Non_Zero_Rows = [Chr1_Chr3_Row1024, Chr1_Chr3_Row1094, Chr1_Chr3_Row1095, Chr1_Chr3_Row1096, Chr1_Chr3_Row1098, Chr1_Chr3_Row1099, Chr1_Chr3_Row1100, Chr1_Chr3_Row1101, Chr1_Chr3_Row1102, Chr1_Chr3_Row1103, Chr1_Chr3_Row1104, Chr1_Chr3_Row1105, Chr1_Chr3_Row1106, Chr1_Chr3_Row1107, Chr1_Chr3_Row1108, Chr1_Chr3_Row1109, Chr1_Chr3_Row1110, Chr1_Chr3_Row1111, Chr1_Chr3_Row1112, Chr1_Chr3_Row1113, Chr1_Chr3_Row1114, Chr1_Chr3_Row1115, Chr1_Chr3_Row1116, Chr1_Chr3_Row1117, Chr1_Chr3_Row1118, Chr1_Chr3_Row1127, Chr1_Chr3_Row1128, Chr1_Chr3_Row1129, Chr1_Chr3_Row1130, Chr1_Chr3_Row1131, Chr1_Chr3_Row1132, Chr1_Chr3_Row1133, Chr1_Chr3_Row1134, Chr1_Chr3_Row1135, Chr1_Chr3_Row1136, Chr1_Chr3_Row1137, Chr1_Chr3_Row1138, Chr1_Chr3_Row1139, Chr1_Chr3_Row1140, Chr1_Chr3_Row1141, Chr1_Chr3_Row1142, Chr1_Chr3_Row1143, Chr1_Chr3_Row1144, Chr1_Chr3_Row1145, Chr1_Chr3_Row1146, Chr1_Chr3_Row1147, Chr1_Chr3_Row1148, Chr1_Chr3_Row1149, Chr1_Chr3_Row1150, Chr1_Chr3_Row1152, Chr1_Chr3_Row1172, Chr1_Chr3_Row1256, Chr1_Chr3_Row1257, Chr1_Chr3_Row1258],

	% The list Freqs has one variable for each non-zero row of the 
	% whole-genome contact map	
	Freqs = [Chr1_Chr3_Freq1024, Chr1_Chr3_Freq1094, Chr1_Chr3_Freq1095, Chr1_Chr3_Freq1096, Chr1_Chr3_Freq1098, Chr1_Chr3_Freq1099, Chr1_Chr3_Freq1100, Chr1_Chr3_Freq1101, Chr1_Chr3_Freq1102, Chr1_Chr3_Freq1103, Chr1_Chr3_Freq1104, Chr1_Chr3_Freq1105, Chr1_Chr3_Freq1106, Chr1_Chr3_Freq1107, Chr1_Chr3_Freq1108, Chr1_Chr3_Freq1109, Chr1_Chr3_Freq1110, Chr1_Chr3_Freq1111, Chr1_Chr3_Freq1112, Chr1_Chr3_Freq1113, Chr1_Chr3_Freq1114, Chr1_Chr3_Freq1115, Chr1_Chr3_Freq1116, Chr1_Chr3_Freq1117, Chr1_Chr3_Freq1118, Chr1_Chr3_Freq1127, Chr1_Chr3_Freq1128, Chr1_Chr3_Freq1129, Chr1_Chr3_Freq1130, Chr1_Chr3_Freq1131, Chr1_Chr3_Freq1132, Chr1_Chr3_Freq1133, Chr1_Chr3_Freq1134, Chr1_Chr3_Freq1135, Chr1_Chr3_Freq1136, Chr1_Chr3_Freq1137, Chr1_Chr3_Freq1138, Chr1_Chr3_Freq1139, Chr1_Chr3_Freq1140, Chr1_Chr3_Freq1141, Chr1_Chr3_Freq1142, Chr1_Chr3_Freq1143, Chr1_Chr3_Freq1144, Chr1_Chr3_Freq1145, Chr1_Chr3_Freq1146, Chr1_Chr3_Freq1147, Chr1_Chr3_Freq1148, Chr1_Chr3_Freq1149, Chr1_Chr3_Freq1150, Chr1_Chr3_Freq1152, Chr1_Chr3_Freq1172, Chr1_Chr3_Freq1256, Chr1_Chr3_Freq1257, Chr1_Chr3_Freq1258],
	
	% Representation of the Genome: 
	% Each Row term can assume a value based on interacting bin 
	% indices where `0' represents an interaction not being 
	% selected and a non-zero value (ranging from 1 to N) 
	% represents which genomic bin is involved in the selected 
	% interaction 
	Chr1_Chr3_Row1024 :: [0, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389],
	Chr1_Chr3_Row1094 :: [0, 165, 225],
	Chr1_Chr3_Row1095 :: [0, 356, 388],
	Chr1_Chr3_Row1096 :: [0, 381, 385, 396],
	Chr1_Chr3_Row1098 :: [0, 369],
	Chr1_Chr3_Row1099 :: [0, 356, 364, 367, 368, 370, 371, 382, 401],
	Chr1_Chr3_Row1100 :: [0, 348, 365, 366, 369, 370, 381],
	Chr1_Chr3_Row1101 :: [0, 365, 366, 369, 370, 381],
	Chr1_Chr3_Row1102 :: [0, 369, 370, 371, 373, 374, 381, 383, 386, 391],
	Chr1_Chr3_Row1103 :: [0, 356, 357, 358, 361, 365, 367, 368, 370, 371, 372, 374, 381, 383, 386, 387, 391],
	Chr1_Chr3_Row1104 :: [0, 358, 363, 364, 366, 367, 368, 369, 370, 371, 372, 373, 374, 381, 382, 383, 385, 388, 389, 391, 394],
	Chr1_Chr3_Row1105 :: [0, 356, 362, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 391, 392, 393],
	Chr1_Chr3_Row1106 :: [0, 349, 353, 356, 358, 363, 364, 365, 366, 367, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 389, 390, 391, 392, 393, 395],
	Chr1_Chr3_Row1107 :: [0, 358, 362, 363, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 393],
	Chr1_Chr3_Row1108 :: [0, 348, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 389, 390, 391, 396],
	Chr1_Chr3_Row1109 :: [0, 356, 358, 359, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 393, 394, 395, 396, 397, 400],
	Chr1_Chr3_Row1110 :: [0, 356, 357, 358, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 396, 397, 399, 401],
	Chr1_Chr3_Row1111 :: [0, 196, 349, 354, 356, 358, 359, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 399, 400, 401],
	Chr1_Chr3_Row1112 :: [0, 353, 354, 358, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 396, 399],
	Chr1_Chr3_Row1113 :: [0, 342, 358, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 396, 399],
	Chr1_Chr3_Row1114 :: [0, 358, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 392],
	Chr1_Chr3_Row1115 :: [0, 358, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 395, 396],
	Chr1_Chr3_Row1116 :: [0, 350, 354, 355, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 396, 397, 398, 399],
	Chr1_Chr3_Row1117 :: [0, 355, 356, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396],
	Chr1_Chr3_Row1118 :: [0, 358, 359, 360, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 394, 399],
	Chr1_Chr3_Row1127 :: [0, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 394, 395, 396],
	Chr1_Chr3_Row1128 :: [0, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 401],
	Chr1_Chr3_Row1129 :: [0, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396],
	Chr1_Chr3_Row1130 :: [0, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 394, 395],
	Chr1_Chr3_Row1131 :: [0, 358, 359, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395],
	Chr1_Chr3_Row1132 :: [0, 356, 357, 358, 359, 360, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 403],
	Chr1_Chr3_Row1133 :: [0, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396],
	Chr1_Chr3_Row1134 :: [0, 356, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395],
	Chr1_Chr3_Row1135 :: [0, 356, 358, 360, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 392],
	Chr1_Chr3_Row1136 :: [0, 349, 352, 356, 358, 359, 360, 361, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 391, 393, 399],
	Chr1_Chr3_Row1137 :: [0, 356, 357, 358, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 394, 395, 396, 397, 399, 400],
	Chr1_Chr3_Row1138 :: [0, 352, 356, 358, 359, 360, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 391, 393, 396],
	Chr1_Chr3_Row1139 :: [0, 355, 356, 360, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 395, 396],
	Chr1_Chr3_Row1140 :: [0, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 387, 388, 389, 391, 394],
	Chr1_Chr3_Row1141 :: [0, 358, 359, 360, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 384, 385, 386, 387, 389, 390, 391, 392, 393, 396],
	Chr1_Chr3_Row1142 :: [0, 345, 351, 354, 356, 360, 362, 364, 365, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 385, 386, 387, 388, 389, 392, 393, 401],
	Chr1_Chr3_Row1143 :: [0, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 380, 381, 382, 383, 385, 386, 387, 390, 393, 394, 396, 401],
	Chr1_Chr3_Row1144 :: [0, 358, 364, 368, 370, 372, 373, 374, 380, 381, 383, 385],
	Chr1_Chr3_Row1145 :: [0, 365, 367, 368, 371, 372, 373, 374, 380],
	Chr1_Chr3_Row1146 :: [0, 367, 368, 372, 373, 374, 381, 382],
	Chr1_Chr3_Row1147 :: [0, 369, 370, 372, 386],
	Chr1_Chr3_Row1148 :: [0, 365, 369, 372, 373, 374, 380, 388],
	Chr1_Chr3_Row1149 :: [0, 364, 365, 373, 382],
	Chr1_Chr3_Row1150 :: [0, 365, 366, 367, 369, 382],
	Chr1_Chr3_Row1152 :: [0, 367],
	Chr1_Chr3_Row1172 :: [0, 338],
	Chr1_Chr3_Row1256 :: [0, 557],
	Chr1_Chr3_Row1257 :: [0, 6, 7, 8, 11, 12, 13, 20, 21, 22, 23, 24, 25, 34, 39, 47, 48, 49, 51, 53, 58, 59, 65, 69, 71, 73, 88, 94, 176, 204, 211, 225, 244, 479, 537, 539, 542, 543, 546, 547, 548, 549, 553, 554, 555, 556, 557],
	Chr1_Chr3_Row1258 :: [0, 8, 10, 12, 20, 30, 33, 36, 56, 58, 215, 225, 537, 541, 542, 543, 546, 550, 553, 554, 555, 556, 557],

	% Each frequency term can assume either the rounded 
	% and scaled integer value (based on the corresponding 
	% interaction frequency from the whole-genome contact 
	% map) or a value of `0' where `0'  represents an 
	% interaction not being selected 
	Chr1_Chr3_Freq1024 :: [0, 1, 2, 3, 7, 4],
	Chr1_Chr3_Freq1094 :: [0, 1],
	Chr1_Chr3_Freq1095 :: [0, 1],
	Chr1_Chr3_Freq1096 :: [0, 1],
	Chr1_Chr3_Freq1098 :: [0, 1],
	Chr1_Chr3_Freq1099 :: [0, 1],
	Chr1_Chr3_Freq1100 :: [0, 1],
	Chr1_Chr3_Freq1101 :: [0, 1],
	Chr1_Chr3_Freq1102 :: [0, 1],
	Chr1_Chr3_Freq1103 :: [0, 1],
	Chr1_Chr3_Freq1104 :: [0, 1],
	Chr1_Chr3_Freq1105 :: [0, 1],
	Chr1_Chr3_Freq1106 :: [0, 1],
	Chr1_Chr3_Freq1107 :: [0, 1],
	Chr1_Chr3_Freq1108 :: [0, 1, 2],
	Chr1_Chr3_Freq1109 :: [0, 1, 2],
	Chr1_Chr3_Freq1110 :: [0, 1, 2],
	Chr1_Chr3_Freq1111 :: [0, 1, 2],
	Chr1_Chr3_Freq1112 :: [0, 1, 2],
	Chr1_Chr3_Freq1113 :: [0, 1, 2],
	Chr1_Chr3_Freq1114 :: [0, 1, 2, 3],
	Chr1_Chr3_Freq1115 :: [0, 1, 2, 3],
	Chr1_Chr3_Freq1116 :: [0, 1, 2, 3, 4],
	Chr1_Chr3_Freq1117 :: [0, 1, 2, 3, 4, 5],
	Chr1_Chr3_Freq1118 :: [0, 1, 2, 3, 4, 5],
	Chr1_Chr3_Freq1127 :: [0, 1, 2, 3, 4, 5, 6, 7],
	Chr1_Chr3_Freq1128 :: [0, 1, 2, 3, 4, 6, 7, 5],
	Chr1_Chr3_Freq1129 :: [0, 1, 2, 3, 5],
	Chr1_Chr3_Freq1130 :: [0, 1, 2, 3, 4],
	Chr1_Chr3_Freq1131 :: [0, 1, 2, 3, 4],
	Chr1_Chr3_Freq1132 :: [0, 1, 2, 3, 4],
	Chr1_Chr3_Freq1133 :: [0, 1, 2, 3],
	Chr1_Chr3_Freq1134 :: [0, 1, 2, 3],
	Chr1_Chr3_Freq1135 :: [0, 1, 2],
	Chr1_Chr3_Freq1136 :: [0, 1, 2],
	Chr1_Chr3_Freq1137 :: [0, 1, 2],
	Chr1_Chr3_Freq1138 :: [0, 1, 2],
	Chr1_Chr3_Freq1139 :: [0, 1],
	Chr1_Chr3_Freq1140 :: [0, 1],
	Chr1_Chr3_Freq1141 :: [0, 1],
	Chr1_Chr3_Freq1142 :: [0, 1],
	Chr1_Chr3_Freq1143 :: [0, 1],
	Chr1_Chr3_Freq1144 :: [0, 1],
	Chr1_Chr3_Freq1145 :: [0, 1],
	Chr1_Chr3_Freq1146 :: [0, 1],
	Chr1_Chr3_Freq1147 :: [0, 1],
	Chr1_Chr3_Freq1148 :: [0, 1],
	Chr1_Chr3_Freq1149 :: [0, 1],
	Chr1_Chr3_Freq1150 :: [0, 1],
	Chr1_Chr3_Freq1152 :: [0, 1],
	Chr1_Chr3_Freq1172 :: [0, 1],
	Chr1_Chr3_Freq1256 :: [0, 1],
	Chr1_Chr3_Freq1257 :: [0, 1, 2, 4],
	Chr1_Chr3_Freq1258 :: [0, 1, 2, 3, 4, 6],

	% Constraints: 
	% Each pair of corresponding (Row<i>, Freq<i>) variables 
	% must assume dependent values based on data from the 
	% whole-genome contact map; A (Row, Freq) pair ground to 
	% (0,0) encodes that nothing is chosen 
	((Chr1_Chr3_Row1024 #= 371) and (Chr1_Chr3_Freq1024 #= 1)) or 
	((Chr1_Chr3_Row1024 #= 372) and (Chr1_Chr3_Freq1024 #= 1)) or 
	((Chr1_Chr3_Row1024 #= 373) and (Chr1_Chr3_Freq1024 #= 2)) or 
	((Chr1_Chr3_Row1024 #= 374) and (Chr1_Chr3_Freq1024 #= 3)) or 
	((Chr1_Chr3_Row1024 #= 380) and (Chr1_Chr3_Freq1024 #= 7)) or 
	((Chr1_Chr3_Row1024 #= 381) and (Chr1_Chr3_Freq1024 #= 4)) or 
	((Chr1_Chr3_Row1024 #= 382) and (Chr1_Chr3_Freq1024 #= 3)) or 
	((Chr1_Chr3_Row1024 #= 383) and (Chr1_Chr3_Freq1024 #= 2)) or 
	((Chr1_Chr3_Row1024 #= 384) and (Chr1_Chr3_Freq1024 #= 1)) or 
	((Chr1_Chr3_Row1024 #= 385) and (Chr1_Chr3_Freq1024 #= 1)) or 
	((Chr1_Chr3_Row1024 #= 386) and (Chr1_Chr3_Freq1024 #= 1)) or 
	((Chr1_Chr3_Row1024 #= 387) and (Chr1_Chr3_Freq1024 #= 1)) or 
	((Chr1_Chr3_Row1024 #= 388) and (Chr1_Chr3_Freq1024 #= 1)) or 
	((Chr1_Chr3_Row1024 #= 389) and (Chr1_Chr3_Freq1024 #= 1)) or 
	((Chr1_Chr3_Freq1024 #= 0) and (Chr1_Chr3_Row1024 #= 0)),

	((Chr1_Chr3_Row1094 #= 165) and (Chr1_Chr3_Freq1094 #= 1)) or 
	((Chr1_Chr3_Row1094 #= 225) and (Chr1_Chr3_Freq1094 #= 1)) or 
	((Chr1_Chr3_Freq1094 #= 0) and (Chr1_Chr3_Row1094 #= 0)),

	((Chr1_Chr3_Row1095 #= 356) and (Chr1_Chr3_Freq1095 #= 1)) or 
	((Chr1_Chr3_Row1095 #= 388) and (Chr1_Chr3_Freq1095 #= 1)) or 
	((Chr1_Chr3_Freq1095 #= 0) and (Chr1_Chr3_Row1095 #= 0)),

	((Chr1_Chr3_Row1096 #= 381) and (Chr1_Chr3_Freq1096 #= 1)) or 
	((Chr1_Chr3_Row1096 #= 385) and (Chr1_Chr3_Freq1096 #= 1)) or 
	((Chr1_Chr3_Row1096 #= 396) and (Chr1_Chr3_Freq1096 #= 1)) or 
	((Chr1_Chr3_Freq1096 #= 0) and (Chr1_Chr3_Row1096 #= 0)),

	((Chr1_Chr3_Row1098 #= 369) and (Chr1_Chr3_Freq1098 #= 1)) or 
	((Chr1_Chr3_Freq1098 #= 0) and (Chr1_Chr3_Row1098 #= 0)),

	((Chr1_Chr3_Row1099 #= 356) and (Chr1_Chr3_Freq1099 #= 1)) or 
	((Chr1_Chr3_Row1099 #= 364) and (Chr1_Chr3_Freq1099 #= 1)) or 
	((Chr1_Chr3_Row1099 #= 367) and (Chr1_Chr3_Freq1099 #= 1)) or 
	((Chr1_Chr3_Row1099 #= 368) and (Chr1_Chr3_Freq1099 #= 1)) or 
	((Chr1_Chr3_Row1099 #= 370) and (Chr1_Chr3_Freq1099 #= 1)) or 
	((Chr1_Chr3_Row1099 #= 371) and (Chr1_Chr3_Freq1099 #= 1)) or 
	((Chr1_Chr3_Row1099 #= 382) and (Chr1_Chr3_Freq1099 #= 1)) or 
	((Chr1_Chr3_Row1099 #= 401) and (Chr1_Chr3_Freq1099 #= 1)) or 
	((Chr1_Chr3_Freq1099 #= 0) and (Chr1_Chr3_Row1099 #= 0)),

	((Chr1_Chr3_Row1100 #= 348) and (Chr1_Chr3_Freq1100 #= 1)) or 
	((Chr1_Chr3_Row1100 #= 365) and (Chr1_Chr3_Freq1100 #= 1)) or 
	((Chr1_Chr3_Row1100 #= 366) and (Chr1_Chr3_Freq1100 #= 1)) or 
	((Chr1_Chr3_Row1100 #= 369) and (Chr1_Chr3_Freq1100 #= 1)) or 
	((Chr1_Chr3_Row1100 #= 370) and (Chr1_Chr3_Freq1100 #= 1)) or 
	((Chr1_Chr3_Row1100 #= 381) and (Chr1_Chr3_Freq1100 #= 1)) or 
	((Chr1_Chr3_Freq1100 #= 0) and (Chr1_Chr3_Row1100 #= 0)),

	((Chr1_Chr3_Row1101 #= 365) and (Chr1_Chr3_Freq1101 #= 1)) or 
	((Chr1_Chr3_Row1101 #= 366) and (Chr1_Chr3_Freq1101 #= 1)) or 
	((Chr1_Chr3_Row1101 #= 369) and (Chr1_Chr3_Freq1101 #= 1)) or 
	((Chr1_Chr3_Row1101 #= 370) and (Chr1_Chr3_Freq1101 #= 1)) or 
	((Chr1_Chr3_Row1101 #= 381) and (Chr1_Chr3_Freq1101 #= 1)) or 
	((Chr1_Chr3_Freq1101 #= 0) and (Chr1_Chr3_Row1101 #= 0)),

	((Chr1_Chr3_Row1102 #= 369) and (Chr1_Chr3_Freq1102 #= 1)) or 
	((Chr1_Chr3_Row1102 #= 370) and (Chr1_Chr3_Freq1102 #= 1)) or 
	((Chr1_Chr3_Row1102 #= 371) and (Chr1_Chr3_Freq1102 #= 1)) or 
	((Chr1_Chr3_Row1102 #= 373) and (Chr1_Chr3_Freq1102 #= 1)) or 
	((Chr1_Chr3_Row1102 #= 374) and (Chr1_Chr3_Freq1102 #= 1)) or 
	((Chr1_Chr3_Row1102 #= 381) and (Chr1_Chr3_Freq1102 #= 1)) or 
	((Chr1_Chr3_Row1102 #= 383) and (Chr1_Chr3_Freq1102 #= 1)) or 
	((Chr1_Chr3_Row1102 #= 386) and (Chr1_Chr3_Freq1102 #= 1)) or 
	((Chr1_Chr3_Row1102 #= 391) and (Chr1_Chr3_Freq1102 #= 1)) or 
	((Chr1_Chr3_Freq1102 #= 0) and (Chr1_Chr3_Row1102 #= 0)),

	((Chr1_Chr3_Row1103 #= 356) and (Chr1_Chr3_Freq1103 #= 1)) or 
	((Chr1_Chr3_Row1103 #= 357) and (Chr1_Chr3_Freq1103 #= 1)) or 
	((Chr1_Chr3_Row1103 #= 358) and (Chr1_Chr3_Freq1103 #= 1)) or 
	((Chr1_Chr3_Row1103 #= 361) and (Chr1_Chr3_Freq1103 #= 1)) or 
	((Chr1_Chr3_Row1103 #= 365) and (Chr1_Chr3_Freq1103 #= 1)) or 
	((Chr1_Chr3_Row1103 #= 367) and (Chr1_Chr3_Freq1103 #= 1)) or 
	((Chr1_Chr3_Row1103 #= 368) and (Chr1_Chr3_Freq1103 #= 1)) or 
	((Chr1_Chr3_Row1103 #= 370) and (Chr1_Chr3_Freq1103 #= 1)) or 
	((Chr1_Chr3_Row1103 #= 371) and (Chr1_Chr3_Freq1103 #= 1)) or 
	((Chr1_Chr3_Row1103 #= 372) and (Chr1_Chr3_Freq1103 #= 1)) or 
	((Chr1_Chr3_Row1103 #= 374) and (Chr1_Chr3_Freq1103 #= 1)) or 
	((Chr1_Chr3_Row1103 #= 381) and (Chr1_Chr3_Freq1103 #= 1)) or 
	((Chr1_Chr3_Row1103 #= 383) and (Chr1_Chr3_Freq1103 #= 1)) or 
	((Chr1_Chr3_Row1103 #= 386) and (Chr1_Chr3_Freq1103 #= 1)) or 
	((Chr1_Chr3_Row1103 #= 387) and (Chr1_Chr3_Freq1103 #= 1)) or 
	((Chr1_Chr3_Row1103 #= 391) and (Chr1_Chr3_Freq1103 #= 1)) or 
	((Chr1_Chr3_Freq1103 #= 0) and (Chr1_Chr3_Row1103 #= 0)),

	((Chr1_Chr3_Row1104 #= 358) and (Chr1_Chr3_Freq1104 #= 1)) or 
	((Chr1_Chr3_Row1104 #= 363) and (Chr1_Chr3_Freq1104 #= 1)) or 
	((Chr1_Chr3_Row1104 #= 364) and (Chr1_Chr3_Freq1104 #= 1)) or 
	((Chr1_Chr3_Row1104 #= 366) and (Chr1_Chr3_Freq1104 #= 1)) or 
	((Chr1_Chr3_Row1104 #= 367) and (Chr1_Chr3_Freq1104 #= 1)) or 
	((Chr1_Chr3_Row1104 #= 368) and (Chr1_Chr3_Freq1104 #= 1)) or 
	((Chr1_Chr3_Row1104 #= 369) and (Chr1_Chr3_Freq1104 #= 1)) or 
	((Chr1_Chr3_Row1104 #= 370) and (Chr1_Chr3_Freq1104 #= 1)) or 
	((Chr1_Chr3_Row1104 #= 371) and (Chr1_Chr3_Freq1104 #= 1)) or 
	((Chr1_Chr3_Row1104 #= 372) and (Chr1_Chr3_Freq1104 #= 1)) or 
	((Chr1_Chr3_Row1104 #= 373) and (Chr1_Chr3_Freq1104 #= 1)) or 
	((Chr1_Chr3_Row1104 #= 374) and (Chr1_Chr3_Freq1104 #= 1)) or 
	((Chr1_Chr3_Row1104 #= 381) and (Chr1_Chr3_Freq1104 #= 1)) or 
	((Chr1_Chr3_Row1104 #= 382) and (Chr1_Chr3_Freq1104 #= 1)) or 
	((Chr1_Chr3_Row1104 #= 383) and (Chr1_Chr3_Freq1104 #= 1)) or 
	((Chr1_Chr3_Row1104 #= 385) and (Chr1_Chr3_Freq1104 #= 1)) or 
	((Chr1_Chr3_Row1104 #= 388) and (Chr1_Chr3_Freq1104 #= 1)) or 
	((Chr1_Chr3_Row1104 #= 389) and (Chr1_Chr3_Freq1104 #= 1)) or 
	((Chr1_Chr3_Row1104 #= 391) and (Chr1_Chr3_Freq1104 #= 1)) or 
	((Chr1_Chr3_Row1104 #= 394) and (Chr1_Chr3_Freq1104 #= 1)) or 
	((Chr1_Chr3_Freq1104 #= 0) and (Chr1_Chr3_Row1104 #= 0)),

	((Chr1_Chr3_Row1105 #= 356) and (Chr1_Chr3_Freq1105 #= 1)) or 
	((Chr1_Chr3_Row1105 #= 362) and (Chr1_Chr3_Freq1105 #= 1)) or 
	((Chr1_Chr3_Row1105 #= 364) and (Chr1_Chr3_Freq1105 #= 1)) or 
	((Chr1_Chr3_Row1105 #= 365) and (Chr1_Chr3_Freq1105 #= 1)) or 
	((Chr1_Chr3_Row1105 #= 366) and (Chr1_Chr3_Freq1105 #= 1)) or 
	((Chr1_Chr3_Row1105 #= 367) and (Chr1_Chr3_Freq1105 #= 1)) or 
	((Chr1_Chr3_Row1105 #= 368) and (Chr1_Chr3_Freq1105 #= 1)) or 
	((Chr1_Chr3_Row1105 #= 369) and (Chr1_Chr3_Freq1105 #= 1)) or 
	((Chr1_Chr3_Row1105 #= 370) and (Chr1_Chr3_Freq1105 #= 1)) or 
	((Chr1_Chr3_Row1105 #= 371) and (Chr1_Chr3_Freq1105 #= 1)) or 
	((Chr1_Chr3_Row1105 #= 372) and (Chr1_Chr3_Freq1105 #= 1)) or 
	((Chr1_Chr3_Row1105 #= 373) and (Chr1_Chr3_Freq1105 #= 1)) or 
	((Chr1_Chr3_Row1105 #= 374) and (Chr1_Chr3_Freq1105 #= 1)) or 
	((Chr1_Chr3_Row1105 #= 380) and (Chr1_Chr3_Freq1105 #= 1)) or 
	((Chr1_Chr3_Row1105 #= 381) and (Chr1_Chr3_Freq1105 #= 1)) or 
	((Chr1_Chr3_Row1105 #= 382) and (Chr1_Chr3_Freq1105 #= 1)) or 
	((Chr1_Chr3_Row1105 #= 383) and (Chr1_Chr3_Freq1105 #= 1)) or 
	((Chr1_Chr3_Row1105 #= 384) and (Chr1_Chr3_Freq1105 #= 1)) or 
	((Chr1_Chr3_Row1105 #= 385) and (Chr1_Chr3_Freq1105 #= 1)) or 
	((Chr1_Chr3_Row1105 #= 386) and (Chr1_Chr3_Freq1105 #= 1)) or 
	((Chr1_Chr3_Row1105 #= 387) and (Chr1_Chr3_Freq1105 #= 1)) or 
	((Chr1_Chr3_Row1105 #= 391) and (Chr1_Chr3_Freq1105 #= 1)) or 
	((Chr1_Chr3_Row1105 #= 392) and (Chr1_Chr3_Freq1105 #= 1)) or 
	((Chr1_Chr3_Row1105 #= 393) and (Chr1_Chr3_Freq1105 #= 1)) or 
	((Chr1_Chr3_Freq1105 #= 0) and (Chr1_Chr3_Row1105 #= 0)),

	((Chr1_Chr3_Row1106 #= 349) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 353) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 356) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 358) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 363) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 364) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 365) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 366) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 367) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 369) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 370) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 371) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 372) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 373) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 374) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 380) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 381) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 382) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 383) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 384) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 385) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 386) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 387) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 389) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 390) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 391) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 392) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 393) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Row1106 #= 395) and (Chr1_Chr3_Freq1106 #= 1)) or 
	((Chr1_Chr3_Freq1106 #= 0) and (Chr1_Chr3_Row1106 #= 0)),

	((Chr1_Chr3_Row1107 #= 358) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 362) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 363) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 365) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 366) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 367) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 368) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 369) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 370) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 371) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 372) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 373) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 374) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 380) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 381) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 382) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 383) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 384) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 385) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 386) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 387) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 388) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 389) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 390) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 391) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Row1107 #= 393) and (Chr1_Chr3_Freq1107 #= 1)) or 
	((Chr1_Chr3_Freq1107 #= 0) and (Chr1_Chr3_Row1107 #= 0)),

	((Chr1_Chr3_Row1108 #= 348) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 358) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 359) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 360) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 361) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 362) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 363) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 364) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 365) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 366) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 367) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 368) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 369) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 370) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 371) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 372) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 373) and (Chr1_Chr3_Freq1108 #= 2)) or 
	((Chr1_Chr3_Row1108 #= 374) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 380) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 381) and (Chr1_Chr3_Freq1108 #= 2)) or 
	((Chr1_Chr3_Row1108 #= 382) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 383) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 384) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 385) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 386) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 387) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 389) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 390) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 391) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Row1108 #= 396) and (Chr1_Chr3_Freq1108 #= 1)) or 
	((Chr1_Chr3_Freq1108 #= 0) and (Chr1_Chr3_Row1108 #= 0)),

	((Chr1_Chr3_Row1109 #= 356) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 358) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 359) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 361) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 362) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 363) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 364) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 365) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 366) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 367) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 368) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 369) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 370) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 371) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 372) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 373) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 374) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 380) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 381) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 382) and (Chr1_Chr3_Freq1109 #= 2)) or 
	((Chr1_Chr3_Row1109 #= 383) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 384) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 385) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 386) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 387) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 388) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 389) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 390) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 391) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 393) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 394) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 395) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 396) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 397) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Row1109 #= 400) and (Chr1_Chr3_Freq1109 #= 1)) or 
	((Chr1_Chr3_Freq1109 #= 0) and (Chr1_Chr3_Row1109 #= 0)),

	((Chr1_Chr3_Row1110 #= 356) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 357) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 358) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 360) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 361) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 362) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 363) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 364) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 365) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 366) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 367) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 368) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 369) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 370) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 371) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 372) and (Chr1_Chr3_Freq1110 #= 2)) or 
	((Chr1_Chr3_Row1110 #= 373) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 374) and (Chr1_Chr3_Freq1110 #= 2)) or 
	((Chr1_Chr3_Row1110 #= 380) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 381) and (Chr1_Chr3_Freq1110 #= 2)) or 
	((Chr1_Chr3_Row1110 #= 382) and (Chr1_Chr3_Freq1110 #= 2)) or 
	((Chr1_Chr3_Row1110 #= 383) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 384) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 385) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 386) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 387) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 388) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 389) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 390) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 391) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 392) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 393) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 396) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 397) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 399) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Row1110 #= 401) and (Chr1_Chr3_Freq1110 #= 1)) or 
	((Chr1_Chr3_Freq1110 #= 0) and (Chr1_Chr3_Row1110 #= 0)),

	((Chr1_Chr3_Row1111 #= 196) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 349) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 354) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 356) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 358) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 359) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 361) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 362) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 363) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 364) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 365) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 366) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 367) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 368) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 369) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 370) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 371) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 372) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 373) and (Chr1_Chr3_Freq1111 #= 2)) or 
	((Chr1_Chr3_Row1111 #= 374) and (Chr1_Chr3_Freq1111 #= 2)) or 
	((Chr1_Chr3_Row1111 #= 380) and (Chr1_Chr3_Freq1111 #= 2)) or 
	((Chr1_Chr3_Row1111 #= 381) and (Chr1_Chr3_Freq1111 #= 2)) or 
	((Chr1_Chr3_Row1111 #= 382) and (Chr1_Chr3_Freq1111 #= 2)) or 
	((Chr1_Chr3_Row1111 #= 383) and (Chr1_Chr3_Freq1111 #= 2)) or 
	((Chr1_Chr3_Row1111 #= 384) and (Chr1_Chr3_Freq1111 #= 2)) or 
	((Chr1_Chr3_Row1111 #= 385) and (Chr1_Chr3_Freq1111 #= 2)) or 
	((Chr1_Chr3_Row1111 #= 386) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 387) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 388) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 389) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 390) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 391) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 392) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 393) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 394) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 395) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 396) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 399) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 400) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Row1111 #= 401) and (Chr1_Chr3_Freq1111 #= 1)) or 
	((Chr1_Chr3_Freq1111 #= 0) and (Chr1_Chr3_Row1111 #= 0)),

	((Chr1_Chr3_Row1112 #= 353) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Row1112 #= 354) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Row1112 #= 358) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Row1112 #= 360) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Row1112 #= 361) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Row1112 #= 362) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Row1112 #= 363) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Row1112 #= 364) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Row1112 #= 365) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Row1112 #= 366) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Row1112 #= 367) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Row1112 #= 368) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Row1112 #= 369) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Row1112 #= 370) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Row1112 #= 371) and (Chr1_Chr3_Freq1112 #= 2)) or 
	((Chr1_Chr3_Row1112 #= 372) and (Chr1_Chr3_Freq1112 #= 2)) or 
	((Chr1_Chr3_Row1112 #= 373) and (Chr1_Chr3_Freq1112 #= 2)) or 
	((Chr1_Chr3_Row1112 #= 374) and (Chr1_Chr3_Freq1112 #= 2)) or 
	((Chr1_Chr3_Row1112 #= 380) and (Chr1_Chr3_Freq1112 #= 2)) or 
	((Chr1_Chr3_Row1112 #= 381) and (Chr1_Chr3_Freq1112 #= 2)) or 
	((Chr1_Chr3_Row1112 #= 382) and (Chr1_Chr3_Freq1112 #= 2)) or 
	((Chr1_Chr3_Row1112 #= 383) and (Chr1_Chr3_Freq1112 #= 2)) or 
	((Chr1_Chr3_Row1112 #= 384) and (Chr1_Chr3_Freq1112 #= 2)) or 
	((Chr1_Chr3_Row1112 #= 385) and (Chr1_Chr3_Freq1112 #= 2)) or 
	((Chr1_Chr3_Row1112 #= 386) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Row1112 #= 387) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Row1112 #= 388) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Row1112 #= 389) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Row1112 #= 390) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Row1112 #= 391) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Row1112 #= 392) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Row1112 #= 393) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Row1112 #= 394) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Row1112 #= 396) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Row1112 #= 399) and (Chr1_Chr3_Freq1112 #= 1)) or 
	((Chr1_Chr3_Freq1112 #= 0) and (Chr1_Chr3_Row1112 #= 0)),

	((Chr1_Chr3_Row1113 #= 342) and (Chr1_Chr3_Freq1113 #= 1)) or 
	((Chr1_Chr3_Row1113 #= 358) and (Chr1_Chr3_Freq1113 #= 1)) or 
	((Chr1_Chr3_Row1113 #= 360) and (Chr1_Chr3_Freq1113 #= 1)) or 
	((Chr1_Chr3_Row1113 #= 361) and (Chr1_Chr3_Freq1113 #= 1)) or 
	((Chr1_Chr3_Row1113 #= 362) and (Chr1_Chr3_Freq1113 #= 1)) or 
	((Chr1_Chr3_Row1113 #= 363) and (Chr1_Chr3_Freq1113 #= 1)) or 
	((Chr1_Chr3_Row1113 #= 364) and (Chr1_Chr3_Freq1113 #= 1)) or 
	((Chr1_Chr3_Row1113 #= 365) and (Chr1_Chr3_Freq1113 #= 1)) or 
	((Chr1_Chr3_Row1113 #= 366) and (Chr1_Chr3_Freq1113 #= 1)) or 
	((Chr1_Chr3_Row1113 #= 367) and (Chr1_Chr3_Freq1113 #= 1)) or 
	((Chr1_Chr3_Row1113 #= 368) and (Chr1_Chr3_Freq1113 #= 2)) or 
	((Chr1_Chr3_Row1113 #= 369) and (Chr1_Chr3_Freq1113 #= 1)) or 
	((Chr1_Chr3_Row1113 #= 370) and (Chr1_Chr3_Freq1113 #= 1)) or 
	((Chr1_Chr3_Row1113 #= 371) and (Chr1_Chr3_Freq1113 #= 2)) or 
	((Chr1_Chr3_Row1113 #= 372) and (Chr1_Chr3_Freq1113 #= 2)) or 
	((Chr1_Chr3_Row1113 #= 373) and (Chr1_Chr3_Freq1113 #= 2)) or 
	((Chr1_Chr3_Row1113 #= 374) and (Chr1_Chr3_Freq1113 #= 2)) or 
	((Chr1_Chr3_Row1113 #= 380) and (Chr1_Chr3_Freq1113 #= 2)) or 
	((Chr1_Chr3_Row1113 #= 381) and (Chr1_Chr3_Freq1113 #= 2)) or 
	((Chr1_Chr3_Row1113 #= 382) and (Chr1_Chr3_Freq1113 #= 2)) or 
	((Chr1_Chr3_Row1113 #= 383) and (Chr1_Chr3_Freq1113 #= 2)) or 
	((Chr1_Chr3_Row1113 #= 384) and (Chr1_Chr3_Freq1113 #= 2)) or 
	((Chr1_Chr3_Row1113 #= 385) and (Chr1_Chr3_Freq1113 #= 2)) or 
	((Chr1_Chr3_Row1113 #= 386) and (Chr1_Chr3_Freq1113 #= 2)) or 
	((Chr1_Chr3_Row1113 #= 387) and (Chr1_Chr3_Freq1113 #= 1)) or 
	((Chr1_Chr3_Row1113 #= 388) and (Chr1_Chr3_Freq1113 #= 1)) or 
	((Chr1_Chr3_Row1113 #= 389) and (Chr1_Chr3_Freq1113 #= 1)) or 
	((Chr1_Chr3_Row1113 #= 390) and (Chr1_Chr3_Freq1113 #= 1)) or 
	((Chr1_Chr3_Row1113 #= 391) and (Chr1_Chr3_Freq1113 #= 1)) or 
	((Chr1_Chr3_Row1113 #= 392) and (Chr1_Chr3_Freq1113 #= 1)) or 
	((Chr1_Chr3_Row1113 #= 393) and (Chr1_Chr3_Freq1113 #= 1)) or 
	((Chr1_Chr3_Row1113 #= 394) and (Chr1_Chr3_Freq1113 #= 1)) or 
	((Chr1_Chr3_Row1113 #= 396) and (Chr1_Chr3_Freq1113 #= 1)) or 
	((Chr1_Chr3_Row1113 #= 399) and (Chr1_Chr3_Freq1113 #= 1)) or 
	((Chr1_Chr3_Freq1113 #= 0) and (Chr1_Chr3_Row1113 #= 0)),

	((Chr1_Chr3_Row1114 #= 358) and (Chr1_Chr3_Freq1114 #= 1)) or 
	((Chr1_Chr3_Row1114 #= 361) and (Chr1_Chr3_Freq1114 #= 1)) or 
	((Chr1_Chr3_Row1114 #= 362) and (Chr1_Chr3_Freq1114 #= 1)) or 
	((Chr1_Chr3_Row1114 #= 363) and (Chr1_Chr3_Freq1114 #= 1)) or 
	((Chr1_Chr3_Row1114 #= 364) and (Chr1_Chr3_Freq1114 #= 1)) or 
	((Chr1_Chr3_Row1114 #= 365) and (Chr1_Chr3_Freq1114 #= 1)) or 
	((Chr1_Chr3_Row1114 #= 366) and (Chr1_Chr3_Freq1114 #= 1)) or 
	((Chr1_Chr3_Row1114 #= 367) and (Chr1_Chr3_Freq1114 #= 1)) or 
	((Chr1_Chr3_Row1114 #= 368) and (Chr1_Chr3_Freq1114 #= 1)) or 
	((Chr1_Chr3_Row1114 #= 369) and (Chr1_Chr3_Freq1114 #= 1)) or 
	((Chr1_Chr3_Row1114 #= 370) and (Chr1_Chr3_Freq1114 #= 1)) or 
	((Chr1_Chr3_Row1114 #= 371) and (Chr1_Chr3_Freq1114 #= 2)) or 
	((Chr1_Chr3_Row1114 #= 372) and (Chr1_Chr3_Freq1114 #= 2)) or 
	((Chr1_Chr3_Row1114 #= 373) and (Chr1_Chr3_Freq1114 #= 2)) or 
	((Chr1_Chr3_Row1114 #= 374) and (Chr1_Chr3_Freq1114 #= 2)) or 
	((Chr1_Chr3_Row1114 #= 380) and (Chr1_Chr3_Freq1114 #= 2)) or 
	((Chr1_Chr3_Row1114 #= 381) and (Chr1_Chr3_Freq1114 #= 3)) or 
	((Chr1_Chr3_Row1114 #= 382) and (Chr1_Chr3_Freq1114 #= 2)) or 
	((Chr1_Chr3_Row1114 #= 383) and (Chr1_Chr3_Freq1114 #= 2)) or 
	((Chr1_Chr3_Row1114 #= 384) and (Chr1_Chr3_Freq1114 #= 1)) or 
	((Chr1_Chr3_Row1114 #= 385) and (Chr1_Chr3_Freq1114 #= 1)) or 
	((Chr1_Chr3_Row1114 #= 386) and (Chr1_Chr3_Freq1114 #= 1)) or 
	((Chr1_Chr3_Row1114 #= 387) and (Chr1_Chr3_Freq1114 #= 1)) or 
	((Chr1_Chr3_Row1114 #= 388) and (Chr1_Chr3_Freq1114 #= 1)) or 
	((Chr1_Chr3_Row1114 #= 389) and (Chr1_Chr3_Freq1114 #= 1)) or 
	((Chr1_Chr3_Row1114 #= 390) and (Chr1_Chr3_Freq1114 #= 1)) or 
	((Chr1_Chr3_Row1114 #= 392) and (Chr1_Chr3_Freq1114 #= 1)) or 
	((Chr1_Chr3_Freq1114 #= 0) and (Chr1_Chr3_Row1114 #= 0)),

	((Chr1_Chr3_Row1115 #= 358) and (Chr1_Chr3_Freq1115 #= 1)) or 
	((Chr1_Chr3_Row1115 #= 360) and (Chr1_Chr3_Freq1115 #= 1)) or 
	((Chr1_Chr3_Row1115 #= 361) and (Chr1_Chr3_Freq1115 #= 1)) or 
	((Chr1_Chr3_Row1115 #= 362) and (Chr1_Chr3_Freq1115 #= 1)) or 
	((Chr1_Chr3_Row1115 #= 363) and (Chr1_Chr3_Freq1115 #= 1)) or 
	((Chr1_Chr3_Row1115 #= 364) and (Chr1_Chr3_Freq1115 #= 1)) or 
	((Chr1_Chr3_Row1115 #= 365) and (Chr1_Chr3_Freq1115 #= 1)) or 
	((Chr1_Chr3_Row1115 #= 366) and (Chr1_Chr3_Freq1115 #= 1)) or 
	((Chr1_Chr3_Row1115 #= 367) and (Chr1_Chr3_Freq1115 #= 1)) or 
	((Chr1_Chr3_Row1115 #= 368) and (Chr1_Chr3_Freq1115 #= 2)) or 
	((Chr1_Chr3_Row1115 #= 369) and (Chr1_Chr3_Freq1115 #= 2)) or 
	((Chr1_Chr3_Row1115 #= 370) and (Chr1_Chr3_Freq1115 #= 2)) or 
	((Chr1_Chr3_Row1115 #= 371) and (Chr1_Chr3_Freq1115 #= 2)) or 
	((Chr1_Chr3_Row1115 #= 372) and (Chr1_Chr3_Freq1115 #= 2)) or 
	((Chr1_Chr3_Row1115 #= 373) and (Chr1_Chr3_Freq1115 #= 3)) or 
	((Chr1_Chr3_Row1115 #= 374) and (Chr1_Chr3_Freq1115 #= 3)) or 
	((Chr1_Chr3_Row1115 #= 380) and (Chr1_Chr3_Freq1115 #= 3)) or 
	((Chr1_Chr3_Row1115 #= 381) and (Chr1_Chr3_Freq1115 #= 3)) or 
	((Chr1_Chr3_Row1115 #= 382) and (Chr1_Chr3_Freq1115 #= 3)) or 
	((Chr1_Chr3_Row1115 #= 383) and (Chr1_Chr3_Freq1115 #= 2)) or 
	((Chr1_Chr3_Row1115 #= 384) and (Chr1_Chr3_Freq1115 #= 2)) or 
	((Chr1_Chr3_Row1115 #= 385) and (Chr1_Chr3_Freq1115 #= 2)) or 
	((Chr1_Chr3_Row1115 #= 386) and (Chr1_Chr3_Freq1115 #= 1)) or 
	((Chr1_Chr3_Row1115 #= 387) and (Chr1_Chr3_Freq1115 #= 2)) or 
	((Chr1_Chr3_Row1115 #= 388) and (Chr1_Chr3_Freq1115 #= 1)) or 
	((Chr1_Chr3_Row1115 #= 389) and (Chr1_Chr3_Freq1115 #= 1)) or 
	((Chr1_Chr3_Row1115 #= 390) and (Chr1_Chr3_Freq1115 #= 1)) or 
	((Chr1_Chr3_Row1115 #= 391) and (Chr1_Chr3_Freq1115 #= 1)) or 
	((Chr1_Chr3_Row1115 #= 392) and (Chr1_Chr3_Freq1115 #= 1)) or 
	((Chr1_Chr3_Row1115 #= 393) and (Chr1_Chr3_Freq1115 #= 1)) or 
	((Chr1_Chr3_Row1115 #= 395) and (Chr1_Chr3_Freq1115 #= 1)) or 
	((Chr1_Chr3_Row1115 #= 396) and (Chr1_Chr3_Freq1115 #= 1)) or 
	((Chr1_Chr3_Freq1115 #= 0) and (Chr1_Chr3_Row1115 #= 0)),

	((Chr1_Chr3_Row1116 #= 350) and (Chr1_Chr3_Freq1116 #= 1)) or 
	((Chr1_Chr3_Row1116 #= 354) and (Chr1_Chr3_Freq1116 #= 1)) or 
	((Chr1_Chr3_Row1116 #= 355) and (Chr1_Chr3_Freq1116 #= 1)) or 
	((Chr1_Chr3_Row1116 #= 358) and (Chr1_Chr3_Freq1116 #= 1)) or 
	((Chr1_Chr3_Row1116 #= 359) and (Chr1_Chr3_Freq1116 #= 1)) or 
	((Chr1_Chr3_Row1116 #= 360) and (Chr1_Chr3_Freq1116 #= 1)) or 
	((Chr1_Chr3_Row1116 #= 361) and (Chr1_Chr3_Freq1116 #= 1)) or 
	((Chr1_Chr3_Row1116 #= 362) and (Chr1_Chr3_Freq1116 #= 1)) or 
	((Chr1_Chr3_Row1116 #= 363) and (Chr1_Chr3_Freq1116 #= 1)) or 
	((Chr1_Chr3_Row1116 #= 364) and (Chr1_Chr3_Freq1116 #= 1)) or 
	((Chr1_Chr3_Row1116 #= 365) and (Chr1_Chr3_Freq1116 #= 1)) or 
	((Chr1_Chr3_Row1116 #= 366) and (Chr1_Chr3_Freq1116 #= 2)) or 
	((Chr1_Chr3_Row1116 #= 367) and (Chr1_Chr3_Freq1116 #= 2)) or 
	((Chr1_Chr3_Row1116 #= 368) and (Chr1_Chr3_Freq1116 #= 2)) or 
	((Chr1_Chr3_Row1116 #= 369) and (Chr1_Chr3_Freq1116 #= 2)) or 
	((Chr1_Chr3_Row1116 #= 370) and (Chr1_Chr3_Freq1116 #= 2)) or 
	((Chr1_Chr3_Row1116 #= 371) and (Chr1_Chr3_Freq1116 #= 3)) or 
	((Chr1_Chr3_Row1116 #= 372) and (Chr1_Chr3_Freq1116 #= 3)) or 
	((Chr1_Chr3_Row1116 #= 373) and (Chr1_Chr3_Freq1116 #= 4)) or 
	((Chr1_Chr3_Row1116 #= 374) and (Chr1_Chr3_Freq1116 #= 3)) or 
	((Chr1_Chr3_Row1116 #= 380) and (Chr1_Chr3_Freq1116 #= 4)) or 
	((Chr1_Chr3_Row1116 #= 381) and (Chr1_Chr3_Freq1116 #= 4)) or 
	((Chr1_Chr3_Row1116 #= 382) and (Chr1_Chr3_Freq1116 #= 3)) or 
	((Chr1_Chr3_Row1116 #= 383) and (Chr1_Chr3_Freq1116 #= 3)) or 
	((Chr1_Chr3_Row1116 #= 384) and (Chr1_Chr3_Freq1116 #= 2)) or 
	((Chr1_Chr3_Row1116 #= 385) and (Chr1_Chr3_Freq1116 #= 2)) or 
	((Chr1_Chr3_Row1116 #= 386) and (Chr1_Chr3_Freq1116 #= 2)) or 
	((Chr1_Chr3_Row1116 #= 387) and (Chr1_Chr3_Freq1116 #= 2)) or 
	((Chr1_Chr3_Row1116 #= 388) and (Chr1_Chr3_Freq1116 #= 1)) or 
	((Chr1_Chr3_Row1116 #= 389) and (Chr1_Chr3_Freq1116 #= 1)) or 
	((Chr1_Chr3_Row1116 #= 390) and (Chr1_Chr3_Freq1116 #= 1)) or 
	((Chr1_Chr3_Row1116 #= 391) and (Chr1_Chr3_Freq1116 #= 1)) or 
	((Chr1_Chr3_Row1116 #= 392) and (Chr1_Chr3_Freq1116 #= 1)) or 
	((Chr1_Chr3_Row1116 #= 393) and (Chr1_Chr3_Freq1116 #= 1)) or 
	((Chr1_Chr3_Row1116 #= 396) and (Chr1_Chr3_Freq1116 #= 1)) or 
	((Chr1_Chr3_Row1116 #= 397) and (Chr1_Chr3_Freq1116 #= 1)) or 
	((Chr1_Chr3_Row1116 #= 398) and (Chr1_Chr3_Freq1116 #= 1)) or 
	((Chr1_Chr3_Row1116 #= 399) and (Chr1_Chr3_Freq1116 #= 1)) or 
	((Chr1_Chr3_Freq1116 #= 0) and (Chr1_Chr3_Row1116 #= 0)),

	((Chr1_Chr3_Row1117 #= 355) and (Chr1_Chr3_Freq1117 #= 1)) or 
	((Chr1_Chr3_Row1117 #= 356) and (Chr1_Chr3_Freq1117 #= 1)) or 
	((Chr1_Chr3_Row1117 #= 358) and (Chr1_Chr3_Freq1117 #= 1)) or 
	((Chr1_Chr3_Row1117 #= 359) and (Chr1_Chr3_Freq1117 #= 1)) or 
	((Chr1_Chr3_Row1117 #= 360) and (Chr1_Chr3_Freq1117 #= 1)) or 
	((Chr1_Chr3_Row1117 #= 361) and (Chr1_Chr3_Freq1117 #= 1)) or 
	((Chr1_Chr3_Row1117 #= 362) and (Chr1_Chr3_Freq1117 #= 1)) or 
	((Chr1_Chr3_Row1117 #= 363) and (Chr1_Chr3_Freq1117 #= 1)) or 
	((Chr1_Chr3_Row1117 #= 364) and (Chr1_Chr3_Freq1117 #= 1)) or 
	((Chr1_Chr3_Row1117 #= 365) and (Chr1_Chr3_Freq1117 #= 1)) or 
	((Chr1_Chr3_Row1117 #= 366) and (Chr1_Chr3_Freq1117 #= 1)) or 
	((Chr1_Chr3_Row1117 #= 367) and (Chr1_Chr3_Freq1117 #= 2)) or 
	((Chr1_Chr3_Row1117 #= 368) and (Chr1_Chr3_Freq1117 #= 2)) or 
	((Chr1_Chr3_Row1117 #= 369) and (Chr1_Chr3_Freq1117 #= 2)) or 
	((Chr1_Chr3_Row1117 #= 370) and (Chr1_Chr3_Freq1117 #= 3)) or 
	((Chr1_Chr3_Row1117 #= 371) and (Chr1_Chr3_Freq1117 #= 3)) or 
	((Chr1_Chr3_Row1117 #= 372) and (Chr1_Chr3_Freq1117 #= 3)) or 
	((Chr1_Chr3_Row1117 #= 373) and (Chr1_Chr3_Freq1117 #= 4)) or 
	((Chr1_Chr3_Row1117 #= 374) and (Chr1_Chr3_Freq1117 #= 4)) or 
	((Chr1_Chr3_Row1117 #= 380) and (Chr1_Chr3_Freq1117 #= 5)) or 
	((Chr1_Chr3_Row1117 #= 381) and (Chr1_Chr3_Freq1117 #= 4)) or 
	((Chr1_Chr3_Row1117 #= 382) and (Chr1_Chr3_Freq1117 #= 4)) or 
	((Chr1_Chr3_Row1117 #= 383) and (Chr1_Chr3_Freq1117 #= 3)) or 
	((Chr1_Chr3_Row1117 #= 384) and (Chr1_Chr3_Freq1117 #= 2)) or 
	((Chr1_Chr3_Row1117 #= 385) and (Chr1_Chr3_Freq1117 #= 3)) or 
	((Chr1_Chr3_Row1117 #= 386) and (Chr1_Chr3_Freq1117 #= 2)) or 
	((Chr1_Chr3_Row1117 #= 387) and (Chr1_Chr3_Freq1117 #= 2)) or 
	((Chr1_Chr3_Row1117 #= 388) and (Chr1_Chr3_Freq1117 #= 2)) or 
	((Chr1_Chr3_Row1117 #= 389) and (Chr1_Chr3_Freq1117 #= 1)) or 
	((Chr1_Chr3_Row1117 #= 390) and (Chr1_Chr3_Freq1117 #= 1)) or 
	((Chr1_Chr3_Row1117 #= 391) and (Chr1_Chr3_Freq1117 #= 1)) or 
	((Chr1_Chr3_Row1117 #= 392) and (Chr1_Chr3_Freq1117 #= 1)) or 
	((Chr1_Chr3_Row1117 #= 393) and (Chr1_Chr3_Freq1117 #= 1)) or 
	((Chr1_Chr3_Row1117 #= 394) and (Chr1_Chr3_Freq1117 #= 1)) or 
	((Chr1_Chr3_Row1117 #= 395) and (Chr1_Chr3_Freq1117 #= 1)) or 
	((Chr1_Chr3_Row1117 #= 396) and (Chr1_Chr3_Freq1117 #= 1)) or 
	((Chr1_Chr3_Freq1117 #= 0) and (Chr1_Chr3_Row1117 #= 0)),

	((Chr1_Chr3_Row1118 #= 358) and (Chr1_Chr3_Freq1118 #= 1)) or 
	((Chr1_Chr3_Row1118 #= 359) and (Chr1_Chr3_Freq1118 #= 1)) or 
	((Chr1_Chr3_Row1118 #= 360) and (Chr1_Chr3_Freq1118 #= 1)) or 
	((Chr1_Chr3_Row1118 #= 362) and (Chr1_Chr3_Freq1118 #= 1)) or 
	((Chr1_Chr3_Row1118 #= 363) and (Chr1_Chr3_Freq1118 #= 1)) or 
	((Chr1_Chr3_Row1118 #= 364) and (Chr1_Chr3_Freq1118 #= 1)) or 
	((Chr1_Chr3_Row1118 #= 365) and (Chr1_Chr3_Freq1118 #= 1)) or 
	((Chr1_Chr3_Row1118 #= 366) and (Chr1_Chr3_Freq1118 #= 1)) or 
	((Chr1_Chr3_Row1118 #= 367) and (Chr1_Chr3_Freq1118 #= 1)) or 
	((Chr1_Chr3_Row1118 #= 368) and (Chr1_Chr3_Freq1118 #= 1)) or 
	((Chr1_Chr3_Row1118 #= 369) and (Chr1_Chr3_Freq1118 #= 2)) or 
	((Chr1_Chr3_Row1118 #= 370) and (Chr1_Chr3_Freq1118 #= 2)) or 
	((Chr1_Chr3_Row1118 #= 371) and (Chr1_Chr3_Freq1118 #= 2)) or 
	((Chr1_Chr3_Row1118 #= 372) and (Chr1_Chr3_Freq1118 #= 3)) or 
	((Chr1_Chr3_Row1118 #= 373) and (Chr1_Chr3_Freq1118 #= 4)) or 
	((Chr1_Chr3_Row1118 #= 374) and (Chr1_Chr3_Freq1118 #= 5)) or 
	((Chr1_Chr3_Row1118 #= 380) and (Chr1_Chr3_Freq1118 #= 5)) or 
	((Chr1_Chr3_Row1118 #= 381) and (Chr1_Chr3_Freq1118 #= 5)) or 
	((Chr1_Chr3_Row1118 #= 382) and (Chr1_Chr3_Freq1118 #= 3)) or 
	((Chr1_Chr3_Row1118 #= 383) and (Chr1_Chr3_Freq1118 #= 4)) or 
	((Chr1_Chr3_Row1118 #= 384) and (Chr1_Chr3_Freq1118 #= 3)) or 
	((Chr1_Chr3_Row1118 #= 385) and (Chr1_Chr3_Freq1118 #= 3)) or 
	((Chr1_Chr3_Row1118 #= 386) and (Chr1_Chr3_Freq1118 #= 2)) or 
	((Chr1_Chr3_Row1118 #= 387) and (Chr1_Chr3_Freq1118 #= 2)) or 
	((Chr1_Chr3_Row1118 #= 388) and (Chr1_Chr3_Freq1118 #= 1)) or 
	((Chr1_Chr3_Row1118 #= 389) and (Chr1_Chr3_Freq1118 #= 1)) or 
	((Chr1_Chr3_Row1118 #= 390) and (Chr1_Chr3_Freq1118 #= 1)) or 
	((Chr1_Chr3_Row1118 #= 391) and (Chr1_Chr3_Freq1118 #= 1)) or 
	((Chr1_Chr3_Row1118 #= 392) and (Chr1_Chr3_Freq1118 #= 1)) or 
	((Chr1_Chr3_Row1118 #= 394) and (Chr1_Chr3_Freq1118 #= 1)) or 
	((Chr1_Chr3_Row1118 #= 399) and (Chr1_Chr3_Freq1118 #= 1)) or 
	((Chr1_Chr3_Freq1118 #= 0) and (Chr1_Chr3_Row1118 #= 0)),

	((Chr1_Chr3_Row1127 #= 358) and (Chr1_Chr3_Freq1127 #= 1)) or 
	((Chr1_Chr3_Row1127 #= 359) and (Chr1_Chr3_Freq1127 #= 1)) or 
	((Chr1_Chr3_Row1127 #= 360) and (Chr1_Chr3_Freq1127 #= 1)) or 
	((Chr1_Chr3_Row1127 #= 361) and (Chr1_Chr3_Freq1127 #= 1)) or 
	((Chr1_Chr3_Row1127 #= 362) and (Chr1_Chr3_Freq1127 #= 1)) or 
	((Chr1_Chr3_Row1127 #= 363) and (Chr1_Chr3_Freq1127 #= 1)) or 
	((Chr1_Chr3_Row1127 #= 364) and (Chr1_Chr3_Freq1127 #= 1)) or 
	((Chr1_Chr3_Row1127 #= 365) and (Chr1_Chr3_Freq1127 #= 2)) or 
	((Chr1_Chr3_Row1127 #= 366) and (Chr1_Chr3_Freq1127 #= 1)) or 
	((Chr1_Chr3_Row1127 #= 367) and (Chr1_Chr3_Freq1127 #= 1)) or 
	((Chr1_Chr3_Row1127 #= 368) and (Chr1_Chr3_Freq1127 #= 2)) or 
	((Chr1_Chr3_Row1127 #= 369) and (Chr1_Chr3_Freq1127 #= 2)) or 
	((Chr1_Chr3_Row1127 #= 370) and (Chr1_Chr3_Freq1127 #= 3)) or 
	((Chr1_Chr3_Row1127 #= 371) and (Chr1_Chr3_Freq1127 #= 2)) or 
	((Chr1_Chr3_Row1127 #= 372) and (Chr1_Chr3_Freq1127 #= 4)) or 
	((Chr1_Chr3_Row1127 #= 373) and (Chr1_Chr3_Freq1127 #= 5)) or 
	((Chr1_Chr3_Row1127 #= 374) and (Chr1_Chr3_Freq1127 #= 6)) or 
	((Chr1_Chr3_Row1127 #= 380) and (Chr1_Chr3_Freq1127 #= 7)) or 
	((Chr1_Chr3_Row1127 #= 381) and (Chr1_Chr3_Freq1127 #= 6)) or 
	((Chr1_Chr3_Row1127 #= 382) and (Chr1_Chr3_Freq1127 #= 6)) or 
	((Chr1_Chr3_Row1127 #= 383) and (Chr1_Chr3_Freq1127 #= 4)) or 
	((Chr1_Chr3_Row1127 #= 384) and (Chr1_Chr3_Freq1127 #= 4)) or 
	((Chr1_Chr3_Row1127 #= 385) and (Chr1_Chr3_Freq1127 #= 3)) or 
	((Chr1_Chr3_Row1127 #= 386) and (Chr1_Chr3_Freq1127 #= 2)) or 
	((Chr1_Chr3_Row1127 #= 387) and (Chr1_Chr3_Freq1127 #= 2)) or 
	((Chr1_Chr3_Row1127 #= 388) and (Chr1_Chr3_Freq1127 #= 2)) or 
	((Chr1_Chr3_Row1127 #= 389) and (Chr1_Chr3_Freq1127 #= 1)) or 
	((Chr1_Chr3_Row1127 #= 390) and (Chr1_Chr3_Freq1127 #= 1)) or 
	((Chr1_Chr3_Row1127 #= 391) and (Chr1_Chr3_Freq1127 #= 1)) or 
	((Chr1_Chr3_Row1127 #= 392) and (Chr1_Chr3_Freq1127 #= 1)) or 
	((Chr1_Chr3_Row1127 #= 394) and (Chr1_Chr3_Freq1127 #= 1)) or 
	((Chr1_Chr3_Row1127 #= 395) and (Chr1_Chr3_Freq1127 #= 1)) or 
	((Chr1_Chr3_Row1127 #= 396) and (Chr1_Chr3_Freq1127 #= 1)) or 
	((Chr1_Chr3_Freq1127 #= 0) and (Chr1_Chr3_Row1127 #= 0)),

	((Chr1_Chr3_Row1128 #= 356) and (Chr1_Chr3_Freq1128 #= 1)) or 
	((Chr1_Chr3_Row1128 #= 357) and (Chr1_Chr3_Freq1128 #= 1)) or 
	((Chr1_Chr3_Row1128 #= 358) and (Chr1_Chr3_Freq1128 #= 1)) or 
	((Chr1_Chr3_Row1128 #= 359) and (Chr1_Chr3_Freq1128 #= 1)) or 
	((Chr1_Chr3_Row1128 #= 360) and (Chr1_Chr3_Freq1128 #= 1)) or 
	((Chr1_Chr3_Row1128 #= 361) and (Chr1_Chr3_Freq1128 #= 1)) or 
	((Chr1_Chr3_Row1128 #= 362) and (Chr1_Chr3_Freq1128 #= 1)) or 
	((Chr1_Chr3_Row1128 #= 363) and (Chr1_Chr3_Freq1128 #= 1)) or 
	((Chr1_Chr3_Row1128 #= 364) and (Chr1_Chr3_Freq1128 #= 1)) or 
	((Chr1_Chr3_Row1128 #= 365) and (Chr1_Chr3_Freq1128 #= 1)) or 
	((Chr1_Chr3_Row1128 #= 366) and (Chr1_Chr3_Freq1128 #= 2)) or 
	((Chr1_Chr3_Row1128 #= 367) and (Chr1_Chr3_Freq1128 #= 2)) or 
	((Chr1_Chr3_Row1128 #= 368) and (Chr1_Chr3_Freq1128 #= 2)) or 
	((Chr1_Chr3_Row1128 #= 369) and (Chr1_Chr3_Freq1128 #= 3)) or 
	((Chr1_Chr3_Row1128 #= 370) and (Chr1_Chr3_Freq1128 #= 3)) or 
	((Chr1_Chr3_Row1128 #= 371) and (Chr1_Chr3_Freq1128 #= 3)) or 
	((Chr1_Chr3_Row1128 #= 372) and (Chr1_Chr3_Freq1128 #= 4)) or 
	((Chr1_Chr3_Row1128 #= 373) and (Chr1_Chr3_Freq1128 #= 6)) or 
	((Chr1_Chr3_Row1128 #= 374) and (Chr1_Chr3_Freq1128 #= 6)) or 
	((Chr1_Chr3_Row1128 #= 380) and (Chr1_Chr3_Freq1128 #= 6)) or 
	((Chr1_Chr3_Row1128 #= 381) and (Chr1_Chr3_Freq1128 #= 7)) or 
	((Chr1_Chr3_Row1128 #= 382) and (Chr1_Chr3_Freq1128 #= 5)) or 
	((Chr1_Chr3_Row1128 #= 383) and (Chr1_Chr3_Freq1128 #= 4)) or 
	((Chr1_Chr3_Row1128 #= 384) and (Chr1_Chr3_Freq1128 #= 3)) or 
	((Chr1_Chr3_Row1128 #= 385) and (Chr1_Chr3_Freq1128 #= 3)) or 
	((Chr1_Chr3_Row1128 #= 386) and (Chr1_Chr3_Freq1128 #= 3)) or 
	((Chr1_Chr3_Row1128 #= 387) and (Chr1_Chr3_Freq1128 #= 2)) or 
	((Chr1_Chr3_Row1128 #= 388) and (Chr1_Chr3_Freq1128 #= 3)) or 
	((Chr1_Chr3_Row1128 #= 389) and (Chr1_Chr3_Freq1128 #= 2)) or 
	((Chr1_Chr3_Row1128 #= 390) and (Chr1_Chr3_Freq1128 #= 1)) or 
	((Chr1_Chr3_Row1128 #= 391) and (Chr1_Chr3_Freq1128 #= 1)) or 
	((Chr1_Chr3_Row1128 #= 392) and (Chr1_Chr3_Freq1128 #= 1)) or 
	((Chr1_Chr3_Row1128 #= 393) and (Chr1_Chr3_Freq1128 #= 1)) or 
	((Chr1_Chr3_Row1128 #= 394) and (Chr1_Chr3_Freq1128 #= 1)) or 
	((Chr1_Chr3_Row1128 #= 395) and (Chr1_Chr3_Freq1128 #= 1)) or 
	((Chr1_Chr3_Row1128 #= 396) and (Chr1_Chr3_Freq1128 #= 1)) or 
	((Chr1_Chr3_Row1128 #= 401) and (Chr1_Chr3_Freq1128 #= 1)) or 
	((Chr1_Chr3_Freq1128 #= 0) and (Chr1_Chr3_Row1128 #= 0)),

	((Chr1_Chr3_Row1129 #= 358) and (Chr1_Chr3_Freq1129 #= 1)) or 
	((Chr1_Chr3_Row1129 #= 359) and (Chr1_Chr3_Freq1129 #= 1)) or 
	((Chr1_Chr3_Row1129 #= 360) and (Chr1_Chr3_Freq1129 #= 1)) or 
	((Chr1_Chr3_Row1129 #= 361) and (Chr1_Chr3_Freq1129 #= 1)) or 
	((Chr1_Chr3_Row1129 #= 362) and (Chr1_Chr3_Freq1129 #= 1)) or 
	((Chr1_Chr3_Row1129 #= 363) and (Chr1_Chr3_Freq1129 #= 1)) or 
	((Chr1_Chr3_Row1129 #= 364) and (Chr1_Chr3_Freq1129 #= 1)) or 
	((Chr1_Chr3_Row1129 #= 365) and (Chr1_Chr3_Freq1129 #= 1)) or 
	((Chr1_Chr3_Row1129 #= 366) and (Chr1_Chr3_Freq1129 #= 2)) or 
	((Chr1_Chr3_Row1129 #= 367) and (Chr1_Chr3_Freq1129 #= 1)) or 
	((Chr1_Chr3_Row1129 #= 368) and (Chr1_Chr3_Freq1129 #= 2)) or 
	((Chr1_Chr3_Row1129 #= 369) and (Chr1_Chr3_Freq1129 #= 2)) or 
	((Chr1_Chr3_Row1129 #= 370) and (Chr1_Chr3_Freq1129 #= 3)) or 
	((Chr1_Chr3_Row1129 #= 371) and (Chr1_Chr3_Freq1129 #= 3)) or 
	((Chr1_Chr3_Row1129 #= 372) and (Chr1_Chr3_Freq1129 #= 3)) or 
	((Chr1_Chr3_Row1129 #= 373) and (Chr1_Chr3_Freq1129 #= 5)) or 
	((Chr1_Chr3_Row1129 #= 374) and (Chr1_Chr3_Freq1129 #= 5)) or 
	((Chr1_Chr3_Row1129 #= 380) and (Chr1_Chr3_Freq1129 #= 5)) or 
	((Chr1_Chr3_Row1129 #= 381) and (Chr1_Chr3_Freq1129 #= 5)) or 
	((Chr1_Chr3_Row1129 #= 382) and (Chr1_Chr3_Freq1129 #= 5)) or 
	((Chr1_Chr3_Row1129 #= 383) and (Chr1_Chr3_Freq1129 #= 3)) or 
	((Chr1_Chr3_Row1129 #= 384) and (Chr1_Chr3_Freq1129 #= 3)) or 
	((Chr1_Chr3_Row1129 #= 385) and (Chr1_Chr3_Freq1129 #= 2)) or 
	((Chr1_Chr3_Row1129 #= 386) and (Chr1_Chr3_Freq1129 #= 3)) or 
	((Chr1_Chr3_Row1129 #= 387) and (Chr1_Chr3_Freq1129 #= 2)) or 
	((Chr1_Chr3_Row1129 #= 388) and (Chr1_Chr3_Freq1129 #= 2)) or 
	((Chr1_Chr3_Row1129 #= 389) and (Chr1_Chr3_Freq1129 #= 1)) or 
	((Chr1_Chr3_Row1129 #= 390) and (Chr1_Chr3_Freq1129 #= 1)) or 
	((Chr1_Chr3_Row1129 #= 391) and (Chr1_Chr3_Freq1129 #= 1)) or 
	((Chr1_Chr3_Row1129 #= 392) and (Chr1_Chr3_Freq1129 #= 1)) or 
	((Chr1_Chr3_Row1129 #= 393) and (Chr1_Chr3_Freq1129 #= 1)) or 
	((Chr1_Chr3_Row1129 #= 394) and (Chr1_Chr3_Freq1129 #= 1)) or 
	((Chr1_Chr3_Row1129 #= 395) and (Chr1_Chr3_Freq1129 #= 1)) or 
	((Chr1_Chr3_Row1129 #= 396) and (Chr1_Chr3_Freq1129 #= 1)) or 
	((Chr1_Chr3_Freq1129 #= 0) and (Chr1_Chr3_Row1129 #= 0)),

	((Chr1_Chr3_Row1130 #= 358) and (Chr1_Chr3_Freq1130 #= 1)) or 
	((Chr1_Chr3_Row1130 #= 359) and (Chr1_Chr3_Freq1130 #= 1)) or 
	((Chr1_Chr3_Row1130 #= 360) and (Chr1_Chr3_Freq1130 #= 1)) or 
	((Chr1_Chr3_Row1130 #= 361) and (Chr1_Chr3_Freq1130 #= 1)) or 
	((Chr1_Chr3_Row1130 #= 362) and (Chr1_Chr3_Freq1130 #= 1)) or 
	((Chr1_Chr3_Row1130 #= 363) and (Chr1_Chr3_Freq1130 #= 1)) or 
	((Chr1_Chr3_Row1130 #= 364) and (Chr1_Chr3_Freq1130 #= 1)) or 
	((Chr1_Chr3_Row1130 #= 365) and (Chr1_Chr3_Freq1130 #= 1)) or 
	((Chr1_Chr3_Row1130 #= 366) and (Chr1_Chr3_Freq1130 #= 2)) or 
	((Chr1_Chr3_Row1130 #= 367) and (Chr1_Chr3_Freq1130 #= 1)) or 
	((Chr1_Chr3_Row1130 #= 368) and (Chr1_Chr3_Freq1130 #= 1)) or 
	((Chr1_Chr3_Row1130 #= 369) and (Chr1_Chr3_Freq1130 #= 2)) or 
	((Chr1_Chr3_Row1130 #= 370) and (Chr1_Chr3_Freq1130 #= 2)) or 
	((Chr1_Chr3_Row1130 #= 371) and (Chr1_Chr3_Freq1130 #= 2)) or 
	((Chr1_Chr3_Row1130 #= 372) and (Chr1_Chr3_Freq1130 #= 3)) or 
	((Chr1_Chr3_Row1130 #= 373) and (Chr1_Chr3_Freq1130 #= 3)) or 
	((Chr1_Chr3_Row1130 #= 374) and (Chr1_Chr3_Freq1130 #= 4)) or 
	((Chr1_Chr3_Row1130 #= 380) and (Chr1_Chr3_Freq1130 #= 4)) or 
	((Chr1_Chr3_Row1130 #= 381) and (Chr1_Chr3_Freq1130 #= 4)) or 
	((Chr1_Chr3_Row1130 #= 382) and (Chr1_Chr3_Freq1130 #= 4)) or 
	((Chr1_Chr3_Row1130 #= 383) and (Chr1_Chr3_Freq1130 #= 3)) or 
	((Chr1_Chr3_Row1130 #= 384) and (Chr1_Chr3_Freq1130 #= 2)) or 
	((Chr1_Chr3_Row1130 #= 385) and (Chr1_Chr3_Freq1130 #= 2)) or 
	((Chr1_Chr3_Row1130 #= 386) and (Chr1_Chr3_Freq1130 #= 2)) or 
	((Chr1_Chr3_Row1130 #= 387) and (Chr1_Chr3_Freq1130 #= 2)) or 
	((Chr1_Chr3_Row1130 #= 388) and (Chr1_Chr3_Freq1130 #= 1)) or 
	((Chr1_Chr3_Row1130 #= 389) and (Chr1_Chr3_Freq1130 #= 1)) or 
	((Chr1_Chr3_Row1130 #= 390) and (Chr1_Chr3_Freq1130 #= 1)) or 
	((Chr1_Chr3_Row1130 #= 391) and (Chr1_Chr3_Freq1130 #= 1)) or 
	((Chr1_Chr3_Row1130 #= 394) and (Chr1_Chr3_Freq1130 #= 1)) or 
	((Chr1_Chr3_Row1130 #= 395) and (Chr1_Chr3_Freq1130 #= 1)) or 
	((Chr1_Chr3_Freq1130 #= 0) and (Chr1_Chr3_Row1130 #= 0)),

	((Chr1_Chr3_Row1131 #= 358) and (Chr1_Chr3_Freq1131 #= 1)) or 
	((Chr1_Chr3_Row1131 #= 359) and (Chr1_Chr3_Freq1131 #= 1)) or 
	((Chr1_Chr3_Row1131 #= 361) and (Chr1_Chr3_Freq1131 #= 1)) or 
	((Chr1_Chr3_Row1131 #= 362) and (Chr1_Chr3_Freq1131 #= 1)) or 
	((Chr1_Chr3_Row1131 #= 363) and (Chr1_Chr3_Freq1131 #= 1)) or 
	((Chr1_Chr3_Row1131 #= 364) and (Chr1_Chr3_Freq1131 #= 1)) or 
	((Chr1_Chr3_Row1131 #= 365) and (Chr1_Chr3_Freq1131 #= 1)) or 
	((Chr1_Chr3_Row1131 #= 366) and (Chr1_Chr3_Freq1131 #= 1)) or 
	((Chr1_Chr3_Row1131 #= 367) and (Chr1_Chr3_Freq1131 #= 1)) or 
	((Chr1_Chr3_Row1131 #= 368) and (Chr1_Chr3_Freq1131 #= 1)) or 
	((Chr1_Chr3_Row1131 #= 369) and (Chr1_Chr3_Freq1131 #= 2)) or 
	((Chr1_Chr3_Row1131 #= 370) and (Chr1_Chr3_Freq1131 #= 2)) or 
	((Chr1_Chr3_Row1131 #= 371) and (Chr1_Chr3_Freq1131 #= 2)) or 
	((Chr1_Chr3_Row1131 #= 372) and (Chr1_Chr3_Freq1131 #= 2)) or 
	((Chr1_Chr3_Row1131 #= 373) and (Chr1_Chr3_Freq1131 #= 3)) or 
	((Chr1_Chr3_Row1131 #= 374) and (Chr1_Chr3_Freq1131 #= 4)) or 
	((Chr1_Chr3_Row1131 #= 380) and (Chr1_Chr3_Freq1131 #= 3)) or 
	((Chr1_Chr3_Row1131 #= 381) and (Chr1_Chr3_Freq1131 #= 3)) or 
	((Chr1_Chr3_Row1131 #= 382) and (Chr1_Chr3_Freq1131 #= 3)) or 
	((Chr1_Chr3_Row1131 #= 383) and (Chr1_Chr3_Freq1131 #= 3)) or 
	((Chr1_Chr3_Row1131 #= 384) and (Chr1_Chr3_Freq1131 #= 2)) or 
	((Chr1_Chr3_Row1131 #= 385) and (Chr1_Chr3_Freq1131 #= 2)) or 
	((Chr1_Chr3_Row1131 #= 386) and (Chr1_Chr3_Freq1131 #= 1)) or 
	((Chr1_Chr3_Row1131 #= 387) and (Chr1_Chr3_Freq1131 #= 1)) or 
	((Chr1_Chr3_Row1131 #= 388) and (Chr1_Chr3_Freq1131 #= 2)) or 
	((Chr1_Chr3_Row1131 #= 389) and (Chr1_Chr3_Freq1131 #= 1)) or 
	((Chr1_Chr3_Row1131 #= 390) and (Chr1_Chr3_Freq1131 #= 1)) or 
	((Chr1_Chr3_Row1131 #= 391) and (Chr1_Chr3_Freq1131 #= 1)) or 
	((Chr1_Chr3_Row1131 #= 392) and (Chr1_Chr3_Freq1131 #= 1)) or 
	((Chr1_Chr3_Row1131 #= 393) and (Chr1_Chr3_Freq1131 #= 1)) or 
	((Chr1_Chr3_Row1131 #= 394) and (Chr1_Chr3_Freq1131 #= 1)) or 
	((Chr1_Chr3_Row1131 #= 395) and (Chr1_Chr3_Freq1131 #= 1)) or 
	((Chr1_Chr3_Freq1131 #= 0) and (Chr1_Chr3_Row1131 #= 0)),

	((Chr1_Chr3_Row1132 #= 356) and (Chr1_Chr3_Freq1132 #= 1)) or 
	((Chr1_Chr3_Row1132 #= 357) and (Chr1_Chr3_Freq1132 #= 1)) or 
	((Chr1_Chr3_Row1132 #= 358) and (Chr1_Chr3_Freq1132 #= 1)) or 
	((Chr1_Chr3_Row1132 #= 359) and (Chr1_Chr3_Freq1132 #= 1)) or 
	((Chr1_Chr3_Row1132 #= 360) and (Chr1_Chr3_Freq1132 #= 1)) or 
	((Chr1_Chr3_Row1132 #= 362) and (Chr1_Chr3_Freq1132 #= 1)) or 
	((Chr1_Chr3_Row1132 #= 363) and (Chr1_Chr3_Freq1132 #= 1)) or 
	((Chr1_Chr3_Row1132 #= 364) and (Chr1_Chr3_Freq1132 #= 1)) or 
	((Chr1_Chr3_Row1132 #= 365) and (Chr1_Chr3_Freq1132 #= 1)) or 
	((Chr1_Chr3_Row1132 #= 366) and (Chr1_Chr3_Freq1132 #= 1)) or 
	((Chr1_Chr3_Row1132 #= 367) and (Chr1_Chr3_Freq1132 #= 2)) or 
	((Chr1_Chr3_Row1132 #= 368) and (Chr1_Chr3_Freq1132 #= 1)) or 
	((Chr1_Chr3_Row1132 #= 369) and (Chr1_Chr3_Freq1132 #= 2)) or 
	((Chr1_Chr3_Row1132 #= 370) and (Chr1_Chr3_Freq1132 #= 2)) or 
	((Chr1_Chr3_Row1132 #= 371) and (Chr1_Chr3_Freq1132 #= 2)) or 
	((Chr1_Chr3_Row1132 #= 372) and (Chr1_Chr3_Freq1132 #= 3)) or 
	((Chr1_Chr3_Row1132 #= 373) and (Chr1_Chr3_Freq1132 #= 3)) or 
	((Chr1_Chr3_Row1132 #= 374) and (Chr1_Chr3_Freq1132 #= 3)) or 
	((Chr1_Chr3_Row1132 #= 380) and (Chr1_Chr3_Freq1132 #= 3)) or 
	((Chr1_Chr3_Row1132 #= 381) and (Chr1_Chr3_Freq1132 #= 4)) or 
	((Chr1_Chr3_Row1132 #= 382) and (Chr1_Chr3_Freq1132 #= 3)) or 
	((Chr1_Chr3_Row1132 #= 383) and (Chr1_Chr3_Freq1132 #= 2)) or 
	((Chr1_Chr3_Row1132 #= 384) and (Chr1_Chr3_Freq1132 #= 3)) or 
	((Chr1_Chr3_Row1132 #= 385) and (Chr1_Chr3_Freq1132 #= 2)) or 
	((Chr1_Chr3_Row1132 #= 386) and (Chr1_Chr3_Freq1132 #= 1)) or 
	((Chr1_Chr3_Row1132 #= 387) and (Chr1_Chr3_Freq1132 #= 1)) or 
	((Chr1_Chr3_Row1132 #= 388) and (Chr1_Chr3_Freq1132 #= 1)) or 
	((Chr1_Chr3_Row1132 #= 389) and (Chr1_Chr3_Freq1132 #= 1)) or 
	((Chr1_Chr3_Row1132 #= 390) and (Chr1_Chr3_Freq1132 #= 1)) or 
	((Chr1_Chr3_Row1132 #= 391) and (Chr1_Chr3_Freq1132 #= 1)) or 
	((Chr1_Chr3_Row1132 #= 392) and (Chr1_Chr3_Freq1132 #= 1)) or 
	((Chr1_Chr3_Row1132 #= 393) and (Chr1_Chr3_Freq1132 #= 1)) or 
	((Chr1_Chr3_Row1132 #= 394) and (Chr1_Chr3_Freq1132 #= 1)) or 
	((Chr1_Chr3_Row1132 #= 395) and (Chr1_Chr3_Freq1132 #= 1)) or 
	((Chr1_Chr3_Row1132 #= 403) and (Chr1_Chr3_Freq1132 #= 1)) or 
	((Chr1_Chr3_Freq1132 #= 0) and (Chr1_Chr3_Row1132 #= 0)),

	((Chr1_Chr3_Row1133 #= 357) and (Chr1_Chr3_Freq1133 #= 1)) or 
	((Chr1_Chr3_Row1133 #= 358) and (Chr1_Chr3_Freq1133 #= 1)) or 
	((Chr1_Chr3_Row1133 #= 359) and (Chr1_Chr3_Freq1133 #= 1)) or 
	((Chr1_Chr3_Row1133 #= 360) and (Chr1_Chr3_Freq1133 #= 1)) or 
	((Chr1_Chr3_Row1133 #= 361) and (Chr1_Chr3_Freq1133 #= 1)) or 
	((Chr1_Chr3_Row1133 #= 362) and (Chr1_Chr3_Freq1133 #= 1)) or 
	((Chr1_Chr3_Row1133 #= 363) and (Chr1_Chr3_Freq1133 #= 1)) or 
	((Chr1_Chr3_Row1133 #= 364) and (Chr1_Chr3_Freq1133 #= 1)) or 
	((Chr1_Chr3_Row1133 #= 365) and (Chr1_Chr3_Freq1133 #= 1)) or 
	((Chr1_Chr3_Row1133 #= 366) and (Chr1_Chr3_Freq1133 #= 1)) or 
	((Chr1_Chr3_Row1133 #= 367) and (Chr1_Chr3_Freq1133 #= 1)) or 
	((Chr1_Chr3_Row1133 #= 368) and (Chr1_Chr3_Freq1133 #= 1)) or 
	((Chr1_Chr3_Row1133 #= 369) and (Chr1_Chr3_Freq1133 #= 2)) or 
	((Chr1_Chr3_Row1133 #= 370) and (Chr1_Chr3_Freq1133 #= 2)) or 
	((Chr1_Chr3_Row1133 #= 371) and (Chr1_Chr3_Freq1133 #= 2)) or 
	((Chr1_Chr3_Row1133 #= 372) and (Chr1_Chr3_Freq1133 #= 2)) or 
	((Chr1_Chr3_Row1133 #= 373) and (Chr1_Chr3_Freq1133 #= 2)) or 
	((Chr1_Chr3_Row1133 #= 374) and (Chr1_Chr3_Freq1133 #= 3)) or 
	((Chr1_Chr3_Row1133 #= 380) and (Chr1_Chr3_Freq1133 #= 3)) or 
	((Chr1_Chr3_Row1133 #= 381) and (Chr1_Chr3_Freq1133 #= 3)) or 
	((Chr1_Chr3_Row1133 #= 382) and (Chr1_Chr3_Freq1133 #= 2)) or 
	((Chr1_Chr3_Row1133 #= 383) and (Chr1_Chr3_Freq1133 #= 3)) or 
	((Chr1_Chr3_Row1133 #= 384) and (Chr1_Chr3_Freq1133 #= 2)) or 
	((Chr1_Chr3_Row1133 #= 385) and (Chr1_Chr3_Freq1133 #= 2)) or 
	((Chr1_Chr3_Row1133 #= 386) and (Chr1_Chr3_Freq1133 #= 2)) or 
	((Chr1_Chr3_Row1133 #= 387) and (Chr1_Chr3_Freq1133 #= 1)) or 
	((Chr1_Chr3_Row1133 #= 388) and (Chr1_Chr3_Freq1133 #= 1)) or 
	((Chr1_Chr3_Row1133 #= 389) and (Chr1_Chr3_Freq1133 #= 1)) or 
	((Chr1_Chr3_Row1133 #= 390) and (Chr1_Chr3_Freq1133 #= 1)) or 
	((Chr1_Chr3_Row1133 #= 391) and (Chr1_Chr3_Freq1133 #= 1)) or 
	((Chr1_Chr3_Row1133 #= 392) and (Chr1_Chr3_Freq1133 #= 1)) or 
	((Chr1_Chr3_Row1133 #= 393) and (Chr1_Chr3_Freq1133 #= 1)) or 
	((Chr1_Chr3_Row1133 #= 394) and (Chr1_Chr3_Freq1133 #= 1)) or 
	((Chr1_Chr3_Row1133 #= 395) and (Chr1_Chr3_Freq1133 #= 1)) or 
	((Chr1_Chr3_Row1133 #= 396) and (Chr1_Chr3_Freq1133 #= 1)) or 
	((Chr1_Chr3_Freq1133 #= 0) and (Chr1_Chr3_Row1133 #= 0)),

	((Chr1_Chr3_Row1134 #= 356) and (Chr1_Chr3_Freq1134 #= 1)) or 
	((Chr1_Chr3_Row1134 #= 358) and (Chr1_Chr3_Freq1134 #= 1)) or 
	((Chr1_Chr3_Row1134 #= 359) and (Chr1_Chr3_Freq1134 #= 1)) or 
	((Chr1_Chr3_Row1134 #= 360) and (Chr1_Chr3_Freq1134 #= 1)) or 
	((Chr1_Chr3_Row1134 #= 361) and (Chr1_Chr3_Freq1134 #= 1)) or 
	((Chr1_Chr3_Row1134 #= 362) and (Chr1_Chr3_Freq1134 #= 1)) or 
	((Chr1_Chr3_Row1134 #= 363) and (Chr1_Chr3_Freq1134 #= 1)) or 
	((Chr1_Chr3_Row1134 #= 364) and (Chr1_Chr3_Freq1134 #= 1)) or 
	((Chr1_Chr3_Row1134 #= 365) and (Chr1_Chr3_Freq1134 #= 1)) or 
	((Chr1_Chr3_Row1134 #= 366) and (Chr1_Chr3_Freq1134 #= 1)) or 
	((Chr1_Chr3_Row1134 #= 367) and (Chr1_Chr3_Freq1134 #= 1)) or 
	((Chr1_Chr3_Row1134 #= 368) and (Chr1_Chr3_Freq1134 #= 1)) or 
	((Chr1_Chr3_Row1134 #= 369) and (Chr1_Chr3_Freq1134 #= 2)) or 
	((Chr1_Chr3_Row1134 #= 370) and (Chr1_Chr3_Freq1134 #= 2)) or 
	((Chr1_Chr3_Row1134 #= 371) and (Chr1_Chr3_Freq1134 #= 2)) or 
	((Chr1_Chr3_Row1134 #= 372) and (Chr1_Chr3_Freq1134 #= 2)) or 
	((Chr1_Chr3_Row1134 #= 373) and (Chr1_Chr3_Freq1134 #= 3)) or 
	((Chr1_Chr3_Row1134 #= 374) and (Chr1_Chr3_Freq1134 #= 3)) or 
	((Chr1_Chr3_Row1134 #= 380) and (Chr1_Chr3_Freq1134 #= 2)) or 
	((Chr1_Chr3_Row1134 #= 381) and (Chr1_Chr3_Freq1134 #= 3)) or 
	((Chr1_Chr3_Row1134 #= 382) and (Chr1_Chr3_Freq1134 #= 2)) or 
	((Chr1_Chr3_Row1134 #= 383) and (Chr1_Chr3_Freq1134 #= 2)) or 
	((Chr1_Chr3_Row1134 #= 384) and (Chr1_Chr3_Freq1134 #= 2)) or 
	((Chr1_Chr3_Row1134 #= 385) and (Chr1_Chr3_Freq1134 #= 2)) or 
	((Chr1_Chr3_Row1134 #= 386) and (Chr1_Chr3_Freq1134 #= 2)) or 
	((Chr1_Chr3_Row1134 #= 387) and (Chr1_Chr3_Freq1134 #= 1)) or 
	((Chr1_Chr3_Row1134 #= 388) and (Chr1_Chr3_Freq1134 #= 1)) or 
	((Chr1_Chr3_Row1134 #= 389) and (Chr1_Chr3_Freq1134 #= 1)) or 
	((Chr1_Chr3_Row1134 #= 390) and (Chr1_Chr3_Freq1134 #= 1)) or 
	((Chr1_Chr3_Row1134 #= 391) and (Chr1_Chr3_Freq1134 #= 1)) or 
	((Chr1_Chr3_Row1134 #= 392) and (Chr1_Chr3_Freq1134 #= 1)) or 
	((Chr1_Chr3_Row1134 #= 393) and (Chr1_Chr3_Freq1134 #= 1)) or 
	((Chr1_Chr3_Row1134 #= 394) and (Chr1_Chr3_Freq1134 #= 1)) or 
	((Chr1_Chr3_Row1134 #= 395) and (Chr1_Chr3_Freq1134 #= 1)) or 
	((Chr1_Chr3_Freq1134 #= 0) and (Chr1_Chr3_Row1134 #= 0)),

	((Chr1_Chr3_Row1135 #= 356) and (Chr1_Chr3_Freq1135 #= 1)) or 
	((Chr1_Chr3_Row1135 #= 358) and (Chr1_Chr3_Freq1135 #= 1)) or 
	((Chr1_Chr3_Row1135 #= 360) and (Chr1_Chr3_Freq1135 #= 1)) or 
	((Chr1_Chr3_Row1135 #= 363) and (Chr1_Chr3_Freq1135 #= 1)) or 
	((Chr1_Chr3_Row1135 #= 364) and (Chr1_Chr3_Freq1135 #= 1)) or 
	((Chr1_Chr3_Row1135 #= 365) and (Chr1_Chr3_Freq1135 #= 1)) or 
	((Chr1_Chr3_Row1135 #= 366) and (Chr1_Chr3_Freq1135 #= 1)) or 
	((Chr1_Chr3_Row1135 #= 367) and (Chr1_Chr3_Freq1135 #= 1)) or 
	((Chr1_Chr3_Row1135 #= 368) and (Chr1_Chr3_Freq1135 #= 1)) or 
	((Chr1_Chr3_Row1135 #= 369) and (Chr1_Chr3_Freq1135 #= 1)) or 
	((Chr1_Chr3_Row1135 #= 370) and (Chr1_Chr3_Freq1135 #= 1)) or 
	((Chr1_Chr3_Row1135 #= 371) and (Chr1_Chr3_Freq1135 #= 2)) or 
	((Chr1_Chr3_Row1135 #= 372) and (Chr1_Chr3_Freq1135 #= 2)) or 
	((Chr1_Chr3_Row1135 #= 373) and (Chr1_Chr3_Freq1135 #= 2)) or 
	((Chr1_Chr3_Row1135 #= 374) and (Chr1_Chr3_Freq1135 #= 2)) or 
	((Chr1_Chr3_Row1135 #= 380) and (Chr1_Chr3_Freq1135 #= 2)) or 
	((Chr1_Chr3_Row1135 #= 381) and (Chr1_Chr3_Freq1135 #= 2)) or 
	((Chr1_Chr3_Row1135 #= 382) and (Chr1_Chr3_Freq1135 #= 2)) or 
	((Chr1_Chr3_Row1135 #= 383) and (Chr1_Chr3_Freq1135 #= 2)) or 
	((Chr1_Chr3_Row1135 #= 384) and (Chr1_Chr3_Freq1135 #= 1)) or 
	((Chr1_Chr3_Row1135 #= 385) and (Chr1_Chr3_Freq1135 #= 1)) or 
	((Chr1_Chr3_Row1135 #= 386) and (Chr1_Chr3_Freq1135 #= 1)) or 
	((Chr1_Chr3_Row1135 #= 387) and (Chr1_Chr3_Freq1135 #= 1)) or 
	((Chr1_Chr3_Row1135 #= 388) and (Chr1_Chr3_Freq1135 #= 1)) or 
	((Chr1_Chr3_Row1135 #= 389) and (Chr1_Chr3_Freq1135 #= 1)) or 
	((Chr1_Chr3_Row1135 #= 390) and (Chr1_Chr3_Freq1135 #= 1)) or 
	((Chr1_Chr3_Row1135 #= 392) and (Chr1_Chr3_Freq1135 #= 1)) or 
	((Chr1_Chr3_Freq1135 #= 0) and (Chr1_Chr3_Row1135 #= 0)),

	((Chr1_Chr3_Row1136 #= 349) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 352) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 356) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 358) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 359) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 360) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 361) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 363) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 364) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 365) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 366) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 367) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 368) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 369) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 370) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 371) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 372) and (Chr1_Chr3_Freq1136 #= 2)) or 
	((Chr1_Chr3_Row1136 #= 373) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 374) and (Chr1_Chr3_Freq1136 #= 2)) or 
	((Chr1_Chr3_Row1136 #= 380) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 381) and (Chr1_Chr3_Freq1136 #= 2)) or 
	((Chr1_Chr3_Row1136 #= 382) and (Chr1_Chr3_Freq1136 #= 2)) or 
	((Chr1_Chr3_Row1136 #= 383) and (Chr1_Chr3_Freq1136 #= 2)) or 
	((Chr1_Chr3_Row1136 #= 384) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 385) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 386) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 387) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 388) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 389) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 391) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 393) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Row1136 #= 399) and (Chr1_Chr3_Freq1136 #= 1)) or 
	((Chr1_Chr3_Freq1136 #= 0) and (Chr1_Chr3_Row1136 #= 0)),

	((Chr1_Chr3_Row1137 #= 356) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 357) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 358) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 360) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 361) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 362) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 363) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 364) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 365) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 366) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 367) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 368) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 369) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 370) and (Chr1_Chr3_Freq1137 #= 2)) or 
	((Chr1_Chr3_Row1137 #= 371) and (Chr1_Chr3_Freq1137 #= 2)) or 
	((Chr1_Chr3_Row1137 #= 372) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 373) and (Chr1_Chr3_Freq1137 #= 2)) or 
	((Chr1_Chr3_Row1137 #= 374) and (Chr1_Chr3_Freq1137 #= 2)) or 
	((Chr1_Chr3_Row1137 #= 380) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 381) and (Chr1_Chr3_Freq1137 #= 2)) or 
	((Chr1_Chr3_Row1137 #= 382) and (Chr1_Chr3_Freq1137 #= 2)) or 
	((Chr1_Chr3_Row1137 #= 383) and (Chr1_Chr3_Freq1137 #= 2)) or 
	((Chr1_Chr3_Row1137 #= 384) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 385) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 386) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 387) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 388) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 389) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 390) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 391) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 392) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 394) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 395) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 396) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 397) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 399) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Row1137 #= 400) and (Chr1_Chr3_Freq1137 #= 1)) or 
	((Chr1_Chr3_Freq1137 #= 0) and (Chr1_Chr3_Row1137 #= 0)),

	((Chr1_Chr3_Row1138 #= 352) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 356) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 358) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 359) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 360) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 362) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 363) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 364) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 365) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 366) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 367) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 368) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 369) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 370) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 371) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 372) and (Chr1_Chr3_Freq1138 #= 2)) or 
	((Chr1_Chr3_Row1138 #= 373) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 374) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 380) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 381) and (Chr1_Chr3_Freq1138 #= 2)) or 
	((Chr1_Chr3_Row1138 #= 382) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 383) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 384) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 385) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 386) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 387) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 388) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 389) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 391) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 393) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Row1138 #= 396) and (Chr1_Chr3_Freq1138 #= 1)) or 
	((Chr1_Chr3_Freq1138 #= 0) and (Chr1_Chr3_Row1138 #= 0)),

	((Chr1_Chr3_Row1139 #= 355) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 356) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 360) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 362) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 363) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 364) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 365) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 366) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 367) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 368) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 369) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 370) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 371) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 372) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 373) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 374) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 380) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 381) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 382) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 383) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 384) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 385) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 386) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 387) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 388) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 389) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 390) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 391) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 392) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 395) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Row1139 #= 396) and (Chr1_Chr3_Freq1139 #= 1)) or 
	((Chr1_Chr3_Freq1139 #= 0) and (Chr1_Chr3_Row1139 #= 0)),

	((Chr1_Chr3_Row1140 #= 359) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 360) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 361) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 362) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 363) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 364) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 365) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 366) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 367) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 368) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 369) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 370) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 371) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 372) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 373) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 374) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 380) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 381) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 382) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 383) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 384) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 385) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 387) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 388) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 389) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 391) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Row1140 #= 394) and (Chr1_Chr3_Freq1140 #= 1)) or 
	((Chr1_Chr3_Freq1140 #= 0) and (Chr1_Chr3_Row1140 #= 0)),

	((Chr1_Chr3_Row1141 #= 358) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 359) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 360) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 362) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 363) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 364) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 365) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 366) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 367) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 368) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 369) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 370) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 371) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 372) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 373) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 374) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 380) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 381) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 382) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 383) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 384) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 385) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 386) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 387) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 389) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 390) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 391) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 392) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 393) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Row1141 #= 396) and (Chr1_Chr3_Freq1141 #= 1)) or 
	((Chr1_Chr3_Freq1141 #= 0) and (Chr1_Chr3_Row1141 #= 0)),

	((Chr1_Chr3_Row1142 #= 345) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 351) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 354) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 356) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 360) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 362) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 364) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 365) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 367) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 368) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 369) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 370) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 371) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 372) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 373) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 374) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 380) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 381) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 382) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 383) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 385) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 386) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 387) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 388) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 389) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 392) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 393) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Row1142 #= 401) and (Chr1_Chr3_Freq1142 #= 1)) or 
	((Chr1_Chr3_Freq1142 #= 0) and (Chr1_Chr3_Row1142 #= 0)),

	((Chr1_Chr3_Row1143 #= 365) and (Chr1_Chr3_Freq1143 #= 1)) or 
	((Chr1_Chr3_Row1143 #= 366) and (Chr1_Chr3_Freq1143 #= 1)) or 
	((Chr1_Chr3_Row1143 #= 367) and (Chr1_Chr3_Freq1143 #= 1)) or 
	((Chr1_Chr3_Row1143 #= 368) and (Chr1_Chr3_Freq1143 #= 1)) or 
	((Chr1_Chr3_Row1143 #= 369) and (Chr1_Chr3_Freq1143 #= 1)) or 
	((Chr1_Chr3_Row1143 #= 370) and (Chr1_Chr3_Freq1143 #= 1)) or 
	((Chr1_Chr3_Row1143 #= 371) and (Chr1_Chr3_Freq1143 #= 1)) or 
	((Chr1_Chr3_Row1143 #= 372) and (Chr1_Chr3_Freq1143 #= 1)) or 
	((Chr1_Chr3_Row1143 #= 373) and (Chr1_Chr3_Freq1143 #= 1)) or 
	((Chr1_Chr3_Row1143 #= 374) and (Chr1_Chr3_Freq1143 #= 1)) or 
	((Chr1_Chr3_Row1143 #= 380) and (Chr1_Chr3_Freq1143 #= 1)) or 
	((Chr1_Chr3_Row1143 #= 381) and (Chr1_Chr3_Freq1143 #= 1)) or 
	((Chr1_Chr3_Row1143 #= 382) and (Chr1_Chr3_Freq1143 #= 1)) or 
	((Chr1_Chr3_Row1143 #= 383) and (Chr1_Chr3_Freq1143 #= 1)) or 
	((Chr1_Chr3_Row1143 #= 385) and (Chr1_Chr3_Freq1143 #= 1)) or 
	((Chr1_Chr3_Row1143 #= 386) and (Chr1_Chr3_Freq1143 #= 1)) or 
	((Chr1_Chr3_Row1143 #= 387) and (Chr1_Chr3_Freq1143 #= 1)) or 
	((Chr1_Chr3_Row1143 #= 390) and (Chr1_Chr3_Freq1143 #= 1)) or 
	((Chr1_Chr3_Row1143 #= 393) and (Chr1_Chr3_Freq1143 #= 1)) or 
	((Chr1_Chr3_Row1143 #= 394) and (Chr1_Chr3_Freq1143 #= 1)) or 
	((Chr1_Chr3_Row1143 #= 396) and (Chr1_Chr3_Freq1143 #= 1)) or 
	((Chr1_Chr3_Row1143 #= 401) and (Chr1_Chr3_Freq1143 #= 1)) or 
	((Chr1_Chr3_Freq1143 #= 0) and (Chr1_Chr3_Row1143 #= 0)),

	((Chr1_Chr3_Row1144 #= 358) and (Chr1_Chr3_Freq1144 #= 1)) or 
	((Chr1_Chr3_Row1144 #= 364) and (Chr1_Chr3_Freq1144 #= 1)) or 
	((Chr1_Chr3_Row1144 #= 368) and (Chr1_Chr3_Freq1144 #= 1)) or 
	((Chr1_Chr3_Row1144 #= 370) and (Chr1_Chr3_Freq1144 #= 1)) or 
	((Chr1_Chr3_Row1144 #= 372) and (Chr1_Chr3_Freq1144 #= 1)) or 
	((Chr1_Chr3_Row1144 #= 373) and (Chr1_Chr3_Freq1144 #= 1)) or 
	((Chr1_Chr3_Row1144 #= 374) and (Chr1_Chr3_Freq1144 #= 1)) or 
	((Chr1_Chr3_Row1144 #= 380) and (Chr1_Chr3_Freq1144 #= 1)) or 
	((Chr1_Chr3_Row1144 #= 381) and (Chr1_Chr3_Freq1144 #= 1)) or 
	((Chr1_Chr3_Row1144 #= 383) and (Chr1_Chr3_Freq1144 #= 1)) or 
	((Chr1_Chr3_Row1144 #= 385) and (Chr1_Chr3_Freq1144 #= 1)) or 
	((Chr1_Chr3_Freq1144 #= 0) and (Chr1_Chr3_Row1144 #= 0)),

	((Chr1_Chr3_Row1145 #= 365) and (Chr1_Chr3_Freq1145 #= 1)) or 
	((Chr1_Chr3_Row1145 #= 367) and (Chr1_Chr3_Freq1145 #= 1)) or 
	((Chr1_Chr3_Row1145 #= 368) and (Chr1_Chr3_Freq1145 #= 1)) or 
	((Chr1_Chr3_Row1145 #= 371) and (Chr1_Chr3_Freq1145 #= 1)) or 
	((Chr1_Chr3_Row1145 #= 372) and (Chr1_Chr3_Freq1145 #= 1)) or 
	((Chr1_Chr3_Row1145 #= 373) and (Chr1_Chr3_Freq1145 #= 1)) or 
	((Chr1_Chr3_Row1145 #= 374) and (Chr1_Chr3_Freq1145 #= 1)) or 
	((Chr1_Chr3_Row1145 #= 380) and (Chr1_Chr3_Freq1145 #= 1)) or 
	((Chr1_Chr3_Freq1145 #= 0) and (Chr1_Chr3_Row1145 #= 0)),

	((Chr1_Chr3_Row1146 #= 367) and (Chr1_Chr3_Freq1146 #= 1)) or 
	((Chr1_Chr3_Row1146 #= 368) and (Chr1_Chr3_Freq1146 #= 1)) or 
	((Chr1_Chr3_Row1146 #= 372) and (Chr1_Chr3_Freq1146 #= 1)) or 
	((Chr1_Chr3_Row1146 #= 373) and (Chr1_Chr3_Freq1146 #= 1)) or 
	((Chr1_Chr3_Row1146 #= 374) and (Chr1_Chr3_Freq1146 #= 1)) or 
	((Chr1_Chr3_Row1146 #= 381) and (Chr1_Chr3_Freq1146 #= 1)) or 
	((Chr1_Chr3_Row1146 #= 382) and (Chr1_Chr3_Freq1146 #= 1)) or 
	((Chr1_Chr3_Freq1146 #= 0) and (Chr1_Chr3_Row1146 #= 0)),

	((Chr1_Chr3_Row1147 #= 369) and (Chr1_Chr3_Freq1147 #= 1)) or 
	((Chr1_Chr3_Row1147 #= 370) and (Chr1_Chr3_Freq1147 #= 1)) or 
	((Chr1_Chr3_Row1147 #= 372) and (Chr1_Chr3_Freq1147 #= 1)) or 
	((Chr1_Chr3_Row1147 #= 386) and (Chr1_Chr3_Freq1147 #= 1)) or 
	((Chr1_Chr3_Freq1147 #= 0) and (Chr1_Chr3_Row1147 #= 0)),

	((Chr1_Chr3_Row1148 #= 365) and (Chr1_Chr3_Freq1148 #= 1)) or 
	((Chr1_Chr3_Row1148 #= 369) and (Chr1_Chr3_Freq1148 #= 1)) or 
	((Chr1_Chr3_Row1148 #= 372) and (Chr1_Chr3_Freq1148 #= 1)) or 
	((Chr1_Chr3_Row1148 #= 373) and (Chr1_Chr3_Freq1148 #= 1)) or 
	((Chr1_Chr3_Row1148 #= 374) and (Chr1_Chr3_Freq1148 #= 1)) or 
	((Chr1_Chr3_Row1148 #= 380) and (Chr1_Chr3_Freq1148 #= 1)) or 
	((Chr1_Chr3_Row1148 #= 388) and (Chr1_Chr3_Freq1148 #= 1)) or 
	((Chr1_Chr3_Freq1148 #= 0) and (Chr1_Chr3_Row1148 #= 0)),

	((Chr1_Chr3_Row1149 #= 364) and (Chr1_Chr3_Freq1149 #= 1)) or 
	((Chr1_Chr3_Row1149 #= 365) and (Chr1_Chr3_Freq1149 #= 1)) or 
	((Chr1_Chr3_Row1149 #= 373) and (Chr1_Chr3_Freq1149 #= 1)) or 
	((Chr1_Chr3_Row1149 #= 382) and (Chr1_Chr3_Freq1149 #= 1)) or 
	((Chr1_Chr3_Freq1149 #= 0) and (Chr1_Chr3_Row1149 #= 0)),

	((Chr1_Chr3_Row1150 #= 365) and (Chr1_Chr3_Freq1150 #= 1)) or 
	((Chr1_Chr3_Row1150 #= 366) and (Chr1_Chr3_Freq1150 #= 1)) or 
	((Chr1_Chr3_Row1150 #= 367) and (Chr1_Chr3_Freq1150 #= 1)) or 
	((Chr1_Chr3_Row1150 #= 369) and (Chr1_Chr3_Freq1150 #= 1)) or 
	((Chr1_Chr3_Row1150 #= 382) and (Chr1_Chr3_Freq1150 #= 1)) or 
	((Chr1_Chr3_Freq1150 #= 0) and (Chr1_Chr3_Row1150 #= 0)),

	((Chr1_Chr3_Row1152 #= 367) and (Chr1_Chr3_Freq1152 #= 1)) or 
	((Chr1_Chr3_Freq1152 #= 0) and (Chr1_Chr3_Row1152 #= 0)),

	((Chr1_Chr3_Row1172 #= 338) and (Chr1_Chr3_Freq1172 #= 1)) or 
	((Chr1_Chr3_Freq1172 #= 0) and (Chr1_Chr3_Row1172 #= 0)),

	((Chr1_Chr3_Row1256 #= 557) and (Chr1_Chr3_Freq1256 #= 1)) or 
	((Chr1_Chr3_Freq1256 #= 0) and (Chr1_Chr3_Row1256 #= 0)),

	((Chr1_Chr3_Row1257 #= 11) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 12) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 13) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 176) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 20) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 204) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 21) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 211) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 22) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 225) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 23) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 24) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 244) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 25) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 34) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 39) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 47) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 479) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 48) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 49) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 51) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 53) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 537) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 539) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 542) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 543) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 546) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 547) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 548) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 549) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 553) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 554) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 555) and (Chr1_Chr3_Freq1257 #= 2)) or 
	((Chr1_Chr3_Row1257 #= 556) and (Chr1_Chr3_Freq1257 #= 4)) or 
	((Chr1_Chr3_Row1257 #= 557) and (Chr1_Chr3_Freq1257 #= 4)) or 
	((Chr1_Chr3_Row1257 #= 58) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 59) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 6) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 65) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 69) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 7) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 71) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 73) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 8) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 88) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Row1257 #= 94) and (Chr1_Chr3_Freq1257 #= 1)) or 
	((Chr1_Chr3_Freq1257 #= 0) and (Chr1_Chr3_Row1257 #= 0)),

	((Chr1_Chr3_Row1258 #= 10) and (Chr1_Chr3_Freq1258 #= 1)) or 
	((Chr1_Chr3_Row1258 #= 12) and (Chr1_Chr3_Freq1258 #= 1)) or 
	((Chr1_Chr3_Row1258 #= 20) and (Chr1_Chr3_Freq1258 #= 1)) or 
	((Chr1_Chr3_Row1258 #= 215) and (Chr1_Chr3_Freq1258 #= 1)) or 
	((Chr1_Chr3_Row1258 #= 225) and (Chr1_Chr3_Freq1258 #= 1)) or 
	((Chr1_Chr3_Row1258 #= 30) and (Chr1_Chr3_Freq1258 #= 1)) or 
	((Chr1_Chr3_Row1258 #= 33) and (Chr1_Chr3_Freq1258 #= 1)) or 
	((Chr1_Chr3_Row1258 #= 36) and (Chr1_Chr3_Freq1258 #= 1)) or 
	((Chr1_Chr3_Row1258 #= 537) and (Chr1_Chr3_Freq1258 #= 1)) or 
	((Chr1_Chr3_Row1258 #= 541) and (Chr1_Chr3_Freq1258 #= 1)) or 
	((Chr1_Chr3_Row1258 #= 542) and (Chr1_Chr3_Freq1258 #= 1)) or 
	((Chr1_Chr3_Row1258 #= 543) and (Chr1_Chr3_Freq1258 #= 1)) or 
	((Chr1_Chr3_Row1258 #= 546) and (Chr1_Chr3_Freq1258 #= 1)) or 
	((Chr1_Chr3_Row1258 #= 550) and (Chr1_Chr3_Freq1258 #= 1)) or 
	((Chr1_Chr3_Row1258 #= 553) and (Chr1_Chr3_Freq1258 #= 1)) or 
	((Chr1_Chr3_Row1258 #= 554) and (Chr1_Chr3_Freq1258 #= 2)) or 
	((Chr1_Chr3_Row1258 #= 555) and (Chr1_Chr3_Freq1258 #= 3)) or 
	((Chr1_Chr3_Row1258 #= 556) and (Chr1_Chr3_Freq1258 #= 4)) or 
	((Chr1_Chr3_Row1258 #= 557) and (Chr1_Chr3_Freq1258 #= 6)) or 
	((Chr1_Chr3_Row1258 #= 56) and (Chr1_Chr3_Freq1258 #= 1)) or 
	((Chr1_Chr3_Row1258 #= 58) and (Chr1_Chr3_Freq1258 #= 1)) or 
	((Chr1_Chr3_Row1258 #= 8) and (Chr1_Chr3_Freq1258 #= 1)) or 
	((Chr1_Chr3_Freq1258 #= 0) and (Chr1_Chr3_Row1258 #= 0)),

	% All of the values assumed by the Row<i> variables must be 
	% all different or zero to ensure each genomic bin is 
	% involved in only 1 interaction; multiple zeros are allowed 
	alldifferent(Non_Zero_Rows),
	atmost(52, Non_Zero_Rows, 0),

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
