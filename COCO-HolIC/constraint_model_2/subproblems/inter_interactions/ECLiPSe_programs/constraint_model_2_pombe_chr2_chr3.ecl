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
	Non_Zero_Rows = [Chr2_Chr3_Row1051, Chr2_Chr3_Row1061, Chr2_Chr3_Row1081, Chr2_Chr3_Row1084, Chr2_Chr3_Row1085, Chr2_Chr3_Row1089, Chr2_Chr3_Row1091, Chr2_Chr3_Row1092, Chr2_Chr3_Row1093, Chr2_Chr3_Row1094, Chr2_Chr3_Row1095, Chr2_Chr3_Row1096, Chr2_Chr3_Row1097, Chr2_Chr3_Row1098, Chr2_Chr3_Row1099, Chr2_Chr3_Row1100, Chr2_Chr3_Row1101, Chr2_Chr3_Row1102, Chr2_Chr3_Row1103, Chr2_Chr3_Row1104, Chr2_Chr3_Row1105, Chr2_Chr3_Row1106, Chr2_Chr3_Row1107, Chr2_Chr3_Row1108, Chr2_Chr3_Row1109, Chr2_Chr3_Row1110, Chr2_Chr3_Row1111, Chr2_Chr3_Row1112, Chr2_Chr3_Row1113, Chr2_Chr3_Row1114, Chr2_Chr3_Row1115, Chr2_Chr3_Row1116, Chr2_Chr3_Row1117, Chr2_Chr3_Row1118, Chr2_Chr3_Row1127, Chr2_Chr3_Row1128, Chr2_Chr3_Row1129, Chr2_Chr3_Row1130, Chr2_Chr3_Row1131, Chr2_Chr3_Row1132, Chr2_Chr3_Row1133, Chr2_Chr3_Row1134, Chr2_Chr3_Row1135, Chr2_Chr3_Row1136, Chr2_Chr3_Row1137, Chr2_Chr3_Row1138, Chr2_Chr3_Row1139, Chr2_Chr3_Row1140, Chr2_Chr3_Row1141, Chr2_Chr3_Row1142, Chr2_Chr3_Row1143, Chr2_Chr3_Row1144, Chr2_Chr3_Row1145, Chr2_Chr3_Row1146, Chr2_Chr3_Row1147, Chr2_Chr3_Row1148, Chr2_Chr3_Row1149, Chr2_Chr3_Row1150, Chr2_Chr3_Row1151, Chr2_Chr3_Row1152, Chr2_Chr3_Row1153, Chr2_Chr3_Row1154, Chr2_Chr3_Row1155, Chr2_Chr3_Row1156, Chr2_Chr3_Row1157, Chr2_Chr3_Row1158, Chr2_Chr3_Row1159, Chr2_Chr3_Row1160, Chr2_Chr3_Row1161, Chr2_Chr3_Row1162, Chr2_Chr3_Row1163, Chr2_Chr3_Row1164, Chr2_Chr3_Row1165, Chr2_Chr3_Row1166, Chr2_Chr3_Row1167, Chr2_Chr3_Row1168, Chr2_Chr3_Row1171, Chr2_Chr3_Row1257, Chr2_Chr3_Row1258],

	% The list Freqs has one variable for each non-zero row of the 
	% whole-genome contact map
	Freqs = [Chr2_Chr3_Freq1051, Chr2_Chr3_Freq1061, Chr2_Chr3_Freq1081, Chr2_Chr3_Freq1084, Chr2_Chr3_Freq1085, Chr2_Chr3_Freq1089, Chr2_Chr3_Freq1091, Chr2_Chr3_Freq1092, Chr2_Chr3_Freq1093, Chr2_Chr3_Freq1094, Chr2_Chr3_Freq1095, Chr2_Chr3_Freq1096, Chr2_Chr3_Freq1097, Chr2_Chr3_Freq1098, Chr2_Chr3_Freq1099, Chr2_Chr3_Freq1100, Chr2_Chr3_Freq1101, Chr2_Chr3_Freq1102, Chr2_Chr3_Freq1103, Chr2_Chr3_Freq1104, Chr2_Chr3_Freq1105, Chr2_Chr3_Freq1106, Chr2_Chr3_Freq1107, Chr2_Chr3_Freq1108, Chr2_Chr3_Freq1109, Chr2_Chr3_Freq1110, Chr2_Chr3_Freq1111, Chr2_Chr3_Freq1112, Chr2_Chr3_Freq1113, Chr2_Chr3_Freq1114, Chr2_Chr3_Freq1115, Chr2_Chr3_Freq1116, Chr2_Chr3_Freq1117, Chr2_Chr3_Freq1118, Chr2_Chr3_Freq1127, Chr2_Chr3_Freq1128, Chr2_Chr3_Freq1129, Chr2_Chr3_Freq1130, Chr2_Chr3_Freq1131, Chr2_Chr3_Freq1132, Chr2_Chr3_Freq1133, Chr2_Chr3_Freq1134, Chr2_Chr3_Freq1135, Chr2_Chr3_Freq1136, Chr2_Chr3_Freq1137, Chr2_Chr3_Freq1138, Chr2_Chr3_Freq1139, Chr2_Chr3_Freq1140, Chr2_Chr3_Freq1141, Chr2_Chr3_Freq1142, Chr2_Chr3_Freq1143, Chr2_Chr3_Freq1144, Chr2_Chr3_Freq1145, Chr2_Chr3_Freq1146, Chr2_Chr3_Freq1147, Chr2_Chr3_Freq1148, Chr2_Chr3_Freq1149, Chr2_Chr3_Freq1150, Chr2_Chr3_Freq1151, Chr2_Chr3_Freq1152, Chr2_Chr3_Freq1153, Chr2_Chr3_Freq1154, Chr2_Chr3_Freq1155, Chr2_Chr3_Freq1156, Chr2_Chr3_Freq1157, Chr2_Chr3_Freq1158, Chr2_Chr3_Freq1159, Chr2_Chr3_Freq1160, Chr2_Chr3_Freq1161, Chr2_Chr3_Freq1162, Chr2_Chr3_Freq1163, Chr2_Chr3_Freq1164, Chr2_Chr3_Freq1165, Chr2_Chr3_Freq1166, Chr2_Chr3_Freq1167, Chr2_Chr3_Freq1168, Chr2_Chr3_Freq1171, Chr2_Chr3_Freq1257, Chr2_Chr3_Freq1258],
	
	% Representation of the Genome: 
	% Each Row term can assume a value based on interacting bin 
	% indices where `0' represents an interaction not being 
	% selected and a non-zero value (ranging from 1 to N) 
	% represents which genomic bin is involved in the selected 
	% interaction 
	Chr2_Chr3_Row1051 :: [0, 690, 706],
	Chr2_Chr3_Row1061 :: [0, 703],
	Chr2_Chr3_Row1081 :: [0, 706],
	Chr2_Chr3_Row1084 :: [0, 704],
	Chr2_Chr3_Row1085 :: [0, 705, 706],
	Chr2_Chr3_Row1089 :: [0, 703, 716, 753],
	Chr2_Chr3_Row1091 :: [0, 703, 707],
	Chr2_Chr3_Row1092 :: [0, 704, 705, 706, 709, 711, 712, 713],
	Chr2_Chr3_Row1093 :: [0, 702, 703, 707, 708, 710, 715, 716, 725],
	Chr2_Chr3_Row1094 :: [0, 616, 703, 704, 705, 706, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 734],
	Chr2_Chr3_Row1095 :: [0, 706, 707, 709, 712, 713, 714, 715, 716, 717, 733],
	Chr2_Chr3_Row1096 :: [0, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 715, 716, 717, 727],
	Chr2_Chr3_Row1097 :: [0, 703, 704, 708, 709, 710, 712, 714, 715, 716, 717, 725, 741],
	Chr2_Chr3_Row1098 :: [0, 699, 703, 704, 705, 706, 707, 708, 710, 711, 712, 714, 715, 716, 717, 725],
	Chr2_Chr3_Row1099 :: [0, 698, 699, 701, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 725, 726, 727, 739],
	Chr2_Chr3_Row1100 :: [0, 694, 699, 700, 702, 704, 708, 709, 710, 712, 714, 715, 716, 725, 726, 732, 733, 741],
	Chr2_Chr3_Row1101 :: [0, 699, 700, 702, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 732, 733],
	Chr2_Chr3_Row1102 :: [0, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 728, 734, 735],
	Chr2_Chr3_Row1103 :: [0, 697, 699, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 734, 737, 738, 739],
	Chr2_Chr3_Row1104 :: [0, 694, 698, 699, 700, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 731, 733, 734, 735, 736, 737, 738, 739],
	Chr2_Chr3_Row1105 :: [0, 694, 698, 699, 700, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 730, 731, 733, 734, 735, 736, 737, 738, 739, 741, 743],
	Chr2_Chr3_Row1106 :: [0, 694, 698, 699, 701, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 740],
	Chr2_Chr3_Row1107 :: [0, 699, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737],
	Chr2_Chr3_Row1108 :: [0, 696, 698, 699, 700, 701, 702, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 740, 744, 752],
	Chr2_Chr3_Row1109 :: [0, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 740, 742, 744, 746],
	Chr2_Chr3_Row1110 :: [0, 696, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 745],
	Chr2_Chr3_Row1111 :: [0, 695, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 744, 746],
	Chr2_Chr3_Row1112 :: [0, 695, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741],
	Chr2_Chr3_Row1113 :: [0, 692, 696, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741],
	Chr2_Chr3_Row1114 :: [0, 694, 695, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 740],
	Chr2_Chr3_Row1115 :: [0, 699, 701, 702, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 740],
	Chr2_Chr3_Row1116 :: [0, 696, 697, 698, 699, 701, 702, 703, 704, 705, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 740, 742, 745, 750],
	Chr2_Chr3_Row1117 :: [0, 694, 695, 696, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 745],
	Chr2_Chr3_Row1118 :: [0, 696, 699, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 740],
	Chr2_Chr3_Row1127 :: [0, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 737, 738],
	Chr2_Chr3_Row1128 :: [0, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738],
	Chr2_Chr3_Row1129 :: [0, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740],
	Chr2_Chr3_Row1130 :: [0, 699, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736],
	Chr2_Chr3_Row1131 :: [0, 701, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736],
	Chr2_Chr3_Row1132 :: [0, 700, 701, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 738, 739],
	Chr2_Chr3_Row1133 :: [0, 695, 697, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 738, 739],
	Chr2_Chr3_Row1134 :: [0, 699, 701, 702, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 732, 734, 735, 737, 740],
	Chr2_Chr3_Row1135 :: [0, 698, 699, 701, 702, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 734, 735, 736],
	Chr2_Chr3_Row1136 :: [0, 698, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 737],
	Chr2_Chr3_Row1137 :: [0, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739],
	Chr2_Chr3_Row1138 :: [0, 695, 697, 698, 699, 700, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 733, 735, 737, 738],
	Chr2_Chr3_Row1139 :: [0, 695, 696, 698, 699, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 734],
	Chr2_Chr3_Row1140 :: [0, 697, 698, 699, 700, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 732, 735],
	Chr2_Chr3_Row1141 :: [0, 692, 693, 694, 695, 696, 697, 698, 699, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 729, 730, 731, 732, 734, 735, 736],
	Chr2_Chr3_Row1142 :: [0, 690, 691, 697, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 730, 731, 732, 733, 734, 738, 739],
	Chr2_Chr3_Row1143 :: [0, 693, 694, 696, 698, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 730, 731, 732, 734, 735],
	Chr2_Chr3_Row1144 :: [0, 694, 695, 697, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 728, 733, 734, 736],
	Chr2_Chr3_Row1145 :: [0, 690, 695, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 725, 726, 728],
	Chr2_Chr3_Row1146 :: [0, 691, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724],
	Chr2_Chr3_Row1147 :: [0, 691, 692, 693, 694, 695, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 727],
	Chr2_Chr3_Row1148 :: [0, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 725, 726, 727, 729, 734],
	Chr2_Chr3_Row1149 :: [0, 686, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 726, 727, 740],
	Chr2_Chr3_Row1150 :: [0, 683, 684, 687, 688, 689, 690, 691, 692, 693, 694, 695, 697, 698, 699, 700, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 724, 725, 728],
	Chr2_Chr3_Row1151 :: [0, 692, 696, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 712, 713, 714, 716, 717],
	Chr2_Chr3_Row1152 :: [0, 685, 688, 689, 690, 691, 692, 693, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 727, 730, 734, 736],
	Chr2_Chr3_Row1153 :: [0, 684, 685, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717],
	Chr2_Chr3_Row1154 :: [0, 684, 685, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 715, 716, 717],
	Chr2_Chr3_Row1155 :: [0, 683, 684, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 712, 713, 715],
	Chr2_Chr3_Row1156 :: [0, 682, 683, 684, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 714, 716],
	Chr2_Chr3_Row1157 :: [0, 678, 683, 684, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711],
	Chr2_Chr3_Row1158 :: [0, 683, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 713, 714, 715],
	Chr2_Chr3_Row1159 :: [0, 684, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711],
	Chr2_Chr3_Row1160 :: [0, 683, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 702, 704, 705, 706, 707, 708],
	Chr2_Chr3_Row1161 :: [0, 688, 689, 690, 691, 692, 693, 694, 699, 700, 702, 703, 705, 706, 707, 708],
	Chr2_Chr3_Row1162 :: [0, 685, 686, 689, 690, 691, 692, 693, 694, 696, 697, 699, 700, 701, 702, 704, 705, 706, 707],
	Chr2_Chr3_Row1163 :: [0, 690, 691, 694, 704, 707, 709],
	Chr2_Chr3_Row1164 :: [0, 677, 689, 690, 699, 700, 702, 704, 707, 709, 755],
	Chr2_Chr3_Row1165 :: [0, 685, 699, 704, 705, 706, 707],
	Chr2_Chr3_Row1166 :: [0, 695, 705, 706],
	Chr2_Chr3_Row1167 :: [0, 690, 704],
	Chr2_Chr3_Row1168 :: [0, 691],
	Chr2_Chr3_Row1171 :: [0, 703],
	Chr2_Chr3_Row1257 :: [0, 562, 563, 568, 571, 582, 584, 591, 592, 593, 594, 596, 599, 601, 603, 604, 605, 665, 681, 690, 691, 694, 695, 704, 823, 877, 898, 908, 917, 922, 923, 924, 925, 926, 927, 951, 964, 967, 968, 969, 971, 974, 975, 976, 978, 979, 980, 981, 982, 983, 984, 985, 987, 988, 990, 992, 994, 995, 996, 997, 998, 999, 1000, 1002, 1003, 1004, 1005, 1006, 1007, 1008],
	Chr2_Chr3_Row1258 :: [0, 564, 566, 580, 585, 591, 594, 682, 924, 968, 970, 972, 979, 980, 981, 982, 983, 984, 985, 986, 987, 992, 994, 998, 999, 1000, 1003, 1004, 1005, 1006, 1007],

	% Each frequency term can assume either the rounded 
	% and scaled integer value (based on the corresponding 
	% interaction frequency from the whole-genome contact 
	% map) or a value of `0' where `0'  represents an 
	% interaction not being selected 
	Chr2_Chr3_Freq1051 :: [0, 1],
	Chr2_Chr3_Freq1061 :: [0, 1],
	Chr2_Chr3_Freq1081 :: [0, 1],
	Chr2_Chr3_Freq1084 :: [0, 1],
	Chr2_Chr3_Freq1085 :: [0, 1],
	Chr2_Chr3_Freq1089 :: [0, 1],
	Chr2_Chr3_Freq1091 :: [0, 1],
	Chr2_Chr3_Freq1092 :: [0, 1],
	Chr2_Chr3_Freq1093 :: [0, 1],
	Chr2_Chr3_Freq1094 :: [0, 1],
	Chr2_Chr3_Freq1095 :: [0, 1],
	Chr2_Chr3_Freq1096 :: [0, 1],
	Chr2_Chr3_Freq1097 :: [0, 1],
	Chr2_Chr3_Freq1098 :: [0, 1],
	Chr2_Chr3_Freq1099 :: [0, 1],
	Chr2_Chr3_Freq1100 :: [0, 1],
	Chr2_Chr3_Freq1101 :: [0, 1],
	Chr2_Chr3_Freq1102 :: [0, 1],
	Chr2_Chr3_Freq1103 :: [0, 1],
	Chr2_Chr3_Freq1104 :: [0, 1],
	Chr2_Chr3_Freq1105 :: [0, 1, 2],
	Chr2_Chr3_Freq1106 :: [0, 1, 2],
	Chr2_Chr3_Freq1107 :: [0, 1, 2],
	Chr2_Chr3_Freq1108 :: [0, 1, 2],
	Chr2_Chr3_Freq1109 :: [0, 1, 2, 3],
	Chr2_Chr3_Freq1110 :: [0, 1, 2],
	Chr2_Chr3_Freq1111 :: [0, 1, 2, 3],
	Chr2_Chr3_Freq1112 :: [0, 1, 2, 3],
	Chr2_Chr3_Freq1113 :: [0, 1, 2, 3],
	Chr2_Chr3_Freq1114 :: [0, 1, 2, 3],
	Chr2_Chr3_Freq1115 :: [0, 1, 2, 3, 4],
	Chr2_Chr3_Freq1116 :: [0, 1, 2, 4, 5, 3],
	Chr2_Chr3_Freq1117 :: [0, 1, 2, 3, 4, 5, 6],
	Chr2_Chr3_Freq1118 :: [0, 1, 2, 3, 4, 5, 6],
	Chr2_Chr3_Freq1127 :: [0, 1, 2, 3, 4, 6],
	Chr2_Chr3_Freq1128 :: [0, 1, 2, 3, 4, 5],
	Chr2_Chr3_Freq1129 :: [0, 1, 2, 3, 4, 5],
	Chr2_Chr3_Freq1130 :: [0, 1, 2, 3, 4],
	Chr2_Chr3_Freq1131 :: [0, 1, 2, 3],
	Chr2_Chr3_Freq1132 :: [0, 1, 2, 3],
	Chr2_Chr3_Freq1133 :: [0, 1, 2, 3],
	Chr2_Chr3_Freq1134 :: [0, 1, 2, 3],
	Chr2_Chr3_Freq1135 :: [0, 1, 2],
	Chr2_Chr3_Freq1136 :: [0, 1, 2],
	Chr2_Chr3_Freq1137 :: [0, 1, 2],
	Chr2_Chr3_Freq1138 :: [0, 1, 2],
	Chr2_Chr3_Freq1139 :: [0, 1, 2],
	Chr2_Chr3_Freq1140 :: [0, 1, 2],
	Chr2_Chr3_Freq1141 :: [0, 1, 2],
	Chr2_Chr3_Freq1142 :: [0, 1],
	Chr2_Chr3_Freq1143 :: [0, 1, 2],
	Chr2_Chr3_Freq1144 :: [0, 1],
	Chr2_Chr3_Freq1145 :: [0, 1],
	Chr2_Chr3_Freq1146 :: [0, 1],
	Chr2_Chr3_Freq1147 :: [0, 1],
	Chr2_Chr3_Freq1148 :: [0, 1],
	Chr2_Chr3_Freq1149 :: [0, 1],
	Chr2_Chr3_Freq1150 :: [0, 1],
	Chr2_Chr3_Freq1151 :: [0, 1],
	Chr2_Chr3_Freq1152 :: [0, 1],
	Chr2_Chr3_Freq1153 :: [0, 1],
	Chr2_Chr3_Freq1154 :: [0, 1],
	Chr2_Chr3_Freq1155 :: [0, 1],
	Chr2_Chr3_Freq1156 :: [0, 1, 2],
	Chr2_Chr3_Freq1157 :: [0, 1, 2],
	Chr2_Chr3_Freq1158 :: [0, 1, 2],
	Chr2_Chr3_Freq1159 :: [0, 1, 2],
	Chr2_Chr3_Freq1160 :: [0, 1],
	Chr2_Chr3_Freq1161 :: [0, 1],
	Chr2_Chr3_Freq1162 :: [0, 1],
	Chr2_Chr3_Freq1163 :: [0, 1],
	Chr2_Chr3_Freq1164 :: [0, 1],
	Chr2_Chr3_Freq1165 :: [0, 1],
	Chr2_Chr3_Freq1166 :: [0, 1],
	Chr2_Chr3_Freq1167 :: [0, 1],
	Chr2_Chr3_Freq1168 :: [0, 1],
	Chr2_Chr3_Freq1171 :: [0, 1],
	Chr2_Chr3_Freq1257 :: [0, 1, 2, 3],
	Chr2_Chr3_Freq1258 :: [0, 1, 2],

	% Constraints: 
	% Each pair of corresponding (Row<i>, Freq<i>) variables 
	% must assume dependent values based on data from the 
	% whole-genome contact map; A (Row, Freq) pair ground to 
	% (0,0) encodes that nothing is chosen 
	((Chr2_Chr3_Row1051 #= 690) and (Chr2_Chr3_Freq1051 #= 1)) or 
	((Chr2_Chr3_Row1051 #= 706) and (Chr2_Chr3_Freq1051 #= 1)) or 
	((Chr2_Chr3_Freq1051 #= 0) and (Chr2_Chr3_Row1051 #= 0)),

	((Chr2_Chr3_Row1061 #= 703) and (Chr2_Chr3_Freq1061 #= 1)) or 
	((Chr2_Chr3_Freq1061 #= 0) and (Chr2_Chr3_Row1061 #= 0)),

	((Chr2_Chr3_Row1081 #= 706) and (Chr2_Chr3_Freq1081 #= 1)) or 
	((Chr2_Chr3_Freq1081 #= 0) and (Chr2_Chr3_Row1081 #= 0)),

	((Chr2_Chr3_Row1084 #= 704) and (Chr2_Chr3_Freq1084 #= 1)) or 
	((Chr2_Chr3_Freq1084 #= 0) and (Chr2_Chr3_Row1084 #= 0)),

	((Chr2_Chr3_Row1085 #= 705) and (Chr2_Chr3_Freq1085 #= 1)) or 
	((Chr2_Chr3_Row1085 #= 706) and (Chr2_Chr3_Freq1085 #= 1)) or 
	((Chr2_Chr3_Freq1085 #= 0) and (Chr2_Chr3_Row1085 #= 0)),

	((Chr2_Chr3_Row1089 #= 703) and (Chr2_Chr3_Freq1089 #= 1)) or 
	((Chr2_Chr3_Row1089 #= 716) and (Chr2_Chr3_Freq1089 #= 1)) or 
	((Chr2_Chr3_Row1089 #= 753) and (Chr2_Chr3_Freq1089 #= 1)) or 
	((Chr2_Chr3_Freq1089 #= 0) and (Chr2_Chr3_Row1089 #= 0)),

	((Chr2_Chr3_Row1091 #= 703) and (Chr2_Chr3_Freq1091 #= 1)) or 
	((Chr2_Chr3_Row1091 #= 707) and (Chr2_Chr3_Freq1091 #= 1)) or 
	((Chr2_Chr3_Freq1091 #= 0) and (Chr2_Chr3_Row1091 #= 0)),

	((Chr2_Chr3_Row1092 #= 704) and (Chr2_Chr3_Freq1092 #= 1)) or 
	((Chr2_Chr3_Row1092 #= 705) and (Chr2_Chr3_Freq1092 #= 1)) or 
	((Chr2_Chr3_Row1092 #= 706) and (Chr2_Chr3_Freq1092 #= 1)) or 
	((Chr2_Chr3_Row1092 #= 709) and (Chr2_Chr3_Freq1092 #= 1)) or 
	((Chr2_Chr3_Row1092 #= 711) and (Chr2_Chr3_Freq1092 #= 1)) or 
	((Chr2_Chr3_Row1092 #= 712) and (Chr2_Chr3_Freq1092 #= 1)) or 
	((Chr2_Chr3_Row1092 #= 713) and (Chr2_Chr3_Freq1092 #= 1)) or 
	((Chr2_Chr3_Freq1092 #= 0) and (Chr2_Chr3_Row1092 #= 0)),

	((Chr2_Chr3_Row1093 #= 702) and (Chr2_Chr3_Freq1093 #= 1)) or 
	((Chr2_Chr3_Row1093 #= 703) and (Chr2_Chr3_Freq1093 #= 1)) or 
	((Chr2_Chr3_Row1093 #= 707) and (Chr2_Chr3_Freq1093 #= 1)) or 
	((Chr2_Chr3_Row1093 #= 708) and (Chr2_Chr3_Freq1093 #= 1)) or 
	((Chr2_Chr3_Row1093 #= 710) and (Chr2_Chr3_Freq1093 #= 1)) or 
	((Chr2_Chr3_Row1093 #= 715) and (Chr2_Chr3_Freq1093 #= 1)) or 
	((Chr2_Chr3_Row1093 #= 716) and (Chr2_Chr3_Freq1093 #= 1)) or 
	((Chr2_Chr3_Row1093 #= 725) and (Chr2_Chr3_Freq1093 #= 1)) or 
	((Chr2_Chr3_Freq1093 #= 0) and (Chr2_Chr3_Row1093 #= 0)),

	((Chr2_Chr3_Row1094 #= 616) and (Chr2_Chr3_Freq1094 #= 1)) or 
	((Chr2_Chr3_Row1094 #= 703) and (Chr2_Chr3_Freq1094 #= 1)) or 
	((Chr2_Chr3_Row1094 #= 704) and (Chr2_Chr3_Freq1094 #= 1)) or 
	((Chr2_Chr3_Row1094 #= 705) and (Chr2_Chr3_Freq1094 #= 1)) or 
	((Chr2_Chr3_Row1094 #= 706) and (Chr2_Chr3_Freq1094 #= 1)) or 
	((Chr2_Chr3_Row1094 #= 709) and (Chr2_Chr3_Freq1094 #= 1)) or 
	((Chr2_Chr3_Row1094 #= 710) and (Chr2_Chr3_Freq1094 #= 1)) or 
	((Chr2_Chr3_Row1094 #= 711) and (Chr2_Chr3_Freq1094 #= 1)) or 
	((Chr2_Chr3_Row1094 #= 712) and (Chr2_Chr3_Freq1094 #= 1)) or 
	((Chr2_Chr3_Row1094 #= 713) and (Chr2_Chr3_Freq1094 #= 1)) or 
	((Chr2_Chr3_Row1094 #= 714) and (Chr2_Chr3_Freq1094 #= 1)) or 
	((Chr2_Chr3_Row1094 #= 715) and (Chr2_Chr3_Freq1094 #= 1)) or 
	((Chr2_Chr3_Row1094 #= 716) and (Chr2_Chr3_Freq1094 #= 1)) or 
	((Chr2_Chr3_Row1094 #= 717) and (Chr2_Chr3_Freq1094 #= 1)) or 
	((Chr2_Chr3_Row1094 #= 724) and (Chr2_Chr3_Freq1094 #= 1)) or 
	((Chr2_Chr3_Row1094 #= 725) and (Chr2_Chr3_Freq1094 #= 1)) or 
	((Chr2_Chr3_Row1094 #= 726) and (Chr2_Chr3_Freq1094 #= 1)) or 
	((Chr2_Chr3_Row1094 #= 727) and (Chr2_Chr3_Freq1094 #= 1)) or 
	((Chr2_Chr3_Row1094 #= 734) and (Chr2_Chr3_Freq1094 #= 1)) or 
	((Chr2_Chr3_Freq1094 #= 0) and (Chr2_Chr3_Row1094 #= 0)),

	((Chr2_Chr3_Row1095 #= 706) and (Chr2_Chr3_Freq1095 #= 1)) or 
	((Chr2_Chr3_Row1095 #= 707) and (Chr2_Chr3_Freq1095 #= 1)) or 
	((Chr2_Chr3_Row1095 #= 709) and (Chr2_Chr3_Freq1095 #= 1)) or 
	((Chr2_Chr3_Row1095 #= 712) and (Chr2_Chr3_Freq1095 #= 1)) or 
	((Chr2_Chr3_Row1095 #= 713) and (Chr2_Chr3_Freq1095 #= 1)) or 
	((Chr2_Chr3_Row1095 #= 714) and (Chr2_Chr3_Freq1095 #= 1)) or 
	((Chr2_Chr3_Row1095 #= 715) and (Chr2_Chr3_Freq1095 #= 1)) or 
	((Chr2_Chr3_Row1095 #= 716) and (Chr2_Chr3_Freq1095 #= 1)) or 
	((Chr2_Chr3_Row1095 #= 717) and (Chr2_Chr3_Freq1095 #= 1)) or 
	((Chr2_Chr3_Row1095 #= 733) and (Chr2_Chr3_Freq1095 #= 1)) or 
	((Chr2_Chr3_Freq1095 #= 0) and (Chr2_Chr3_Row1095 #= 0)),

	((Chr2_Chr3_Row1096 #= 702) and (Chr2_Chr3_Freq1096 #= 1)) or 
	((Chr2_Chr3_Row1096 #= 703) and (Chr2_Chr3_Freq1096 #= 1)) or 
	((Chr2_Chr3_Row1096 #= 704) and (Chr2_Chr3_Freq1096 #= 1)) or 
	((Chr2_Chr3_Row1096 #= 705) and (Chr2_Chr3_Freq1096 #= 1)) or 
	((Chr2_Chr3_Row1096 #= 706) and (Chr2_Chr3_Freq1096 #= 1)) or 
	((Chr2_Chr3_Row1096 #= 707) and (Chr2_Chr3_Freq1096 #= 1)) or 
	((Chr2_Chr3_Row1096 #= 708) and (Chr2_Chr3_Freq1096 #= 1)) or 
	((Chr2_Chr3_Row1096 #= 709) and (Chr2_Chr3_Freq1096 #= 1)) or 
	((Chr2_Chr3_Row1096 #= 710) and (Chr2_Chr3_Freq1096 #= 1)) or 
	((Chr2_Chr3_Row1096 #= 711) and (Chr2_Chr3_Freq1096 #= 1)) or 
	((Chr2_Chr3_Row1096 #= 712) and (Chr2_Chr3_Freq1096 #= 1)) or 
	((Chr2_Chr3_Row1096 #= 713) and (Chr2_Chr3_Freq1096 #= 1)) or 
	((Chr2_Chr3_Row1096 #= 715) and (Chr2_Chr3_Freq1096 #= 1)) or 
	((Chr2_Chr3_Row1096 #= 716) and (Chr2_Chr3_Freq1096 #= 1)) or 
	((Chr2_Chr3_Row1096 #= 717) and (Chr2_Chr3_Freq1096 #= 1)) or 
	((Chr2_Chr3_Row1096 #= 727) and (Chr2_Chr3_Freq1096 #= 1)) or 
	((Chr2_Chr3_Freq1096 #= 0) and (Chr2_Chr3_Row1096 #= 0)),

	((Chr2_Chr3_Row1097 #= 703) and (Chr2_Chr3_Freq1097 #= 1)) or 
	((Chr2_Chr3_Row1097 #= 704) and (Chr2_Chr3_Freq1097 #= 1)) or 
	((Chr2_Chr3_Row1097 #= 708) and (Chr2_Chr3_Freq1097 #= 1)) or 
	((Chr2_Chr3_Row1097 #= 709) and (Chr2_Chr3_Freq1097 #= 1)) or 
	((Chr2_Chr3_Row1097 #= 710) and (Chr2_Chr3_Freq1097 #= 1)) or 
	((Chr2_Chr3_Row1097 #= 712) and (Chr2_Chr3_Freq1097 #= 1)) or 
	((Chr2_Chr3_Row1097 #= 714) and (Chr2_Chr3_Freq1097 #= 1)) or 
	((Chr2_Chr3_Row1097 #= 715) and (Chr2_Chr3_Freq1097 #= 1)) or 
	((Chr2_Chr3_Row1097 #= 716) and (Chr2_Chr3_Freq1097 #= 1)) or 
	((Chr2_Chr3_Row1097 #= 717) and (Chr2_Chr3_Freq1097 #= 1)) or 
	((Chr2_Chr3_Row1097 #= 725) and (Chr2_Chr3_Freq1097 #= 1)) or 
	((Chr2_Chr3_Row1097 #= 741) and (Chr2_Chr3_Freq1097 #= 1)) or 
	((Chr2_Chr3_Freq1097 #= 0) and (Chr2_Chr3_Row1097 #= 0)),

	((Chr2_Chr3_Row1098 #= 699) and (Chr2_Chr3_Freq1098 #= 1)) or 
	((Chr2_Chr3_Row1098 #= 703) and (Chr2_Chr3_Freq1098 #= 1)) or 
	((Chr2_Chr3_Row1098 #= 704) and (Chr2_Chr3_Freq1098 #= 1)) or 
	((Chr2_Chr3_Row1098 #= 705) and (Chr2_Chr3_Freq1098 #= 1)) or 
	((Chr2_Chr3_Row1098 #= 706) and (Chr2_Chr3_Freq1098 #= 1)) or 
	((Chr2_Chr3_Row1098 #= 707) and (Chr2_Chr3_Freq1098 #= 1)) or 
	((Chr2_Chr3_Row1098 #= 708) and (Chr2_Chr3_Freq1098 #= 1)) or 
	((Chr2_Chr3_Row1098 #= 710) and (Chr2_Chr3_Freq1098 #= 1)) or 
	((Chr2_Chr3_Row1098 #= 711) and (Chr2_Chr3_Freq1098 #= 1)) or 
	((Chr2_Chr3_Row1098 #= 712) and (Chr2_Chr3_Freq1098 #= 1)) or 
	((Chr2_Chr3_Row1098 #= 714) and (Chr2_Chr3_Freq1098 #= 1)) or 
	((Chr2_Chr3_Row1098 #= 715) and (Chr2_Chr3_Freq1098 #= 1)) or 
	((Chr2_Chr3_Row1098 #= 716) and (Chr2_Chr3_Freq1098 #= 1)) or 
	((Chr2_Chr3_Row1098 #= 717) and (Chr2_Chr3_Freq1098 #= 1)) or 
	((Chr2_Chr3_Row1098 #= 725) and (Chr2_Chr3_Freq1098 #= 1)) or 
	((Chr2_Chr3_Freq1098 #= 0) and (Chr2_Chr3_Row1098 #= 0)),

	((Chr2_Chr3_Row1099 #= 698) and (Chr2_Chr3_Freq1099 #= 1)) or 
	((Chr2_Chr3_Row1099 #= 699) and (Chr2_Chr3_Freq1099 #= 1)) or 
	((Chr2_Chr3_Row1099 #= 701) and (Chr2_Chr3_Freq1099 #= 1)) or 
	((Chr2_Chr3_Row1099 #= 703) and (Chr2_Chr3_Freq1099 #= 1)) or 
	((Chr2_Chr3_Row1099 #= 704) and (Chr2_Chr3_Freq1099 #= 1)) or 
	((Chr2_Chr3_Row1099 #= 705) and (Chr2_Chr3_Freq1099 #= 1)) or 
	((Chr2_Chr3_Row1099 #= 706) and (Chr2_Chr3_Freq1099 #= 1)) or 
	((Chr2_Chr3_Row1099 #= 707) and (Chr2_Chr3_Freq1099 #= 1)) or 
	((Chr2_Chr3_Row1099 #= 708) and (Chr2_Chr3_Freq1099 #= 1)) or 
	((Chr2_Chr3_Row1099 #= 709) and (Chr2_Chr3_Freq1099 #= 1)) or 
	((Chr2_Chr3_Row1099 #= 710) and (Chr2_Chr3_Freq1099 #= 1)) or 
	((Chr2_Chr3_Row1099 #= 711) and (Chr2_Chr3_Freq1099 #= 1)) or 
	((Chr2_Chr3_Row1099 #= 712) and (Chr2_Chr3_Freq1099 #= 1)) or 
	((Chr2_Chr3_Row1099 #= 713) and (Chr2_Chr3_Freq1099 #= 1)) or 
	((Chr2_Chr3_Row1099 #= 714) and (Chr2_Chr3_Freq1099 #= 1)) or 
	((Chr2_Chr3_Row1099 #= 715) and (Chr2_Chr3_Freq1099 #= 1)) or 
	((Chr2_Chr3_Row1099 #= 716) and (Chr2_Chr3_Freq1099 #= 1)) or 
	((Chr2_Chr3_Row1099 #= 717) and (Chr2_Chr3_Freq1099 #= 1)) or 
	((Chr2_Chr3_Row1099 #= 725) and (Chr2_Chr3_Freq1099 #= 1)) or 
	((Chr2_Chr3_Row1099 #= 726) and (Chr2_Chr3_Freq1099 #= 1)) or 
	((Chr2_Chr3_Row1099 #= 727) and (Chr2_Chr3_Freq1099 #= 1)) or 
	((Chr2_Chr3_Row1099 #= 739) and (Chr2_Chr3_Freq1099 #= 1)) or 
	((Chr2_Chr3_Freq1099 #= 0) and (Chr2_Chr3_Row1099 #= 0)),

	((Chr2_Chr3_Row1100 #= 694) and (Chr2_Chr3_Freq1100 #= 1)) or 
	((Chr2_Chr3_Row1100 #= 699) and (Chr2_Chr3_Freq1100 #= 1)) or 
	((Chr2_Chr3_Row1100 #= 700) and (Chr2_Chr3_Freq1100 #= 1)) or 
	((Chr2_Chr3_Row1100 #= 702) and (Chr2_Chr3_Freq1100 #= 1)) or 
	((Chr2_Chr3_Row1100 #= 704) and (Chr2_Chr3_Freq1100 #= 1)) or 
	((Chr2_Chr3_Row1100 #= 708) and (Chr2_Chr3_Freq1100 #= 1)) or 
	((Chr2_Chr3_Row1100 #= 709) and (Chr2_Chr3_Freq1100 #= 1)) or 
	((Chr2_Chr3_Row1100 #= 710) and (Chr2_Chr3_Freq1100 #= 1)) or 
	((Chr2_Chr3_Row1100 #= 712) and (Chr2_Chr3_Freq1100 #= 1)) or 
	((Chr2_Chr3_Row1100 #= 714) and (Chr2_Chr3_Freq1100 #= 1)) or 
	((Chr2_Chr3_Row1100 #= 715) and (Chr2_Chr3_Freq1100 #= 1)) or 
	((Chr2_Chr3_Row1100 #= 716) and (Chr2_Chr3_Freq1100 #= 1)) or 
	((Chr2_Chr3_Row1100 #= 725) and (Chr2_Chr3_Freq1100 #= 1)) or 
	((Chr2_Chr3_Row1100 #= 726) and (Chr2_Chr3_Freq1100 #= 1)) or 
	((Chr2_Chr3_Row1100 #= 732) and (Chr2_Chr3_Freq1100 #= 1)) or 
	((Chr2_Chr3_Row1100 #= 733) and (Chr2_Chr3_Freq1100 #= 1)) or 
	((Chr2_Chr3_Row1100 #= 741) and (Chr2_Chr3_Freq1100 #= 1)) or 
	((Chr2_Chr3_Freq1100 #= 0) and (Chr2_Chr3_Row1100 #= 0)),

	((Chr2_Chr3_Row1101 #= 699) and (Chr2_Chr3_Freq1101 #= 1)) or 
	((Chr2_Chr3_Row1101 #= 700) and (Chr2_Chr3_Freq1101 #= 1)) or 
	((Chr2_Chr3_Row1101 #= 702) and (Chr2_Chr3_Freq1101 #= 1)) or 
	((Chr2_Chr3_Row1101 #= 704) and (Chr2_Chr3_Freq1101 #= 1)) or 
	((Chr2_Chr3_Row1101 #= 705) and (Chr2_Chr3_Freq1101 #= 1)) or 
	((Chr2_Chr3_Row1101 #= 706) and (Chr2_Chr3_Freq1101 #= 1)) or 
	((Chr2_Chr3_Row1101 #= 707) and (Chr2_Chr3_Freq1101 #= 1)) or 
	((Chr2_Chr3_Row1101 #= 708) and (Chr2_Chr3_Freq1101 #= 1)) or 
	((Chr2_Chr3_Row1101 #= 709) and (Chr2_Chr3_Freq1101 #= 1)) or 
	((Chr2_Chr3_Row1101 #= 710) and (Chr2_Chr3_Freq1101 #= 1)) or 
	((Chr2_Chr3_Row1101 #= 711) and (Chr2_Chr3_Freq1101 #= 1)) or 
	((Chr2_Chr3_Row1101 #= 712) and (Chr2_Chr3_Freq1101 #= 1)) or 
	((Chr2_Chr3_Row1101 #= 713) and (Chr2_Chr3_Freq1101 #= 1)) or 
	((Chr2_Chr3_Row1101 #= 714) and (Chr2_Chr3_Freq1101 #= 1)) or 
	((Chr2_Chr3_Row1101 #= 715) and (Chr2_Chr3_Freq1101 #= 1)) or 
	((Chr2_Chr3_Row1101 #= 716) and (Chr2_Chr3_Freq1101 #= 1)) or 
	((Chr2_Chr3_Row1101 #= 717) and (Chr2_Chr3_Freq1101 #= 1)) or 
	((Chr2_Chr3_Row1101 #= 724) and (Chr2_Chr3_Freq1101 #= 1)) or 
	((Chr2_Chr3_Row1101 #= 725) and (Chr2_Chr3_Freq1101 #= 1)) or 
	((Chr2_Chr3_Row1101 #= 726) and (Chr2_Chr3_Freq1101 #= 1)) or 
	((Chr2_Chr3_Row1101 #= 732) and (Chr2_Chr3_Freq1101 #= 1)) or 
	((Chr2_Chr3_Row1101 #= 733) and (Chr2_Chr3_Freq1101 #= 1)) or 
	((Chr2_Chr3_Freq1101 #= 0) and (Chr2_Chr3_Row1101 #= 0)),

	((Chr2_Chr3_Row1102 #= 699) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Row1102 #= 700) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Row1102 #= 701) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Row1102 #= 702) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Row1102 #= 703) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Row1102 #= 704) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Row1102 #= 705) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Row1102 #= 706) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Row1102 #= 707) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Row1102 #= 708) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Row1102 #= 709) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Row1102 #= 710) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Row1102 #= 711) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Row1102 #= 712) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Row1102 #= 713) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Row1102 #= 714) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Row1102 #= 715) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Row1102 #= 716) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Row1102 #= 717) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Row1102 #= 724) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Row1102 #= 725) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Row1102 #= 726) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Row1102 #= 728) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Row1102 #= 734) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Row1102 #= 735) and (Chr2_Chr3_Freq1102 #= 1)) or 
	((Chr2_Chr3_Freq1102 #= 0) and (Chr2_Chr3_Row1102 #= 0)),

	((Chr2_Chr3_Row1103 #= 697) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Row1103 #= 699) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Row1103 #= 702) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Row1103 #= 703) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Row1103 #= 704) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Row1103 #= 705) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Row1103 #= 706) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Row1103 #= 707) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Row1103 #= 708) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Row1103 #= 709) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Row1103 #= 710) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Row1103 #= 711) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Row1103 #= 712) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Row1103 #= 713) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Row1103 #= 714) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Row1103 #= 715) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Row1103 #= 716) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Row1103 #= 717) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Row1103 #= 724) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Row1103 #= 725) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Row1103 #= 726) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Row1103 #= 734) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Row1103 #= 737) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Row1103 #= 738) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Row1103 #= 739) and (Chr2_Chr3_Freq1103 #= 1)) or 
	((Chr2_Chr3_Freq1103 #= 0) and (Chr2_Chr3_Row1103 #= 0)),

	((Chr2_Chr3_Row1104 #= 694) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 698) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 699) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 700) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 702) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 703) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 704) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 705) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 706) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 707) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 708) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 709) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 710) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 711) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 712) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 713) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 714) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 715) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 716) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 717) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 724) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 725) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 726) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 727) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 728) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 729) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 731) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 733) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 734) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 735) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 736) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 737) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 738) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Row1104 #= 739) and (Chr2_Chr3_Freq1104 #= 1)) or 
	((Chr2_Chr3_Freq1104 #= 0) and (Chr2_Chr3_Row1104 #= 0)),

	((Chr2_Chr3_Row1105 #= 694) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 698) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 699) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 700) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 702) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 703) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 704) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 705) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 706) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 707) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 708) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 709) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 710) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 711) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 712) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 713) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 714) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 715) and (Chr2_Chr3_Freq1105 #= 2)) or 
	((Chr2_Chr3_Row1105 #= 716) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 717) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 724) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 725) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 726) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 727) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 728) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 730) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 731) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 733) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 734) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 735) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 736) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 737) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 738) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 739) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 741) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Row1105 #= 743) and (Chr2_Chr3_Freq1105 #= 1)) or 
	((Chr2_Chr3_Freq1105 #= 0) and (Chr2_Chr3_Row1105 #= 0)),

	((Chr2_Chr3_Row1106 #= 694) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 698) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 699) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 701) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 703) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 704) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 705) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 706) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 707) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 708) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 709) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 710) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 711) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 712) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 713) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 714) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 715) and (Chr2_Chr3_Freq1106 #= 2)) or 
	((Chr2_Chr3_Row1106 #= 716) and (Chr2_Chr3_Freq1106 #= 2)) or 
	((Chr2_Chr3_Row1106 #= 717) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 724) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 725) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 726) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 727) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 728) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 729) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 730) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 731) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 732) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 733) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 734) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 735) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 736) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Row1106 #= 740) and (Chr2_Chr3_Freq1106 #= 1)) or 
	((Chr2_Chr3_Freq1106 #= 0) and (Chr2_Chr3_Row1106 #= 0)),

	((Chr2_Chr3_Row1107 #= 699) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 701) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 702) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 703) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 704) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 705) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 706) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 707) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 708) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 709) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 710) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 711) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 712) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 713) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 714) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 715) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 716) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 717) and (Chr2_Chr3_Freq1107 #= 2)) or 
	((Chr2_Chr3_Row1107 #= 724) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 725) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 726) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 727) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 728) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 729) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 730) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 731) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 732) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 733) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 734) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 735) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 736) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Row1107 #= 737) and (Chr2_Chr3_Freq1107 #= 1)) or 
	((Chr2_Chr3_Freq1107 #= 0) and (Chr2_Chr3_Row1107 #= 0)),

	((Chr2_Chr3_Row1108 #= 696) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 698) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 699) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 700) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 701) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 702) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 704) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 705) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 706) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 707) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 708) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 709) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 710) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 711) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 712) and (Chr2_Chr3_Freq1108 #= 2)) or 
	((Chr2_Chr3_Row1108 #= 713) and (Chr2_Chr3_Freq1108 #= 2)) or 
	((Chr2_Chr3_Row1108 #= 714) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 715) and (Chr2_Chr3_Freq1108 #= 2)) or 
	((Chr2_Chr3_Row1108 #= 716) and (Chr2_Chr3_Freq1108 #= 2)) or 
	((Chr2_Chr3_Row1108 #= 717) and (Chr2_Chr3_Freq1108 #= 2)) or 
	((Chr2_Chr3_Row1108 #= 724) and (Chr2_Chr3_Freq1108 #= 2)) or 
	((Chr2_Chr3_Row1108 #= 725) and (Chr2_Chr3_Freq1108 #= 2)) or 
	((Chr2_Chr3_Row1108 #= 726) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 727) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 728) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 729) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 730) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 731) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 732) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 733) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 734) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 735) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 736) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 740) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 744) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Row1108 #= 752) and (Chr2_Chr3_Freq1108 #= 1)) or 
	((Chr2_Chr3_Freq1108 #= 0) and (Chr2_Chr3_Row1108 #= 0)),

	((Chr2_Chr3_Row1109 #= 695) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 696) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 697) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 698) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 699) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 700) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 701) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 702) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 703) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 704) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 705) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 706) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 707) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 708) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 709) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 710) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 711) and (Chr2_Chr3_Freq1109 #= 2)) or 
	((Chr2_Chr3_Row1109 #= 712) and (Chr2_Chr3_Freq1109 #= 2)) or 
	((Chr2_Chr3_Row1109 #= 713) and (Chr2_Chr3_Freq1109 #= 2)) or 
	((Chr2_Chr3_Row1109 #= 714) and (Chr2_Chr3_Freq1109 #= 2)) or 
	((Chr2_Chr3_Row1109 #= 715) and (Chr2_Chr3_Freq1109 #= 2)) or 
	((Chr2_Chr3_Row1109 #= 716) and (Chr2_Chr3_Freq1109 #= 2)) or 
	((Chr2_Chr3_Row1109 #= 717) and (Chr2_Chr3_Freq1109 #= 3)) or 
	((Chr2_Chr3_Row1109 #= 724) and (Chr2_Chr3_Freq1109 #= 2)) or 
	((Chr2_Chr3_Row1109 #= 725) and (Chr2_Chr3_Freq1109 #= 2)) or 
	((Chr2_Chr3_Row1109 #= 726) and (Chr2_Chr3_Freq1109 #= 2)) or 
	((Chr2_Chr3_Row1109 #= 727) and (Chr2_Chr3_Freq1109 #= 2)) or 
	((Chr2_Chr3_Row1109 #= 728) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 729) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 730) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 731) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 732) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 733) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 734) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 735) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 736) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 740) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 742) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 744) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Row1109 #= 746) and (Chr2_Chr3_Freq1109 #= 1)) or 
	((Chr2_Chr3_Freq1109 #= 0) and (Chr2_Chr3_Row1109 #= 0)),

	((Chr2_Chr3_Row1110 #= 696) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 698) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 699) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 700) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 701) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 702) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 703) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 704) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 705) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 706) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 707) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 708) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 709) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 710) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 711) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 712) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 713) and (Chr2_Chr3_Freq1110 #= 2)) or 
	((Chr2_Chr3_Row1110 #= 714) and (Chr2_Chr3_Freq1110 #= 2)) or 
	((Chr2_Chr3_Row1110 #= 715) and (Chr2_Chr3_Freq1110 #= 2)) or 
	((Chr2_Chr3_Row1110 #= 716) and (Chr2_Chr3_Freq1110 #= 2)) or 
	((Chr2_Chr3_Row1110 #= 717) and (Chr2_Chr3_Freq1110 #= 2)) or 
	((Chr2_Chr3_Row1110 #= 724) and (Chr2_Chr3_Freq1110 #= 2)) or 
	((Chr2_Chr3_Row1110 #= 725) and (Chr2_Chr3_Freq1110 #= 2)) or 
	((Chr2_Chr3_Row1110 #= 726) and (Chr2_Chr3_Freq1110 #= 2)) or 
	((Chr2_Chr3_Row1110 #= 727) and (Chr2_Chr3_Freq1110 #= 2)) or 
	((Chr2_Chr3_Row1110 #= 728) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 729) and (Chr2_Chr3_Freq1110 #= 2)) or 
	((Chr2_Chr3_Row1110 #= 730) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 731) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 732) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 733) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 734) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 735) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 736) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 737) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 738) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 739) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 740) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 741) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Row1110 #= 745) and (Chr2_Chr3_Freq1110 #= 1)) or 
	((Chr2_Chr3_Freq1110 #= 0) and (Chr2_Chr3_Row1110 #= 0)),

	((Chr2_Chr3_Row1111 #= 695) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 697) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 698) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 699) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 700) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 701) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 702) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 703) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 704) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 705) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 706) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 707) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 708) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 709) and (Chr2_Chr3_Freq1111 #= 2)) or 
	((Chr2_Chr3_Row1111 #= 710) and (Chr2_Chr3_Freq1111 #= 2)) or 
	((Chr2_Chr3_Row1111 #= 711) and (Chr2_Chr3_Freq1111 #= 2)) or 
	((Chr2_Chr3_Row1111 #= 712) and (Chr2_Chr3_Freq1111 #= 2)) or 
	((Chr2_Chr3_Row1111 #= 713) and (Chr2_Chr3_Freq1111 #= 2)) or 
	((Chr2_Chr3_Row1111 #= 714) and (Chr2_Chr3_Freq1111 #= 2)) or 
	((Chr2_Chr3_Row1111 #= 715) and (Chr2_Chr3_Freq1111 #= 2)) or 
	((Chr2_Chr3_Row1111 #= 716) and (Chr2_Chr3_Freq1111 #= 3)) or 
	((Chr2_Chr3_Row1111 #= 717) and (Chr2_Chr3_Freq1111 #= 3)) or 
	((Chr2_Chr3_Row1111 #= 724) and (Chr2_Chr3_Freq1111 #= 2)) or 
	((Chr2_Chr3_Row1111 #= 725) and (Chr2_Chr3_Freq1111 #= 2)) or 
	((Chr2_Chr3_Row1111 #= 726) and (Chr2_Chr3_Freq1111 #= 2)) or 
	((Chr2_Chr3_Row1111 #= 727) and (Chr2_Chr3_Freq1111 #= 2)) or 
	((Chr2_Chr3_Row1111 #= 728) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 729) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 730) and (Chr2_Chr3_Freq1111 #= 2)) or 
	((Chr2_Chr3_Row1111 #= 731) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 732) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 733) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 734) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 735) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 736) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 737) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 738) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 739) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 740) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 744) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Row1111 #= 746) and (Chr2_Chr3_Freq1111 #= 1)) or 
	((Chr2_Chr3_Freq1111 #= 0) and (Chr2_Chr3_Row1111 #= 0)),

	((Chr2_Chr3_Row1112 #= 695) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 697) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 698) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 699) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 700) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 701) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 702) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 703) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 704) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 705) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 706) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 707) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 708) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 709) and (Chr2_Chr3_Freq1112 #= 2)) or 
	((Chr2_Chr3_Row1112 #= 710) and (Chr2_Chr3_Freq1112 #= 2)) or 
	((Chr2_Chr3_Row1112 #= 711) and (Chr2_Chr3_Freq1112 #= 2)) or 
	((Chr2_Chr3_Row1112 #= 712) and (Chr2_Chr3_Freq1112 #= 2)) or 
	((Chr2_Chr3_Row1112 #= 713) and (Chr2_Chr3_Freq1112 #= 2)) or 
	((Chr2_Chr3_Row1112 #= 714) and (Chr2_Chr3_Freq1112 #= 2)) or 
	((Chr2_Chr3_Row1112 #= 715) and (Chr2_Chr3_Freq1112 #= 3)) or 
	((Chr2_Chr3_Row1112 #= 716) and (Chr2_Chr3_Freq1112 #= 2)) or 
	((Chr2_Chr3_Row1112 #= 717) and (Chr2_Chr3_Freq1112 #= 3)) or 
	((Chr2_Chr3_Row1112 #= 724) and (Chr2_Chr3_Freq1112 #= 3)) or 
	((Chr2_Chr3_Row1112 #= 725) and (Chr2_Chr3_Freq1112 #= 3)) or 
	((Chr2_Chr3_Row1112 #= 726) and (Chr2_Chr3_Freq1112 #= 2)) or 
	((Chr2_Chr3_Row1112 #= 727) and (Chr2_Chr3_Freq1112 #= 2)) or 
	((Chr2_Chr3_Row1112 #= 728) and (Chr2_Chr3_Freq1112 #= 2)) or 
	((Chr2_Chr3_Row1112 #= 729) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 730) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 731) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 732) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 733) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 734) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 735) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 736) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 737) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 738) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 739) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 740) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Row1112 #= 741) and (Chr2_Chr3_Freq1112 #= 1)) or 
	((Chr2_Chr3_Freq1112 #= 0) and (Chr2_Chr3_Row1112 #= 0)),

	((Chr2_Chr3_Row1113 #= 692) and (Chr2_Chr3_Freq1113 #= 1)) or 
	((Chr2_Chr3_Row1113 #= 696) and (Chr2_Chr3_Freq1113 #= 1)) or 
	((Chr2_Chr3_Row1113 #= 698) and (Chr2_Chr3_Freq1113 #= 1)) or 
	((Chr2_Chr3_Row1113 #= 699) and (Chr2_Chr3_Freq1113 #= 1)) or 
	((Chr2_Chr3_Row1113 #= 700) and (Chr2_Chr3_Freq1113 #= 1)) or 
	((Chr2_Chr3_Row1113 #= 701) and (Chr2_Chr3_Freq1113 #= 1)) or 
	((Chr2_Chr3_Row1113 #= 702) and (Chr2_Chr3_Freq1113 #= 1)) or 
	((Chr2_Chr3_Row1113 #= 703) and (Chr2_Chr3_Freq1113 #= 1)) or 
	((Chr2_Chr3_Row1113 #= 704) and (Chr2_Chr3_Freq1113 #= 1)) or 
	((Chr2_Chr3_Row1113 #= 705) and (Chr2_Chr3_Freq1113 #= 1)) or 
	((Chr2_Chr3_Row1113 #= 706) and (Chr2_Chr3_Freq1113 #= 1)) or 
	((Chr2_Chr3_Row1113 #= 707) and (Chr2_Chr3_Freq1113 #= 2)) or 
	((Chr2_Chr3_Row1113 #= 708) and (Chr2_Chr3_Freq1113 #= 1)) or 
	((Chr2_Chr3_Row1113 #= 709) and (Chr2_Chr3_Freq1113 #= 2)) or 
	((Chr2_Chr3_Row1113 #= 710) and (Chr2_Chr3_Freq1113 #= 2)) or 
	((Chr2_Chr3_Row1113 #= 711) and (Chr2_Chr3_Freq1113 #= 2)) or 
	((Chr2_Chr3_Row1113 #= 712) and (Chr2_Chr3_Freq1113 #= 2)) or 
	((Chr2_Chr3_Row1113 #= 713) and (Chr2_Chr3_Freq1113 #= 2)) or 
	((Chr2_Chr3_Row1113 #= 714) and (Chr2_Chr3_Freq1113 #= 2)) or 
	((Chr2_Chr3_Row1113 #= 715) and (Chr2_Chr3_Freq1113 #= 3)) or 
	((Chr2_Chr3_Row1113 #= 716) and (Chr2_Chr3_Freq1113 #= 3)) or 
	((Chr2_Chr3_Row1113 #= 717) and (Chr2_Chr3_Freq1113 #= 3)) or 
	((Chr2_Chr3_Row1113 #= 724) and (Chr2_Chr3_Freq1113 #= 3)) or 
	((Chr2_Chr3_Row1113 #= 725) and (Chr2_Chr3_Freq1113 #= 3)) or 
	((Chr2_Chr3_Row1113 #= 726) and (Chr2_Chr3_Freq1113 #= 2)) or 
	((Chr2_Chr3_Row1113 #= 727) and (Chr2_Chr3_Freq1113 #= 2)) or 
	((Chr2_Chr3_Row1113 #= 728) and (Chr2_Chr3_Freq1113 #= 1)) or 
	((Chr2_Chr3_Row1113 #= 729) and (Chr2_Chr3_Freq1113 #= 2)) or 
	((Chr2_Chr3_Row1113 #= 730) and (Chr2_Chr3_Freq1113 #= 2)) or 
	((Chr2_Chr3_Row1113 #= 731) and (Chr2_Chr3_Freq1113 #= 1)) or 
	((Chr2_Chr3_Row1113 #= 732) and (Chr2_Chr3_Freq1113 #= 1)) or 
	((Chr2_Chr3_Row1113 #= 733) and (Chr2_Chr3_Freq1113 #= 1)) or 
	((Chr2_Chr3_Row1113 #= 734) and (Chr2_Chr3_Freq1113 #= 1)) or 
	((Chr2_Chr3_Row1113 #= 735) and (Chr2_Chr3_Freq1113 #= 1)) or 
	((Chr2_Chr3_Row1113 #= 736) and (Chr2_Chr3_Freq1113 #= 1)) or 
	((Chr2_Chr3_Row1113 #= 737) and (Chr2_Chr3_Freq1113 #= 1)) or 
	((Chr2_Chr3_Row1113 #= 738) and (Chr2_Chr3_Freq1113 #= 1)) or 
	((Chr2_Chr3_Row1113 #= 739) and (Chr2_Chr3_Freq1113 #= 1)) or 
	((Chr2_Chr3_Row1113 #= 740) and (Chr2_Chr3_Freq1113 #= 1)) or 
	((Chr2_Chr3_Row1113 #= 741) and (Chr2_Chr3_Freq1113 #= 1)) or 
	((Chr2_Chr3_Freq1113 #= 0) and (Chr2_Chr3_Row1113 #= 0)),

	((Chr2_Chr3_Row1114 #= 694) and (Chr2_Chr3_Freq1114 #= 1)) or 
	((Chr2_Chr3_Row1114 #= 695) and (Chr2_Chr3_Freq1114 #= 1)) or 
	((Chr2_Chr3_Row1114 #= 698) and (Chr2_Chr3_Freq1114 #= 1)) or 
	((Chr2_Chr3_Row1114 #= 699) and (Chr2_Chr3_Freq1114 #= 1)) or 
	((Chr2_Chr3_Row1114 #= 700) and (Chr2_Chr3_Freq1114 #= 1)) or 
	((Chr2_Chr3_Row1114 #= 701) and (Chr2_Chr3_Freq1114 #= 1)) or 
	((Chr2_Chr3_Row1114 #= 702) and (Chr2_Chr3_Freq1114 #= 1)) or 
	((Chr2_Chr3_Row1114 #= 703) and (Chr2_Chr3_Freq1114 #= 1)) or 
	((Chr2_Chr3_Row1114 #= 704) and (Chr2_Chr3_Freq1114 #= 1)) or 
	((Chr2_Chr3_Row1114 #= 705) and (Chr2_Chr3_Freq1114 #= 1)) or 
	((Chr2_Chr3_Row1114 #= 706) and (Chr2_Chr3_Freq1114 #= 1)) or 
	((Chr2_Chr3_Row1114 #= 707) and (Chr2_Chr3_Freq1114 #= 1)) or 
	((Chr2_Chr3_Row1114 #= 708) and (Chr2_Chr3_Freq1114 #= 1)) or 
	((Chr2_Chr3_Row1114 #= 709) and (Chr2_Chr3_Freq1114 #= 2)) or 
	((Chr2_Chr3_Row1114 #= 710) and (Chr2_Chr3_Freq1114 #= 1)) or 
	((Chr2_Chr3_Row1114 #= 711) and (Chr2_Chr3_Freq1114 #= 2)) or 
	((Chr2_Chr3_Row1114 #= 712) and (Chr2_Chr3_Freq1114 #= 2)) or 
	((Chr2_Chr3_Row1114 #= 713) and (Chr2_Chr3_Freq1114 #= 2)) or 
	((Chr2_Chr3_Row1114 #= 714) and (Chr2_Chr3_Freq1114 #= 3)) or 
	((Chr2_Chr3_Row1114 #= 715) and (Chr2_Chr3_Freq1114 #= 3)) or 
	((Chr2_Chr3_Row1114 #= 716) and (Chr2_Chr3_Freq1114 #= 3)) or 
	((Chr2_Chr3_Row1114 #= 717) and (Chr2_Chr3_Freq1114 #= 3)) or 
	((Chr2_Chr3_Row1114 #= 724) and (Chr2_Chr3_Freq1114 #= 3)) or 
	((Chr2_Chr3_Row1114 #= 725) and (Chr2_Chr3_Freq1114 #= 3)) or 
	((Chr2_Chr3_Row1114 #= 726) and (Chr2_Chr3_Freq1114 #= 2)) or 
	((Chr2_Chr3_Row1114 #= 727) and (Chr2_Chr3_Freq1114 #= 2)) or 
	((Chr2_Chr3_Row1114 #= 728) and (Chr2_Chr3_Freq1114 #= 2)) or 
	((Chr2_Chr3_Row1114 #= 729) and (Chr2_Chr3_Freq1114 #= 2)) or 
	((Chr2_Chr3_Row1114 #= 730) and (Chr2_Chr3_Freq1114 #= 1)) or 
	((Chr2_Chr3_Row1114 #= 731) and (Chr2_Chr3_Freq1114 #= 1)) or 
	((Chr2_Chr3_Row1114 #= 732) and (Chr2_Chr3_Freq1114 #= 1)) or 
	((Chr2_Chr3_Row1114 #= 733) and (Chr2_Chr3_Freq1114 #= 1)) or 
	((Chr2_Chr3_Row1114 #= 734) and (Chr2_Chr3_Freq1114 #= 1)) or 
	((Chr2_Chr3_Row1114 #= 735) and (Chr2_Chr3_Freq1114 #= 1)) or 
	((Chr2_Chr3_Row1114 #= 736) and (Chr2_Chr3_Freq1114 #= 1)) or 
	((Chr2_Chr3_Row1114 #= 737) and (Chr2_Chr3_Freq1114 #= 1)) or 
	((Chr2_Chr3_Row1114 #= 740) and (Chr2_Chr3_Freq1114 #= 1)) or 
	((Chr2_Chr3_Freq1114 #= 0) and (Chr2_Chr3_Row1114 #= 0)),

	((Chr2_Chr3_Row1115 #= 699) and (Chr2_Chr3_Freq1115 #= 1)) or 
	((Chr2_Chr3_Row1115 #= 701) and (Chr2_Chr3_Freq1115 #= 1)) or 
	((Chr2_Chr3_Row1115 #= 702) and (Chr2_Chr3_Freq1115 #= 1)) or 
	((Chr2_Chr3_Row1115 #= 704) and (Chr2_Chr3_Freq1115 #= 1)) or 
	((Chr2_Chr3_Row1115 #= 705) and (Chr2_Chr3_Freq1115 #= 1)) or 
	((Chr2_Chr3_Row1115 #= 706) and (Chr2_Chr3_Freq1115 #= 1)) or 
	((Chr2_Chr3_Row1115 #= 707) and (Chr2_Chr3_Freq1115 #= 1)) or 
	((Chr2_Chr3_Row1115 #= 708) and (Chr2_Chr3_Freq1115 #= 1)) or 
	((Chr2_Chr3_Row1115 #= 709) and (Chr2_Chr3_Freq1115 #= 2)) or 
	((Chr2_Chr3_Row1115 #= 710) and (Chr2_Chr3_Freq1115 #= 2)) or 
	((Chr2_Chr3_Row1115 #= 711) and (Chr2_Chr3_Freq1115 #= 2)) or 
	((Chr2_Chr3_Row1115 #= 712) and (Chr2_Chr3_Freq1115 #= 2)) or 
	((Chr2_Chr3_Row1115 #= 713) and (Chr2_Chr3_Freq1115 #= 3)) or 
	((Chr2_Chr3_Row1115 #= 714) and (Chr2_Chr3_Freq1115 #= 3)) or 
	((Chr2_Chr3_Row1115 #= 715) and (Chr2_Chr3_Freq1115 #= 3)) or 
	((Chr2_Chr3_Row1115 #= 716) and (Chr2_Chr3_Freq1115 #= 4)) or 
	((Chr2_Chr3_Row1115 #= 717) and (Chr2_Chr3_Freq1115 #= 4)) or 
	((Chr2_Chr3_Row1115 #= 724) and (Chr2_Chr3_Freq1115 #= 3)) or 
	((Chr2_Chr3_Row1115 #= 725) and (Chr2_Chr3_Freq1115 #= 4)) or 
	((Chr2_Chr3_Row1115 #= 726) and (Chr2_Chr3_Freq1115 #= 3)) or 
	((Chr2_Chr3_Row1115 #= 727) and (Chr2_Chr3_Freq1115 #= 3)) or 
	((Chr2_Chr3_Row1115 #= 728) and (Chr2_Chr3_Freq1115 #= 2)) or 
	((Chr2_Chr3_Row1115 #= 729) and (Chr2_Chr3_Freq1115 #= 2)) or 
	((Chr2_Chr3_Row1115 #= 730) and (Chr2_Chr3_Freq1115 #= 1)) or 
	((Chr2_Chr3_Row1115 #= 731) and (Chr2_Chr3_Freq1115 #= 2)) or 
	((Chr2_Chr3_Row1115 #= 732) and (Chr2_Chr3_Freq1115 #= 1)) or 
	((Chr2_Chr3_Row1115 #= 733) and (Chr2_Chr3_Freq1115 #= 1)) or 
	((Chr2_Chr3_Row1115 #= 734) and (Chr2_Chr3_Freq1115 #= 1)) or 
	((Chr2_Chr3_Row1115 #= 735) and (Chr2_Chr3_Freq1115 #= 1)) or 
	((Chr2_Chr3_Row1115 #= 736) and (Chr2_Chr3_Freq1115 #= 1)) or 
	((Chr2_Chr3_Row1115 #= 737) and (Chr2_Chr3_Freq1115 #= 1)) or 
	((Chr2_Chr3_Row1115 #= 738) and (Chr2_Chr3_Freq1115 #= 1)) or 
	((Chr2_Chr3_Row1115 #= 740) and (Chr2_Chr3_Freq1115 #= 1)) or 
	((Chr2_Chr3_Freq1115 #= 0) and (Chr2_Chr3_Row1115 #= 0)),

	((Chr2_Chr3_Row1116 #= 696) and (Chr2_Chr3_Freq1116 #= 1)) or 
	((Chr2_Chr3_Row1116 #= 697) and (Chr2_Chr3_Freq1116 #= 1)) or 
	((Chr2_Chr3_Row1116 #= 698) and (Chr2_Chr3_Freq1116 #= 1)) or 
	((Chr2_Chr3_Row1116 #= 699) and (Chr2_Chr3_Freq1116 #= 1)) or 
	((Chr2_Chr3_Row1116 #= 701) and (Chr2_Chr3_Freq1116 #= 1)) or 
	((Chr2_Chr3_Row1116 #= 702) and (Chr2_Chr3_Freq1116 #= 1)) or 
	((Chr2_Chr3_Row1116 #= 703) and (Chr2_Chr3_Freq1116 #= 1)) or 
	((Chr2_Chr3_Row1116 #= 704) and (Chr2_Chr3_Freq1116 #= 1)) or 
	((Chr2_Chr3_Row1116 #= 705) and (Chr2_Chr3_Freq1116 #= 1)) or 
	((Chr2_Chr3_Row1116 #= 707) and (Chr2_Chr3_Freq1116 #= 2)) or 
	((Chr2_Chr3_Row1116 #= 708) and (Chr2_Chr3_Freq1116 #= 2)) or 
	((Chr2_Chr3_Row1116 #= 709) and (Chr2_Chr3_Freq1116 #= 2)) or 
	((Chr2_Chr3_Row1116 #= 710) and (Chr2_Chr3_Freq1116 #= 2)) or 
	((Chr2_Chr3_Row1116 #= 711) and (Chr2_Chr3_Freq1116 #= 2)) or 
	((Chr2_Chr3_Row1116 #= 712) and (Chr2_Chr3_Freq1116 #= 2)) or 
	((Chr2_Chr3_Row1116 #= 713) and (Chr2_Chr3_Freq1116 #= 4)) or 
	((Chr2_Chr3_Row1116 #= 714) and (Chr2_Chr3_Freq1116 #= 4)) or 
	((Chr2_Chr3_Row1116 #= 715) and (Chr2_Chr3_Freq1116 #= 4)) or 
	((Chr2_Chr3_Row1116 #= 716) and (Chr2_Chr3_Freq1116 #= 5)) or 
	((Chr2_Chr3_Row1116 #= 717) and (Chr2_Chr3_Freq1116 #= 5)) or 
	((Chr2_Chr3_Row1116 #= 724) and (Chr2_Chr3_Freq1116 #= 4)) or 
	((Chr2_Chr3_Row1116 #= 725) and (Chr2_Chr3_Freq1116 #= 5)) or 
	((Chr2_Chr3_Row1116 #= 726) and (Chr2_Chr3_Freq1116 #= 4)) or 
	((Chr2_Chr3_Row1116 #= 727) and (Chr2_Chr3_Freq1116 #= 3)) or 
	((Chr2_Chr3_Row1116 #= 728) and (Chr2_Chr3_Freq1116 #= 2)) or 
	((Chr2_Chr3_Row1116 #= 729) and (Chr2_Chr3_Freq1116 #= 2)) or 
	((Chr2_Chr3_Row1116 #= 730) and (Chr2_Chr3_Freq1116 #= 2)) or 
	((Chr2_Chr3_Row1116 #= 731) and (Chr2_Chr3_Freq1116 #= 2)) or 
	((Chr2_Chr3_Row1116 #= 732) and (Chr2_Chr3_Freq1116 #= 1)) or 
	((Chr2_Chr3_Row1116 #= 733) and (Chr2_Chr3_Freq1116 #= 1)) or 
	((Chr2_Chr3_Row1116 #= 734) and (Chr2_Chr3_Freq1116 #= 1)) or 
	((Chr2_Chr3_Row1116 #= 735) and (Chr2_Chr3_Freq1116 #= 2)) or 
	((Chr2_Chr3_Row1116 #= 736) and (Chr2_Chr3_Freq1116 #= 2)) or 
	((Chr2_Chr3_Row1116 #= 737) and (Chr2_Chr3_Freq1116 #= 1)) or 
	((Chr2_Chr3_Row1116 #= 740) and (Chr2_Chr3_Freq1116 #= 1)) or 
	((Chr2_Chr3_Row1116 #= 742) and (Chr2_Chr3_Freq1116 #= 1)) or 
	((Chr2_Chr3_Row1116 #= 745) and (Chr2_Chr3_Freq1116 #= 1)) or 
	((Chr2_Chr3_Row1116 #= 750) and (Chr2_Chr3_Freq1116 #= 1)) or 
	((Chr2_Chr3_Freq1116 #= 0) and (Chr2_Chr3_Row1116 #= 0)),

	((Chr2_Chr3_Row1117 #= 694) and (Chr2_Chr3_Freq1117 #= 1)) or 
	((Chr2_Chr3_Row1117 #= 695) and (Chr2_Chr3_Freq1117 #= 1)) or 
	((Chr2_Chr3_Row1117 #= 696) and (Chr2_Chr3_Freq1117 #= 1)) or 
	((Chr2_Chr3_Row1117 #= 699) and (Chr2_Chr3_Freq1117 #= 1)) or 
	((Chr2_Chr3_Row1117 #= 700) and (Chr2_Chr3_Freq1117 #= 1)) or 
	((Chr2_Chr3_Row1117 #= 701) and (Chr2_Chr3_Freq1117 #= 1)) or 
	((Chr2_Chr3_Row1117 #= 702) and (Chr2_Chr3_Freq1117 #= 1)) or 
	((Chr2_Chr3_Row1117 #= 703) and (Chr2_Chr3_Freq1117 #= 1)) or 
	((Chr2_Chr3_Row1117 #= 704) and (Chr2_Chr3_Freq1117 #= 1)) or 
	((Chr2_Chr3_Row1117 #= 705) and (Chr2_Chr3_Freq1117 #= 1)) or 
	((Chr2_Chr3_Row1117 #= 706) and (Chr2_Chr3_Freq1117 #= 1)) or 
	((Chr2_Chr3_Row1117 #= 707) and (Chr2_Chr3_Freq1117 #= 2)) or 
	((Chr2_Chr3_Row1117 #= 708) and (Chr2_Chr3_Freq1117 #= 1)) or 
	((Chr2_Chr3_Row1117 #= 709) and (Chr2_Chr3_Freq1117 #= 2)) or 
	((Chr2_Chr3_Row1117 #= 710) and (Chr2_Chr3_Freq1117 #= 2)) or 
	((Chr2_Chr3_Row1117 #= 711) and (Chr2_Chr3_Freq1117 #= 2)) or 
	((Chr2_Chr3_Row1117 #= 712) and (Chr2_Chr3_Freq1117 #= 3)) or 
	((Chr2_Chr3_Row1117 #= 713) and (Chr2_Chr3_Freq1117 #= 3)) or 
	((Chr2_Chr3_Row1117 #= 714) and (Chr2_Chr3_Freq1117 #= 4)) or 
	((Chr2_Chr3_Row1117 #= 715) and (Chr2_Chr3_Freq1117 #= 4)) or 
	((Chr2_Chr3_Row1117 #= 716) and (Chr2_Chr3_Freq1117 #= 5)) or 
	((Chr2_Chr3_Row1117 #= 717) and (Chr2_Chr3_Freq1117 #= 5)) or 
	((Chr2_Chr3_Row1117 #= 724) and (Chr2_Chr3_Freq1117 #= 5)) or 
	((Chr2_Chr3_Row1117 #= 725) and (Chr2_Chr3_Freq1117 #= 6)) or 
	((Chr2_Chr3_Row1117 #= 726) and (Chr2_Chr3_Freq1117 #= 4)) or 
	((Chr2_Chr3_Row1117 #= 727) and (Chr2_Chr3_Freq1117 #= 3)) or 
	((Chr2_Chr3_Row1117 #= 728) and (Chr2_Chr3_Freq1117 #= 3)) or 
	((Chr2_Chr3_Row1117 #= 729) and (Chr2_Chr3_Freq1117 #= 3)) or 
	((Chr2_Chr3_Row1117 #= 730) and (Chr2_Chr3_Freq1117 #= 2)) or 
	((Chr2_Chr3_Row1117 #= 731) and (Chr2_Chr3_Freq1117 #= 2)) or 
	((Chr2_Chr3_Row1117 #= 732) and (Chr2_Chr3_Freq1117 #= 1)) or 
	((Chr2_Chr3_Row1117 #= 733) and (Chr2_Chr3_Freq1117 #= 1)) or 
	((Chr2_Chr3_Row1117 #= 734) and (Chr2_Chr3_Freq1117 #= 1)) or 
	((Chr2_Chr3_Row1117 #= 735) and (Chr2_Chr3_Freq1117 #= 1)) or 
	((Chr2_Chr3_Row1117 #= 736) and (Chr2_Chr3_Freq1117 #= 1)) or 
	((Chr2_Chr3_Row1117 #= 737) and (Chr2_Chr3_Freq1117 #= 1)) or 
	((Chr2_Chr3_Row1117 #= 738) and (Chr2_Chr3_Freq1117 #= 1)) or 
	((Chr2_Chr3_Row1117 #= 739) and (Chr2_Chr3_Freq1117 #= 1)) or 
	((Chr2_Chr3_Row1117 #= 740) and (Chr2_Chr3_Freq1117 #= 1)) or 
	((Chr2_Chr3_Row1117 #= 745) and (Chr2_Chr3_Freq1117 #= 1)) or 
	((Chr2_Chr3_Freq1117 #= 0) and (Chr2_Chr3_Row1117 #= 0)),

	((Chr2_Chr3_Row1118 #= 696) and (Chr2_Chr3_Freq1118 #= 1)) or 
	((Chr2_Chr3_Row1118 #= 699) and (Chr2_Chr3_Freq1118 #= 1)) or 
	((Chr2_Chr3_Row1118 #= 701) and (Chr2_Chr3_Freq1118 #= 1)) or 
	((Chr2_Chr3_Row1118 #= 702) and (Chr2_Chr3_Freq1118 #= 1)) or 
	((Chr2_Chr3_Row1118 #= 703) and (Chr2_Chr3_Freq1118 #= 1)) or 
	((Chr2_Chr3_Row1118 #= 704) and (Chr2_Chr3_Freq1118 #= 1)) or 
	((Chr2_Chr3_Row1118 #= 705) and (Chr2_Chr3_Freq1118 #= 1)) or 
	((Chr2_Chr3_Row1118 #= 706) and (Chr2_Chr3_Freq1118 #= 1)) or 
	((Chr2_Chr3_Row1118 #= 707) and (Chr2_Chr3_Freq1118 #= 1)) or 
	((Chr2_Chr3_Row1118 #= 708) and (Chr2_Chr3_Freq1118 #= 1)) or 
	((Chr2_Chr3_Row1118 #= 709) and (Chr2_Chr3_Freq1118 #= 1)) or 
	((Chr2_Chr3_Row1118 #= 710) and (Chr2_Chr3_Freq1118 #= 2)) or 
	((Chr2_Chr3_Row1118 #= 711) and (Chr2_Chr3_Freq1118 #= 2)) or 
	((Chr2_Chr3_Row1118 #= 712) and (Chr2_Chr3_Freq1118 #= 2)) or 
	((Chr2_Chr3_Row1118 #= 713) and (Chr2_Chr3_Freq1118 #= 3)) or 
	((Chr2_Chr3_Row1118 #= 714) and (Chr2_Chr3_Freq1118 #= 3)) or 
	((Chr2_Chr3_Row1118 #= 715) and (Chr2_Chr3_Freq1118 #= 4)) or 
	((Chr2_Chr3_Row1118 #= 716) and (Chr2_Chr3_Freq1118 #= 5)) or 
	((Chr2_Chr3_Row1118 #= 717) and (Chr2_Chr3_Freq1118 #= 6)) or 
	((Chr2_Chr3_Row1118 #= 724) and (Chr2_Chr3_Freq1118 #= 5)) or 
	((Chr2_Chr3_Row1118 #= 725) and (Chr2_Chr3_Freq1118 #= 5)) or 
	((Chr2_Chr3_Row1118 #= 726) and (Chr2_Chr3_Freq1118 #= 5)) or 
	((Chr2_Chr3_Row1118 #= 727) and (Chr2_Chr3_Freq1118 #= 4)) or 
	((Chr2_Chr3_Row1118 #= 728) and (Chr2_Chr3_Freq1118 #= 3)) or 
	((Chr2_Chr3_Row1118 #= 729) and (Chr2_Chr3_Freq1118 #= 3)) or 
	((Chr2_Chr3_Row1118 #= 730) and (Chr2_Chr3_Freq1118 #= 2)) or 
	((Chr2_Chr3_Row1118 #= 731) and (Chr2_Chr3_Freq1118 #= 2)) or 
	((Chr2_Chr3_Row1118 #= 732) and (Chr2_Chr3_Freq1118 #= 1)) or 
	((Chr2_Chr3_Row1118 #= 733) and (Chr2_Chr3_Freq1118 #= 1)) or 
	((Chr2_Chr3_Row1118 #= 734) and (Chr2_Chr3_Freq1118 #= 1)) or 
	((Chr2_Chr3_Row1118 #= 735) and (Chr2_Chr3_Freq1118 #= 1)) or 
	((Chr2_Chr3_Row1118 #= 736) and (Chr2_Chr3_Freq1118 #= 1)) or 
	((Chr2_Chr3_Row1118 #= 737) and (Chr2_Chr3_Freq1118 #= 1)) or 
	((Chr2_Chr3_Row1118 #= 738) and (Chr2_Chr3_Freq1118 #= 1)) or 
	((Chr2_Chr3_Row1118 #= 740) and (Chr2_Chr3_Freq1118 #= 1)) or 
	((Chr2_Chr3_Freq1118 #= 0) and (Chr2_Chr3_Row1118 #= 0)),
	 
	((Chr2_Chr3_Row1127 #= 701) and (Chr2_Chr3_Freq1127 #= 1)) or 
	((Chr2_Chr3_Row1127 #= 702) and (Chr2_Chr3_Freq1127 #= 1)) or 
	((Chr2_Chr3_Row1127 #= 703) and (Chr2_Chr3_Freq1127 #= 1)) or 
	((Chr2_Chr3_Row1127 #= 704) and (Chr2_Chr3_Freq1127 #= 1)) or 
	((Chr2_Chr3_Row1127 #= 705) and (Chr2_Chr3_Freq1127 #= 1)) or 
	((Chr2_Chr3_Row1127 #= 706) and (Chr2_Chr3_Freq1127 #= 1)) or 
	((Chr2_Chr3_Row1127 #= 707) and (Chr2_Chr3_Freq1127 #= 1)) or 
	((Chr2_Chr3_Row1127 #= 708) and (Chr2_Chr3_Freq1127 #= 1)) or 
	((Chr2_Chr3_Row1127 #= 709) and (Chr2_Chr3_Freq1127 #= 1)) or 
	((Chr2_Chr3_Row1127 #= 710) and (Chr2_Chr3_Freq1127 #= 2)) or 
	((Chr2_Chr3_Row1127 #= 711) and (Chr2_Chr3_Freq1127 #= 2)) or 
	((Chr2_Chr3_Row1127 #= 712) and (Chr2_Chr3_Freq1127 #= 2)) or 
	((Chr2_Chr3_Row1127 #= 713) and (Chr2_Chr3_Freq1127 #= 2)) or 
	((Chr2_Chr3_Row1127 #= 714) and (Chr2_Chr3_Freq1127 #= 3)) or 
	((Chr2_Chr3_Row1127 #= 715) and (Chr2_Chr3_Freq1127 #= 4)) or 
	((Chr2_Chr3_Row1127 #= 716) and (Chr2_Chr3_Freq1127 #= 4)) or 
	((Chr2_Chr3_Row1127 #= 717) and (Chr2_Chr3_Freq1127 #= 4)) or 
	((Chr2_Chr3_Row1127 #= 724) and (Chr2_Chr3_Freq1127 #= 6)) or 
	((Chr2_Chr3_Row1127 #= 725) and (Chr2_Chr3_Freq1127 #= 6)) or 
	((Chr2_Chr3_Row1127 #= 726) and (Chr2_Chr3_Freq1127 #= 4)) or 
	((Chr2_Chr3_Row1127 #= 727) and (Chr2_Chr3_Freq1127 #= 3)) or 
	((Chr2_Chr3_Row1127 #= 728) and (Chr2_Chr3_Freq1127 #= 3)) or 
	((Chr2_Chr3_Row1127 #= 729) and (Chr2_Chr3_Freq1127 #= 2)) or 
	((Chr2_Chr3_Row1127 #= 730) and (Chr2_Chr3_Freq1127 #= 2)) or 
	((Chr2_Chr3_Row1127 #= 731) and (Chr2_Chr3_Freq1127 #= 1)) or 
	((Chr2_Chr3_Row1127 #= 732) and (Chr2_Chr3_Freq1127 #= 1)) or 
	((Chr2_Chr3_Row1127 #= 733) and (Chr2_Chr3_Freq1127 #= 1)) or 
	((Chr2_Chr3_Row1127 #= 734) and (Chr2_Chr3_Freq1127 #= 1)) or 
	((Chr2_Chr3_Row1127 #= 735) and (Chr2_Chr3_Freq1127 #= 1)) or 
	((Chr2_Chr3_Row1127 #= 737) and (Chr2_Chr3_Freq1127 #= 1)) or 
	((Chr2_Chr3_Row1127 #= 738) and (Chr2_Chr3_Freq1127 #= 1)) or 
	((Chr2_Chr3_Freq1127 #= 0) and (Chr2_Chr3_Row1127 #= 0)),

	((Chr2_Chr3_Row1128 #= 700) and (Chr2_Chr3_Freq1128 #= 1)) or 
	((Chr2_Chr3_Row1128 #= 701) and (Chr2_Chr3_Freq1128 #= 1)) or 
	((Chr2_Chr3_Row1128 #= 702) and (Chr2_Chr3_Freq1128 #= 1)) or 
	((Chr2_Chr3_Row1128 #= 703) and (Chr2_Chr3_Freq1128 #= 1)) or 
	((Chr2_Chr3_Row1128 #= 704) and (Chr2_Chr3_Freq1128 #= 1)) or 
	((Chr2_Chr3_Row1128 #= 705) and (Chr2_Chr3_Freq1128 #= 1)) or 
	((Chr2_Chr3_Row1128 #= 706) and (Chr2_Chr3_Freq1128 #= 1)) or 
	((Chr2_Chr3_Row1128 #= 707) and (Chr2_Chr3_Freq1128 #= 1)) or 
	((Chr2_Chr3_Row1128 #= 708) and (Chr2_Chr3_Freq1128 #= 2)) or 
	((Chr2_Chr3_Row1128 #= 709) and (Chr2_Chr3_Freq1128 #= 2)) or 
	((Chr2_Chr3_Row1128 #= 710) and (Chr2_Chr3_Freq1128 #= 2)) or 
	((Chr2_Chr3_Row1128 #= 711) and (Chr2_Chr3_Freq1128 #= 2)) or 
	((Chr2_Chr3_Row1128 #= 712) and (Chr2_Chr3_Freq1128 #= 2)) or 
	((Chr2_Chr3_Row1128 #= 713) and (Chr2_Chr3_Freq1128 #= 3)) or 
	((Chr2_Chr3_Row1128 #= 714) and (Chr2_Chr3_Freq1128 #= 3)) or 
	((Chr2_Chr3_Row1128 #= 715) and (Chr2_Chr3_Freq1128 #= 4)) or 
	((Chr2_Chr3_Row1128 #= 716) and (Chr2_Chr3_Freq1128 #= 4)) or 
	((Chr2_Chr3_Row1128 #= 717) and (Chr2_Chr3_Freq1128 #= 5)) or 
	((Chr2_Chr3_Row1128 #= 724) and (Chr2_Chr3_Freq1128 #= 5)) or 
	((Chr2_Chr3_Row1128 #= 725) and (Chr2_Chr3_Freq1128 #= 5)) or 
	((Chr2_Chr3_Row1128 #= 726) and (Chr2_Chr3_Freq1128 #= 4)) or 
	((Chr2_Chr3_Row1128 #= 727) and (Chr2_Chr3_Freq1128 #= 4)) or 
	((Chr2_Chr3_Row1128 #= 728) and (Chr2_Chr3_Freq1128 #= 3)) or 
	((Chr2_Chr3_Row1128 #= 729) and (Chr2_Chr3_Freq1128 #= 3)) or 
	((Chr2_Chr3_Row1128 #= 730) and (Chr2_Chr3_Freq1128 #= 2)) or 
	((Chr2_Chr3_Row1128 #= 731) and (Chr2_Chr3_Freq1128 #= 2)) or 
	((Chr2_Chr3_Row1128 #= 732) and (Chr2_Chr3_Freq1128 #= 1)) or 
	((Chr2_Chr3_Row1128 #= 733) and (Chr2_Chr3_Freq1128 #= 1)) or 
	((Chr2_Chr3_Row1128 #= 734) and (Chr2_Chr3_Freq1128 #= 2)) or 
	((Chr2_Chr3_Row1128 #= 735) and (Chr2_Chr3_Freq1128 #= 1)) or 
	((Chr2_Chr3_Row1128 #= 736) and (Chr2_Chr3_Freq1128 #= 1)) or 
	((Chr2_Chr3_Row1128 #= 737) and (Chr2_Chr3_Freq1128 #= 1)) or 
	((Chr2_Chr3_Row1128 #= 738) and (Chr2_Chr3_Freq1128 #= 1)) or 
	((Chr2_Chr3_Freq1128 #= 0) and (Chr2_Chr3_Row1128 #= 0)),

	((Chr2_Chr3_Row1129 #= 700) and (Chr2_Chr3_Freq1129 #= 1)) or 
	((Chr2_Chr3_Row1129 #= 701) and (Chr2_Chr3_Freq1129 #= 1)) or 
	((Chr2_Chr3_Row1129 #= 702) and (Chr2_Chr3_Freq1129 #= 1)) or 
	((Chr2_Chr3_Row1129 #= 703) and (Chr2_Chr3_Freq1129 #= 1)) or 
	((Chr2_Chr3_Row1129 #= 704) and (Chr2_Chr3_Freq1129 #= 1)) or 
	((Chr2_Chr3_Row1129 #= 705) and (Chr2_Chr3_Freq1129 #= 1)) or 
	((Chr2_Chr3_Row1129 #= 706) and (Chr2_Chr3_Freq1129 #= 1)) or 
	((Chr2_Chr3_Row1129 #= 707) and (Chr2_Chr3_Freq1129 #= 1)) or 
	((Chr2_Chr3_Row1129 #= 708) and (Chr2_Chr3_Freq1129 #= 2)) or 
	((Chr2_Chr3_Row1129 #= 709) and (Chr2_Chr3_Freq1129 #= 2)) or 
	((Chr2_Chr3_Row1129 #= 710) and (Chr2_Chr3_Freq1129 #= 1)) or 
	((Chr2_Chr3_Row1129 #= 711) and (Chr2_Chr3_Freq1129 #= 2)) or 
	((Chr2_Chr3_Row1129 #= 712) and (Chr2_Chr3_Freq1129 #= 2)) or 
	((Chr2_Chr3_Row1129 #= 713) and (Chr2_Chr3_Freq1129 #= 2)) or 
	((Chr2_Chr3_Row1129 #= 714) and (Chr2_Chr3_Freq1129 #= 2)) or 
	((Chr2_Chr3_Row1129 #= 715) and (Chr2_Chr3_Freq1129 #= 3)) or 
	((Chr2_Chr3_Row1129 #= 716) and (Chr2_Chr3_Freq1129 #= 4)) or 
	((Chr2_Chr3_Row1129 #= 717) and (Chr2_Chr3_Freq1129 #= 5)) or 
	((Chr2_Chr3_Row1129 #= 724) and (Chr2_Chr3_Freq1129 #= 5)) or 
	((Chr2_Chr3_Row1129 #= 725) and (Chr2_Chr3_Freq1129 #= 5)) or 
	((Chr2_Chr3_Row1129 #= 726) and (Chr2_Chr3_Freq1129 #= 3)) or 
	((Chr2_Chr3_Row1129 #= 727) and (Chr2_Chr3_Freq1129 #= 3)) or 
	((Chr2_Chr3_Row1129 #= 728) and (Chr2_Chr3_Freq1129 #= 2)) or 
	((Chr2_Chr3_Row1129 #= 729) and (Chr2_Chr3_Freq1129 #= 2)) or 
	((Chr2_Chr3_Row1129 #= 730) and (Chr2_Chr3_Freq1129 #= 2)) or 
	((Chr2_Chr3_Row1129 #= 731) and (Chr2_Chr3_Freq1129 #= 1)) or 
	((Chr2_Chr3_Row1129 #= 732) and (Chr2_Chr3_Freq1129 #= 1)) or 
	((Chr2_Chr3_Row1129 #= 733) and (Chr2_Chr3_Freq1129 #= 1)) or 
	((Chr2_Chr3_Row1129 #= 734) and (Chr2_Chr3_Freq1129 #= 1)) or 
	((Chr2_Chr3_Row1129 #= 735) and (Chr2_Chr3_Freq1129 #= 1)) or 
	((Chr2_Chr3_Row1129 #= 736) and (Chr2_Chr3_Freq1129 #= 1)) or 
	((Chr2_Chr3_Row1129 #= 737) and (Chr2_Chr3_Freq1129 #= 1)) or 
	((Chr2_Chr3_Row1129 #= 738) and (Chr2_Chr3_Freq1129 #= 1)) or 
	((Chr2_Chr3_Row1129 #= 739) and (Chr2_Chr3_Freq1129 #= 1)) or 
	((Chr2_Chr3_Row1129 #= 740) and (Chr2_Chr3_Freq1129 #= 1)) or 
	((Chr2_Chr3_Freq1129 #= 0) and (Chr2_Chr3_Row1129 #= 0)),

	((Chr2_Chr3_Row1130 #= 699) and (Chr2_Chr3_Freq1130 #= 1)) or 
	((Chr2_Chr3_Row1130 #= 701) and (Chr2_Chr3_Freq1130 #= 1)) or 
	((Chr2_Chr3_Row1130 #= 702) and (Chr2_Chr3_Freq1130 #= 1)) or 
	((Chr2_Chr3_Row1130 #= 703) and (Chr2_Chr3_Freq1130 #= 1)) or 
	((Chr2_Chr3_Row1130 #= 704) and (Chr2_Chr3_Freq1130 #= 1)) or 
	((Chr2_Chr3_Row1130 #= 705) and (Chr2_Chr3_Freq1130 #= 1)) or 
	((Chr2_Chr3_Row1130 #= 706) and (Chr2_Chr3_Freq1130 #= 1)) or 
	((Chr2_Chr3_Row1130 #= 707) and (Chr2_Chr3_Freq1130 #= 1)) or 
	((Chr2_Chr3_Row1130 #= 708) and (Chr2_Chr3_Freq1130 #= 1)) or 
	((Chr2_Chr3_Row1130 #= 709) and (Chr2_Chr3_Freq1130 #= 1)) or 
	((Chr2_Chr3_Row1130 #= 710) and (Chr2_Chr3_Freq1130 #= 1)) or 
	((Chr2_Chr3_Row1130 #= 711) and (Chr2_Chr3_Freq1130 #= 2)) or 
	((Chr2_Chr3_Row1130 #= 712) and (Chr2_Chr3_Freq1130 #= 2)) or 
	((Chr2_Chr3_Row1130 #= 713) and (Chr2_Chr3_Freq1130 #= 2)) or 
	((Chr2_Chr3_Row1130 #= 714) and (Chr2_Chr3_Freq1130 #= 2)) or 
	((Chr2_Chr3_Row1130 #= 715) and (Chr2_Chr3_Freq1130 #= 3)) or 
	((Chr2_Chr3_Row1130 #= 716) and (Chr2_Chr3_Freq1130 #= 3)) or 
	((Chr2_Chr3_Row1130 #= 717) and (Chr2_Chr3_Freq1130 #= 4)) or 
	((Chr2_Chr3_Row1130 #= 724) and (Chr2_Chr3_Freq1130 #= 4)) or 
	((Chr2_Chr3_Row1130 #= 725) and (Chr2_Chr3_Freq1130 #= 3)) or 
	((Chr2_Chr3_Row1130 #= 726) and (Chr2_Chr3_Freq1130 #= 3)) or 
	((Chr2_Chr3_Row1130 #= 727) and (Chr2_Chr3_Freq1130 #= 2)) or 
	((Chr2_Chr3_Row1130 #= 728) and (Chr2_Chr3_Freq1130 #= 2)) or 
	((Chr2_Chr3_Row1130 #= 729) and (Chr2_Chr3_Freq1130 #= 2)) or 
	((Chr2_Chr3_Row1130 #= 730) and (Chr2_Chr3_Freq1130 #= 2)) or 
	((Chr2_Chr3_Row1130 #= 731) and (Chr2_Chr3_Freq1130 #= 1)) or 
	((Chr2_Chr3_Row1130 #= 732) and (Chr2_Chr3_Freq1130 #= 1)) or 
	((Chr2_Chr3_Row1130 #= 733) and (Chr2_Chr3_Freq1130 #= 1)) or 
	((Chr2_Chr3_Row1130 #= 734) and (Chr2_Chr3_Freq1130 #= 1)) or 
	((Chr2_Chr3_Row1130 #= 735) and (Chr2_Chr3_Freq1130 #= 1)) or 
	((Chr2_Chr3_Row1130 #= 736) and (Chr2_Chr3_Freq1130 #= 1)) or 
	((Chr2_Chr3_Freq1130 #= 0) and (Chr2_Chr3_Row1130 #= 0)),

	((Chr2_Chr3_Row1131 #= 701) and (Chr2_Chr3_Freq1131 #= 1)) or 
	((Chr2_Chr3_Row1131 #= 704) and (Chr2_Chr3_Freq1131 #= 1)) or 
	((Chr2_Chr3_Row1131 #= 705) and (Chr2_Chr3_Freq1131 #= 1)) or 
	((Chr2_Chr3_Row1131 #= 706) and (Chr2_Chr3_Freq1131 #= 1)) or 
	((Chr2_Chr3_Row1131 #= 707) and (Chr2_Chr3_Freq1131 #= 1)) or 
	((Chr2_Chr3_Row1131 #= 708) and (Chr2_Chr3_Freq1131 #= 1)) or 
	((Chr2_Chr3_Row1131 #= 709) and (Chr2_Chr3_Freq1131 #= 1)) or 
	((Chr2_Chr3_Row1131 #= 710) and (Chr2_Chr3_Freq1131 #= 1)) or 
	((Chr2_Chr3_Row1131 #= 711) and (Chr2_Chr3_Freq1131 #= 1)) or 
	((Chr2_Chr3_Row1131 #= 712) and (Chr2_Chr3_Freq1131 #= 2)) or 
	((Chr2_Chr3_Row1131 #= 713) and (Chr2_Chr3_Freq1131 #= 2)) or 
	((Chr2_Chr3_Row1131 #= 714) and (Chr2_Chr3_Freq1131 #= 2)) or 
	((Chr2_Chr3_Row1131 #= 715) and (Chr2_Chr3_Freq1131 #= 2)) or 
	((Chr2_Chr3_Row1131 #= 716) and (Chr2_Chr3_Freq1131 #= 3)) or 
	((Chr2_Chr3_Row1131 #= 717) and (Chr2_Chr3_Freq1131 #= 3)) or 
	((Chr2_Chr3_Row1131 #= 724) and (Chr2_Chr3_Freq1131 #= 3)) or 
	((Chr2_Chr3_Row1131 #= 725) and (Chr2_Chr3_Freq1131 #= 3)) or 
	((Chr2_Chr3_Row1131 #= 726) and (Chr2_Chr3_Freq1131 #= 3)) or 
	((Chr2_Chr3_Row1131 #= 727) and (Chr2_Chr3_Freq1131 #= 2)) or 
	((Chr2_Chr3_Row1131 #= 728) and (Chr2_Chr3_Freq1131 #= 1)) or 
	((Chr2_Chr3_Row1131 #= 729) and (Chr2_Chr3_Freq1131 #= 1)) or 
	((Chr2_Chr3_Row1131 #= 730) and (Chr2_Chr3_Freq1131 #= 1)) or 
	((Chr2_Chr3_Row1131 #= 731) and (Chr2_Chr3_Freq1131 #= 1)) or 
	((Chr2_Chr3_Row1131 #= 732) and (Chr2_Chr3_Freq1131 #= 1)) or 
	((Chr2_Chr3_Row1131 #= 733) and (Chr2_Chr3_Freq1131 #= 1)) or 
	((Chr2_Chr3_Row1131 #= 734) and (Chr2_Chr3_Freq1131 #= 1)) or 
	((Chr2_Chr3_Row1131 #= 735) and (Chr2_Chr3_Freq1131 #= 1)) or 
	((Chr2_Chr3_Row1131 #= 736) and (Chr2_Chr3_Freq1131 #= 1)) or 
	((Chr2_Chr3_Freq1131 #= 0) and (Chr2_Chr3_Row1131 #= 0)),

	((Chr2_Chr3_Row1132 #= 700) and (Chr2_Chr3_Freq1132 #= 1)) or 
	((Chr2_Chr3_Row1132 #= 701) and (Chr2_Chr3_Freq1132 #= 1)) or 
	((Chr2_Chr3_Row1132 #= 703) and (Chr2_Chr3_Freq1132 #= 1)) or 
	((Chr2_Chr3_Row1132 #= 704) and (Chr2_Chr3_Freq1132 #= 1)) or 
	((Chr2_Chr3_Row1132 #= 705) and (Chr2_Chr3_Freq1132 #= 1)) or 
	((Chr2_Chr3_Row1132 #= 706) and (Chr2_Chr3_Freq1132 #= 1)) or 
	((Chr2_Chr3_Row1132 #= 707) and (Chr2_Chr3_Freq1132 #= 1)) or 
	((Chr2_Chr3_Row1132 #= 708) and (Chr2_Chr3_Freq1132 #= 1)) or 
	((Chr2_Chr3_Row1132 #= 709) and (Chr2_Chr3_Freq1132 #= 1)) or 
	((Chr2_Chr3_Row1132 #= 710) and (Chr2_Chr3_Freq1132 #= 2)) or 
	((Chr2_Chr3_Row1132 #= 711) and (Chr2_Chr3_Freq1132 #= 2)) or 
	((Chr2_Chr3_Row1132 #= 712) and (Chr2_Chr3_Freq1132 #= 1)) or 
	((Chr2_Chr3_Row1132 #= 713) and (Chr2_Chr3_Freq1132 #= 2)) or 
	((Chr2_Chr3_Row1132 #= 714) and (Chr2_Chr3_Freq1132 #= 2)) or 
	((Chr2_Chr3_Row1132 #= 715) and (Chr2_Chr3_Freq1132 #= 3)) or 
	((Chr2_Chr3_Row1132 #= 716) and (Chr2_Chr3_Freq1132 #= 3)) or 
	((Chr2_Chr3_Row1132 #= 717) and (Chr2_Chr3_Freq1132 #= 3)) or 
	((Chr2_Chr3_Row1132 #= 724) and (Chr2_Chr3_Freq1132 #= 3)) or 
	((Chr2_Chr3_Row1132 #= 725) and (Chr2_Chr3_Freq1132 #= 3)) or 
	((Chr2_Chr3_Row1132 #= 726) and (Chr2_Chr3_Freq1132 #= 3)) or 
	((Chr2_Chr3_Row1132 #= 727) and (Chr2_Chr3_Freq1132 #= 2)) or 
	((Chr2_Chr3_Row1132 #= 728) and (Chr2_Chr3_Freq1132 #= 2)) or 
	((Chr2_Chr3_Row1132 #= 729) and (Chr2_Chr3_Freq1132 #= 2)) or 
	((Chr2_Chr3_Row1132 #= 730) and (Chr2_Chr3_Freq1132 #= 1)) or 
	((Chr2_Chr3_Row1132 #= 731) and (Chr2_Chr3_Freq1132 #= 1)) or 
	((Chr2_Chr3_Row1132 #= 732) and (Chr2_Chr3_Freq1132 #= 1)) or 
	((Chr2_Chr3_Row1132 #= 733) and (Chr2_Chr3_Freq1132 #= 1)) or 
	((Chr2_Chr3_Row1132 #= 734) and (Chr2_Chr3_Freq1132 #= 1)) or 
	((Chr2_Chr3_Row1132 #= 735) and (Chr2_Chr3_Freq1132 #= 1)) or 
	((Chr2_Chr3_Row1132 #= 736) and (Chr2_Chr3_Freq1132 #= 1)) or 
	((Chr2_Chr3_Row1132 #= 738) and (Chr2_Chr3_Freq1132 #= 1)) or 
	((Chr2_Chr3_Row1132 #= 739) and (Chr2_Chr3_Freq1132 #= 1)) or 
	((Chr2_Chr3_Freq1132 #= 0) and (Chr2_Chr3_Row1132 #= 0)),

	((Chr2_Chr3_Row1133 #= 695) and (Chr2_Chr3_Freq1133 #= 1)) or 
	((Chr2_Chr3_Row1133 #= 697) and (Chr2_Chr3_Freq1133 #= 1)) or 
	((Chr2_Chr3_Row1133 #= 699) and (Chr2_Chr3_Freq1133 #= 1)) or 
	((Chr2_Chr3_Row1133 #= 700) and (Chr2_Chr3_Freq1133 #= 1)) or 
	((Chr2_Chr3_Row1133 #= 701) and (Chr2_Chr3_Freq1133 #= 1)) or 
	((Chr2_Chr3_Row1133 #= 702) and (Chr2_Chr3_Freq1133 #= 1)) or 
	((Chr2_Chr3_Row1133 #= 703) and (Chr2_Chr3_Freq1133 #= 1)) or 
	((Chr2_Chr3_Row1133 #= 704) and (Chr2_Chr3_Freq1133 #= 1)) or 
	((Chr2_Chr3_Row1133 #= 705) and (Chr2_Chr3_Freq1133 #= 1)) or 
	((Chr2_Chr3_Row1133 #= 706) and (Chr2_Chr3_Freq1133 #= 1)) or 
	((Chr2_Chr3_Row1133 #= 707) and (Chr2_Chr3_Freq1133 #= 1)) or 
	((Chr2_Chr3_Row1133 #= 708) and (Chr2_Chr3_Freq1133 #= 1)) or 
	((Chr2_Chr3_Row1133 #= 709) and (Chr2_Chr3_Freq1133 #= 1)) or 
	((Chr2_Chr3_Row1133 #= 710) and (Chr2_Chr3_Freq1133 #= 2)) or 
	((Chr2_Chr3_Row1133 #= 711) and (Chr2_Chr3_Freq1133 #= 2)) or 
	((Chr2_Chr3_Row1133 #= 712) and (Chr2_Chr3_Freq1133 #= 2)) or 
	((Chr2_Chr3_Row1133 #= 713) and (Chr2_Chr3_Freq1133 #= 2)) or 
	((Chr2_Chr3_Row1133 #= 714) and (Chr2_Chr3_Freq1133 #= 2)) or 
	((Chr2_Chr3_Row1133 #= 715) and (Chr2_Chr3_Freq1133 #= 3)) or 
	((Chr2_Chr3_Row1133 #= 716) and (Chr2_Chr3_Freq1133 #= 3)) or 
	((Chr2_Chr3_Row1133 #= 717) and (Chr2_Chr3_Freq1133 #= 3)) or 
	((Chr2_Chr3_Row1133 #= 724) and (Chr2_Chr3_Freq1133 #= 3)) or 
	((Chr2_Chr3_Row1133 #= 725) and (Chr2_Chr3_Freq1133 #= 3)) or 
	((Chr2_Chr3_Row1133 #= 726) and (Chr2_Chr3_Freq1133 #= 2)) or 
	((Chr2_Chr3_Row1133 #= 727) and (Chr2_Chr3_Freq1133 #= 2)) or 
	((Chr2_Chr3_Row1133 #= 728) and (Chr2_Chr3_Freq1133 #= 2)) or 
	((Chr2_Chr3_Row1133 #= 729) and (Chr2_Chr3_Freq1133 #= 1)) or 
	((Chr2_Chr3_Row1133 #= 730) and (Chr2_Chr3_Freq1133 #= 1)) or 
	((Chr2_Chr3_Row1133 #= 731) and (Chr2_Chr3_Freq1133 #= 1)) or 
	((Chr2_Chr3_Row1133 #= 732) and (Chr2_Chr3_Freq1133 #= 1)) or 
	((Chr2_Chr3_Row1133 #= 733) and (Chr2_Chr3_Freq1133 #= 1)) or 
	((Chr2_Chr3_Row1133 #= 734) and (Chr2_Chr3_Freq1133 #= 1)) or 
	((Chr2_Chr3_Row1133 #= 735) and (Chr2_Chr3_Freq1133 #= 1)) or 
	((Chr2_Chr3_Row1133 #= 736) and (Chr2_Chr3_Freq1133 #= 1)) or 
	((Chr2_Chr3_Row1133 #= 738) and (Chr2_Chr3_Freq1133 #= 1)) or 
	((Chr2_Chr3_Row1133 #= 739) and (Chr2_Chr3_Freq1133 #= 1)) or 
	((Chr2_Chr3_Freq1133 #= 0) and (Chr2_Chr3_Row1133 #= 0)),

	((Chr2_Chr3_Row1134 #= 699) and (Chr2_Chr3_Freq1134 #= 1)) or 
	((Chr2_Chr3_Row1134 #= 701) and (Chr2_Chr3_Freq1134 #= 1)) or 
	((Chr2_Chr3_Row1134 #= 702) and (Chr2_Chr3_Freq1134 #= 1)) or 
	((Chr2_Chr3_Row1134 #= 704) and (Chr2_Chr3_Freq1134 #= 1)) or 
	((Chr2_Chr3_Row1134 #= 705) and (Chr2_Chr3_Freq1134 #= 1)) or 
	((Chr2_Chr3_Row1134 #= 706) and (Chr2_Chr3_Freq1134 #= 1)) or 
	((Chr2_Chr3_Row1134 #= 707) and (Chr2_Chr3_Freq1134 #= 1)) or 
	((Chr2_Chr3_Row1134 #= 708) and (Chr2_Chr3_Freq1134 #= 1)) or 
	((Chr2_Chr3_Row1134 #= 709) and (Chr2_Chr3_Freq1134 #= 1)) or 
	((Chr2_Chr3_Row1134 #= 710) and (Chr2_Chr3_Freq1134 #= 1)) or 
	((Chr2_Chr3_Row1134 #= 711) and (Chr2_Chr3_Freq1134 #= 1)) or 
	((Chr2_Chr3_Row1134 #= 712) and (Chr2_Chr3_Freq1134 #= 2)) or 
	((Chr2_Chr3_Row1134 #= 713) and (Chr2_Chr3_Freq1134 #= 2)) or 
	((Chr2_Chr3_Row1134 #= 714) and (Chr2_Chr3_Freq1134 #= 2)) or 
	((Chr2_Chr3_Row1134 #= 715) and (Chr2_Chr3_Freq1134 #= 2)) or 
	((Chr2_Chr3_Row1134 #= 716) and (Chr2_Chr3_Freq1134 #= 2)) or 
	((Chr2_Chr3_Row1134 #= 717) and (Chr2_Chr3_Freq1134 #= 3)) or 
	((Chr2_Chr3_Row1134 #= 724) and (Chr2_Chr3_Freq1134 #= 2)) or 
	((Chr2_Chr3_Row1134 #= 725) and (Chr2_Chr3_Freq1134 #= 2)) or 
	((Chr2_Chr3_Row1134 #= 726) and (Chr2_Chr3_Freq1134 #= 2)) or 
	((Chr2_Chr3_Row1134 #= 727) and (Chr2_Chr3_Freq1134 #= 2)) or 
	((Chr2_Chr3_Row1134 #= 728) and (Chr2_Chr3_Freq1134 #= 1)) or 
	((Chr2_Chr3_Row1134 #= 729) and (Chr2_Chr3_Freq1134 #= 1)) or 
	((Chr2_Chr3_Row1134 #= 730) and (Chr2_Chr3_Freq1134 #= 1)) or 
	((Chr2_Chr3_Row1134 #= 731) and (Chr2_Chr3_Freq1134 #= 1)) or 
	((Chr2_Chr3_Row1134 #= 732) and (Chr2_Chr3_Freq1134 #= 1)) or 
	((Chr2_Chr3_Row1134 #= 734) and (Chr2_Chr3_Freq1134 #= 1)) or 
	((Chr2_Chr3_Row1134 #= 735) and (Chr2_Chr3_Freq1134 #= 1)) or 
	((Chr2_Chr3_Row1134 #= 737) and (Chr2_Chr3_Freq1134 #= 1)) or 
	((Chr2_Chr3_Row1134 #= 740) and (Chr2_Chr3_Freq1134 #= 1)) or 
	((Chr2_Chr3_Freq1134 #= 0) and (Chr2_Chr3_Row1134 #= 0)),

	((Chr2_Chr3_Row1135 #= 698) and (Chr2_Chr3_Freq1135 #= 1)) or 
	((Chr2_Chr3_Row1135 #= 699) and (Chr2_Chr3_Freq1135 #= 1)) or 
	((Chr2_Chr3_Row1135 #= 701) and (Chr2_Chr3_Freq1135 #= 1)) or 
	((Chr2_Chr3_Row1135 #= 702) and (Chr2_Chr3_Freq1135 #= 1)) or 
	((Chr2_Chr3_Row1135 #= 704) and (Chr2_Chr3_Freq1135 #= 1)) or 
	((Chr2_Chr3_Row1135 #= 705) and (Chr2_Chr3_Freq1135 #= 1)) or 
	((Chr2_Chr3_Row1135 #= 706) and (Chr2_Chr3_Freq1135 #= 1)) or 
	((Chr2_Chr3_Row1135 #= 707) and (Chr2_Chr3_Freq1135 #= 1)) or 
	((Chr2_Chr3_Row1135 #= 708) and (Chr2_Chr3_Freq1135 #= 1)) or 
	((Chr2_Chr3_Row1135 #= 709) and (Chr2_Chr3_Freq1135 #= 1)) or 
	((Chr2_Chr3_Row1135 #= 710) and (Chr2_Chr3_Freq1135 #= 2)) or 
	((Chr2_Chr3_Row1135 #= 711) and (Chr2_Chr3_Freq1135 #= 2)) or 
	((Chr2_Chr3_Row1135 #= 712) and (Chr2_Chr3_Freq1135 #= 1)) or 
	((Chr2_Chr3_Row1135 #= 713) and (Chr2_Chr3_Freq1135 #= 2)) or 
	((Chr2_Chr3_Row1135 #= 714) and (Chr2_Chr3_Freq1135 #= 1)) or 
	((Chr2_Chr3_Row1135 #= 715) and (Chr2_Chr3_Freq1135 #= 2)) or 
	((Chr2_Chr3_Row1135 #= 716) and (Chr2_Chr3_Freq1135 #= 2)) or 
	((Chr2_Chr3_Row1135 #= 717) and (Chr2_Chr3_Freq1135 #= 2)) or 
	((Chr2_Chr3_Row1135 #= 724) and (Chr2_Chr3_Freq1135 #= 2)) or 
	((Chr2_Chr3_Row1135 #= 725) and (Chr2_Chr3_Freq1135 #= 2)) or 
	((Chr2_Chr3_Row1135 #= 726) and (Chr2_Chr3_Freq1135 #= 2)) or 
	((Chr2_Chr3_Row1135 #= 727) and (Chr2_Chr3_Freq1135 #= 1)) or 
	((Chr2_Chr3_Row1135 #= 728) and (Chr2_Chr3_Freq1135 #= 1)) or 
	((Chr2_Chr3_Row1135 #= 729) and (Chr2_Chr3_Freq1135 #= 1)) or 
	((Chr2_Chr3_Row1135 #= 730) and (Chr2_Chr3_Freq1135 #= 1)) or 
	((Chr2_Chr3_Row1135 #= 731) and (Chr2_Chr3_Freq1135 #= 1)) or 
	((Chr2_Chr3_Row1135 #= 734) and (Chr2_Chr3_Freq1135 #= 1)) or 
	((Chr2_Chr3_Row1135 #= 735) and (Chr2_Chr3_Freq1135 #= 1)) or 
	((Chr2_Chr3_Row1135 #= 736) and (Chr2_Chr3_Freq1135 #= 1)) or 
	((Chr2_Chr3_Freq1135 #= 0) and (Chr2_Chr3_Row1135 #= 0)),

	((Chr2_Chr3_Row1136 #= 698) and (Chr2_Chr3_Freq1136 #= 1)) or 
	((Chr2_Chr3_Row1136 #= 701) and (Chr2_Chr3_Freq1136 #= 1)) or 
	((Chr2_Chr3_Row1136 #= 702) and (Chr2_Chr3_Freq1136 #= 1)) or 
	((Chr2_Chr3_Row1136 #= 703) and (Chr2_Chr3_Freq1136 #= 1)) or 
	((Chr2_Chr3_Row1136 #= 704) and (Chr2_Chr3_Freq1136 #= 1)) or 
	((Chr2_Chr3_Row1136 #= 705) and (Chr2_Chr3_Freq1136 #= 1)) or 
	((Chr2_Chr3_Row1136 #= 706) and (Chr2_Chr3_Freq1136 #= 1)) or 
	((Chr2_Chr3_Row1136 #= 707) and (Chr2_Chr3_Freq1136 #= 1)) or 
	((Chr2_Chr3_Row1136 #= 708) and (Chr2_Chr3_Freq1136 #= 1)) or 
	((Chr2_Chr3_Row1136 #= 709) and (Chr2_Chr3_Freq1136 #= 1)) or 
	((Chr2_Chr3_Row1136 #= 710) and (Chr2_Chr3_Freq1136 #= 1)) or 
	((Chr2_Chr3_Row1136 #= 711) and (Chr2_Chr3_Freq1136 #= 1)) or 
	((Chr2_Chr3_Row1136 #= 712) and (Chr2_Chr3_Freq1136 #= 2)) or 
	((Chr2_Chr3_Row1136 #= 713) and (Chr2_Chr3_Freq1136 #= 2)) or 
	((Chr2_Chr3_Row1136 #= 714) and (Chr2_Chr3_Freq1136 #= 2)) or 
	((Chr2_Chr3_Row1136 #= 715) and (Chr2_Chr3_Freq1136 #= 2)) or 
	((Chr2_Chr3_Row1136 #= 716) and (Chr2_Chr3_Freq1136 #= 2)) or 
	((Chr2_Chr3_Row1136 #= 717) and (Chr2_Chr3_Freq1136 #= 2)) or 
	((Chr2_Chr3_Row1136 #= 724) and (Chr2_Chr3_Freq1136 #= 1)) or 
	((Chr2_Chr3_Row1136 #= 725) and (Chr2_Chr3_Freq1136 #= 1)) or 
	((Chr2_Chr3_Row1136 #= 726) and (Chr2_Chr3_Freq1136 #= 2)) or 
	((Chr2_Chr3_Row1136 #= 727) and (Chr2_Chr3_Freq1136 #= 1)) or 
	((Chr2_Chr3_Row1136 #= 728) and (Chr2_Chr3_Freq1136 #= 1)) or 
	((Chr2_Chr3_Row1136 #= 729) and (Chr2_Chr3_Freq1136 #= 1)) or 
	((Chr2_Chr3_Row1136 #= 730) and (Chr2_Chr3_Freq1136 #= 1)) or 
	((Chr2_Chr3_Row1136 #= 731) and (Chr2_Chr3_Freq1136 #= 1)) or 
	((Chr2_Chr3_Row1136 #= 732) and (Chr2_Chr3_Freq1136 #= 1)) or 
	((Chr2_Chr3_Row1136 #= 733) and (Chr2_Chr3_Freq1136 #= 1)) or 
	((Chr2_Chr3_Row1136 #= 734) and (Chr2_Chr3_Freq1136 #= 1)) or 
	((Chr2_Chr3_Row1136 #= 737) and (Chr2_Chr3_Freq1136 #= 1)) or 
	((Chr2_Chr3_Freq1136 #= 0) and (Chr2_Chr3_Row1136 #= 0)),

	((Chr2_Chr3_Row1137 #= 698) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 699) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 700) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 701) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 702) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 703) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 704) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 705) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 706) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 707) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 708) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 709) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 710) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 711) and (Chr2_Chr3_Freq1137 #= 2)) or 
	((Chr2_Chr3_Row1137 #= 712) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 713) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 714) and (Chr2_Chr3_Freq1137 #= 2)) or 
	((Chr2_Chr3_Row1137 #= 715) and (Chr2_Chr3_Freq1137 #= 2)) or 
	((Chr2_Chr3_Row1137 #= 716) and (Chr2_Chr3_Freq1137 #= 2)) or 
	((Chr2_Chr3_Row1137 #= 717) and (Chr2_Chr3_Freq1137 #= 2)) or 
	((Chr2_Chr3_Row1137 #= 724) and (Chr2_Chr3_Freq1137 #= 2)) or 
	((Chr2_Chr3_Row1137 #= 725) and (Chr2_Chr3_Freq1137 #= 2)) or 
	((Chr2_Chr3_Row1137 #= 726) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 727) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 728) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 729) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 730) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 731) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 732) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 733) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 734) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 735) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 736) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 737) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 738) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Row1137 #= 739) and (Chr2_Chr3_Freq1137 #= 1)) or 
	((Chr2_Chr3_Freq1137 #= 0) and (Chr2_Chr3_Row1137 #= 0)),

	((Chr2_Chr3_Row1138 #= 695) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 697) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 698) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 699) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 700) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 702) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 703) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 704) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 705) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 706) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 707) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 708) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 709) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 710) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 711) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 712) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 713) and (Chr2_Chr3_Freq1138 #= 2)) or 
	((Chr2_Chr3_Row1138 #= 714) and (Chr2_Chr3_Freq1138 #= 2)) or 
	((Chr2_Chr3_Row1138 #= 715) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 716) and (Chr2_Chr3_Freq1138 #= 2)) or 
	((Chr2_Chr3_Row1138 #= 717) and (Chr2_Chr3_Freq1138 #= 2)) or 
	((Chr2_Chr3_Row1138 #= 724) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 725) and (Chr2_Chr3_Freq1138 #= 2)) or 
	((Chr2_Chr3_Row1138 #= 726) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 727) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 728) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 729) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 730) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 731) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 733) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 735) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 737) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Row1138 #= 738) and (Chr2_Chr3_Freq1138 #= 1)) or 
	((Chr2_Chr3_Freq1138 #= 0) and (Chr2_Chr3_Row1138 #= 0)),

	((Chr2_Chr3_Row1139 #= 695) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 696) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 698) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 699) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 701) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 702) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 703) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 704) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 705) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 706) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 707) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 708) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 709) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 710) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 711) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 712) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 713) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 714) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 715) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 716) and (Chr2_Chr3_Freq1139 #= 2)) or 
	((Chr2_Chr3_Row1139 #= 717) and (Chr2_Chr3_Freq1139 #= 2)) or 
	((Chr2_Chr3_Row1139 #= 724) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 725) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 726) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 727) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 728) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 729) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 730) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 731) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Row1139 #= 734) and (Chr2_Chr3_Freq1139 #= 1)) or 
	((Chr2_Chr3_Freq1139 #= 0) and (Chr2_Chr3_Row1139 #= 0)),

	((Chr2_Chr3_Row1140 #= 697) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Row1140 #= 698) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Row1140 #= 699) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Row1140 #= 700) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Row1140 #= 702) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Row1140 #= 703) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Row1140 #= 704) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Row1140 #= 705) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Row1140 #= 706) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Row1140 #= 707) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Row1140 #= 708) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Row1140 #= 709) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Row1140 #= 710) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Row1140 #= 711) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Row1140 #= 712) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Row1140 #= 713) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Row1140 #= 714) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Row1140 #= 715) and (Chr2_Chr3_Freq1140 #= 2)) or 
	((Chr2_Chr3_Row1140 #= 716) and (Chr2_Chr3_Freq1140 #= 2)) or 
	((Chr2_Chr3_Row1140 #= 717) and (Chr2_Chr3_Freq1140 #= 2)) or 
	((Chr2_Chr3_Row1140 #= 724) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Row1140 #= 725) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Row1140 #= 726) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Row1140 #= 727) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Row1140 #= 728) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Row1140 #= 729) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Row1140 #= 732) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Row1140 #= 735) and (Chr2_Chr3_Freq1140 #= 1)) or 
	((Chr2_Chr3_Freq1140 #= 0) and (Chr2_Chr3_Row1140 #= 0)),

	((Chr2_Chr3_Row1141 #= 692) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 693) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 694) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 695) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 696) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 697) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 698) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 699) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 701) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 702) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 703) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 704) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 705) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 706) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 707) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 708) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 709) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 710) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 711) and (Chr2_Chr3_Freq1141 #= 2)) or 
	((Chr2_Chr3_Row1141 #= 712) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 713) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 714) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 715) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 716) and (Chr2_Chr3_Freq1141 #= 2)) or 
	((Chr2_Chr3_Row1141 #= 717) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 724) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 725) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 726) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 727) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 728) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 729) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 730) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 731) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 732) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 734) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 735) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Row1141 #= 736) and (Chr2_Chr3_Freq1141 #= 1)) or 
	((Chr2_Chr3_Freq1141 #= 0) and (Chr2_Chr3_Row1141 #= 0)),

	((Chr2_Chr3_Row1142 #= 690) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 691) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 697) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 699) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 700) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 701) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 702) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 703) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 704) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 705) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 706) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 707) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 708) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 709) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 710) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 711) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 712) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 713) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 714) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 715) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 716) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 717) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 724) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 725) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 726) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 727) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 728) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 730) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 731) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 732) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 733) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 734) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 738) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Row1142 #= 739) and (Chr2_Chr3_Freq1142 #= 1)) or 
	((Chr2_Chr3_Freq1142 #= 0) and (Chr2_Chr3_Row1142 #= 0)),

	((Chr2_Chr3_Row1143 #= 693) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 694) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 696) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 698) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 700) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 701) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 702) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 703) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 704) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 705) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 706) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 707) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 708) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 709) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 710) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 711) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 712) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 713) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 714) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 715) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 716) and (Chr2_Chr3_Freq1143 #= 2)) or 
	((Chr2_Chr3_Row1143 #= 717) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 724) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 725) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 726) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 727) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 728) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 730) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 731) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 732) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 734) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Row1143 #= 735) and (Chr2_Chr3_Freq1143 #= 1)) or 
	((Chr2_Chr3_Freq1143 #= 0) and (Chr2_Chr3_Row1143 #= 0)),

	((Chr2_Chr3_Row1144 #= 694) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 695) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 697) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 699) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 700) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 701) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 702) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 703) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 704) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 705) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 706) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 707) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 708) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 709) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 710) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 711) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 712) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 713) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 714) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 715) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 716) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 717) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 724) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 725) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 726) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 727) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 728) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 733) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 734) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Row1144 #= 736) and (Chr2_Chr3_Freq1144 #= 1)) or 
	((Chr2_Chr3_Freq1144 #= 0) and (Chr2_Chr3_Row1144 #= 0)),

	((Chr2_Chr3_Row1145 #= 690) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 695) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 697) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 698) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 699) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 700) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 701) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 702) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 703) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 704) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 705) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 706) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 707) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 708) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 709) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 710) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 711) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 712) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 713) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 714) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 715) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 716) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 717) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 725) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 726) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Row1145 #= 728) and (Chr2_Chr3_Freq1145 #= 1)) or 
	((Chr2_Chr3_Freq1145 #= 0) and (Chr2_Chr3_Row1145 #= 0)),

	((Chr2_Chr3_Row1146 #= 691) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 693) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 694) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 695) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 696) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 697) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 698) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 699) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 700) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 701) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 702) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 703) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 704) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 705) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 706) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 707) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 708) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 709) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 710) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 711) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 712) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 713) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 714) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 715) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 716) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 717) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Row1146 #= 724) and (Chr2_Chr3_Freq1146 #= 1)) or 
	((Chr2_Chr3_Freq1146 #= 0) and (Chr2_Chr3_Row1146 #= 0)),

	((Chr2_Chr3_Row1147 #= 691) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 692) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 693) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 694) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 695) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 697) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 698) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 699) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 700) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 701) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 702) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 703) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 704) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 705) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 706) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 707) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 708) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 709) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 710) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 711) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 712) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 713) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 714) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 715) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 716) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 717) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Row1147 #= 727) and (Chr2_Chr3_Freq1147 #= 1)) or 
	((Chr2_Chr3_Freq1147 #= 0) and (Chr2_Chr3_Row1147 #= 0)),

	((Chr2_Chr3_Row1148 #= 690) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 691) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 692) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 693) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 694) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 695) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 696) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 697) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 698) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 699) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 700) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 701) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 702) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 703) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 704) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 705) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 706) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 707) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 708) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 709) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 710) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 711) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 712) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 713) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 714) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 715) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 716) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 717) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 725) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 726) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 727) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 729) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Row1148 #= 734) and (Chr2_Chr3_Freq1148 #= 1)) or 
	((Chr2_Chr3_Freq1148 #= 0) and (Chr2_Chr3_Row1148 #= 0)),

	((Chr2_Chr3_Row1149 #= 686) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 690) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 691) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 692) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 693) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 694) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 695) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 696) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 697) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 698) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 699) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 700) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 701) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 702) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 703) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 704) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 705) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 706) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 707) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 708) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 709) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 710) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 711) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 712) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 713) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 714) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 715) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 716) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 717) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 724) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 725) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 726) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 727) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Row1149 #= 740) and (Chr2_Chr3_Freq1149 #= 1)) or 
	((Chr2_Chr3_Freq1149 #= 0) and (Chr2_Chr3_Row1149 #= 0)),

	((Chr2_Chr3_Row1150 #= 683) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 684) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 687) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 688) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 689) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 690) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 691) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 692) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 693) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 694) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 695) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 697) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 698) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 699) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 700) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 702) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 703) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 704) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 705) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 706) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 707) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 708) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 709) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 710) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 711) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 712) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 713) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 714) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 715) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 716) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 717) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 724) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 725) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Row1150 #= 728) and (Chr2_Chr3_Freq1150 #= 1)) or 
	((Chr2_Chr3_Freq1150 #= 0) and (Chr2_Chr3_Row1150 #= 0)),

	((Chr2_Chr3_Row1151 #= 692) and (Chr2_Chr3_Freq1151 #= 1)) or 
	((Chr2_Chr3_Row1151 #= 696) and (Chr2_Chr3_Freq1151 #= 1)) or 
	((Chr2_Chr3_Row1151 #= 698) and (Chr2_Chr3_Freq1151 #= 1)) or 
	((Chr2_Chr3_Row1151 #= 699) and (Chr2_Chr3_Freq1151 #= 1)) or 
	((Chr2_Chr3_Row1151 #= 700) and (Chr2_Chr3_Freq1151 #= 1)) or 
	((Chr2_Chr3_Row1151 #= 701) and (Chr2_Chr3_Freq1151 #= 1)) or 
	((Chr2_Chr3_Row1151 #= 702) and (Chr2_Chr3_Freq1151 #= 1)) or 
	((Chr2_Chr3_Row1151 #= 703) and (Chr2_Chr3_Freq1151 #= 1)) or 
	((Chr2_Chr3_Row1151 #= 704) and (Chr2_Chr3_Freq1151 #= 1)) or 
	((Chr2_Chr3_Row1151 #= 705) and (Chr2_Chr3_Freq1151 #= 1)) or 
	((Chr2_Chr3_Row1151 #= 706) and (Chr2_Chr3_Freq1151 #= 1)) or 
	((Chr2_Chr3_Row1151 #= 707) and (Chr2_Chr3_Freq1151 #= 1)) or 
	((Chr2_Chr3_Row1151 #= 708) and (Chr2_Chr3_Freq1151 #= 1)) or 
	((Chr2_Chr3_Row1151 #= 709) and (Chr2_Chr3_Freq1151 #= 1)) or 
	((Chr2_Chr3_Row1151 #= 710) and (Chr2_Chr3_Freq1151 #= 1)) or 
	((Chr2_Chr3_Row1151 #= 712) and (Chr2_Chr3_Freq1151 #= 1)) or 
	((Chr2_Chr3_Row1151 #= 713) and (Chr2_Chr3_Freq1151 #= 1)) or 
	((Chr2_Chr3_Row1151 #= 714) and (Chr2_Chr3_Freq1151 #= 1)) or 
	((Chr2_Chr3_Row1151 #= 716) and (Chr2_Chr3_Freq1151 #= 1)) or 
	((Chr2_Chr3_Row1151 #= 717) and (Chr2_Chr3_Freq1151 #= 1)) or 
	((Chr2_Chr3_Freq1151 #= 0) and (Chr2_Chr3_Row1151 #= 0)),

	((Chr2_Chr3_Row1152 #= 685) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 688) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 689) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 690) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 691) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 692) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 693) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 695) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 696) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 697) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 698) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 699) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 700) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 701) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 702) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 703) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 704) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 705) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 706) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 707) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 708) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 709) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 710) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 711) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 712) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 713) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 714) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 715) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 716) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 717) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 727) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 730) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 734) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Row1152 #= 736) and (Chr2_Chr3_Freq1152 #= 1)) or 
	((Chr2_Chr3_Freq1152 #= 0) and (Chr2_Chr3_Row1152 #= 0)),

	((Chr2_Chr3_Row1153 #= 684) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 685) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 688) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 689) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 690) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 691) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 692) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 693) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 694) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 695) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 696) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 697) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 698) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 699) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 700) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 701) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 702) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 703) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 704) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 705) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 706) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 707) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 708) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 709) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 710) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 711) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 712) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 713) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 714) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 715) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 716) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Row1153 #= 717) and (Chr2_Chr3_Freq1153 #= 1)) or 
	((Chr2_Chr3_Freq1153 #= 0) and (Chr2_Chr3_Row1153 #= 0)),

	((Chr2_Chr3_Row1154 #= 684) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 685) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 688) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 689) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 690) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 691) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 692) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 693) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 694) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 695) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 696) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 697) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 698) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 699) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 700) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 701) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 702) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 703) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 704) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 705) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 706) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 707) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 708) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 709) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 710) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 711) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 712) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 715) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 716) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Row1154 #= 717) and (Chr2_Chr3_Freq1154 #= 1)) or 
	((Chr2_Chr3_Freq1154 #= 0) and (Chr2_Chr3_Row1154 #= 0)),

	((Chr2_Chr3_Row1155 #= 683) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 684) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 685) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 686) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 687) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 688) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 689) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 690) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 691) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 692) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 693) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 694) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 695) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 696) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 697) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 698) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 699) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 700) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 701) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 702) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 703) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 704) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 705) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 706) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 707) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 708) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 709) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 710) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 712) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 713) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Row1155 #= 715) and (Chr2_Chr3_Freq1155 #= 1)) or 
	((Chr2_Chr3_Freq1155 #= 0) and (Chr2_Chr3_Row1155 #= 0)),

	((Chr2_Chr3_Row1156 #= 682) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 683) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 684) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 685) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 686) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 687) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 688) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 689) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 690) and (Chr2_Chr3_Freq1156 #= 2)) or 
	((Chr2_Chr3_Row1156 #= 691) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 692) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 693) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 694) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 695) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 696) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 697) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 698) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 699) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 700) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 701) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 702) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 703) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 704) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 705) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 706) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 707) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 708) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 709) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 710) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 711) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 712) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 714) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Row1156 #= 716) and (Chr2_Chr3_Freq1156 #= 1)) or 
	((Chr2_Chr3_Freq1156 #= 0) and (Chr2_Chr3_Row1156 #= 0)),

	((Chr2_Chr3_Row1157 #= 678) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 683) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 684) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 685) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 686) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 687) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 688) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 689) and (Chr2_Chr3_Freq1157 #= 2)) or 
	((Chr2_Chr3_Row1157 #= 690) and (Chr2_Chr3_Freq1157 #= 2)) or 
	((Chr2_Chr3_Row1157 #= 691) and (Chr2_Chr3_Freq1157 #= 2)) or 
	((Chr2_Chr3_Row1157 #= 692) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 693) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 694) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 695) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 696) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 697) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 698) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 699) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 700) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 701) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 702) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 703) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 704) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 705) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 706) and (Chr2_Chr3_Freq1157 #= 2)) or 
	((Chr2_Chr3_Row1157 #= 707) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 708) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 709) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 710) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Row1157 #= 711) and (Chr2_Chr3_Freq1157 #= 1)) or 
	((Chr2_Chr3_Freq1157 #= 0) and (Chr2_Chr3_Row1157 #= 0)),

	((Chr2_Chr3_Row1158 #= 683) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Row1158 #= 685) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Row1158 #= 686) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Row1158 #= 687) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Row1158 #= 688) and (Chr2_Chr3_Freq1158 #= 2)) or 
	((Chr2_Chr3_Row1158 #= 689) and (Chr2_Chr3_Freq1158 #= 2)) or 
	((Chr2_Chr3_Row1158 #= 690) and (Chr2_Chr3_Freq1158 #= 2)) or 
	((Chr2_Chr3_Row1158 #= 691) and (Chr2_Chr3_Freq1158 #= 2)) or 
	((Chr2_Chr3_Row1158 #= 692) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Row1158 #= 693) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Row1158 #= 694) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Row1158 #= 695) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Row1158 #= 696) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Row1158 #= 697) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Row1158 #= 698) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Row1158 #= 699) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Row1158 #= 700) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Row1158 #= 701) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Row1158 #= 702) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Row1158 #= 703) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Row1158 #= 704) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Row1158 #= 705) and (Chr2_Chr3_Freq1158 #= 2)) or 
	((Chr2_Chr3_Row1158 #= 706) and (Chr2_Chr3_Freq1158 #= 2)) or 
	((Chr2_Chr3_Row1158 #= 707) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Row1158 #= 708) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Row1158 #= 709) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Row1158 #= 710) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Row1158 #= 711) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Row1158 #= 713) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Row1158 #= 714) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Row1158 #= 715) and (Chr2_Chr3_Freq1158 #= 1)) or 
	((Chr2_Chr3_Freq1158 #= 0) and (Chr2_Chr3_Row1158 #= 0)),

	((Chr2_Chr3_Row1159 #= 684) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Row1159 #= 685) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Row1159 #= 686) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Row1159 #= 687) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Row1159 #= 688) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Row1159 #= 689) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Row1159 #= 690) and (Chr2_Chr3_Freq1159 #= 2)) or 
	((Chr2_Chr3_Row1159 #= 691) and (Chr2_Chr3_Freq1159 #= 2)) or 
	((Chr2_Chr3_Row1159 #= 692) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Row1159 #= 693) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Row1159 #= 694) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Row1159 #= 695) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Row1159 #= 696) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Row1159 #= 697) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Row1159 #= 698) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Row1159 #= 699) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Row1159 #= 700) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Row1159 #= 701) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Row1159 #= 702) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Row1159 #= 703) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Row1159 #= 704) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Row1159 #= 705) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Row1159 #= 706) and (Chr2_Chr3_Freq1159 #= 2)) or 
	((Chr2_Chr3_Row1159 #= 707) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Row1159 #= 708) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Row1159 #= 709) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Row1159 #= 710) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Row1159 #= 711) and (Chr2_Chr3_Freq1159 #= 1)) or 
	((Chr2_Chr3_Freq1159 #= 0) and (Chr2_Chr3_Row1159 #= 0)),

	((Chr2_Chr3_Row1160 #= 683) and (Chr2_Chr3_Freq1160 #= 1)) or 
	((Chr2_Chr3_Row1160 #= 685) and (Chr2_Chr3_Freq1160 #= 1)) or 
	((Chr2_Chr3_Row1160 #= 686) and (Chr2_Chr3_Freq1160 #= 1)) or 
	((Chr2_Chr3_Row1160 #= 687) and (Chr2_Chr3_Freq1160 #= 1)) or 
	((Chr2_Chr3_Row1160 #= 688) and (Chr2_Chr3_Freq1160 #= 1)) or 
	((Chr2_Chr3_Row1160 #= 689) and (Chr2_Chr3_Freq1160 #= 1)) or 
	((Chr2_Chr3_Row1160 #= 690) and (Chr2_Chr3_Freq1160 #= 1)) or 
	((Chr2_Chr3_Row1160 #= 691) and (Chr2_Chr3_Freq1160 #= 1)) or 
	((Chr2_Chr3_Row1160 #= 692) and (Chr2_Chr3_Freq1160 #= 1)) or 
	((Chr2_Chr3_Row1160 #= 693) and (Chr2_Chr3_Freq1160 #= 1)) or 
	((Chr2_Chr3_Row1160 #= 694) and (Chr2_Chr3_Freq1160 #= 1)) or 
	((Chr2_Chr3_Row1160 #= 695) and (Chr2_Chr3_Freq1160 #= 1)) or 
	((Chr2_Chr3_Row1160 #= 696) and (Chr2_Chr3_Freq1160 #= 1)) or 
	((Chr2_Chr3_Row1160 #= 697) and (Chr2_Chr3_Freq1160 #= 1)) or 
	((Chr2_Chr3_Row1160 #= 698) and (Chr2_Chr3_Freq1160 #= 1)) or 
	((Chr2_Chr3_Row1160 #= 699) and (Chr2_Chr3_Freq1160 #= 1)) or 
	((Chr2_Chr3_Row1160 #= 700) and (Chr2_Chr3_Freq1160 #= 1)) or 
	((Chr2_Chr3_Row1160 #= 702) and (Chr2_Chr3_Freq1160 #= 1)) or 
	((Chr2_Chr3_Row1160 #= 704) and (Chr2_Chr3_Freq1160 #= 1)) or 
	((Chr2_Chr3_Row1160 #= 705) and (Chr2_Chr3_Freq1160 #= 1)) or 
	((Chr2_Chr3_Row1160 #= 706) and (Chr2_Chr3_Freq1160 #= 1)) or 
	((Chr2_Chr3_Row1160 #= 707) and (Chr2_Chr3_Freq1160 #= 1)) or 
	((Chr2_Chr3_Row1160 #= 708) and (Chr2_Chr3_Freq1160 #= 1)) or 
	((Chr2_Chr3_Freq1160 #= 0) and (Chr2_Chr3_Row1160 #= 0)),

	((Chr2_Chr3_Row1161 #= 688) and (Chr2_Chr3_Freq1161 #= 1)) or 
	((Chr2_Chr3_Row1161 #= 689) and (Chr2_Chr3_Freq1161 #= 1)) or 
	((Chr2_Chr3_Row1161 #= 690) and (Chr2_Chr3_Freq1161 #= 1)) or 
	((Chr2_Chr3_Row1161 #= 691) and (Chr2_Chr3_Freq1161 #= 1)) or 
	((Chr2_Chr3_Row1161 #= 692) and (Chr2_Chr3_Freq1161 #= 1)) or 
	((Chr2_Chr3_Row1161 #= 693) and (Chr2_Chr3_Freq1161 #= 1)) or 
	((Chr2_Chr3_Row1161 #= 694) and (Chr2_Chr3_Freq1161 #= 1)) or 
	((Chr2_Chr3_Row1161 #= 699) and (Chr2_Chr3_Freq1161 #= 1)) or 
	((Chr2_Chr3_Row1161 #= 700) and (Chr2_Chr3_Freq1161 #= 1)) or 
	((Chr2_Chr3_Row1161 #= 702) and (Chr2_Chr3_Freq1161 #= 1)) or 
	((Chr2_Chr3_Row1161 #= 703) and (Chr2_Chr3_Freq1161 #= 1)) or 
	((Chr2_Chr3_Row1161 #= 705) and (Chr2_Chr3_Freq1161 #= 1)) or 
	((Chr2_Chr3_Row1161 #= 706) and (Chr2_Chr3_Freq1161 #= 1)) or 
	((Chr2_Chr3_Row1161 #= 707) and (Chr2_Chr3_Freq1161 #= 1)) or 
	((Chr2_Chr3_Row1161 #= 708) and (Chr2_Chr3_Freq1161 #= 1)) or 
	((Chr2_Chr3_Freq1161 #= 0) and (Chr2_Chr3_Row1161 #= 0)),

	((Chr2_Chr3_Row1162 #= 685) and (Chr2_Chr3_Freq1162 #= 1)) or 
	((Chr2_Chr3_Row1162 #= 686) and (Chr2_Chr3_Freq1162 #= 1)) or 
	((Chr2_Chr3_Row1162 #= 689) and (Chr2_Chr3_Freq1162 #= 1)) or 
	((Chr2_Chr3_Row1162 #= 690) and (Chr2_Chr3_Freq1162 #= 1)) or 
	((Chr2_Chr3_Row1162 #= 691) and (Chr2_Chr3_Freq1162 #= 1)) or 
	((Chr2_Chr3_Row1162 #= 692) and (Chr2_Chr3_Freq1162 #= 1)) or 
	((Chr2_Chr3_Row1162 #= 693) and (Chr2_Chr3_Freq1162 #= 1)) or 
	((Chr2_Chr3_Row1162 #= 694) and (Chr2_Chr3_Freq1162 #= 1)) or 
	((Chr2_Chr3_Row1162 #= 696) and (Chr2_Chr3_Freq1162 #= 1)) or 
	((Chr2_Chr3_Row1162 #= 697) and (Chr2_Chr3_Freq1162 #= 1)) or 
	((Chr2_Chr3_Row1162 #= 699) and (Chr2_Chr3_Freq1162 #= 1)) or 
	((Chr2_Chr3_Row1162 #= 700) and (Chr2_Chr3_Freq1162 #= 1)) or 
	((Chr2_Chr3_Row1162 #= 701) and (Chr2_Chr3_Freq1162 #= 1)) or 
	((Chr2_Chr3_Row1162 #= 702) and (Chr2_Chr3_Freq1162 #= 1)) or 
	((Chr2_Chr3_Row1162 #= 704) and (Chr2_Chr3_Freq1162 #= 1)) or 
	((Chr2_Chr3_Row1162 #= 705) and (Chr2_Chr3_Freq1162 #= 1)) or 
	((Chr2_Chr3_Row1162 #= 706) and (Chr2_Chr3_Freq1162 #= 1)) or 
	((Chr2_Chr3_Row1162 #= 707) and (Chr2_Chr3_Freq1162 #= 1)) or 
	((Chr2_Chr3_Freq1162 #= 0) and (Chr2_Chr3_Row1162 #= 0)),

	((Chr2_Chr3_Row1163 #= 690) and (Chr2_Chr3_Freq1163 #= 1)) or 
	((Chr2_Chr3_Row1163 #= 691) and (Chr2_Chr3_Freq1163 #= 1)) or 
	((Chr2_Chr3_Row1163 #= 694) and (Chr2_Chr3_Freq1163 #= 1)) or 
	((Chr2_Chr3_Row1163 #= 704) and (Chr2_Chr3_Freq1163 #= 1)) or 
	((Chr2_Chr3_Row1163 #= 707) and (Chr2_Chr3_Freq1163 #= 1)) or 
	((Chr2_Chr3_Row1163 #= 709) and (Chr2_Chr3_Freq1163 #= 1)) or 
	((Chr2_Chr3_Freq1163 #= 0) and (Chr2_Chr3_Row1163 #= 0)),

	((Chr2_Chr3_Row1164 #= 677) and (Chr2_Chr3_Freq1164 #= 1)) or 
	((Chr2_Chr3_Row1164 #= 689) and (Chr2_Chr3_Freq1164 #= 1)) or 
	((Chr2_Chr3_Row1164 #= 690) and (Chr2_Chr3_Freq1164 #= 1)) or 
	((Chr2_Chr3_Row1164 #= 699) and (Chr2_Chr3_Freq1164 #= 1)) or 
	((Chr2_Chr3_Row1164 #= 700) and (Chr2_Chr3_Freq1164 #= 1)) or 
	((Chr2_Chr3_Row1164 #= 702) and (Chr2_Chr3_Freq1164 #= 1)) or 
	((Chr2_Chr3_Row1164 #= 704) and (Chr2_Chr3_Freq1164 #= 1)) or 
	((Chr2_Chr3_Row1164 #= 707) and (Chr2_Chr3_Freq1164 #= 1)) or 
	((Chr2_Chr3_Row1164 #= 709) and (Chr2_Chr3_Freq1164 #= 1)) or 
	((Chr2_Chr3_Row1164 #= 755) and (Chr2_Chr3_Freq1164 #= 1)) or 
	((Chr2_Chr3_Freq1164 #= 0) and (Chr2_Chr3_Row1164 #= 0)),

	((Chr2_Chr3_Row1165 #= 685) and (Chr2_Chr3_Freq1165 #= 1)) or 
	((Chr2_Chr3_Row1165 #= 699) and (Chr2_Chr3_Freq1165 #= 1)) or 
	((Chr2_Chr3_Row1165 #= 704) and (Chr2_Chr3_Freq1165 #= 1)) or 
	((Chr2_Chr3_Row1165 #= 705) and (Chr2_Chr3_Freq1165 #= 1)) or 
	((Chr2_Chr3_Row1165 #= 706) and (Chr2_Chr3_Freq1165 #= 1)) or 
	((Chr2_Chr3_Row1165 #= 707) and (Chr2_Chr3_Freq1165 #= 1)) or 
	((Chr2_Chr3_Freq1165 #= 0) and (Chr2_Chr3_Row1165 #= 0)),

	((Chr2_Chr3_Row1166 #= 695) and (Chr2_Chr3_Freq1166 #= 1)) or 
	((Chr2_Chr3_Row1166 #= 705) and (Chr2_Chr3_Freq1166 #= 1)) or 
	((Chr2_Chr3_Row1166 #= 706) and (Chr2_Chr3_Freq1166 #= 1)) or 
	((Chr2_Chr3_Freq1166 #= 0) and (Chr2_Chr3_Row1166 #= 0)),

	((Chr2_Chr3_Row1167 #= 690) and (Chr2_Chr3_Freq1167 #= 1)) or 
	((Chr2_Chr3_Row1167 #= 704) and (Chr2_Chr3_Freq1167 #= 1)) or 
	((Chr2_Chr3_Freq1167 #= 0) and (Chr2_Chr3_Row1167 #= 0)),

	((Chr2_Chr3_Row1168 #= 691) and (Chr2_Chr3_Freq1168 #= 1)) or 
	((Chr2_Chr3_Freq1168 #= 0) and (Chr2_Chr3_Row1168 #= 0)),
	 
	((Chr2_Chr3_Row1171 #= 703) and (Chr2_Chr3_Freq1171 #= 1)) or 
	((Chr2_Chr3_Freq1171 #= 0) and (Chr2_Chr3_Row1171 #= 0)),

	((Chr2_Chr3_Row1257 #= 1000) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 1002) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 1003) and (Chr2_Chr3_Freq1257 #= 3)) or 
	((Chr2_Chr3_Row1257 #= 1004) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 1005) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 1006) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 1007) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 1008) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 562) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 563) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 568) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 571) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 582) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 584) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 591) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 592) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 593) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 594) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 596) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 599) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 601) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 603) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 604) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 605) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 665) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 681) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 690) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 691) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 694) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 695) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 704) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 823) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 877) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 898) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 908) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 917) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 922) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 923) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 924) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 925) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 926) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 927) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 951) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 964) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 967) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 968) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 969) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 971) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 974) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 975) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 976) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 978) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 979) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 980) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 981) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 982) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 983) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 984) and (Chr2_Chr3_Freq1257 #= 2)) or 
	((Chr2_Chr3_Row1257 #= 985) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 987) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 988) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 990) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 992) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 994) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 995) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 996) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 997) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 998) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Row1257 #= 999) and (Chr2_Chr3_Freq1257 #= 1)) or 
	((Chr2_Chr3_Freq1257 #= 0) and (Chr2_Chr3_Row1257 #= 0)),

	((Chr2_Chr3_Row1258 #= 1000) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 1003) and (Chr2_Chr3_Freq1258 #= 2)) or 
	((Chr2_Chr3_Row1258 #= 1004) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 1005) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 1006) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 1007) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 564) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 566) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 580) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 585) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 591) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 594) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 682) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 924) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 968) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 970) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 972) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 979) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 980) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 981) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 982) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 983) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 984) and (Chr2_Chr3_Freq1258 #= 2)) or 
	((Chr2_Chr3_Row1258 #= 985) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 986) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 987) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 992) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 994) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 998) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Row1258 #= 999) and (Chr2_Chr3_Freq1258 #= 1)) or 
	((Chr2_Chr3_Freq1258 #= 0) and (Chr2_Chr3_Row1258 #= 0)),

	% All of the values assumed by the Row<i> variables must be 
	% all different or zero to ensure each genomic bin is 
	% involved in only 1 interaction; multiple zeros are allowed 
	alldifferent_except(Non_Zero_Rows),
	atmost(77, Non_Zero_Rows, 0),

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
