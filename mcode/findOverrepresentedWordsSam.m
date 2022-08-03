function [nmersSorted, freqDiffSorted] = findOverrepresentedWordsSam(seq, seq0, W)
%% findOverrepresentedWordsSam
% findOverrepresentedWords Function of the Matlab with Slight
% Modifications

% Inputs:
%   - seq      :  Input RW Seq. on Real Indices of Cell Types
%   - seq0     :  Input RW Seq. on Shuffled Indices of Cell Types
%   - W        :  Length of the Supposed Motif

% Outputs:
%   - nmersSorted        :  Sorted Frequent nmers
%   - freqDiffSorted     :  Frq. Difference of the Real and Shuffled nmers


% Copyright 2007 The MathWorks, Inc.

%=== find and count words of length W
numCandidCells=100; % Just first 100 Dominant nmers are considerred 

[nmers0Cell,nmers0Count] = nmercount(seq0, W,1);
nmers0Cell=num2cell(nmers0Cell(1:numCandidCells, :),2);
nmers0Count=nmers0Count(1:numCandidCells);

[nmersCell,nmersCount] = nmercount(seq, W,1);
nmersCell=num2cell(nmersCell(1:numCandidCells, :),2);
nmersCount=nmersCount(1:numCandidCells);

%=== compute frequency of words 
f = nmersCount/(length(seq) - W + 1);
f0 = nmers0Count/(length(seq0) - W + 1);


%=== determine words common to both set 
[nmersInt, i1, i2] = intersect(nmersCell,nmers0Cell, 'stable');
freqDiffInt = (f(i1) - f0(i2))';

%=== determine words specific to one set only
[~, i3, i4] = setxor(nmersCell,nmers0Cell);
c0 = nmersCell(i3);
d0 = nmers0Cell(i4);
nmersXOr = [c0; d0]; 
freqDiffXOr = [f(i3) ;-f0(i4)];

%=== define all words and their difference in frequency (margin)
nmersAll = [nmersInt; nmersXOr];
freqDiff = [freqDiffInt(:); freqDiffXOr(:)];

%=== sort according to descending difference in frequency
[freqDiffSorted, freqDiffSortedIndex] = sort(freqDiff, 'descend'); 
nmersSorted = nmersAll(freqDiffSortedIndex);

