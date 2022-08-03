function motifAnalyze(randomWalkSeq, s_randomWalkSeq, RBPContract, RBPDiscard)
%% motifAnalyze
% Find Possible Motifs in the input randomWalkSeq against shuffled s_randomWalkSeq

% Inputs:
%   - randomWalkSeq       :  Input RW Seq. on Real Indices of Cell Types
%   - s_randomWalkSeq     :  Input RW Seq. on Shuffled Indices of Cell Types
%   - OContract           : Input Seq. is 'O' Contracted if True

% === compute length and number of sequences in each set
L1 = size(s_randomWalkSeq, 2);
L2 = size(randomWalkSeq, 2);

N1 = size(s_randomWalkSeq, 1);
N2 = size(randomWalkSeq, 1);

chars=['A', 'B','C', 'D','E', 'F','G', 'H','I', 'J','K', 'L','M', 'N', 'O'];
% CellNumbers=[133,141,253,121,344,150,254,257,235,173,314,157,23,51,1293];

% Convert Cell Labels to Chars
s2=chars(randomWalkSeq);
s2=s2.';
s2=s2(:).';



% === add separator between sequences
seq2 = seqinsertgaps(s2, 1:L2:(L2*N2)+N2, 1);
% Motif Word Length
if (RBPDiscard)
    W=3;
else
    if (RBPContract)
        W =8;
    else
        W=25;
    end
end

seq1t=[];
numShuffle=size(s_randomWalkSeq, 3);
motifs10Cell=cell(10* size(s_randomWalkSeq, 3), 1);
margins10Matrix=zeros(10*size(s_randomWalkSeq, 3),1);

for it=1:size(s_randomWalkSeq, 3)
    
    s1=chars(s_randomWalkSeq(:,:, it));
    s1=s1.';
    s1=s1(:).';
    seq1 = seqinsertgaps(s1, 1:L1:(L1*N1)+N1, 1);
    [words, freqDiff] = findOverrepresentedWordsSam(seq2, seq1,W);
    motifs10Cell((it-1)*10+1:it*10, :)=words(1:10);
    margins10Matrix((it-1)*10+1:it*10)=freqDiff(1:10);
    seq1t=[seq1t, seq1];
    
end

[motif,ia,ic]=unique(motifs10Cell);
margin=zeros(length(motif), 1);
for ii =1:length(motif)
    margin(ii)=sum(margins10Matrix(ic==ii));
end

[margin, marginSortedIndex]=sort(margin, 'descend');
motif=motif(marginSortedIndex);


% === consider the top 10 motifs
motif = words(1:10)
margin = freqDiff(1:10)

CPMotif = motif{1};
CPMargin = margin(1);

%% Assessing the Statistical Significance of Margins


% === generate random sequences 
ctrlIter=numShuffle;
ctrlMargin=zeros(1,ctrlIter);
ctrlMotif=cell(ctrlIter,1);

seq1Length=length(seq1t)/numShuffle;
rseqIdx=(1:seq1Length);
rfreqDiffMaxT=0;

for iter=1:ctrlIter


    rfreqDiffMax=0;
    for itt1=1:10
        FrameIdx=floor(100*rand(1, 2))+1;
        s1Idx=(FrameIdx(1)-1)*seq1Length+rseqIdx;
        s2Idx=(FrameIdx(2)-1)*seq1Length+rseqIdx;
        rseq1=seq1t(s1Idx);
        rseq2=seq1t(s2Idx);
        
        % === compute margins for control set
        [rwords, rfreqDiff] = findOverrepresentedWordsSam(rseq2, rseq1,W);
        if rfreqDiff(1)>rfreqDiffMax
            rfreqDiffMax=rfreqDiff(1);
            rwordMax=rwords(1);
        end
        if rfreqDiffMax>rfreqDiffMaxT
            rfreqDiffMaxT=rfreqDiffMax;
            rseq2M=rseq2;
        end
    end
    ctrlMargin(1, iter)=rfreqDiffMax;
    ctrlMotif(iter)=rwordMax;
    
end
% === estimate parameters of distribution
nCtrl = length(ctrlMargin);
buckets = ceil(nCtrl/10);
parmhat = gevfit(ctrlMargin);
k = parmhat(1);     % shape parameter
sigma = parmhat(2); % scale parameter
mu = parmhat(3);    % location parameter

% === compute probability density function
x = linspace(min(ctrlMargin), max([ctrlMargin CPMargin]));
y = gevpdf(x, k, sigma, mu);

% === scale probability density function
[~, c] = hist(ctrlMargin,buckets);
binWidth = c(2) - c(1);
scaleFactor = nCtrl * binWidth;

% === overlay
figure()
hold on;
hist(ctrlMargin, buckets);
h = findobj(gca,'Type','patch');
h.FaceColor = [.9 .9 .9];
plot(x, scaleFactor * y, 'r');
stem(CPMargin, 1, 'b');
xlabel('Margin');
ylabel('Number of sequences');
legend('Ctrl Margins', 'EVD pdf', 'RWM Margin');

%% Determining the Cell Position Motif Presence Among RW Seqs.
CPCount = zeros(2,1);

% === determine positions where CP motif occurs
loc1 = strfind(rseq2M, CPMotif);
if ~isempty(loc1)
loc1i=(loc1>(0:1001:length(rseq2M)-1001).' & loc1<(1001:1001:length(rseq2M)).');
CPCount(1)=sum(sum(loc1i, 2)>0);
else
    CPCount(1)=0;
end
loc2 = strfind(seq2, CPMotif);
loc2i=(loc2>(0:1001:length(rseq2M)-1001).' & loc2<(1001:1001:length(rseq2M)).');
CPCount(2)=sum(sum(loc2i, 2)>0);

% === find proportions of RWs with CP Motif
NumRWs = [N1; N2] ;
CPProp = CPCount ./ NumRWs;

% === plot
figure()
bar(CPProp, 0.5);
ylabel('Proportion of RWs containing Top Motif');
xlabel('RW Cell Set');
title('Presence of Top Motif');
ylim([0 1])
ax = gca;
ax.XTickLabel = {'Shuffled', 'Real'};
grid on