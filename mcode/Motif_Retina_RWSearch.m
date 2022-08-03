clc
close all
clearvars

%% Hyperparameters 
% Hyperparameters for Creating Graph 
shuffleNumber=100;          % Number of shuffled RWs
edgeWeighUniform=true;      % Edge Weighs are considerred unform if set true
RBPDiscard=false;           % Discard RBP Cells from downward analysis.

% Hyperparameters for Random Walk
lazyFactor=0.15;        % Lazy RW Factor, Could be Set to Zero
randomWalkLen=100000;   % Length of the RW, Should be Long Enough to Make Sure Mixing TIme has Passed
numSeqs=300;            % Number of independent RW Seqs
seqLength=1000;         % Seq Length from the End of the RW Seq that is Chosen for Downward Analysis
RBPContract=true;       % Contract Consecutive 'O' Sequences(15 SUbCellType)

%% Generate Graph


% Loading input Retina Data
filename='..\data\Retina1.csv';
T = readtable(filename);
sectNumber=T.SectionNumber;
cellSubtypeVec=T.Subtype(sectNumber==1);
cellSubtypeVec0=cellSubtypeVec;
xcoords=T.Remapped_X(sectNumber==1);
ycoords=T.Remapped_Y(sectNumber==1);
if RBPDiscard
    xcoords=xcoords(~(cellSubtypeVec==15));
    ycoords=ycoords(~(cellSubtypeVec==15));
    cellSubtypeVec=cellSubtypeVec(~(cellSubtypeVec==15));
end

zcoords=[xcoords, ycoords];



Dmat = pdist2(zcoords,zcoords, 'euclidean'); % Distance Matrix, could be used to weight the edges

% Delauney Triangle Computation
Delauney_Triangle=delaunayTriangulation(zcoords);

% Plotting Output Resulted Graph
figure
triplot(Delauney_Triangle)
title('Delauney Triangle of the Input Coordinates Data')

% Extract Edge Data from Connected Triangle Information

connectTriangle=Delauney_Triangle.ConnectivityList;
edges=[connectTriangle(:,1:2); connectTriangle(:,2:3); connectTriangle(:,[1,3])];

% Discard Repeated Edges
edges=sort(edges, 2);
[edges_unique,ia,ic]=unique(edges, 'rows', 'stable');

% Computing Weights with W=exp(-d^2/dave^2), with d the distance between
% respective vertices, and dave the mean distance between all vertices
% [David Harel and Yehuda Koren , On Clustering Using Random Walks]

edge_dist_idx=(edges_unique(:,2)-1)*length(Dmat)+edges_unique(:,1);
edge_dist=Dmat(edge_dist_idx);
Weigths=exp(-edge_dist.^2/mean(edge_dist)^2);

% Discarding Edges with Weights less than 1e-6
edges_unique=edges_unique(Weigths>1e-6, :);
Weigths=Weigths(Weigths>1e-6);

% Weights are Set to All 1's if Unform Weithing is Assumed
if (edgeWeighUniform)
    Weigths=ones(length(Weigths),1);
end
% Define graph object
G=graph(edges_unique(:,1),edges_unique(:,2), Weigths);

% Plot the Output Pruned Graph
figure
h=plot(G,'XData',xcoords,'YData',ycoords);
title('Output Pruned Graph')
grid on
hold on

vertexVec=(1:length(cellSubtypeVec));
ncolor=[1, 0.1, 0.7];
% Highligt RBP cell types 
if ~RBPDiscard
    highlight(h,vertexVec(cellSubtypeVec==15),'NodeColor',ncolor,'EdgeColor',ncolor)
end

%      
% Obtain Shuffled Cell Subtype Vectors
s_cellSubtypeVec=zeros(length(cellSubtypeVec), shuffleNumber);
for ii =1:shuffleNumber 
    s_cellSubtypeVec(:,ii)=cellSubtypeVec(randperm(length(cellSubtypeVec)));
end

%% Random Walk

% Initialize Seq Vectors
randomWalkSeq=zeros(numSeqs, seqLength);                    % RW on Real Data
s_randomWalkSeq=zeros(numSeqs, seqLength, shuffleNumber);   % RW on Shuffled Data
selectVerticeVecMat=zeros(numSeqs,seqLength);               % Matrix of Selected Edges Used to Verify that RWs cover the Entire Retina Length



for seqID = 1: numSeqs
    
    % Initial Vertice of the RW
    initialVertice=edges_unique(randi(length(edges_unique)));
    selectVerticeVec=zeros(randomWalkLen+1, 1);
    selectVerticeVec(1)=initialVertice;
    
    % F^-1(u) Has the Distribution of f(x), Where F is CDF of the Edges' Weights at
    % Each Respective Nodes, f(x) is the pfd of Weights of Edges,
    % and u is a Random Variable with Uniform Distribution, randVec is
    % Later Used to Decide About the Next Node in RW Seq. 
    randVec=rand(randomWalkLen, 1);
    for randomWalki=1:randomWalkLen
        
        % Edges Incident to Selected Vertice
        startVertex=(edges_unique==initialVertice);
        startVertex2=(startVertex(:, 1) | startVertex(:, 2));
        
        % Vertice Position in Possible Edges
        edgeFlag=startVertex(startVertex2,:);
        % Candidate Edges
        candidEdges=edges_unique(startVertex2,:);
        % Adding the Self Loop
        candidEdgese=[candidEdges;[candidEdges(edgeFlag(1,:)),candidEdges(edgeFlag(1,:))]];
        
        % Select the Next Vertice with RW on Weighted Edges
        lazyP=lazyFactor*sum(Weigths(startVertex2));
        Pvec=[Weigths(startVertex2)*(1-lazyFactor);lazyP];
        randP=sum(Pvec)*randVec(randomWalki);
        selectIdx=sum((randP-cumsum(Pvec))>0)+1;
        selectEdge=candidEdgese(selectIdx, :);
        % Add Another Edge to Account for Lazy RW
        edgeFlage=[edgeFlag;0,1];
        
        % The Next Vertice is the Other End of the Selected Edge
        initialVertice=selectEdge(~edgeFlage(selectIdx, :));
        
        % Record Obtained Vertices
        selectVerticeVec(randomWalki+1)=initialVertice;

    end
    % Contract 'O' Long Seqs if Set True
    if (RBPContract && ~RBPDiscard) 
         randomWalkSeq(seqID, :)= seqLengthChoose(cellSubtypeVec(selectVerticeVec), seqLength) ;        
    else
         randomWalkSeq(seqID, :)=cellSubtypeVec(selectVerticeVec(end-seqLength+1:end));    
    end
        
        selectVerticeVecMat(seqID, :)=selectVerticeVec(end-seqLength+1:end);

    for shufflei=1:shuffleNumber
        % The Shuffled cellSubtype
        s_cellSubtypeVecIter=s_cellSubtypeVec(:,shufflei);
        % Contract 'O' Long Seqs if Set True
        if (RBPContract && ~RBPDiscard) 
            s_randomWalkSeq(seqID, :, shufflei)= seqLengthChoose(s_cellSubtypeVecIter(selectVerticeVec), seqLength) ; 
        else
            s_randomWalkSeq(seqID, :, shufflei)= s_cellSubtypeVecIter(selectVerticeVec(end-seqLength+1:end));
        end

    end

end

% Highlight the Selected RW Seqs in the Graph Data
for seqID = 1: numSeqs
    ncolor=rand(1,3);
     highlight(h,selectVerticeVecMat(seqID, :),'NodeColor',ncolor,'EdgeColor',ncolor)
end

%% Identify Motifs
motifAnalyze(randomWalkSeq, s_randomWalkSeq, RBPContract,  RBPDiscard)

