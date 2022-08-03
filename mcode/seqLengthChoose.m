function  randomWalkSeqo= seqLengthChoose(randomWalkSeqi, seqLength) 
%% seqLengthChoose
% Contract '15' or 'O' sequences in the Input Stream

% Inputs:
%   - randomWalkSeqi:  Input RW Seq.
%   - seqLength     :  Required Length of the Output Seq.
% Outputs:
%   -randomWalkSeqo : Output 'O' Contracted Seq. with seqLength Length

randomWalkSeqo=zeros(1,seqLength);
    io=length(randomWalkSeqi);
    ioc=1;
    while ioc<seqLength+1
        if ~(all(randomWalkSeqi(io-1:io)==[15, 15]))
            randomWalkSeqo(ioc)=randomWalkSeqi(io);
            ioc=ioc+1;
            
        end
        io=io-1;
    end
         