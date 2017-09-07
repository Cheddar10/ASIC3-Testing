function errors = biterrors(bits, bitseq, offset)
%computes errors between bits and target sequency given a phase offset

nseq = length(bitseq);
dubseq = repmat(bitseq,2,1);
bitrepeat = repmat(dubseq(offset:offset+nseq-1),ceil(length(bits)/nseq),1);
bitrepeat = bitrepeat(1:length(bitseq));
bitrepeat = repmat(bitrepeat,ceil(length(bits)/length(bitseq)),1);
bitrepeat = bitrepeat(1:length(bits));
errors = bitrepeat == bits;
end