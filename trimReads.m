function [fastq_seq_new] = trimReads(fastq_seq, trimIdx,trimLongReadsIdx)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
header = {fastq_seq(:).Header}';
seq = {fastq_seq(:).Sequence}';
qual = {fastq_seq(:).Quality}';

seq = seq(trimLongReadsIdx);
trimIdx = trimIdx(trimLongReadsIdx);
qual = qual(trimLongReadsIdx);

headerNew = header(trimLongReadsIdx);
seqNew = cellfun(@(x,y) x(y), seq, trimIdx,'UniformOutput',false);
qualNew = cellfun(@(x,y) x(y), qual, trimIdx, 'UniformOutput',false);

fastq_seq_new = struct('Header',headerNew,'Sequence',seqNew,'Quality',qualNew);


end

