function [fastq_seq_new] = trimReads(fastq_seq, trimIdx,trimLongReadsIdx)

% This function will trim reads using the index from 'trimReadIndex.m'. It
% returns a new fastq file with the trimmed sequences.
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

