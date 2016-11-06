%TO DO: When two stretches of trimmed reads have same length, I choose the
%one with higher average quality. May want to change this if one of them
%covers the CRISPR site and still has avg > 30.

%TO DO: Maybe provide an output % of reads that are length greater than 50.

% [fileName, pathName] = uigetfile('E:\Zon Lab\Sequencing\MGH SEQUENCING\*.fastq');
[fileName, pathName] = uigetfile('/Users/jonathanhenninger2/Desktop/MGH CRISPR Sequencing/*.fastq');


%%  Read in fastq file
fastq_seq = fastqread(fullfile(pathName,fileName));

%% Convert quality scores
quality = {fastq_seq(:).Quality}';
quality = cellfun(@(x) double(x)-33,quality,'UniformOutput',false);

%% Trim reads so that quality for every sequence is >= 30.
trimIdx = cellfun(@(x) trimReadIndex(x,30),quality,'UniformOutput',false);
trimReadLengths = cellfun('length',trimIdx);
trimLongReadsIdx = find(trimReadLengths>=50);

[fastq_seq_new] = trimReads(fastq_seq, trimIdx,trimLongReadsIdx);



%% Separate paired end sequencing into 2 files

% Vectorized
idx = ~cellfun('isempty',strfind({fastq_seq_new.Header},'1:N'));
pairedEndTest1 = fastq_seq_new(idx);
pairedEndTest2 = fastq_seq_new(~idx);

% Take reverse complement of the 2nd paired end reads
p2 = {pairedEndTest2(:).Sequence}';
p3 = cellfun(@(x) seqrcomplement(x), p2,'UniformOutput',false);
[pairedEndTest2.Sequence] = p3{:};

%% Find paired reads
[pairedIdx, debug_idx] = findPairedReads({pairedEndTest1(:).Header},{pairedEndTest2(:).Header});


% 
% if ~exist(fullfile(pathName,'Parsed Output'),'dir')
%     mkdir(pathName,'Parsed Output');
% end
% 
% if ~exist(fullfile(pathName,'Parsed Output','Paired1'),'dir')
%     mkdir(fullfile(pathName,'Parsed Output'),'Paired1');
% end
% 
% if ~exist(fullfile(pathName,'Parsed Output','Paired2'),'dir')
%     mkdir(fullfile(pathName,'Parsed Output'),'Paired2');
% end
% 
% fastqwrite(fullfile(pathName,'Parsed Output','Paired1',[fileName,'(1).fastq']),pairedEndTest1);
% fastqwrite(fullfile(pathName,'Parsed Output','Paired2',[fileName,'(2).fastq']),pairedEndTest2);



