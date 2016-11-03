%TO DO: When trimming reads, consider the case where there is one low
%quality read in the exact middle of the sequence. Need to choose what side
%of the read to use (may want to use whichever end has higher overall
%quality).
%NOTE: This can also happen if you have trimmed reads of equal length.
%Should probably compute average quality to choose which one.

[fileName, pathName] = uigetfile('E:\Zon Lab\Sequencing\MGH SEQUENCING\*.fastq');


% Read in fastq file
fastq_seq = fastqread(fullfile(pathName,fileName));

%Convert quality scores
quality = {fastq_seq(:).Quality}';
quality = cellfun(@(x) double(x)-33,quality,'UniformOutput',false);

trimIndex = cellfun(@(x) trimReads(x,30),quality,'UniformOutput',false);


%% Separate paired end sequencing into 2 files

% % Vectorized
% idx = ~cellfun('isempty',strfind({fastq_seq.Header},'1:N'));
% pairedEndTest1 = fastq_seq(idx);
% pairedEndTest2 = fastq_seq(~idx);
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



function seqIdx = trimReads(scores,threshold)

A = scores>=threshold;
x1 = find(A == 0);

numOfZeros = size(x1,1);
if numOfZeros == 0
    seqIdx = NaN;
elseif numOfZeros > 1
    x2 = diff(x1);
    [m, idx] = max(x2);
    x3 = x1(idx)+1;
    seqIdx = x3:(x3+m-2);
elseif numOfZeros == 1
    if x1 <= size(scores,1)/2
        seqIdx = (x1+1):size(scores,1);
    else
        seqIdx = 1:(x1-1);
    end
    
end

end