[fileName, pathName] = uigetfile('E:\Zon Lab\Sequencing\MGH SEQUENCING\*.fastq');


% Read in fastq file
fastq_seq = fastqread(fullfile(pathName,fileName));

%% Separate paired end sequencing into 2 files

% Vectorized
idx = ~cellfun('isempty',strfind({fastq_seq.Header},'1:N'));
pairedEndTest1 = fastq_seq(idx);
pairedEndTest2 = fastq_seq(~idx);

if ~exist(fullfile(pathName,'Parsed Output'),'dir')
    mkdir(pathName,'Parsed Output');
end

if ~exist(fullfile(pathName,'Parsed Output','Paired1'),'dir')
    mkdir(fullfile(pathName,'Parsed Output'),'Paired1');
end

if ~exist(fullfile(pathName,'Parsed Output','Paired2'),'dir')
    mkdir(fullfile(pathName,'Parsed Output'),'Paired2');
end

fastqwrite(fullfile(pathName,'Parsed Output','Paired1',[fileName,'(1).fastq']),pairedEndTest1);
fastqwrite(fullfile(pathName,'Parsed Output','Paired2',[fileName,'(2).fastq']),pairedEndTest2);
