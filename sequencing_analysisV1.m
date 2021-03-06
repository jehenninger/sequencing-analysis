%TO DO: When two stretches of trimmed reads have same length, I choose the
%one with higher average quality. May want to change this if one of them
%covers the CRISPR site and still has avg > 30.

%TO DO: Maybe provide an output % of reads that are length greater than 50.

%TO DO: Maybe want to change some of these cellfun's to for loops in order
%to add progress bars. Probably won't change performance much.

%TO DO: Talk to Song about the genome file. It seems like my reads aren't
%lining up with the Exons that he has provided.





%% ALL MOVED TO crisprSeqAnalysis.m
% % [fileName, pathName] = uigetfile('E:\Zon Lab\Sequencing\MGH SEQUENCING\*.fastq');
% % [fileName, pathName] = uigetfile('D:\Jon-HDD\Google Drive\Sequencing\MGH CRISPR SEQUENCING\*.fastq');
% [fileName, pathName] = uigetfile('/Users/jonathanhenninger2/Desktop/seq_test/*.fastq');
% 
% [~,sampleName] = fileparts(fileName);
% 
% %%  Read in fastq file
% fastq_seq = fastqread(fullfile(pathName,fileName));
% 
% %% Split paired end reads assuming odds and evens
% paired1 = fastq_seq(1:2:end); 
% paired2 = fastq_seq(2:2:end);
% 
% %% Convert quality scores
% quality1 = {paired1(:).Quality}';
% quality1 = cellfun(@(x) double(x)-33,quality1,'UniformOutput',false);
% 
% quality2 = {paired2(:).Quality}';
% quality2 = cellfun(@(x) double(x)-33,quality2,'UniformOutput',false);
% 
% %% Trim reads so that quality for every sequence is >= 30.
% trimIdx1 = cellfun(@(x) trimReadIndex(x,30),quality1,'UniformOutput',false);
% trimReadLengths1 = cellfun('length',trimIdx1);
% trimLongReadsIdx1 = find(trimReadLengths1>=50);
% 
% [paired1New] = trimReads(paired1, trimIdx1,trimLongReadsIdx1);
% 
% trimIdx2 = cellfun(@(x) trimReadIndex(x,30),quality2,'UniformOutput',false);
% trimReadLengths2 = cellfun('length',trimIdx2);
% trimLongReadsIdx2 = find(trimReadLengths2>=50);
% 
% [paired2New] = trimReads(paired2, trimIdx2,trimLongReadsIdx2);
% 
% 
% 
% %% Separate paired end sequencing into 2 files
% % 
% % % Vectorized
% % idx = ~cellfun('isempty',strfind({fastq_seq_new.Header},'1:N'));
% % pairedEndTest1 = fastq_seq_new(idx);
% % pairedEndTest2 = fastq_seq_new(~idx);
% 
% % % Take reverse complement of the 2nd paired end reads
% % % NEEDED IF YOU WANT TO JOIN SEQUENCES INTO ONE READ
% % p2 = {paired2New(:).Sequence}';
% % p3 = cellfun(@(x) seqrcomplement(x), p2,'UniformOutput',false);
% % [paired2New.Sequence] = p3{:};
% % clear p2 p3
% 
% % %% Find paired reads
% % pairedIdx = findPairedReads({pairedEndTest1(:).Header},{pairedEndTest2(:).Header});
% 
% 
% %% Merge matched paired reads
% [~, pairs1Idx, pairs2Idx] = intersect(trimLongReadsIdx1, trimLongReadsIdx2);
% 
% [~, unpaired1Idx] = setdiff(trimLongReadsIdx1,trimLongReadsIdx2,'stable');
% [~, unpaired2Idx] = setdiff(trimLongReadsIdx2,trimLongReadsIdx1,'stable');
% 
% paired1NewFilt = paired1New(pairs1Idx);
% unpaired1NewFilt = paired1New(unpaired1Idx);
% 
% paired2NewFilt = paired2New(pairs2Idx);
% unpaired2NewFilt = paired2New(unpaired2Idx);
% 
% %%If we want to join the actual sequences together
% % pairedJoined = cellfun(@(x,y) joinseq(x,y),{paired1NewFilt(:).Sequence}',{paired2NewFilt(:).Sequence}');
% 
% % Join paired sequences into one fastq
% numOfCellsNeeded = numel(pairs1Idx) + numel(pairs2Idx);
% paired = struct('Header',cell(numOfCellsNeeded,1),...
%     'Sequence',numOfCellsNeeded,...
%     'Quality',numOfCellsNeeded);
% paired(2:2:end) = paired2NewFilt(:);
% paired(1:2:end) = paired1NewFilt(:);
% paired = paired';
% 
% % Join unpaired sequences into one fastq
% numOfCellsNeeded = numel(unpaired1NewFilt) + numel(unpaired2NewFilt);
% unpaired = struct('Header',cell(numOfCellsNeeded,1),...
%     'Sequence',numOfCellsNeeded,...
%     'Quality',numOfCellsNeeded);
% unpaired(1:numel(unpaired1NewFilt)) = unpaired1NewFilt(:);
% unpaired((numel(unpaired1NewFilt)+1):end) = unpaired2NewFilt(:);
% unpaired = unpaired';
% 
% %% Write fastq files
% disp('Writing fastq files...');
% if ~exist(fullfile(pathName,'Fastq_trimmed'),'dir')
%     mkdir(pathName,'Fastq_trimmed');
% end
% 
% % fastqwrite(fullfile(pathName, 'Fastq trimmed',[sampleName,'_paired.fastq']),paired);
% % fastqwrite(fullfile(pathName, 'Fastq trimmed',[sampleName,'_unpaired.fastq']),unpaired);
% pair1File = fullfile(pathName,'Fastq_trimmed',[sampleName,'_paired_1.fastq']);
% pair2File = fullfile(pathName,'Fastq_trimmed',[sampleName,'_paired_2.fastq']);
% unpairedFile = fullfile(pathName,'Fastq_trimmed',[sampleName,'_unpaired.fastq']);
% 
% fastqwrite(pair1File,paired1NewFilt(:)');
% fastqwrite(pair2File,paired2NewFilt(:)');
% 
% fastqwrite(unpairedFile,unpaired);
% 
% % %% Align reads with Bowtie
% % 
% % [genomeFile, genomePath] = uigetfile(...
% %     '/Users/jonathanhenninger2/Desktop/MGH CRISPR SEQUENCING/custom_genome/*.fa');
% % genomeFile = genomeFile(1:(end-3));
% % 
% % if ~exist(fullfile(pathName,'BAM Output'),'dir')
% %     mkdir(pathName,'BAM Output');
% % end
% % 
% % bowTieFilePaired = fullfile(pathName,'BAM Output',[sampleName,'_paired.bam']);
% % bowtie(fullfile(genomePath,genomeFile),...
% %     {pair1File,pair2File},...
% %     bowTieFilePaired,...
% %     'Paired',true,...
% %     'BowtieOptions','-l 15',...
% %     'BowtieOptions','-e 10000');
% % 
% % bowTieFileUnPaired = fullfile(pathName,'BAM Output',[sampleName,'_unpaired.bam']);
% % bowtie(fullfile(genomePath,genomeFile),...
% %     {unpairedFile},...
% %     bowTieFileUnPaired,...
% %     'BowtieOptions','-l 15',...
% %     'BowtieOptions','-e 10000');
% % 
% % bmPaired = BioMap(bowTieFilePaired);
% % bmUnpaired = BioMap(bowTieFileUnPaired);
% 
% % 
% % if ~exist(fullfile(pathName,'Parsed Output'),'dir')
% %     mkdir(pathName,'Parsed Output');
% % end
% % 
% % if ~exist(fullfile(pathName,'Parsed Output','Paired1'),'dir')
% %     mkdir(fullfile(pathName,'Parsed Output'),'Paired1');
% % end
% % 
% % if ~exist(fullfile(pathName,'Parsed Output','Paired2'),'dir')
% %     mkdir(fullfile(pathName,'Parsed Output'),'Paired2');
% % end
% % 
% % fastqwrite(fullfile(pathName,'Parsed Output','Paired1',[fileName,'(1).fastq']),pairedEndTest1);
% % fastqwrite(fullfile(pathName,'Parsed Output','Paired2',[fileName,'(2).fastq']),pairedEndTest2);
% 
% 
% 
