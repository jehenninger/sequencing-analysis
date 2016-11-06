function [string_idx,debug_idx] = findPairedReads(p1Header,p2Header)

% p1Header and p2Header are the Header fields from the separated paired end
% reads

%TO DO: Not every read is going to have a match. How to know which ones do?
p1Header = p1Header';
p2Header = p2Header';

sizeTest = numel(p1Header)-numel(p2Header);

p1Strings = cellfun(@(x) x(1:(end-8)),p1Header,...
    'UniformOutput',false);
p2Strings = cellfun(@(x) x(1:(end-8)),p2Header,...
    'UniformOutput',false);

if sizeTest >= 0
    %     padding = cell(sizeTest,1);
    %     padding(:) = {'a'};
    %     p2Strings = [p2Strings; padding];
    hWait = waitbar(0,'Finding paired read matches...');
    steps = numel(p2Strings);
    
    %     string_idx = zeros(numel(p2Strings),1);
    string_idx = cell(numel(p2Strings),1);
%     debug_idx = zeros(numel(p2Strings),1);
    for ii = 1:numel(p2Strings)
        %         try
        string_idx{ii} = find(strcmp(p2Strings{ii},p1Strings)== 1);
        if isempty(string_idx{ii})
            string_idx{ii} = -1;
        end
        %         catch
        %             debug_idx(ii) = ii;
        %         end
        waitbar(ii/steps);
    end
    close(hWait);
    
    string_idx = cell2mat(string_idx);
    
    %         string_idx = cellfun(@(x,y) find(strcmp(x,y)),p1Strings,p2Strings,...
    %             'UniformOutput',false);
else
    %     p1Strings = cellfun(@(x) x(1:(end-8)),p1Header,...
    %         'UniformOutput',false);
    %     string_idx = cellfun(@(x,y) find(cell2mat(strfind(x,y)) == 1),p2Header,p1Strings,...
    %         'UniformOutput',false);
    
end


end

