function trimIdx = trimReadIndex(scores,threshold)
%%==== This function will find the best trimmed sequence containing reads
%%==== with quality score > threshold. Scores is a cell array of quality
%%==== scores. Threshold is the quality threshold (30 is
%%==== usually good).

A = scores >= threshold;

B = SplitVec(A);
sizeOfVecs = cellfun('length',B);
[maxValue, maxIndex] = max(sizeOfVecs);

maxTest = find(sizeOfVecs == maxValue);

firstIndex = zeros(numel(B),1);
lastIndex = zeros(numel(B),1);
for ii = 1:numel(B)
    if ii == 1
        firstIndex(ii) = 1;
        lastIndex(ii) = sizeOfVecs(ii);
        count = sizeOfVecs(ii);
    else
        firstIndex(ii) = firstIndex(ii-1)+sizeOfVecs(ii-1);
        lastIndex(ii) = count+ sizeOfVecs(ii);
        count = count + sizeOfVecs(ii);
    end
end

if numel(maxTest) == 1
    trimIdx = firstIndex(maxIndex) : lastIndex(maxIndex);
else
    trimIdxTest = cell(numel(maxTest),1);
    meanQual = zeros(numel(maxTest),1);
    for kk = 1:numel(maxTest)
        trimIdxTest{kk} = firstIndex(maxTest(kk)) : lastIndex(maxTest(kk));
        meanQual(kk) = mean(scores(trimIdxTest{kk}));
        
    end
    [~, highestQualIndex] = max(meanQual);
    trimIdx = firstIndex(maxTest(highestQualIndex)) : lastIndex(maxTest(highestQualIndex));
end







% %%Old way
% x1 = find(A == 0);
%
% numOfZeros = numel(x1);
% if isempty(x1) == 1
%     seqIdx = NaN;
% elseif numOfZeros > 1
%     x2 = diff(x1);
%     [m, idx] = max(x2);
%     %Test here to see if there are multiple maxes in the matrix. Then
%     %compute mean quality score of each stretch and pick one that is best.
%     x3 = x1(idx)+1;
%     seqIdx = x3:(x3+m-2);
% elseif numOfZeros == 1
%     if x1 <= size(scores,2)/2
%         seqIdx = (x1+1):size(scores,2);
%     else
%         seqIdx = 1:(x1-1);
%     end
%
% end

end