addpath('../functions/'); 
timeNum = [1576, 1604, 1624]; 

for i = 1:length(timeNum)
    fprintf('Time Slice number %d:\n', timeNum(i)); 
    
    load(strcat('labeled_', num2str(timeNum(i)), '.mat'));
    load(strcat('zhang_', num2str(timeNum(i)), '.mat'));
    idx1 = simIdx(bwZ, bw); 
    fprintf('Similarity Index of the method proposed by Zhang is %1.4f\n', idx1); 
    
    load(strcat('proposed_', num2str(timeNum(i)), '_noColor.mat'));
    idx2 = simIdx(bwP, bw); 
    fprintf('Similarity Index of the proposed method w.o. Color is %1.4f\n', idx2); 
    
    load(strcat('proposed_', num2str(timeNum(i)), '_Color.mat'));
    idx3 = simIdx(bwP, bw); 
    fprintf('Similarity Index of the proposed method w. Color is %1.4f\n', idx3); 
    fprintf('\n\n'); 
end



% 
% load labeled_1576.mat % bw
% load zhang_1576.mat % bwZ
% addpath('..\functions');
% 
% idxZhang = simIdx(bwZ, bw)
% 
% load proposed_1576_noColor.mat % bwP
% idxProposed_noColor = simIdx(bwP, bw)
% 
% load proposed_1576_Color.mat % bwP
% idxProposed = simIdx(bwP, bw)