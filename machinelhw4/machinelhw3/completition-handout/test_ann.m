function [y] = test_ann(XTest)
    %% BEGIN SOLUTION
    % You are free to modify any code in this section. Just make sure to 
    % load your trained weights from weights.mat
    load('weights.mat');
    [nTest,f] = size(XTest);
    h = logsig([XTest,ones(nTest,1)]*Wih);
    [~,y] = max([h,ones(nTest,1)]*Who,[],2);
    y = y-1;
    %% END SOLUTION
end