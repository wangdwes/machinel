load('digits.mat')
load('weights.mat')

% This file just verify if your output matrix has correct dimensions
% Generate some random data and labels
XTest = rand(randi([8000,11000]),784);
nTest = size(XTest,1);
yTest = randi(9,nTest,1);
% check test accuracy
y = test_ann(XTest);
accu_test = sum((y-yTest) == 0)/nTest;
% final score: (just an example)
% disp(accu_test*100)
if accu_test >=0 && accu_test <=1
    disp('Your code is in good shape')
end