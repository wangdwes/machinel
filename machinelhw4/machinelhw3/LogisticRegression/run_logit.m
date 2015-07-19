function errdiff = run_logit()
% 10-601 HW3 starter code.
% PART 1: LOGISTIC REGRESSION
% Do not modify this file.

addpath ./helpers
load ../data/mnist.mat

lambda = 1; % Don't change this.

% filter 4s and 7s in the training and test set.

[X_train, Y_train] = filter4s7s(X_train, Y_train);
[X_test, Y_test] = filter4s7s(X_test, Y_test);

% train your logistic regression code.
theta = trainLR(X_train, Y_train, lambda); % you need to implement this function
% make predictions.

partition = PartitionCrossSet(prod(size(Y_test)), 10);
errorsLR = zeros(10, 1);

for label = 1: 10
  yPred = predictLR(X_test(find(partition == label), :), theta);
  errorsLR(label) = mean(yPred != Y_test(find(partition == label)));
endfor

n_classes = 2;
opt.hidden_sizes = 64;
opt.lambda = 0.1;   
opt.MaxIter = 400; % max iterations for minimization function.    
opt.beta = 0.0;
opt.p = 0.01;
    
% add 1 to labels so that they are in 1:10, instead of 0:9
theta = nnTrainClassification(X_train, Y_train + 1, opt);
    
% visualize the weights we learned.
% Ws = unflattenParameters(theta, [28*28; opt.hidden_sizes; n_classes]);
% W1 = Ws{1};
% show_centroids((W1 - min(min(W1)))./max(max(W1 - min(min(W1)))),28,28);
    
errorsNN = zeros(10, 1);

for label = 1: 10
  yPred = nnPredictClassification(X_test(find(partition == label), :), theta, n_classes, opt);
  errorsNN(label) = mean((yPred - 1) != Y_test(find(partition == label))');
endfor


errorsLR
errorsNN
errdiff = errorsLR - errorsNN;

end

function [X_f, Y_f] = filter4s7s(X, Y)
    X_f = X(Y == 4 | Y == 7, :);
    Y_f = Y(Y == 4 | Y == 7);
    Y_f = Y_f == 7;
end
