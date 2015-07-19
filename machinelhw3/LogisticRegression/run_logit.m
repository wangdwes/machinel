function [] = run_logit()
% 10-601 HW3 starter code.
% PART 1: LOGISTIC REGRESSION
% Do not modify this file.

load ../data/mnist.mat

lambda = 1; % Don't change this.

% filter 4s and 7s in the training and test set.

[X_train, Y_train] = filter4s7s(X_train, Y_train);
[X_test, Y_test] = filter4s7s(X_test, Y_test);

% train your logistic regression code.
theta = trainLR(X_train, Y_train, lambda); % you need to implement this function
% make predictions.
Y_preds = predictLR(X_test, theta); % you need to implement this function

% show the test accuracy and display misclassified digits.
corrects = Y_preds == Y_test;
accuracy = 100*mean(corrects);
fprintf('Logistic Regression Accuracy: %.3f%%\n', accuracy);

% show the misclassifications.
if sum(not(corrects)) < 1600
    show_centroids(X_test(not(corrects),:),28,28);
end

end

function [X_f, Y_f] = filter4s7s(X, Y)
    X_f = X(Y == 4 | Y == 7, :);
    Y_f = Y(Y == 4 | Y == 7);
    Y_f = Y_f == 7;
end
