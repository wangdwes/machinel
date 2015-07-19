
% One use case of the implementation. Here we compared the training error
% and test error of non-kernel and Gaussion kernel SVM.

clear
clc
load('ps4-svm.mat')

C = 0.5;

% Linear Kernel
% kernel = @(x,z) x'*z;
% model_linear = svm_train(x_train, y_train, C, kernel);
% pred_linear = svm_classify(model_linear, x_test);

xMerged = vertcat(x_train, x_test);
yMerged = vertcat(y_train, y_test);

numberOfInstances = prod(size(yMerged));

while (1)

tic
rand('seed', rand());

partition = PartitionHeldOut(numberOfInstances, 2);
xTrain = xMerged(find(!partition));
yTrain = yMerged(find(!partition));

modelLinear = svm_train(xTrain, yTrain, C, @(x, z)x'*z);
predLinear = svm_classify(modelLinear, xMerged(find(partition)));

[accuracy, lowerInterval, upperInterval] = ConstructInterval(predLinear, yMerged(find(partition)), 0.95);
fprintf('%.1f\\%% & ', 100 * accuracy);
fprintf('[%.3f,%.3f] & ', lowerInterval, upperInterval);

[accuracy, lowerInterval, upperInterval] = ConstructInterval(predLinear, yMerged(find(partition)), 0.99);
fprintf('[%.3f,%.3f] & ', lowerInterval, upperInterval);

numberOfInstances = prod(size(yMerged));
toc
partition = PartitionHeldOut(numberOfInstances, 10);
xTrain = xMerged(find(!partition));
yTrain = yMerged(find(!partition));

modelLinear = svm_train(xTrain, yTrain, C, @(x, z)x'*z);
predLinear = svm_classify(modelLinear, xMerged(find(partition)));

[accuracy, lowerInterval, upperInterval] = ConstructInterval(predLinear, yMerged(find(partition)), 0.95);
fprintf('%.1f\\%% & ', 100 * accuracy);
fprintf('[%.3f,%.3f] & ', lowerInterval, upperInterval);

[accuracy, lowerInterval, upperInterval] = ConstructInterval(predLinear, yMerged(find(partition)), 0.99);
fprintf('[%.3f,%.3f] \\\\\n', lowerInterval, upperInterval);

endwhile


% kernel = @(x,z) x'*z;
% model_linear = svm_train([x_train;x_test], [y_train;y_test], C, kernel);
% pred_linear = svm_classify(model_linear, x_test);

% Gussian kernel. We can also specify the polynomial kernal here.
% kernel = @(x,z) gaussian_kernel(x, z);
% model_gaussian = svm_train(x_train, y_train, C, kernel);
% pred_gaussian = svm_classify(model_gaussian, x_test);

% polynomial kernel.
% kernel = @(x,z) polynomial_kernel(x, z);
% model_polynomial = svm_train(x_train, y_train, C, kernel);
% pred_polynomial = svm_classify(model_polynomial, x_test);

% fprintf('C = %.3f%5\n', C);
% fprintf('Linear Accuracy: %.3f%%\n', 100 * mean(y_test == pred_linear));
% fprintf('95% Confidence Interval = [%.3f%%, %.3f%%]',  )

% fprintf('Gaussian Accuracy: %.3f%%\n', 100 * mean(y_test == pred_gaussian));
% fprintf('Polynomial Accuracy: %.3f%%\n', 100 * mean(y_test == pred_polynomial));
