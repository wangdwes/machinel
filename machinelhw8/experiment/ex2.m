% extract the images from the file specified by file, 
% with proper labels set by parsing the filename. 
function [labels, instances] = extract(file)
  fid = fopen(file); instances = labels = [];
  while !feof(fid)
    path = fgetl(fid); 
    instances = [instances; double(imread(path)(:)')];
    labels = [labels; any(strfind(path, 'sunglasses'))];
  end
end

% scale the features as required by the libsvm svmtrain method. 
function [scaled] = scale(instances, lb, ub)
  scaled = (instances - repmat(lb, size(instances, 1), 1)) ./ ...
    repmat(ub - lb, size(instances, 1), 1);
end

% this is the beginning of this script.
% include the libsvm interface for matlab/octave. 
addpath('libsvm-3.20/matlab/');

% train the support vector machine. 
[train_labels, train_instances] = extract('all_train.list');
ub = max(train_instances); lb = min(train_instances); 
scaled_train_instances = scale(train_instances, lb, ub);
[test1_labels, test1_instances] = extract('all_test1.list');
scaled_test1_instances = scale(test1_instances, lb, ub);
[test2_labels, test2_instances] = extract('all_test2.list');
scaled_test2_instances = scale(test2_instances, lb, ub);

% generate some baseline classification accuracies...
fprintf('Baseline:\n'); tic;
model = svmtrain(train_labels, scaled_train_instances, '-q -c 100 -g 0.01');
predicted1_label = svmpredict(test1_labels, scaled_test1_instances, model);
predicted2_label = svmpredict(test2_labels, scaled_test2_instances, model); toc;

% now apply principal component analysis to see if the performance gets better...
% apparently we're supposed to use the same basis vectors... or the scores don't make sense.
[all_labels, all_instances] = extract('all.list');
[~, ~, d] = svd(all_instances);
reduced_train_instances = center(scaled_train_instances) * d;
ub = max(reduced_train_instances); lb = min(reduced_train_instances);
reduced_train_instances = scale(reduced_train_instances, lb, ub);
reduced_test1_instances = scale(center(scaled_test1_instances) * d, lb, ub);
reduced_test2_instances = scale(center(scaled_test2_instances) * d, lb, ub);

fprintf('Reduced to 50 features...\n'); tic;
model_pca50 = svmtrain(train_labels, reduced_train_instances(:, 1: 50), '-q -c 500 -g 0.07');
svmpredict(test1_labels, reduced_test1_instances(:, 1: 50), model_pca50);
svmpredict(test2_labels, reduced_test2_instances(:, 1: 50), model_pca50); toc;

fprintf('Reduced to 150 features...\n'); tic;
model_pca150 = svmtrain(train_labels, reduced_train_instances(:, 1: 150), '-q -c 100 -g 0.16');
svmpredict(test1_labels, reduced_test1_instances(:, 1: 150), model_pca150);
svmpredict(test2_labels, reduced_test2_instances(:, 1: 150), model_pca150); toc;

fprintf('Reduced to 200 features...\n'); tic;
model_pca200 = svmtrain(train_labels, reduced_train_instances(:, 1: 200), '-q -c 100 -g 0.07');
svmpredict(test1_labels, reduced_test1_instances(:, 1: 200), model_pca200);
svmpredict(test2_labels, reduced_test2_instances(:, 1: 200), model_pca200); toc;

