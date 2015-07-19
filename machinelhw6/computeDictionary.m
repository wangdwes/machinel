clear all; addpath('image');addpath('image/inst');
load('data.mat'); % load image paths and filter bank  

k = 169;      % Number of clusters. 
ratio = 0.03; % What percentage of total points should be considered?
data = [];    % Initialize something so vertcat doesn't get baffled.

for index = 1: length(imagePaths) tic; 

    % Extract the feature points, as instructed, then randomly pick a few data points. 
    features = extractFilterResponses(imread(imagePaths{index}), filterBank);
    data = [data; features(randperm(size(features, 1)) < size(features, 1) * ratio, :)];

    % Save the user from endless waiting. 
    fprintf('Image %d/%d processed: ', index, length(imagePaths)); toc;

end 

% Friendly reminder to the user such that he/she can grab some food. 
fprintf('Running K-means with %d clusters and %d data points...\n', k, size(data, 1)); tic;
[dictionary, ~] = k_means(data, k); save('dictionary.mat', 'dictionary'); 
fprintf('K-means completed and a dictionary has been saved. '); toc;

