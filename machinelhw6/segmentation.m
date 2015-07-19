clear all; addpath('image');addpath('image/inst');
load('data.mat'); % load filter bank
load('dictionary.mat'); % load dictionary

% Read the image and extract the features...
img = imread('myBasilica.jpg'); 
img_size = [size(img, 1), size(img, 2)];

% Extracting features...
fprintf('Extracting features on %d pixels...\n', prod(img_size)); tic;
warning("off"); features = extractFilterResponses(img, filterBank); toc; 

% Segmenting...
fprintf('Segmenting...\n'); tic; 
distances = bsxfun(@plus, sumsq(dictionary, 2), sumsq(features', 1)) - 2 * dictionary * features';
[~, belongings] = min(distances); segments = reshape(belongings, img_size); toc;

cmap = rainbow(size(dictionary, 1));
imagesc(segments); colormap(cmap);
imwrite(segments, cmap, 'myBasilicaSegmented.jpg');
