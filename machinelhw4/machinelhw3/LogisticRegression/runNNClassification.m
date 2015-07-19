function [] = runNNClassification()
% fits a neural network to a non-linearly seperable 2D binary dataset and
% then plots the estimated probabilities.

addpath ./helpers
load ../data/non_linear.mat

%X = [randn(100, 2);randn(100, 2) + 3.0];
%y = [ones(100, 1); zeros(100,1)];

lambda = 0.1;
hidden_size = 3;

lambda_overfit = 1e-4;
hidden_size_overfit = 20;

lambda_underfit = 1.0;
hidden_size_underfit = 1;

point_size = 40;
if isOctave(), point_size = 4; end

subplot(1, 3, 1);
hold on;
plotProbabilities(X, y, lambda_underfit, hidden_size_underfit);
scatter(X(:, 1), X(:, 2), point_size, y, 'filled');
pbaspect([1, 1, 1]);
if ~isOctave()
    title('\lambda = 1.0, hidden units = 1, Underfit', 'FontSize', 17);
end
hold off;

subplot(1, 3, 2);
hold on;
plotProbabilities(X, y, lambda, hidden_size);
scatter(X(:, 1), X(:, 2), point_size, y, 'filled');
pbaspect([1, 1, 1]);
if ~isOctave()
    title('\lambda = 0.1, hidden units = 4, Good fit', 'FontSize', 17);
end
hold off;

subplot(1, 3, 3);
hold on;
plotProbabilities(X, y, lambda_overfit, hidden_size_overfit);
scatter(X(:, 1), X(:, 2), point_size, y, 'filled');
pbaspect([1, 1, 1]);
if ~isOctave()
    title('\lambda = 1\times10^{-4}, hidden units = 20, Overfit',...
          'FontSize', 17);
end
hold off;

end

function [] = plotProbabilities(X, y, lambda, hidden_size)
% plot a heatmap of predicted probabilities.
    opt.beta = 0.0;
    opt.hidden_sizes = hidden_size;
    opt.lambda = lambda;
    theta = nnTrainClassification(X, y+1, opt);

    x1s = linspace(min(X(:, 1)), max(X(:, 1)), 300);
    x2s = linspace(min(X(:, 2)), max(X(:, 2)), 300);
    [x1, x2] = meshgrid(x1s, x2s);
    
    x1s = size(x1);

    X_ = [x1(:), x2(:)];

    a3 = nnComputeActivations(theta, X_, 2, opt);
    a3 = a3(2,:)./sum(a3);
    
    X1 = reshape(X_(:, 1), x1s);
    X2 = reshape(X_(:, 2), x1s);
    y_ = reshape(a3, x1s);
    
    % training accuracy
    preds = nnPredictClassification(X, theta, 2, opt);
    corrects = preds-1 == (y');
    accuracy = 100*mean(corrects);
    fprintf('Train accuracy: %.3f%%\n', accuracy);
    
    % plot probabilities
    contourf(X1, X2, y_, 50, 'linecolor', 'none');
    if ~isOctave(), alpha(0.1); end
    
end