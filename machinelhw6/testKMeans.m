function testKMeans()
% runs a simple k-means demo in 2D.
    X = [];
    K = 2; % change this to add more clusters.
    for i = 1:K
        X = [X; bsxfun(@plus, 0.1*randn(60, 2), randn(1,2))];
    end
%    [y, assigns] = k_means(X, K, X(1:K,:));
    [y, assigns] = k_means(X, K);

    figure;
    hold on;
    for i = 1:K
        scatter(X(assigns == i, 1), X(assigns == i, 2), [], rand(1, 3));
    end
    plot(y(:, 1), y(:, 2), 'bx');
    hold off;
end
