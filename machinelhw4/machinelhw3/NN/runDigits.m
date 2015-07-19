function errors = runDigits()
% Trains a single layer NN on a subset of the MNIST set and evaluates the
% NN.
    addpath ./helpers

    load ../data/mnist.mat

    % don't change the parameters in this file.
    % we get ~94.8% accuracy with these parameters. 
    n_classes = 10;
    opt.hidden_sizes = 64;
    opt.lambda = 0.1;
    
    opt.MaxIter = 400; % max iterations for minimization function.
    
    opt.beta = 0.0;
    opt.p = 0.01;
    
    % add 1 to labels so that they are in 1:10, instead of 0:9
    theta = nnTrainClassification(X_train, Y_train+1, opt);
    
    % visualize the weights we learned.
    Ws = unflattenParameters(theta, [28*28; opt.hidden_sizes; n_classes]);
    W1 = Ws{1};
    show_centroids((W1 - min(min(W1)))./max(max(W1 - min(min(W1)))),28,28);
    
    % show the accuracy.
    trainpreds = nnPredictClassification(X_train, theta, n_classes, opt);

    partition = PartitionCrossSet(prod(size(Y_test), 10);
    errors = zeros(10, 1);

    for label = 1: 10
      yPred = nnPredictClassification(X_test(find(partition == label)), theta, n_classes, opt);
      errors(label) = mean(yPred != Y_test(find(partition == label))');
    endfor

end
