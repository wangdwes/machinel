function [] = runDigits()
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

    preds = nnPredictClassification(X_test, theta, n_classes, opt);
    
    % dont forget to subtract 1 from the labels we added to earlier!
    traincorrects = trainpreds-1 == (Y_train');
    trainaccuracy = 100*mean(traincorrects);
    fprintf('Train accuracy: %.3f%%\n', trainaccuracy);

    corrects = preds-1 == (Y_test');
    accuracy = 100*mean(corrects);
    fprintf('Test accuracy: %.3f%%\n', accuracy);
    
    % show the misclassifications.
    if sum(not(corrects)) < 1600
        show_centroids(X_test(not(corrects),:),28,28);
    end

end