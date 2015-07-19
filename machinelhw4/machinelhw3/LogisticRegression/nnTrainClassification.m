function theta = nnTrainClassification(X, y, opt)
%  Wrapper function for nnTrain for classification problems.
%  Converts the m x 1 vector of labels into an k x m indicator matrix,
%  where the ith col will have k - 1 zeros and 1 one in the jth row
%  indicating that the ith observation is in class j. This indicator matrix
%  is then passed to nnTrain.
%
%  e.g. if we get y = [1; 3; 1; 2] then we would get an indicator matrix of
%       [1 0 1 0;
%        0 0 0 1;
%        0 1 0 0]
%
%  X           - m x n design matrix
%  y           - m x 1 labels (must be discrete).
%                IMPORTANT: If there are k classes then the class labels
%                should be in 1:k.
%  hidden_size - number of hidden units in the neural network.
%  output_size - number of units in output layer. Should be the number of
%                classes.
%
%  preds - m x 1 vector of predictions
% 

    Y_train = full(sparse(y, 1:length(y), 1));
    theta = nnTrain(X, Y_train, opt);
    
end
