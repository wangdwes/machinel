function [] = ploterr(c)

  load('ps4-svm.mat');
  
  xMerged = vertcat(x_train, x_test);
  yMerged = vertcat(y_train, y_test);
  yPredict = zeros(size(yMerged));

  aveTrainError = [];
  aveTestError = [];
  numberOfInstances = prod(size(yMerged));
  partition = PartitionCrossSet(numberOfInstances, 4);

  for index = 1: prod(size(c)) 

    trainError = [];
    testError = [];

    for label = 1: 4

      isHeldOut = (partition == label); 

      model = svm_train(xMerged(!isHeldOut, :), yMerged(!isHeldOut), c(index));
      trainError = [trainError, mean(svm_classify(model, xMerged(!isHeldOut, :)) != yMerged(!isHeldOut))];
      testError = [testError, mean(svm_classify(model, xMerged(isHeldOut, :)) != yMerged(isHeldOut))];

    endfor

    aveTrainError = [aveTrainError, mean(trainError)];
    aveTestError = [aveTestError, mean(testError)];

  endfor

  plot(c, aveTrainError,'linewidth', 2, 'r'); hold on;
  plot(c, aveTestError,'linewidth', 2, 'g');

endfunction
