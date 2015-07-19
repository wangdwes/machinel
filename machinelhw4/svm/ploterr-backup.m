function [] = ploterr(c)

  load('ps4-svm.mat');
  
  xMerged = vertcat(x_train, x_test);
  yMerged = vertcat(y_train, y_test);
  yPredict = zeros(size(yMerged));

  trainError = zeros(4, 1);
  aveTrainError = zeros(size(c));
  aveTestError = zeros(size(c));

  numberOfInstances = prod(size(yMerged));
  partition = PartitionCrossSet(numberOfInstances, 4);

  for index = 1: prod(size(c)) 

    label = 1;
    do 

      isHeldOut = (partition == label);

      model = svm_train(xMerged(find(!isHeldOut),:), yMerged(find(!isHeldOut)), c(index));
      yPredict(find(isHeldOut)) = svm_classify(model, xMerged(find(isHeldOut),:));
      trainError(label) = mean(svm_classify(model, xMerged(find(!isHeldOut),:)) != yMerged(find(!isHeldOut)) );

    until (!any(partition == ++label))

    aveTrainError(index) = mean(trainError);
    aveTestError(index) = mean(yPredict != yMerged);

  endfor
 
  plot(c, aveTrainError,'linewidth', 2, 'r'); hold on;
  plot(c, aveTestError,'linewidth', 2, 'g');

endfunction
