function yPredict = TrainCrossSet(xTrain, yTrain, label)

  yPredict = ones(size(yTrain)); index = 1;
  do 
    isHeldOut = (label == index);
    yPredict(find(isHeldOut)) = TrainHeldOut(xTrain, yTrain, isHeldOut);
  until (!any(label == ++index))

end
