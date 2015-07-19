function yPredict = TrainHeldOut(xTrain, yTrain, label)
  model = nb_train(xTrain(find(!label), :), yTrain(find(!label), :));
  yPredict = nb_test(model, xTrain(find(label), :));
end
