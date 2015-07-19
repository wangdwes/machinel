function [accuracy, lowerInterval, upperInterval] = ConstructInterval(yPredict, yTest, confLevel)

  % https://people.richland.edu/james/lecture/m170/ch08-int.html
  % z-score with confience levels of 50%, 80%, 90%, 95%, 98%, 99%.
  zScores = [0.674, 1.282, 1.645, 1.960, 2.326, 2.576];
  confLevels = [0.50, 0.80, 0.90, 0.95, 0.98, 0.99];
  
  accuracy = mean(yPredict == yTest);
  zScore = zScores(sum(confLevel >= confLevels));
  lowerInterval = accuracy - zScore * sqrt(accuracy * (1 - accuracy) / prod(size(yTest)));
  upperInterval = accuracy + accuracy - lowerInterval;
  
end
