10-601 Fall 2014 Assignment 3: Logistic Regression and Neural Networks.

Primary TAs for this assignment:

Harry Gifford - hgifford@andrew.cmu.edu
Jin Sun       - jins@andrew.cmu.edu

Please use Piazza if you have any questions.

Please see handout and individual files for more details.

- Directory structure:

(*) indicates you should modify this file.

./data/            - contains all of the training and test data you will use. Do not touch anything in this folder.

./run_logit.m      - train and test the logistic regression classifier on a subset of MNIST [1]. Simply run this script in octave/matlab and it will print out the accuracy as well as display the misclassified digits. Our test accuracy is 98.202%.

(*)./costLR.m      - logistic cost function.
(*)./minimize.m    - find a local minimum of a function given a starting location.
(*)./trainLR.m     - function used to train the logistic regression classifier.
(*)./predictLR.m   - function to make predictions, given a trained set of weights.

./show_centroids.m - display a subset of digits.

Again, you should modify (probably in this order):
  costLR.m
  minimize.m
  trainLR.m
  predictLR.m

You should NOT modify:
  show_centroids.m
  run_logit.m

- Sources:

[1] Y. LeCun, L. Bottou, Y. Bengio, and P. Haffner. "Gradient-based learning applied to document recognition." Proceedings of the IEEE, 86(11):2278-2324, November 1998.

[2] An Analysis of Single-Layer Networks in Unsupervised Feature Learning, Adam Coates, Honglak Lee, and Andrew Y. Ng.

[3] Hinton, G. E. & Salakhutdinov, R. R. (2006), 'Reducing the dimensionality of data with neural networks', Science 313 (5786) , 504-507 .