clear
clc
load('ps4-svm.mat')

C = [0, 0.1, 0.3, 0.5, 1, 2, 5, 8, 10];
%C = 2;

new_x = [x_test;x_train];
new_y = [y_test;y_train];

crossSetLabel = PartitionCrossSet(size(new_x,1),4);
trainErr = [];
testErr =[];

for i = 1:length(C)
    trainErrk = [];
    testErrk = [];

    for k=1:max(crossSetLabel)
        testInstanceLabel = zeros(size(new_y,1),1);
        testInstanceLabel(crossSetLabel==k,:)=1;

        partitionx = new_x(testInstanceLabel==0,:);
        partitiony = new_y(testInstanceLabel==0,:);

        model = svm_train(partitionx, partitiony, C(i));
        predictedTrain = svm_classify(model, new_x(testInstanceLabel==0,:));
        trainErrk = [trainErrk 1-mean(predictedTrain==partitiony)];

        predicted =  svm_classify(model, new_x(testInstanceLabel==1,:));
        testErrk = [testErrk 1-mean(predicted==new_y(testInstanceLabel==1))];
    end
    trainErr = [trainErr mean(trainErrk)];
    testErr = [testErr mean(testErrk)];
end

plot(C,testErr,'-r');
hold on
plot(C,trainErr,'-b');
legend('Test Error','Train Error');
title('2.3 Training Error vs Testing Error');
ylabel('Error');
xlabel('C');
