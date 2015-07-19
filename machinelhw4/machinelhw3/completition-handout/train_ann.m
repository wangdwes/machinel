function [] = train_ann()
    % Write code for training a neural net
    %% BEGIN SOLUTION
    % You are free to modify any code in this section. Just make sure to 
    % store your weights
    load('digits.mat');
    [nTrain,f] = size(XTrain);
    % initialize weights
    Wih = rand(f+1,10);
    Who = rand(11,10);
    % store initial weights
    save('ini_weights.mat','Wih','Who');
    % put your training code here
    % save your final weights
    save('weights.mat','Wih','Who');
    %% END SOLUTION
end