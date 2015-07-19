% Function   : train_svm
% 
% Purpose    : Train a svm based on the dataset.
% 
% Parameters : trainning data feature, training data
% labels, test data feature, parameter, A kernal function
% 
% Return     : SVM model.

function model = svm_train(x_train, y_train, C, kernel)
    
    % default linear kernel
    if nargin < 4
        kernel = @(x,y) x'*y;
    end
    
    n = size(x_train, 1);
    H = zeros(n, n);

    for i = 1:n
        for j = 1:n
            H(i,j) = y_train(i)*y_train(j)*kernel(x_train(i,:)', x_train(j,:)');
        end
    end

    f = ones(n, 1);
    
    Aeq = y_train';
    beq = 0;
%    
 %   a = quadprog(H,-f,[],[],Aeq,beq, zeros(n,1), C*ones(n,1), [], optimset('Algorithm', 'interior-point-convex', 'Display', 'off'));

    %comment out the above line and use the following statement instead if you are using octave
    a = qp([], H, -f, Aeq, beq, zeros(n,1), C*ones(n,1),
           optimset('Algorithm', 'interior-point-convex', 'Display', 'off'));

    b = 0;
    c = 0;
    for i=1:n
        if a(i)> 0 && a(i)<C
            b = b + y_train(i);
            c = c + 1;
            for j=1:n
                if (a(i) > 0)
                    b = b-a(j)*y_train(j)*kernel(x_train(j,:)', x_train(i,:)');
                end;
            end  
        end
    end
    b = b / c;
    model.a = a;
    model.b = b;
    model.x_train = x_train;
    model.y_train = y_train;
    model.kernel = kernel;    
end
