function predicted =  svm_classify(model, x_test)
    
    n = size(model.a,1);
    m = size(x_test, 1);
    predicted = zeros(m,1);
    x_train = model.x_train;
    y_train = model.y_train;
    a = model.a;
    
    for i=1:m
        s = 0;
        for j = 1:n
            if (a(i) > 0)
                s = s + a(j)*y_train(j)*model.kernel(x_train(j,:)', x_test(i,:)');
            end;
        end
        s = s + model.b;
        if s>0
            predicted(i) = 1;
        else
            predicted(i) = -1;
        end
    end    
end

