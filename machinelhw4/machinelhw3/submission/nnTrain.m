function theta = nnTrain(X, y, opt)
%NNTRAIN  Trains a neural network with opt.hidden_sizes hidden units.
%
% function theta = nnTrain(X, y, opt)
%
%  X             - m x n design matrix
%  y             - k x m labels for each x_i
%  opt           - Struct containing NN parameters:
%                    lambda and hidden_sizes.
%
%  theta - flattened vector of ALL parameters in the neural network.
% 
    addpath ./helpers
    
    % Set default parameters in case opt doesn't define them all.
    if ~isfield(opt, 'lambda')
        opt.lambda = 0;
    end
    if ~isfield(opt, 'hidden_sizes')
        opt.hidden_sizes = 64;
    end
    
    all_layer_sizes = [size(X, 2); opt.hidden_sizes; size(y, 1)];
    
    if isfield(opt,'init_theta')
        init_theta = opt.init_theta;
    else
        init_theta = initializeParameters(all_layer_sizes);
    end
    
    % For this part of the assignment we will use a minimization
    % tool, minFunc. The interface is similar to the gradient descent
    % optimizer you wrote, but it uses some fancy tricks to try and
    % approximate the Hessian in order to reduce the number of times J must
    % be called.
    addpath ./minFunc
    addpath ./minFunc/autoDif
    addpath ./minFunc/compiled
    
    % you may have to run mexAll in the minFunc directory to compile the
    % minimizer.
    
    % check that the gradients look reasonable.
    toyopt.hidden_sizes = 4;
    toyopt.lambda = 0.0;
    toyw = initializeParameters([3; 4; 2]);
    toyopt.init_theta = toyw;
    toyX = rand(10, 3);
    toyY = rand(1, 10) < 0.5;
    toyY = [toyY; ~toyY];
    checkGradient(@(x) costNN(toyX, toyY, x, toyopt), toyw);
    
    theta = minFunc(@(x) costNN(X, y, x, opt), init_theta, opt);
end

function checkGradient(costfn, theta)
% CHECKGRADIENT  Checks your gradient against finite differences gradient
%               computed from function evaluation. Prints a warning if
%               L2 norm difference (sum of squares of differences) is
%               too large.
%
%   J should be a function of the form
%        [cost : scalar, grad : mx1 vector] = J(theta : mx1 vector)
%
    tol = 1e-6;
    n_dims = size(theta, 1);
    dx = 1e-5;
    id = eye(n_dims)*dx;
    grad_hat = zeros(size(theta));
    [~, grad] = costfn(theta); % get true analytical gradient.
    for i = 1:n_dims % now compute finite differences gradient.
        xp1 = costfn(theta + id(:, i));
        xm1 = costfn(theta - id(:, i));
        grad_hat(i) = (xp1 - xm1)/(2*dx); 
    end
    
    diff = norm(grad_hat - grad);
    
    if diff > tol
        fprintf(['Warning: Your gradients differ by %f. Your gradient'...
                 ' or cost function may be incorrect.\n'], diff);
    end
    assert(diff <= tol);
end
