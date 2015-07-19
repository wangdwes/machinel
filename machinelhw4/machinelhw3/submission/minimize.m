function theta = minimize(f, init_theta)
% MINIMIZE  Find a local minima of a function f, starting at init_theta.
%          f - function to be minimized. f should be of the form:
%              [cost, grad(, hess)] = f(theta)
    tol = 1e-5; % you should stop optimization when the absolute difference
                % in cost between two iterations is less than tol.
                
    maxIter = 1000; % you should alternatively break after maxIter.
    alpha = 0.1;
    alpha_decay = 0.998;
    % Write your solution below. You should use either gradient descent or
    % Newton-Raphson to find the (local) minimum of the function.
    % Our solution is ~10 lines

    %% BEGIN SOLUTION (GRADIENT DESCENT) --- toc - tic > 1 second 

    cost = 0; theta = init_theta; 
    do  % iterate until no improvement on f... 
        prevcost = cost; [cost, grad] = f(theta);
        alpha *= alpha_decay; theta -= alpha * grad;
    until (!(maxIter-- && abs(prevcost - cost) >= tol))

    %% BEGIN SOLUTION (NEWTON'S METHOD)
  
    %% END SOLUTION
end
