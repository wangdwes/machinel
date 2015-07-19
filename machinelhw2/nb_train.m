function model = nb_train(x, y)
  cumcount = [!y'*x; y'*x]; % concatenate & count; 
  model.prior = [!y'*!y; y'*y] / rows(y); % compute priors; 
  model.condpr = bsxfun(@rdivide, cumcount + 1, sum(cumcount, 2) + columns(x));
endfunction
