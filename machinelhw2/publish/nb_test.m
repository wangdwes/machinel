function Pred_nb = nb_test(model, x)
  likelihood = bsxfun(@plus, log(model.condpr)*x', model.prior);
  Pred_nb = diff(likelihood)' > 0;
endfunction
