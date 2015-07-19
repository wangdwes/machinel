function kernel = gaussian_kernel(x, z)
  kernel = exp(-0.5 * sumsq(x - z));
end 
