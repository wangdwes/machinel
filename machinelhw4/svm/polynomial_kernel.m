function kernel = polynomial_kernel(x, z)
  kernel = (x' * z + 1) ^ 2;
end
