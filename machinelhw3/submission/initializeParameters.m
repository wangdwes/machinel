function theta = initializeParameters(layer_sizes)
% initializes parameters for neural network randomly.

    n_layers = length(layer_sizes)-1;
    bs = cell(n_layers, 1);
    Ws = cell(n_layers, 1);
    for i = 1:n_layers
        isize = layer_sizes(i);
        i1size = layer_sizes(i+1);
        bs{i} = zeros(i1size, 1);
        r = sqrt(6)/sqrt(isize+i1size+1);
        Ws{i} = 2*r*rand(i1size, isize) - r;
    end
    
    theta = flattenParameters(Ws, bs);
end