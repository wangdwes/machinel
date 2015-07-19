% read the images and save them in a matrix.

fid = fopen('all.list'); faces = [];
while !feof(fid)
  f = imread(fgetl(fid));
  faces = [faces; double(f(:)')];
end

% note that we have far more features than samples, it is more computationally feasible
% to only find out the singular value decomposition for the covariance matrix instead of the original.
[u, s, d] = svd(faces); 

% extract the first five eigenfaces, scale them into the range of [0, 255]
% and plot them in a nice row with their indices as titles. 
for index = 1: 5
  ef = d(:, index); scaled = uint8((ef - min(ef)) ./ (max(ef) - min(ef)) * 256); 
  subplot (1, 5, index, 'align'); imshow (reshape(scaled, size(f))); title(num2str(index));
end

% project kawamura into the new space - though we don't need all the scores but for now - 
% let's compute everything and extract some of them in the loop. then scale and plot it, 
% don't forget to add the mean and compute the reconstruction error and display it.

reduced_kawamura = (faces(1, :) - mean(faces(1, :))) * d;
% reduced_kawamura = center(faces(1, :)) * d; 
index = 0; 
for n = 50: 50: 600
  kawamura = reduced_kawamura(:, 1: n) * d(:, 1: n)' + mean(faces(1, :));
  scaled = uint8((kawamura - min(kawamura)) ./ (max(kawamura) - min(kawamura)) * 256); 
  subplot (3, 4, ++index, 'align'); imshow (reshape(scaled, size(f)));
  title (sprintf('%d (%.3f)', n, sum(sumsq(kawamura - faces(1, :)))));
end


