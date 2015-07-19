
function [centers, belongings] = k_means(data, k, start_points)

  datanum = size(data, 1); 
  belongings = zeros(1, datanum);
  if (nargin == 3) centers = start_points; 
  else % centers = resize(data, k, size(data, 2))(randperm(k), :); % .... or k-means++ here. 
    centers(1, :) = data(randi(datanum), :);
    for index = 2: k 
%      % Note that when computing the distance of a point to itself, numerical error could occur - max with zero. 
      nearest_distances = min(bsxfun(@plus, sumsq(centers, 2), sumsq(data', 1)) - 2 * centers * data', [], 1);
      centers = [centers; data(discrete_rnd(linspace(1, datanum, datanum), max(nearest_distances, 0), 1), :)]; 
    end  
  end

  do 
    old_belongings = belongings; 

    % See the report for detailed explanation on computing the distances. 
    % Then, find belongings and derive a kronecker delta matrix. 
    distances = bsxfun(@plus, sumsq(centers, 2), sumsq(data', 1)) - 2 * centers * data';
    [~, belongings] = min(distances); kronecker = bsxfun(@eq, linspace(1, k, k)', belongings);

    % Cannot take the mean if there is no element in the operand - division by zero. 
    % This is a walkaround. Pick those clusters with non zeros and find the mean, then merge. 
    nonzeros = find(sum(kronecker, 2) > 0); 
    centers(nonzeros, :) = bsxfun(@rdivide, (kronecker * data)(nonzeros, :), sum(kronecker, 2)(nonzeros)); 

  % Iterate until belongings no longer changes. 
  until (all(belongings == old_belongings))

  % Transpose the matrix to meet the spec.
  belongings = belongings';

end
