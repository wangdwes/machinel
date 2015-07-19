
function [centroids, classes] = kmeans(observ, k, initial_points)

  obsnum = size(observ, 1); classes = zeros(obsnum, 1);
  prevclasses = ones(obsnum, 1); centroids = observ(randperm(obsnum, k), :);

  while any(prevclasses ~= classes), prevclasses = classes; 
    [~, classes] = min(pdist2(observ, centroids), [], 2); 
    for uniqidx = unique(classes'),
      centroids(uniqidx, :) = mean(observ(classes == uniqidx, :), 1); 
    end 
  end

end
