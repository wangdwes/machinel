function yCrossSetLabel = PartitionCrossSet(numberOfInstances, k)
  yCrossSetLabel = ceil(randperm(numberOfInstances) * k / numberOfInstances)';
end

