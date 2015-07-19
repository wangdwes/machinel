function testInstanceLabel = PartitionHeldOut(numberOfInstances, k)
  testInstanceLabel = (randperm(numberOfInstances) <= numberOfInstances / k)';
end
