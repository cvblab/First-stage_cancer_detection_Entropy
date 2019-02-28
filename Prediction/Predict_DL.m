function [labels, probabilities] = Predict_DL(modelNumber, modelsDir, inputI)

modelfile = [modelsDir '\model' num2str(modelNumber) '.json'];
files = dir([modelsDir '\model' num2str(modelNumber) '*.hdf5']);
weights = [modelsDir '\' files(end).name];
net = importKerasNetwork(modelfile,'WeightFile',weights, ...
      'OutputLayerType','classification');

[labels, probabilities] = classify(net, inputI, 'MiniBatchSize', 4, 'ExecutionEnvironment', 'gpu');