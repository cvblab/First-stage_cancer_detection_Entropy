function [net, targets_train] = Train_NN_mc(inputs_train,targets_train, classes)

%% Inputs and targets
inputs_train = inputs_train';
targetsTrain = targets_train';
num = unique(targetsTrain);
targets_train = zeros(classes, size(targetsTrain,2));
for i = 1:length(num)
   idx = find(targetsTrain == num(i));
   targets_train(i,idx) = 1;
end

%% Definition of a new neural network
newHiddenSizes = 15;
net = patternnet(newHiddenSizes);
% Parameters
net.PerformFcn = 'crossentropy'; % Loss function
net.trainFcn = 'traingdx'; % Gradient descent with a momentum and adaptative learning rate backpropagation
net.trainParam.lr = 1e-2; 
net.trainParam.lr_inc = 1.5;
net.trainParam.lr_dec = 0.5;
net.trainParam.mc = 0.95;
net.trainParam.epochs = 1000; % Maximum number of epochs
net.trainParam.max_fail = 20; % Stop criterion
% Indices
net.divideFcn = 'divideind';
trInd = 1:size(inputs_train,2);
vlInd = randi([1 size(inputs_train,2)],1,round(0.1*size(inputs_train,2)));
trInd(vlInd) = [];
net.divideParam.trainInd = trInd;
net.divideParam.valInd = vlInd;
% Network configuration
net = configure(net,inputs_train,targets_train);
% Network training
net.trainParam.showWindow = false;
[net,tr] = train(net,inputs_train,targets_train);


