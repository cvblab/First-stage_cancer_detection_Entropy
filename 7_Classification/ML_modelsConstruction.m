%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% BUILDING OF THE CLASSIFICATION MODELS %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc; 
set(groot,'ShowHiddenHandles','on'); delete(get(groot,'Children'));
rng('default')

%% Type of classifier
classifier = 1;         % 1 = quadratic SVM Classifier
                        % 2 = Neural Network Classifier (MLP)
version = 'Masks'; 

%% Paths
addpath('..\5_StatisticalAnalysis');
addpath('..\6_DataPartitioning');

%% Generation of Sets
load sets.mat
load featuresSelected.mat;

featTest1 = SETS(1).TEST;       featTrain1 = SETS(1).TRAIN;
featTest2 = SETS(2).TEST;       featTrain2 = SETS(2).TRAIN;
featTest3 = SETS(3).TEST;       featTrain3 = SETS(3).TRAIN;
featTest4 = SETS(4).TEST;       featTrain4 = SETS(4).TRAIN;
featTest5 = SETS(5).TEST;       featTrain5 = SETS(5).TRAIN;

tests = {featTest1; featTest2; featTest3; featTest4; featTest5};
trains = {featTrain1; featTrain2; featTrain3; featTrain4; featTrain5};

%% Definition of waitbar
leng = 5; % number of folds
ini = num2str(leng);
step = 1;
h = waitbar(0,ini,'Name','Building the MODEL...', 'CreateCancelBtn',...
    'setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0)

%% Building Models of Machine Learning
for i = 1:5
    % Sets train
    setTrain = trains{i};
    featTrain = setTrain(:,1:end-1);
    respTrain = setTrain(:,end);

    % Models construction
    if classifier == 1
        SVMclassifier{i} = Train_q_SVM_mc(setTrain, namesFeat);
    else
        [net,~] = Train_NN_mc(featTrain,respTrain, 3);
        netTrain{i} = net;
    end
    
    waitbar(step/leng,h,leng-step);
    step = step+1;
end
delete(h)

%% Save classifiers
if classifier == 1
    save qSVMclassifier_mc_v1.mat SVMclassifier mu sigma;
else
    save NNclassifier_mc_v1.mat netTrain mu sigma;
end
