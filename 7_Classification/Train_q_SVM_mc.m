function [finalClassification, C, gamma] = Train_q_SVM_mc(trainingData, namesFeat)

inputTable = array2table(trainingData, 'VariableNames', [namesFeat; 'response']);
predictorNames = namesFeat;
predictors = inputTable(:, predictorNames);
response = inputTable.response;

%% Cross validation
c = cvpartition(size(predictors,1),'Kfold',10);

%% Grade of polynomial
t = templateSVM('KernelFunction','polynomial',...
                'PolynomialOrder', 2,...
                'Standardize', true); 

%% Options of optimization
opts = struct('Optimizer','bayesopt','ShowPlots',true,'CVPartition',c,...
    'AcquisitionFunctionName','expected-improvement-plus', 'Verbose',1,...
    'SaveIntermediateResults', true, 'MaxObjectiveEvaluations',20, 'UseParallel', true);

params = hyperparameters('fitcecoc',predictors,response,'svm');
paramsToOptimize = params(2:3);
paramsToOptimize(1).Range = [1e-2,1e2];
paramsToOptimize(2).Range = [1e-2,1e2];

%% Multiclass SVM classifier
classificationSVM = fitcecoc(predictors,response,'Coding','onevsall','Learners',t,...
'OptimizeHyperparameters',paramsToOptimize,...
'HyperparameterOptimizationOptions', opts);

C = classificationSVM.HyperparameterOptimizationResults.XAtMinObjective.BoxConstraint;
gamma = classificationSVM.HyperparameterOptimizationResults.XAtMinObjective.KernelScale;

%% Final Classifier
t2 = templateSVM('KernelFunction','polynomial',...
                'PolynomialOrder', 2,...
                'Standardize', true, ...
                'BoxConstraint', C, ...
                'KernelScale', gamma);
            
finalClassification = fitcecoc(predictors,response, 'Coding','onevsall', 'Learners', t2);
