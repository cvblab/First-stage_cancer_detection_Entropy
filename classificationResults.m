%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% EVALUATION OF RESULTS (TEST) %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; %clc; 
set(groot,'ShowHiddenHandles','on'); delete(get(groot,'Children'));
rng('default')

%% INPUTS
typePrediction = 2;
typeStudy = 1;
classifier = 1;
% Type of prediction
    % 1 = Machine Learning
    % 2 = Deep Learning
% Type of study
    % 1 = Artifacts Vs Glands
    % 2 = Healthy Vs Pathologist 
% classifier
    % 1 = quadratic SVM Classifier
    % 2 = Neural Network Classifier
    
%% Deep Learning
if typePrediction == 2
    addpath('Deep_Learning_Results');
end

%% Paths
addpath('4_FeatureExtraction'); 
addpath('5_StatisticalAnalysis');
addpath('6_DataPartitioning');
addpath('7_Classification');
addpath('Prediction');

    
%% Sets
if typePrediction == 1
    load featuresSelected.mat;
    load sets.mat
    
    featTest1 = SETS(1).TEST;       featTrain1 = SETS(1).TRAIN;
    featTest2 = SETS(2).TEST;       featTrain2 = SETS(2).TRAIN;
    featTest3 = SETS(3).TEST;       featTrain3 = SETS(3).TRAIN;
    featTest4 = SETS(4).TEST;       featTrain4 = SETS(4).TRAIN;
    featTest5 = SETS(5).TEST;       featTrain5 = SETS(5).TRAIN;
    tests = {featTest1; featTest2; featTest3; featTest4; featTest5};
    trains = {featTrain1; featTrain2; featTrain3; featTrain4; featTrain5};

    %% Classification Models
    if classifier == 1
        load qSVMclassifier_mc.mat;
    else
        load NNclassifier_mc.mat
    end

    %% Prediction
    figure
    for i = 1:5
        % setTest
        setTest = tests{i};
        featTest = setTest(:,1:end-1);
        respTest = setTest(:,end);     
        % Prediction phase
        if classifier == 1 
            [label_predictions, score] = predict(SVMclassifier{i},featTest);
        else
            score = sim(netTrain{i},featTest');
            label_predictions = obtainLabelPredictionsNeuralNetwork(score);
        end

        %% Results CLASSPERF
        if typeStudy == 1 % Artifacts Vs Glands
            CP = classperf(respTest,label_predictions,'Positive',0,'Negative',[1,2]);
        elseif typeStudy == 2 % Healthy Vs Pathological
            CP = classperf(respTest,label_predictions,'Positive',1,'Negative',2);
        end
        Sensitivity(i) = CP.Sensitivity;
        Specificity(i) = CP.Specificity;
        Ppv(i) = CP.PositivePredictiveValue;
        Npv(i) = CP.NegativePredictiveValue;
        F1Score(i) = 2.*Ppv(i).*Sensitivity(i)./(Ppv(i)+Sensitivity(i));
        Accuracy(i) = CP.LastCorrectRate;  
        confusionmatrix = CP.CountingMatrix;
        if typeStudy == 1
            TP = confusionmatrix(1,1);
            TN = confusionmatrix(2,2)+confusionmatrix(2,3)+confusionmatrix(3,2)+confusionmatrix(3,3);
            FP = confusionmatrix(2,1)+confusionmatrix(3,1);
            FN = confusionmatrix(1,2)+confusionmatrix(1,3);
            accu(i) = (TP+TN)/(TP+TN+FP+FN);
        else
            TP = confusionmatrix(2,2);
            TN = confusionmatrix(3,3);
            FP = confusionmatrix(3,2);
            FN = confusionmatrix(2,3);
            accu(i) = (TP+TN)/(TP+TN+FP+FN);        
        end
        
        %% Confusion matrix
        outputs = zeros(size(label_predictions,1),3);
        outputs(find(label_predictions==0),1) = 1;
        outputs(find(label_predictions==1),2) = 1;
        outputs(find(label_predictions==2),3) = 1;
        targets = zeros(size(respTest,1),3);
        targets(find(respTest==0),1) = 1;
        targets(find(respTest==1),2) = 1;
        targets(find(respTest==2),3) = 1;
%         figure
%         plotconfusion(targets',outputs');
        
        %% ROC curve
        XVals = 0:0.0001:1;
        if typeStudy == 1 % Artifacts Vs Glands
            if classifier == 1 | classifier == 2
                [at(:,i),bt(:,i),TT,AUC_test(i),opt] = perfcurve(respTest,score(:,1),0,'NegClass',[1,2],'XVals',XVals, 'UseNearest', 'off');
            else
                [at(:,i),bt(:,i),TT,AUC_test(i),opt] = perfcurve(respTest,score(1,:),0,'NegClass',[1,2],'XVals',XVals, 'UseNearest', 'off'); % ROC curve
            end
        elseif typeStudy == 2
            if classifier == 1 | classifier == 2
                [at(:,i),bt(:,i),TT,AUC_test(i),opt] = perfcurve(respTest,score(:,2),1,'NegClass',2, 'XVals',XVals, 'UseNearest', 'off');
            else
                [at(:,i),bt(:,i),TT,AUC_test(i),opt] = perfcurve(respTest,score(2,:),1,'NegClass',2, 'XVals',XVals, 'UseNearest', 'off');
            end
        end
        plot(at(:,i),bt(:,i),'LineWidth', 2);
        grid on; hold on;
        title('ROC curve', 'FontSize', 14); 
        xlabel('False positive rate', 'FontSize',15); ylabel('True positive rate', 'FontSize',15);
    end
    legend('ROC Fold-1', 'ROC Fold-2', 'ROC Fold-3', 'ROC Fold-4', 'ROC Fold-5', 'Location', 'southeast');

    % Save ROCS
    if classifier == 1 % quadratic SVM
        if typeStudy == 1
            save RocQSVMag.mat at bt;
        elseif typeStudy == 2
            save RocQSVMgg.mat at bt;
        end
    else % multilayerPerceptron
        if typeStudy == 1
            save RocNNag.mat at bt;
        elseif typeStudy == 2
            save RocNNgg.mat at bt;
        end
    end

else
    addpath('Deep_Learning_Results');
    
    %% Train & Prediction Deep Learning
    figure
    for i = 1:5
        i_str = num2str(i);
        predictions = read_csv_data_3classes(['predictionsTest' i_str 'Results.csv']);
        groundTruth = predictions.gt;

        %% Scores & Labels
        score_artifact = predictions.prob_artefact;
        score_grade3 = predictions.prob_grado3;
        score_healthy = predictions.prob_healthy;

        label_predictions = predictions.pred_label;

        %% ROC curve
        if typeStudy == 1
            [x,y,T,AUC_test(i), OPTROCPT] = perfcurve(groundTruth, score_artifact, 0, 'NegClass',[1,2]);
            id = find(x == OPTROCPT(1) & y == OPTROCPT(2));
            th = T(id);
            plot(x,y, 'LineWidth', 2);
            grid on; hold on;
            title('ROC curve', 'FontSize', 14); 
            xlabel('False positive rate', 'FontSize',15); ylabel('True positive rate', 'FontSize',15);
        elseif typeStudy == 2
            [x,y,T,AUC_test(i), OPTROCPT] = perfcurve(groundTruth, score_grade3, 1, 'NegClass', 2);
            id = find(x == OPTROCPT(1) & y == OPTROCPT(2));
            th = T(id);
            plot(x,y, 'LineWidth', 2);
            grid on; hold on;
            title('ROC curve', 'FontSize', 14); 
            xlabel('False positive rate', 'FontSize',15); ylabel('True positive rate', 'FontSize',15);       
        end
        
        %% Classperf
        if typeStudy == 1
            CP = classperf(groundTruth,label_predictions,'Positive',0,'Negative',[1,2]);
        elseif typeStudy == 2
            CP = classperf(groundTruth,label_predictions,'Positive',1,'Negative',2);
        end
        Sensitivity(i) = CP.Sensitivity;
        Specificity(i) = CP.Specificity;
        Ppv(i) = CP.PositivePredictiveValue;
        Npv(i) = CP.NegativePredictiveValue;
        F1Score(i) = 2.*Ppv(i).*Sensitivity(i)./(Ppv(i)+Sensitivity(i));
        Accuracy(i) = CP.CorrectRate;        
        confusionmatrix = CP.CountingMatrix;
        if typeStudy == 1
            TP = confusionmatrix(1,1);
            TN = confusionmatrix(2,2)+confusionmatrix(2,3)+confusionmatrix(3,2)+confusionmatrix(3,3);
            FP = confusionmatrix(2,1)+confusionmatrix(3,1);
            FN = confusionmatrix(1,2)+confusionmatrix(1,3);
            accu(i) = (TP+TN)/(TP+TN+FP+FN);
        else
            TP = confusionmatrix(2,2);
            TN = confusionmatrix(3,3);
            FP = confusionmatrix(3,2);
            FN = confusionmatrix(2,3);
            accu(i) = (TP+TN)/(TP+TN+FP+FN);        
        end
    end
    legend('ROC Fold-1', 'ROC Fold-2', 'ROC Fold-3', 'ROC Fold-4', 'ROC Fold-5', 'Location', 'southeast');
end


%% FINAL RESULTS
S = mean(Sensitivity);      s = std(Sensitivity);
E = mean(Specificity);      e = std(Specificity);
PPV = mean(Ppv);            ppv = std(Ppv);
NPV = mean(Npv);            npv = std(Npv);
F1_S = mean(F1Score);       f1_s = std(F1Score);
ACCURACY = mean(accu);      accuracy = std(accu);   
AUC = mean(AUC_test);       auc = std(AUC_test);
ACC = mean(Accuracy);       acc = std(Accuracy);

MEAN = [S; E; PPV; NPV;  F1_S; ACCURACY; AUC; ACC];
STD = [s; e; ppv; npv;  f1_s; accuracy; auc; acc];
list_parameters = {'Sensitivity (S)'; 'Specificity (E)'; 'Positive Predictive Value (PPV)';...
'Negative Predictive Value (NPV)'; 'F1Score'; 'Accuracy'; 'AUC'; 'Multiclass Accuracy'};
indicators_parameters = table(MEAN,STD, 'RowNames', list_parameters)

