%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% MAIN. END-TO-END %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all; 
set(groot,'ShowHiddenHandles','on'); delete(get(groot,'Children'));

rng('default');

%% INPUTS
% Type of prediction
classifier = 1;             % 1 = quadratic SVM - Machine learning
                            % 2 = Neural Network - Machine learning
                            % 3 = Deep Learning 
                            
%% Directories and name of image
%% Grade 3
addpath('samplesFromNewPatients\pathologicalPattern');
file = '18B0005477G_Block_Region0_3_3.jpg';
% file = '18B0005477G_Block_Region0_4_5.jpg';
% file = '18B0003894H_Block_Region2_2_2.jpg';   

%% Both
% addpath('samplesFromNewPatients\combinedPattern')
% file = '18B0003894H_Block_Region2_2_4.jpg';
% file = '18B0003894H_Block_Region8_2_6.jpg';
% file = '18B0004349A_Block_Region9_0_5.jpg';

%% Healthy
% addpath('samplesFromNewPatients\benignPattern');
% file = '18B0006623J_Block_Region1_1_7.jpg';
% file = '18B0006623A_Block_Region0_1_3.jpg';
% file = '18B0006623A_Block_Region0_0_0.jpg';

%% Paths
addpath('1_DataBase_Preparation');
addpath('2_Clustering');
addpath('3_Segmentation');
addpath('4_FeatureExtraction');
addpath('5_StatisticalAnalysis');
addpath('7_Classification');
addpath('Prediction');

%% Load Classifiers Models
load featuresSelected.mat;
if classifier == 1
    load qSVMclassifier_mc.mat
    trainedClass = SVMclassifier;
elseif classifier == 3
    load NNclassifier_mc.mat
    trainedClass = netTrain;   
end

%% RGB Image
img = imread(file);
img = imresize(img,0.5);

%% CLUSTERING to obtain masks and candidates to lumen
disp('Perfoming CLUSTERING stage...');
[~,~,~,cyan,s,H,magenta] = Chanel_color(img);
[mask_lumen, mask_black] = Lumen_mask(img, s);
[nucleous, nuclei_mask] = Nuclei_mask(img);
[cytoplasm_mask, stroma_mask, over2] = Masks(img, mask_black, nucleous, cyan);
% Resize
img = imresize(img,2);
mask_lumen = imresize(mask_lumen, 2); mask_black = imresize(mask_black,2);
nuclei_mask = imresize(nuclei_mask,2);
cytoplasm_mask = imresize(cytoplasm_mask,2);
stroma_mask = imresize(stroma_mask,2);
over2 = imresize(over2,2);
s = imresize(s,2); cyan = imresize(cyan,2); 
H = imresize(H,2); magenta = imresize(magenta,2);

%% SEGMENTATION of all candidates
disp('Perfoming SEGMENTATION stage...');
whole = tic;
[img_seg] = Segmentation(img,over2,nuclei_mask,mask_lumen,mask_black,stroma_mask);
time = toc(whole);
disp(['Time of segmentation = ' num2str(time)]);

%% Masked BBox  of the gland candidates
disp('STORING gland candidates separately...');
[candidatesBBox, candidatesMask, fileNormal] = obtainBBoxCandidates(file, img_seg, img, mask_black);
[imgMask, bb] = MaskBBox(fileNormal, candidatesMask, img);


%% MACHINE LEARNING
if classifier == 1 | classifier == 2 | classifier == 3 % Machine Learning

    % Matching features
    disp('EXTRACTING FEATURES of TEST...');
    featTest_mc = Extract_Features_Test_mc(imgMask, bb, img, mask_black, nuclei_mask, cytoplasm_mask, cyan);
    fT = bsxfun(@minus,featTest_mc, mu);
    featTest = bsxfun(@rdivide, fT, sigma);

    %% PREDICTION
    disp('Performing the PREDICTION phase...');
    if classifier == 1 | classifier == 2
        for i = 1:5
            [lab(:,i), score(:,:,i)] = predict(trainedClass{i},featTest);
        end
        labels = mode(lab,2);
    else
        for i = 1:5
            score = sim(trainedClass{i},featTest');
            label_predictions(i,:) = obtainLabelPredictionsNeuralNetwork(score);
        end
        labels = mean(label_predictions);
        labels_DL = [];
    end

%% DEEP LEARNING
else % Deep Learning
    modelsDir = 'Deep_Learning_Models';
    inputI = [];
    for i = 1:length(candidatesMask)
        candidateMask = candidatesMask{i};
        candidate = imresize(candidateMask, [224, 224]);
        inputI(:,:,:,i) = candidate;
    end
    for modelNumber = 1:5
        [labels_pred, probabilities] = Predict_DL(modelNumber, modelsDir, inputI);
    end
    labels = double(labels_pred);
end

%% VISUALISATION
magenta = imgaussfilt(magenta,2);
th = graythresh(magenta);
magenta = im2bw(magenta,th);
tissueWhite = imfill(magenta, 'holes');
tissue = uint8(tissueWhite.*double(img));

%% Highlight prediction glands
figure, imshow(tissue)
idStr = strfind(file,'_');
file(idStr) = '-';
title(file);
for jj = 1:length(candidatesMask)
    [x,y] = find(imgMask{jj} == 1);
    mask = zeros(size(img,1),size(img,2));
    for j = 1:length(x)
        mask(x(j)+bb(jj,2),y(j)+bb(jj,1)) = 1;
    end   
    if labels(jj) == 0
        [x,y] = find(mask == 1); 
        %hold on, visboundaries(mask, 'Color', 'k','LineWidth', 2);
    elseif labels(jj) == 1
        [x,y] = find(mask == 1);
        hold on, visboundaries(mask, 'Color', 'r', 'LineWidth', 2);
    elseif labels(jj) == 2
        [x,y] = find(mask == 1);
        hold on, visboundaries(mask, 'Color', 'g', 'LineWidth', 2);
    end
end
