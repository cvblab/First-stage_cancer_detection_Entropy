%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% DATA PARTITIONING - MACHINE LEARNING APPROACH %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng('default');
clear; close all; clc; 
set(groot,'ShowHiddenHandles','on'); delete(get(groot,'Children'));

style = 'Masks';

addpath('..\1_DataBase_Preparation');
addpath('..\5_StatisticalAnalysis');
load featuresSelected.mat; 

%% Folds Partition
[~, namesPerFold] = Folds_construction(style);

%% Load Names
n1_Artifacts = namesPerFold{3}{1};     n1_Gleason3 = namesPerFold{2}{1};      n1_Healthy = namesPerFold{1}{1};
n2_Artifacts = namesPerFold{3}{2};     n2_Gleason3 = namesPerFold{2}{2};      n2_Healthy = namesPerFold{1}{2};
n3_Artifacts = namesPerFold{3}{3};     n3_Gleason3 = namesPerFold{2}{3};      n3_Healthy = namesPerFold{1}{3};
n4_Artifacts = namesPerFold{3}{4};     n4_Gleason3 = namesPerFold{2}{4};      n4_Healthy = namesPerFold{1}{4};
n5_Artifacts = namesPerFold{3}{5};     n5_Gleason3 = namesPerFold{2}{5};      n5_Healthy = namesPerFold{1}{5};

%% Definition of waitbar
leng = 5;
ini = num2str(leng);
step = 1;
h = waitbar(0,ini,'Name','Organizing SETS...', 'CreateCancelBtn',...
    'setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0)

%% SET 1
setTest1 = [n1_Artifacts n1_Gleason3 n1_Healthy];
indexTest1 = [];
for i = 1:length(finalMatrix)
    file = files{i};
    fileT = [];
    for j = 1:length(setTest1)
        fileT = setTest1{j};
        if isequal(file,fileT)
            indexTest1 = [indexTest1 i];
            break;
        end
    end
end
featTest1 = finalMatrix(indexTest1,:);
featTrain1 = finalMatrix;
featTrain1(indexTest1,:) = [];
waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
step = step+1;
            
%% SET 2
setTest2 = [n2_Artifacts, n2_Gleason3, n2_Healthy];
indexTest2 = [];
for i = 1:length(finalMatrix)
    file = files{i};
    fileT = [];
    for j = 1:length(setTest2)
        fileT = setTest2{j};
        if isequal(file,fileT)
            indexTest2 = [indexTest2 i];
            break;
        end
    end
end
featTest2 = finalMatrix(indexTest2,:);
featTrain2 = finalMatrix;
featTrain2(indexTest2,:) = [];
waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
step = step+1;

%% SET 3
setTest3 = [n3_Artifacts, n3_Gleason3, n3_Healthy];
indexTest3 = [];
for i = 1:length(finalMatrix)
    file = files{i};
    fileT = [];
    for j = 1:length(setTest3)
        fileT = setTest3{j};
        if isequal(file,fileT)
            indexTest3 = [indexTest3 i];
            break;
        end
    end
end
featTest3 = finalMatrix(indexTest3,:);
featTrain3 = finalMatrix;
featTrain3(indexTest3,:) = [];
waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
step = step+1;

%% SET 4
setTest4 = [n4_Artifacts, n4_Gleason3, n4_Healthy];
indexTest4 = [];
for i = 1:length(finalMatrix)
    file = files{i};
    fileT = [];
    for j = 1:length(setTest4)
        fileT = setTest4{j};
        if isequal(file,fileT)
            indexTest4 = [indexTest4 i];
            break;
        end
    end
end
featTest4 = finalMatrix(indexTest4,:);
featTrain4 = finalMatrix;
featTrain4(indexTest4,:) = [];
waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
step = step+1;

%% SET 5
setTest5 = [n5_Artifacts, n5_Gleason3, n5_Healthy];
indexTest5 = [];
for i = 1:length(finalMatrix)
    file = files{i};
    fileT = [];
    for j = 1:length(setTest5)
        fileT = setTest5{j};
        if isequal(file,fileT)
            indexTest5 = [indexTest5 i];
            break;
        end
    end
end
featTest5 = finalMatrix(indexTest5,:);
featTrain5 = finalMatrix;
featTrain5(indexTest5,:) = [];
waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
step = step+1;
delete(h)

%% SETS ORGANIZATION
TEST = {featTest1 featTest2 featTest3 featTest4 featTest5};
TRAIN = {featTrain1 featTrain2 featTrain3 featTrain4 featTrain5};
SETS = struct('TEST',TEST, 'TRAIN', TRAIN);

save sets_v1.mat SETS;

