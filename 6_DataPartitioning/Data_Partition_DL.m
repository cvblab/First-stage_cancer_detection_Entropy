%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% DATA PARTITIONING - DEEP LEARNING APPROACH %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc; rng('default');
set(groot,'ShowHiddenHandles','on'); delete(get(groot,'Children'));

addpath('..\1_DataBase_Preparation');

%% Load
style = 'Masks';

[Folds, Names] = Folds_construction(style);

n1_Artifacts = Names{3}{1};     n1_Gleason3 = Names{2}{1};      n1_Healthy = Names{1}{1};
n2_Artifacts = Names{3}{2};     n2_Gleason3 = Names{2}{2};      n2_Healthy = Names{1}{2};
n3_Artifacts = Names{3}{3};     n3_Gleason3 = Names{2}{3};      n3_Healthy = Names{1}{3};
n4_Artifacts = Names{3}{4};     n4_Gleason3 = Names{2}{4};      n4_Healthy = Names{1}{4};
n5_Artifacts = Names{3}{5};     n5_Gleason3 = Names{2}{5};      n5_Healthy = Names{1}{5};

f1_Artifacts = Folds{3}{1};     f1_Gleason3 = Folds{2}{1};      f1_Healthy = Folds{1}{1};
f2_Artifacts = Folds{3}{2};     f2_Gleason3 = Folds{2}{2};      f2_Healthy = Folds{1}{2};
f3_Artifacts = Folds{3}{3};     f3_Gleason3 = Folds{2}{3};      f3_Healthy = Folds{1}{3};
f4_Artifacts = Folds{3}{4};     f4_Gleason3 = Folds{2}{4};      f4_Healthy = Folds{1}{4};
f5_Artifacts = Folds{3}{5};     f5_Gleason3 = Folds{2}{5};      f5_Healthy = Folds{1}{5};

%% Conjuntos
F_artifacts = {f1_Artifacts, f2_Artifacts, f3_Artifacts, f4_Artifacts, f5_Artifacts};
N_artifacts = {n1_Artifacts, n2_Artifacts, n3_Artifacts, n4_Artifacts, n5_Artifacts};
F_grade3 = {f1_Gleason3, f2_Gleason3, f3_Gleason3, f4_Gleason3, f5_Gleason3};
N_grade3 = {n1_Gleason3, n2_Gleason3, n3_Gleason3, n4_Gleason3, n5_Gleason3};
F_healthy = {f1_Healthy, f2_Healthy, f3_Healthy, f4_Healthy, f5_Healthy};
N_healthy = {n1_Healthy, n2_Healthy, n3_Healthy, n4_Healthy, n5_Healthy};

%% Definition of waitbar
leng = (3200+3195+3000)*5;   ini = num2str(leng);
step = 1;
h = waitbar(0,ini,'Name','Organising candidates...', 'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0)

%% Modelos
for ii = 1:5
    % Images                  % Names                
    Fx_Art = F_artifacts;     Nx_Art = N_artifacts;
    Fx_Gra = F_grade3;        Nx_Gra = N_grade3;
    Fx_Hea = F_healthy;       Nx_Hea = N_healthy;
    
    Fx_Art(ii) = [];          Nx_Art(ii) = [];
    Fx_Gra(ii) = [];          Nx_Gra(ii) = [];
    Fx_Hea(ii) = [];          Nx_Hea(ii) = [];
    
    indVal = find(randperm(4) == 1);
    indTrain = [1,2,3,4];
    indTrain(indTrain == indVal) = [];
    i_str = num2str(ii);
    % Train
    directoryTrainArtifacts = ['..\Deep_Learning_' style '\train' i_str '\artifacts'];
    directoryTrainGrade3 = ['..\Deep_Learning_' style '\train' i_str '\grade3'];
    directoryTrainHealthy = ['..\Deep_Learning_' style '\train' i_str '\healthy'];
    idTrAr = Fx_Art(indTrain);        inTrAr = Nx_Art(indTrain);
    idTrGr = Fx_Gra(indTrain);        inTrGr = Nx_Gra(indTrain);
    idTrHe = Fx_Hea(indTrain);        inTrHe = Nx_Hea(indTrain); 
    foldTrainArtifacts = [idTrAr{1},idTrAr{2}, idTrAr{3}];
    nameTrainArtifacts = [inTrAr{1},inTrAr{2}, inTrAr{3}];
    foldTrainGrade3 = [idTrGr{1}, idTrGr{2}, idTrGr{3}];
    nameTrainGrade3 = [inTrGr{1},inTrGr{2}, inTrGr{3}];
    foldTrainHealthy = [idTrHe{1}, idTrHe{2}, idTrHe{3}];
    nameTrainHealthy = [inTrHe{1},inTrHe{2}, inTrHe{3}];
    for i = 1:length(foldTrainArtifacts)
        ima = foldTrainArtifacts{i};
        name = nameTrainArtifacts{i};
        save_dir = fullfile(directoryTrainArtifacts,name);
        imwrite(ima,save_dir);
        
        waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
        step = step+1;
    end
    for j = 1:length(foldTrainGrade3)
        ima = foldTrainGrade3{j};
        name = nameTrainGrade3{j};
        save_dir = fullfile(directoryTrainGrade3,name);
        imwrite(ima,save_dir);
        
        waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
        step = step+1;
    end
    for k = 1:length(foldTrainHealthy)
        ima = foldTrainHealthy{k};
        name = nameTrainHealthy{k};
        save_dir = fullfile(directoryTrainHealthy,name);
        imwrite(ima,save_dir);
        
        waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
        step = step+1;
    end
    
    % Val
    directoryValArtifacts = ['..\Deep_Learning_' style '\val' i_str '\artifacts'];
    directoryValGrade3 = ['..\Deep_Learning_' style '\val' i_str '\grade3'];
    directoryValHealthy = ['..\Deep_Learning_' style '\val' i_str '\healthy'];
    idVlAr = Fx_Art(indVal);          inVlAr = Nx_Art(indVal);
    idVlGr = Fx_Gra(indVal);          inVlGr = Nx_Gra(indVal);
    idVlHe = Fx_Hea(indVal);          inVlHe = Nx_Hea(indVal);  
    foldValArtifacts = idVlAr{1};      nameValArtifacts = inVlAr{1};     
    foldValGrade3 = idVlGr{1};         nameValGrade3 = inVlGr{1};
    foldValHealthy = idVlHe{1};        nameValHealthy = inVlHe{1};
    for i = 1:length(foldValArtifacts)
        ima = foldValArtifacts{i};
        name = nameValArtifacts{i};
        save_dir = fullfile(directoryValArtifacts,name);
        imwrite(ima,save_dir);
        
        waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
        step = step+1;
    end
    for j = 1:length(foldValGrade3)
        ima = foldValGrade3{j};
        name = nameValGrade3{j};
        save_dir = fullfile(directoryValGrade3,name);
        imwrite(ima,save_dir);
        
        waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
        step = step+1;
    end
    for k = 1:length(foldValHealthy)
        ima = foldValHealthy{k};
        name = nameValHealthy{k};
        save_dir = fullfile(directoryValHealthy,name);
        imwrite(ima,save_dir);
        
        waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
        step = step+1;
    end
    % Test
    directoryTestArtifacts = ['..\Deep_Learning_' style '\test' i_str '\artifacts'];
    directoryTestGrade3 = ['..\Deep_Learning_' style '\test' i_str '\grade3'];
    directoryTestHealthy = ['..\Deep_Learning_' style '\test' i_str '\healthy'];
    foldTestArtifacts = F_artifacts(ii);            nameTestArtifacts = N_artifacts(ii);
    foldTestArtifacts = foldTestArtifacts{1};       nameTestArtifacts = nameTestArtifacts{1};
    foldTestGrade3 = F_grade3(ii);                  nameTestGrade3 = N_grade3(ii);
    foldTestGrade3 = foldTestGrade3{1};             nameTestGrade3 = nameTestGrade3{1};
    foldTestHealthy = F_healthy(ii);                nameTestHealthy = N_healthy(ii);
    foldTestHealthy = foldTestHealthy{1};           nameTestHealthy = nameTestHealthy{1};
    for i = 1:length(foldTestArtifacts)
        ima = foldTestArtifacts{i};
        name = nameTestArtifacts{i};
        save_dir = fullfile(directoryTestArtifacts,name);
        imwrite(ima,save_dir);
        
        waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
        step = step+1;
    end
    for j = 1:length(foldTestGrade3)
        ima = foldTestGrade3{j};
        name = nameTestGrade3{j};
        save_dir = fullfile(directoryTestGrade3,name);
        imwrite(ima,save_dir);
        
        waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
        step = step+1;
    end
    for k = 1:length(foldTestHealthy)
        ima = foldTestHealthy{k};
        name = nameTestHealthy{k};
        save_dir = fullfile(directoryTestHealthy,name);
        imwrite(ima,save_dir);
        
        waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
        step = step+1;
    end
end
delete(h)
