function [name_hhcc_parts] = Name_HHCC_struct(folder, response)
 

%% 5 Sessions
load IMG_HHCC_HEALTHY.mat;
load IMG_HHCC_G3.mat;
load IMG_HHCC_G4.mat;
load IMG_HHCC_G5.mat;

%% Name of Slides
fileVect = dir([folder '\*.jpg']);
fileName = {fileVect.name};
for i = 1:length(fileName)
    file = fileName{i};
    name{i} = file(1:10); % All names
end
names = unique(name');
for i = 1:length(names)
    nam(i,:) = names{i}; % Only some names
end

%% Names images and Clinical Historials
% Healthy.
imgHeal = IMG_HHCC_H.IMG;           hhccHeal = IMG_HHCC_H.HHCC; % Totales
imgHeale1e2 = imgHeal(1:17,:);      hhccHeale1e2 = hhccHeal(1:17,:); % Entregas 1 y 2

% Grade 3
img3 = IMG_HHCC.IMG;                hhcc3 = IMG_HHCC.HHCC;     % Totales
img3e1e2 = img3(1:18,:);            hhcc3e1e2 = hhcc3(1:18,:); % Entregas 1 y 2

% Grade 4. 
img4 = IMG_HHCC_G4.IMG;             hhcc4 = IMG_HHCC_G4.HHCC;  % Totales
img4e1e2 = img4(1:15,:);            hhcc4e1e2 = hhcc4(1:15,:); % Entregas 1 y 2

% Grade 5.
img5 = IMG_HHCC_G5.IMG;             hhcc5 = IMG_HHCC_G5.HHCC;  % Totales
img5e1e2 = img5(1:10,:);            hhcc5e1e2 = hhcc5(1:10,:); % Entregas 1 y 2

if response == 1    % Healthy
    nameImages= imgHeale1e2;    nameHHCCs = hhccHeale1e2;
elseif response == 2 % Grades 3, 4, 5. Sessions 1 y 2
    nameImages= [img3e1e2;img4e1e2;img5e1e2];    nameHHCCs = [hhcc3e1e2;hhcc4e1e2;hhcc5e1e2];
else
    disp('ERROR');
end


%% Verification
pos = [];
pos2 = [];
for i = 1:size(nameImages,1)
    for j = 1:size(nam,1)
        resp = isequal(nam(j,:), nameImages(i,:));
        if resp == 1
            pos = [pos; j];
            pos2 = [pos2; i];
        end
    end
end

%% Construction of STRUCT
hc = zeros(size(pos));
for k = 1:length(pos)
    hc(pos(k)) = nameHHCCs(pos2(k));
end

name_hhcc_parts = struct('name',nam,'hhcc',hc);

