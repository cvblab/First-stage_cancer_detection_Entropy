%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% DATA PARTITIONING PER GLANDS FROM THE PATIENTS %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Folds, Names] = Folds_Partition

for k = 1:3

    if k == 1
        read_folder = '..\FINAL_BoundingBox\healthyMasks'; 
    elseif k == 2
        read_folder = '..\FINAL_BoundingBox\grade3Masks'; 
    else
        read_folder = '..\FINAL_BoundingBox\artifactsMasks'; 
    end
    response = k;

    %% Names of all images
    fileVect = dir([read_folder, '\*.jpg']);
    fileName = {fileVect.name};

    if response <3
        %% Catch the HHCC of each fold and the image's names corresponding to their HHCC
        [historial_clinico, fold1, fold2, fold3, fold4, fold5] = Folds_construction(read_folder, response);
        f1 = [];    I1 = 1;     name1 = [];
        f2 = [];    I2 = 1;     name2 = [];
        f3 = [];    I3 = 1;     name3 = [];
        f4 = [];    I4 = 1;     name4 = [];
        f5 = [];    I5 = 1;     name5 = [];

        %% Definition of waitbar
        leng = length(fileName);
        ini = num2str(leng);
        step = 1;
        h = waitbar(0,ini,'Name','Organizing images in folds...', 'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
        setappdata(h,'canceling',0)

        %% Assign each image to their corresponding fold
        folds = {fold1;fold2;fold3;fold4;fold5};
        fs = {f1; f2; f3; f4; f5};
        Is = {I1; I2; I3; I4; I5};
        names = {name1; name2; name3; name4; name5};
        for i = 1:length(fileName)
            file = fileName{i}(1:10);
            p = 0;
            for m = 1:5
                foldStudy = folds{m};
                for n = 1:length(foldStudy)
                    if file == foldStudy{n}
                        p = 1;
                        ima = imread([read_folder '\' fileName{i}]);
                        fs{m}{Is{m}} = ima;
                        names{m}{Is{m}} = fileName{i};
                        Is{m} = Is{m}+1;
                        break;
                    else
                        continue;
                    end
                end
                if p == 1
                    break;
                else
                    continue;
                end
            end 

        waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
        step = step+1;
        end
        delete(h);

        if response == 1
            f1_Healthy = fs{1};    n1_Healthy = names{1}; 
            f2_Healthy = fs{2};    n2_Healthy = names{2};
            f3_Healthy = fs{3};    n3_Healthy = names{3};
            f4_Healthy = fs{4};    n4_Healthy = names{4};
            f5_Healthy = fs{5};    n5_Healthy = names{5};

            Folds{k} = {f1_Healthy; f2_Healthy; f3_Healthy; f4_Healthy; f5_Healthy};
            Names{k} = {n1_Healthy; n2_Healthy; n3_Healthy; n4_Healthy; n5_Healthy};

        elseif response == 2
            f1_Gleason3 = fs{1};   n1_Gleason3 = names{1};
            f2_Gleason3 = fs{2};   n2_Gleason3 = names{2};
            f3_Gleason3 = fs{3};   n3_Gleason3 = names{3};
            f4_Gleason3 = fs{4};   n4_Gleason3 = names{4};
            f5_Gleason3 = fs{5};   n5_Gleason3 = names{5};

            Folds{k} = {f1_Gleason3; f2_Gleason3; f3_Gleason3; f4_Gleason3; f5_Gleason3};
            Names{k} = {n1_Gleason3; n2_Gleason3; n3_Gleason3; n4_Gleason3; n5_Gleason3};

        end

    elseif response == 3

        %% Definition of waitbar
        leng = length(fileName);
        ini = num2str(leng);
        step = 1;
        h = waitbar(0,ini,'Name','Organizing images in folds...', 'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
        setappdata(h,'canceling',0)

        %% Artifacts
        rng('default');
        indx = crossvalind('Kfold',length(fileName),5);
        % Artifacts1
        id1 = find(indx==1);
        for i = 1:length(id1)
            file = fileName{id1(i)};
            ima = imread([read_folder '\' file]);
            f1_Artifacts{i} = ima;
            n1_Artifacts{i} = file;
            waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
            step = step+1;
        end
        % Artifacts2
        id2 = find(indx==2);
        for i = 1:length(id2)
            file = fileName{id2(i)};
            ima = imread([read_folder '\' file]);
            f2_Artifacts{i} = ima;
            n2_Artifacts{i} = file;
            waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
            step = step+1;
        end
        % Artifacts3
        id3 = find(indx==3);
        for i = 1:length(id3)
            file = fileName{id3(i)};
            ima = imread([read_folder '\' file]);
            f3_Artifacts{i} = ima;
            n3_Artifacts{i} = file;
            waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
            step = step+1;
        end
        % Artifacts4
        id4 = find(indx==4);
        for i = 1:length(id4)
            file = fileName{id4(i)};
            ima = imread([read_folder '\' file]);
            f4_Artifacts{i} = ima;
            n4_Artifacts{i} = file;
            waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
            step = step+1;
        end
        % Artifacts5
        id5 = find(indx==5);
        for i = 1:length(id5)
            file = fileName{id5(i)};
            ima = imread([read_folder '\' file]);
            f5_Artifacts{i} = ima;
            n5_Artifacts{i} = file;
            waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
            step = step+1;
        end
        delete(h);

        Folds{k} = {f1_Artifacts; f2_Artifacts; f3_Artifacts; f4_Artifacts; f5_Artifacts};
        Names{k} = {n1_Artifacts; n2_Artifacts; n3_Artifacts; n4_Artifacts; n5_Artifacts};

    end
end

