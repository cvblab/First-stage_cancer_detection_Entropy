%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% DATA PARTITIONING PER PATIENT AT THE SLIDE LEVEL %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [glandsPerFold, namesPerFold] = Folds_construction(style)

%% Folders
for d = 1:3
    if d == 1
        folder = ['..\FINAL_BoundingBox\healthy' style]; 
    elseif d == 2
        folder = ['..\FINAL_BoundingBox\grade3' style]; 
    else
        folder = ['..\FINAL_BoundingBox\artifacts' style]; 
    end
    response = d;
    fileVect = dir([folder '\*.jpg']);
    fileName = {fileVect.name};
    
    %% Definition of waitbar
    leng = length(fileName);
    ini = num2str(leng);
    step = 1;
    h = waitbar(0,ini,'Name','Organizing images in folds...', 'CreateCancelBtn',...
        'setappdata(gcbf,''canceling'',1)');
    setappdata(h,'canceling',0)
    
    %% Parameters
    num_folds = 5;
    
    if response < 3
        %% Clinical history
        [name_hhcc_parts] = Name_HHCC_struct(folder,response);
        hhcc = name_hhcc_parts.hhcc;
        uniqueHHCC = unique(hhcc);

        %% Name and number of glands candidates
        for i = 1:length(fileName)
            file = fileName{i};
            all_names{i} = file(1:10); % Todos los nombres de las imágenes 
        end
        table = tabulate(all_names); % Nombre y número de muestras  de imágenes diferentes
        nameSamples = table(:,1);

        %% Slides per patient
        for k = 1:length(nameSamples)
            S(k) = struct('patients',hhcc(k), 'slides',nameSamples{k});
        end

        %% Glands per Slide
        for m = 1:length(fileName)
            file = fileName{m};
            fileNameinitial = file(1:10);
            for n = 1:length(S)
                if isequal(fileNameinitial,S(n).slides)
                    st{n,m} = file;
                    break;
                end
            end
        end
        for p = 1:length(nameSamples)
            glandPerSlide{p} = st(p,~cellfun(@isempty,st(p,:)));  
        end
        glandPerSlide = glandPerSlide';

        %% Algorithm (PATIENTS per FOLD)
        GMN = length(all_names)/num_folds; % número medio de glándulas por cada fold
        uni = uniqueHHCC; % pacientes únicos
        j = 1;
        while j <= 5
            cach = [];
            x = [];
            for i = 1:length(uni)
                id = find(uni(i) == hhcc);
                nG = [];
                for k = 1:length(id)
                    nG(k) = length(glandPerSlide{id(k)});
                end
                numGlands = sum(nG); 
                cach = [cach numGlands];
                sumGlands = sum(cach);     
                if numGlands>GMN 
                    x = [x uni(i)]; 
                end        
                if sumGlands<GMN 
                    hc(j,i) = uni(i); 
                else
                    if j == num_folds
                        break;
                    end
                    cach(end-length(numGlands)+1:end) = []; 
                end
            end
            if ~isempty(x)
                for k = 1:length(x)   
                    hc(j+k,j+k) = x(k);
                end
            else
                k = 0;
            end
            uni(ismember(uni,hc)) = []; 
            j = j+k+1;
        end
        FP = [];
        for ii = 1:size(hc,1)
            idx = find(hc(ii,:) ~= 0);
            FP{ii} = hc(ii,idx);    
        end

        %% Algorithm (SLIDES per FOLD)
        k = 1;
        FS = [];
        for f = 1:num_folds
            patients = FP{f};
            for p = 1:length(patients)
                id = find(hhcc == patients(p));
                for i = 1:length(id)
                    FS{k,p} = nameSamples{id(i)};
                    p = p+1;
                end
            end
            k = k+1;
        end

        %% Algorithm (GLANDS per FOLD)
        FN = [];
        for f = 1:num_folds
            ss = 1;
            fs = FS(f,~cellfun(@isempty,FS(f,:)));
            for s = 1:length(fs)
                vectCell = strfind(nameSamples,fs(s));
                id = find(~cellfun(@isempty,vectCell));
                for b = 1:length(glandPerSlide{id})
                    gland = imread([folder '\' glandPerSlide{id}{b}]);
                    FG{f,ss} = gland;
                    FN{f,ss} = glandPerSlide{id}{b};
                    ss = ss+1;
                    waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
                    step = step+1;
                end
            end
        end

        name1 = FN(1,~cellfun(@isempty,FN(1,:)));
        name2 = FN(2,~cellfun(@isempty,FN(2,:)));
        name3 = FN(3,~cellfun(@isempty,FN(3,:)));
        name4 = FN(4,~cellfun(@isempty,FN(4,:)));
        name5 = FN(5,~cellfun(@isempty,FN(5,:)));
        fold1 = FG(1,~cellfun(@isempty,FG(1,:))); 
        fold2 = FG(2,~cellfun(@isempty,FG(2,:))); 
        fold3 = FG(3,~cellfun(@isempty,FG(3,:))); 
        fold4 = FG(4,~cellfun(@isempty,FG(4,:))); 
        fold5 = FG(5,~cellfun(@isempty,FG(5,:))); 
             
    else %% Artifacts
        
        rng('default');
        indx = crossvalind('Kfold',length(fileName),5);
        name = fileName;   
        name1 = name(indx==1);
        name2 = name(indx==2);
        name3 = name(indx==3);
        name4 = name(indx==4);
        name5 = name(indx==5);
        
        id1 = find(indx==1);
        for i = 1:length(id1)
            file = fileName{id1(i)};
            ima = imread([folder '\' file]);
            fold1{i} = ima;
            waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
            step = step+1;
        end
        id2 = find(indx==2);
        for i = 1:length(id2)
            file = fileName{id2(i)};
            ima = imread([folder '\' file]);
            fold2{i} = ima;
            waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
            step = step+1;
        end
        id3 = find(indx==3);
        for i = 1:length(id3)
            file = fileName{id3(i)};
            ima = imread([folder '\' file]);
            fold3{i} = ima;
            waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
            step = step+1;
        end
        id4 = find(indx==4);
        for i = 1:length(id4)
            file = fileName{id4(i)};
            ima = imread([folder '\' file]);
            fold4{i} = ima;
            waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
            step = step+1;
        end
        id5 = find(indx==5);
        for i = 1:length(id5)
            file = fileName{id5(i)};
            ima = imread([folder '\' file]);
            fold5{i} = ima;
            waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
            step = step+1;
        end
    end
        
    glandsPerFold{d} = {fold1; fold2; fold3; fold4; fold5};
    namesPerFold{d} = {name1; name2; name3; name4; name5};
    delete(h);
    clearvars -except glandsPerFold namesPerFold style;
end