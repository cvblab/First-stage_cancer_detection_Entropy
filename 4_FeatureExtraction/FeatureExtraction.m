%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% FEATURE EXTRACTION STAGE %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;
set(groot,'ShowHiddenHandles','on'); delete(get(groot,'Children'));

addpath('..\2_Clustering');
rng('default');

%% Input parameter
sizeOriginalImage = 1024;
Neigh = 8; % LBPs
Radius = 1; % LBPs

%% Feature Extraction
folderImages = 'C:\Users\jogarpa7\Desktop\TFM\All_Images_1024x1024x3';
type = {'artifacts'; 'grade3'; 'healthy'};

for k = 1:3
    folderNormal = ['C:\Users\jogarpa7\Desktop\TFM\FINAL_BoundingBox\' [type{k} 'Normal']];
    fileVectNormal = dir([folderNormal '\*.jpg']);
    fileNameNormal = {fileVectNormal.name};
    folderMask = ['C:\Users\jogarpa7\Desktop\TFM\FINAL_BoundingBox\' [type{k} 'Masks']];

    %% Definition of waitbar
    leng = length(fileNameNormal);
    ini = num2str(leng);
    step = 1;
    h = waitbar(0,ini,'Name','Extracting Features...', 'CreateCancelBtn',...
        'setappdata(gcbf,''canceling'',1)');
    setappdata(h,'canceling',0)

    %% EXTRACTION OF FEATURES
    featMatrix = [];
    for i = 1:length(fileNameNormal)
        % bbox NORMAL
        fileNormal = fileNameNormal{i};
        imgRGBnormal = imread([folderNormal '\' fileNormal]);

        % bbox MASK
        imgRGBmask = imread([folderMask '\' fileNormal]);  
        r = double(imgRGBmask(:,:,1));
        g = double(imgRGBmask(:,:,2));
        b = double(imgRGBmask(:,:,3));
        RGBs = r+g+b; RGBs = RGBs./max(unique(RGBs));
        th = 0.1;
        bw = im2bw(RGBs,th);
        imgMask = bwareaopen(bw,20,4);
        imgMask = imfill(logical(imgMask), 'holes');

        [~, numimg] = bwlabel(imgMask);
        if numimg>1
            disp(['Error. + de 1 gland in ' num2str(i) ' and class: ' num2str(k)]);
            continue;
        end
        if max(unique(imgMask)) == 0
            disp(['Error. Black Gland in ' num2str(i) ' and class: ' num2str(k)]);
            continue;
        end

    %% %%%%%%%%%%%%%%%%%%%
    %%% SHAPE FEATURES %%% (20 FEATURES)
    %%%%%%%%%%%%%%%%%%%%%%
    shapeFeatures = [];
    shapeNames = {};

        % GLAND (10)
        propGland = regionprops(imgMask, 'Area', 'BoundingBox', 'Centroid', 'ConvexArea', ...
            'Eccentricity', 'EquivDiameter', 'Extent', 'Orientation', 'Perimeter', 'Solidity');
        areaGland = cat(1,propGland.Area);
        convexAreaGland = cat(1,propGland.ConvexArea);
        eccentGland = cat(1,propGland.Eccentricity);
        equivDiamGland = cat(1,propGland.EquivDiameter);
        extentGland = cat(1,propGland.Extent);
        orientationGland = cat(1, propGland.Orientation);
        perimeterGland = sum(cat(1, propGland.Perimeter));
        solidityGland = sum(cat(1, propGland.Solidity));
            rG=sqrt(areaGland/pi);
        if areaGland ~= 0
            roundnessGland = rG*perimeterGland/areaGland;
            compactnessGland = perimeterGland/sqrt(areaGland);
        else
            roundnessGland = 0;
            compactnessGland = 0;
        end
        % No features
        ind = find(fileNormal == 'x', 1, 'last');
        indx = find(fileNormal == 'y', 1, 'last');
        bb = cat(1,propGland.BoundingBox);
        bb(1) = str2num(fileNormal(ind+1:indx-2)); 
        bb(2) = str2num(fileNormal(indx+1:end-4));
        bb(3) = bb(1)+size(imgMask,2);
        bb(4) = bb(2)+size(imgMask,1);
        if bb(1) < 1
            imgMask(:,1) = [];
        end
        if bb(2) < 1
            imgMask(1,:) = [];
        end
        if bb(3) >sizeOriginalImage
            imgMask(:,end) = [];
        end
        if bb(4) > sizeOriginalImage
            imgMask(end,:) = [];
        end
        imgMask = bwareaopen(imgMask,20);
        bb(bb<1) = 1;
        bb(bb>sizeOriginalImage) = sizeOriginalImage; % size of original image
        cent = cat(1,propGland.Centroid);

        % LUMEN (10)
        % Complet image
        ind = find(fileNormal == '_', 2, 'last');
        fileComplet = fileNormal(1:ind(1)-1);
        completImage = imread([folderImages '\' fileComplet '.jpg']);
        completImage = imresize(completImage, 0.5);
        [~,~,~,cyan,s,~,~] = Chanel_color(completImage);
        [~, lumens] = Lumen_mask(completImage, s);
        lumens = imresize(lumens,2);
        lumens = imfill(lumens,'holes');
        % Bbox lumen
        lumeni_mask = lumens(bb(2)+1:bb(4),bb(1)+1:bb(3));
        lumeni_mask = bwareaopen(lumeni_mask, 20);
        if max(unique(lumeni_mask)) == 0
            disp(['Error. Black LUMEN in ' num2str(i) ' and class: ' num2str(k)]);
            continue;
        end
        if ~isequal(size(lumeni_mask),size(imgMask))
            disp(['Error. DIFFERENT Bbox LUMEN and Bbox GLAND in ' num2str(i) ' and class: ' num2str(k)]);
            continue;
        else
            lumeni = lumeni_mask;
        end
        lumenMask = lumeni.*imgMask;

        [lum_img, num_lum] = bwlabel(lumenMask);
        if num_lum>1
            p1 = round(cent(2));
            p2 = round(cent(1));
            dist = [];
            for ll = 1:num_lum
                elem = lum_img == ll;
                pr = regionprops(elem, 'Centroid');
                centr = cat(1,pr.Centroid);
                c1 = round(centr(2));
                c2 = round(centr(1));
                dist(ll) = pdist2(centr,[p2,p1], 'Euclidean');
            end
            id_dist = find(dist == min(dist));
            lumenMask = lum_img == id_dist;
        end

        % Features Lumen
        propLumen = regionprops(lumenMask, 'Area', 'BoundingBox', 'Centroid', 'ConvexArea',...
            'Eccentricity', 'EquivDiameter', 'Extent', 'Orientation', 'Perimeter', 'Solidity');           
        areaLumen = cat(1,propLumen.Area); 
        convexAreaLumen = cat(1,propLumen.ConvexArea);
        eccentLumen = cat(1,propLumen.Eccentricity);
        equivDiamLumen = cat(1,propLumen.EquivDiameter);
        extentLumen = cat(1,propLumen.Extent);
        orientationLumen = cat(1, propLumen.Orientation);
        perimeterLumen = cat(1, propLumen.Perimeter);
        solidityLumen = cat(1, propLumen.Solidity);
            rG=sqrt(areaLumen/pi);
        if ~isempty(areaLumen)
            roundnessLumen = rG*perimeterLumen/areaLumen;
            compactnessLumen = perimeterLumen/sqrt(areaLumen);
        else
            roundnessLumen = 0;
            compactnessLumen = 0;
        end


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% CONTEXTUAL FEATURES related to shape variables %%% (12 FEATURES)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Nuclei (8)
        [~, nuclei_mask] = Nuclei_mask(completImage);
        nuclei_mask = imresize(nuclei_mask,2);
        nuc_in_BB = nuclei_mask(bb(2)+1:bb(4),bb(1)+1:bb(3)); % Núcleos en la BB
        if ~isequal(size(nuc_in_BB),size(imgMask))
            disp(['Error. DIFFERENT Bbox NUCLEI and Bbox GLAND in ' num2str(i) ' and class: ' num2str(k)]);
            continue;
        else
            nucMask = nuc_in_BB;
        end
        nuc_in_mask = nucMask.*imgMask; % Nuclei in the gland mask

        propNuclei = regionprops(imgMask, 'Area');
        areaGland = cat(1, propNuclei.Area);
        sizeBB = [size(imgMask,1),size(imgMask,2)];       

        % Features Nuclei
        [~, numNucBB] = bwlabel(nuc_in_BB); % número de núcleos en la BB
        numNucBBprop = numNucBB/(sizeBB(1)*sizeBB(2)); % proporcion de núcleos en la BB
        [~, numNucGland] = bwlabel(nuc_in_mask); % número de núcleos en la MASK
        numNucGlandProp = numNucGland/areaGland; % proporción de núcleos en la MASK

        pixNucBB = sum(nuc_in_BB(:)); % pixeles de núcleos en la BB
        pixNucBBprop = pixNucBB/(sizeBB(1)*sizeBB(2)); % proporción de píxeles en la BB
        pixNucGland = sum(nuc_in_mask(:)); %  píxeles de núcleos en la MASK
        pixNucGlandProp = pixNucGland/areaGland;

        % Cytoplasm(4)
        lumi = imresize(lumens,0.5);    nucli = imresize(nuclei_mask, 0.5);
        [cytoplasm_mask, stroma_mask] = Masks(completImage, lumi, nucli, cyan);
        cytoplasm_mask = imresize(cytoplasm_mask, 2);
        cyto_in_BB = cytoplasm_mask(bb(2)+1:bb(4), bb(1)+1:bb(3)); % Cyto in BB
        if ~isequal(size(cyto_in_BB),size(imgMask))
            disp(['Error. DIFFERENT Bbox CYTO and Bbox GLAND in ' num2str(i) ' and class: ' num2str(k)]);
            continue;
        else
            cytoplasmMask = cyto_in_BB;
        end
        cyto_in_Gland = cytoplasmMask.*imgMask; % Cyto in Gland

        % Features 
        pixCytoBB = sum(cyto_in_BB(:)); % píxeles de cytoplasma en la BB
        pixCytoBBprop = pixCytoBB/(sizeBB(1)*sizeBB(2)); % proporción de píxeles en la BB
        pixCytoGland = sum(cyto_in_Gland(:)); % píxeles en la MASK
        pixCytoGlandProp = pixCytoGland/areaGland;

    % FINAL Shape FEATURES
    shapeFeatures = [areaGland, convexAreaGland, eccentGland, equivDiamGland,...
        extentGland, orientationGland, perimeterGland, solidityGland, roundnessGland, compactnessGland, ...
        areaLumen, convexAreaLumen, eccentLumen, equivDiamLumen, extentLumen, orientationLumen, ... 
        perimeterLumen, solidityLumen, roundnessLumen, compactnessLumen, ...
        numNucBB, numNucBBprop, numNucGland, pixNucBB, pixNucBBprop, pixNucGland, numNucGlandProp, pixNucGlandProp, ...
        pixCytoBB, pixCytoBBprop, pixCytoGland, pixCytoGlandProp];

    shapeNames = {'areaGland', 'convexAreaGland', 'eccentGland', 'equivDiamGland',...
        'extentGland', 'orientationGland', 'perimeterGland', 'solidityGland', 'roundnessGland', 'compactnessGland', ...
        'areaLumen', 'convexAreaLumen', 'eccentLumen', 'equivDiamLumen', 'extentLumen', 'orientationLumen', ... 
        'perimeterLumen', 'solidityLumen', 'roundnessLumen', 'compactnessLumen', ...
        'numNucBB', 'numNucBBprop', 'numNucGland', 'pixNucBB', 'pixNucBBprop', 'pixNucGland', 'numNucGlandProp', 'pixNucGlandProp', ...
        'pixCytoBB', 'pixCytoBBprop', 'pixCytoGland', 'pixCytoGlandProp'};


    %% %%%%%%%%%%%%%%%%%%%%%%%%
    %%% CONTEXTUAL FEATURES %%% (8 FEATURES)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    contextualFeatures = [];
    contextualNames = {};               
            ring = logical(imgMask-lumenMask);
            nucRing = logical(ring.*nucMask);
            cytoRing = logical(ring.*cytoplasmMask);
            lumRing = logical(ring.*lumenMask);
            white = logical(nucRing+cytoRing+lumRing);
        pixNCL = length(find(white==1)); % Número de píxeles 
        if max(unique(ring))>0
            pixNCLprop = length(find(white==1))./length(find(ring == 1));
        else
            pixNCLprop = 1; % Glándula sana
        end
        LGrelation = areaLumen/areaGland; % Relación tamaño lumen-glándula
            frontiersLumen = bwboundaries(lumenMask);
            if isempty(frontiersLumen)
                continue;
            end
            fLumen = frontiersLumen{1};
            x = fLumen(:,1);
            y = fLumen(:,2);
            centLumen = cat(1,propLumen.Centroid);
            for j = 1:length(x)
                dist(j)=sqrt((centLumen(2)-x(j))^2+(centLumen(1)-y(j))^2);
            end
        meanDistLumen= mean(dist); % Distancia media del centro del lumen a los píxeles del borde
        varDistLumen = var(dist); % Desvest del centro del lumen a los píxeles del borde
        if numNucBB == 0
            relNumNucBBgland = 0;
            relPixNucBBgland = 0;               
        else
            relNumNucBBgland = numNucGland/numNucBB;
            relPixNucBBgland = pixNucGland/pixNucBB;
        end
        if  pixCytoBB == 0
            relPixCytoBBgland = 0;
        else
            relPixCytoBBgland = pixCytoGland/pixCytoBB;
        end

    % FINAL Contextual FEATURES
    contextualFeatures = [LGrelation, meanDistLumen, varDistLumen, relNumNucBBgland, relPixNucBBgland,...
        relPixCytoBBgland, pixNCL, pixNCLprop];

    contextualNames = {'LGrelation', 'meanDistLumen', 'varDistLumen', 'relNumNucBBgland', 'relPixNucBBgland',...
                    'relPixCytoBBgland', 'pixNCL', 'pixNCLprop'};

%% %%%%%%%%%%%%%%%%%%%%%
%%% TEXTURE FEATURES %%% (186)
%%%%%%%%%%%%%%%%%%%%%%%%
textureFeatures = [];
textureNames = {};
    % GLCM    
        % CYAN. Phase 0 (21)
        cyan_mask = imresize(cyan, 2);
        cyanMask = cyan_mask(bb(2):bb(4),bb(1):bb(3));
        [Homogeneity_c0, Contrast_c0, Energy_c0,...
            Mean_c0, Desvest_c0, Entropy_c0, Correlation_c0] = grayLevelCoocurrenceMatrix(cyanMask, [2 0]);
        % Phase 45 (21)
        [Homogeneity_c45, Contrast_c45, Energy_c45,...
            Mean_c45, Desvest_c45, Entropy_c45, Correlation_c45] = grayLevelCoocurrenceMatrix(cyanMask, [-2 2]);

        % HEMATOXYLIN. Phase 0 (21)                
        [H,E] = colour_deconvolution(completImage, 'H&E');
        H = imresize(H,2);
        hematoMask = H(bb(2):bb(4),bb(1):bb(3));
        [Homogeneity_h0, Contrast_h0, Energy_h0,...
            Mean_h0, Desvest_h0, Entropy_h0, Correlation_h0] = grayLevelCoocurrenceMatrix(hematoMask, [2 0]);      
        % Phase 45 (21)
        [Homogeneity_h45, Contrast_h45, Energy_h45,...
            Mean_h45, Desvest_h45, Entropy_h45, Correlation_h45] = grayLevelCoocurrenceMatrix(hematoMask, [-2 2]); 

        % EOSIN. Phase 0 (21)
        E = imresize(E,2);
        eosinMask = E(bb(2):bb(4),bb(1):bb(3));
        [Homogeneity_e0, Contrast_e0, Energy_e0,...
            Mean_e0, Desvest_e0, Entropy_e0, Correlation_e0] = grayLevelCoocurrenceMatrix(eosinMask, [2 0]);  
        % Phase 45 (21)
        [Homogeneity_e45, Contrast_e45, Energy_e45,...
            Mean_e45, Desvest_e45, Entropy_e45, Correlation_e45] = grayLevelCoocurrenceMatrix(eosinMask, [-2 2]); 

    % LBPs 
            % Cyan (20)
            maping = getmapping(Neigh,'riu2');
            new_imgMask = imgMask(Radius+1:end,Radius+1:end);
            LBPnoMaskCyan_mask = lbp(cyanMask,Radius,Neigh,maping,'i');
            VARnoMaskCyan_mask = cont(cyanMask,Radius,Neigh);
            if ~isequal(size(LBPnoMaskCyan_mask),size(new_imgMask))
                disp(['Error. DIFFERENT LBP CYAN and Bbox GLAND in ' num2str(i) ' and class: ' num2str(k)]);
                continue;
            else
                LBPnoMaskCyan = LBPnoMaskCyan_mask;
                VARnoMaskCyan = VARnoMaskCyan_mask;
            end
            LBPmaskCyan = double(LBPnoMaskCyan).*new_imgMask; % Mask LBP
            VARmaskCyan = double(VARnoMaskCyan).*new_imgMask; % Mask lbp VAR
            new_imgMask = logical(new_imgMask);
        [hist_LBPcyan, hist_LBP_VARcyan] = histogramOfLBP(LBPnoMaskCyan, VARnoMaskCyan, new_imgMask);

            % Hematoxilyn (20)
            LBPnoMaskHemato_mask = lbp(hematoMask,Radius,Neigh,maping,'i');
            VARnoMaskHemato_mask = cont(hematoMask,Radius,Neigh);
            if ~isequal(size(LBPnoMaskHemato_mask),size(new_imgMask))
                disp(['Error. DIFFERENT LBP HEMATO and Bbox GLAND in ' num2str(i) ' and class: ' num2str(k)]);
                continue;
            else
                LBPnoMaskHemato = LBPnoMaskHemato_mask;
                VARnoMaskHemato = VARnoMaskHemato_mask;
            end
            LBPmaskHemato = double(LBPnoMaskHemato).*new_imgMask; % Mask LBP
            VARmaskHemato = double(VARnoMaskHemato).*new_imgMask; % Mask lbp VAR
            new_imgMask = logical(new_imgMask);
        [hist_LBPhemato, hist_LBP_VARhemato] = histogramOfLBP(LBPnoMaskHemato, VARnoMaskHemato, new_imgMask);

            % Eosin (20)
            LBPnoMaskEosin_mask = lbp(eosinMask,Radius,Neigh,maping,'i');
            VARnoMaskEosin_mask = cont(eosinMask,Radius,Neigh);
            if ~isequal(size(LBPnoMaskEosin_mask),size(new_imgMask))
                disp(['Error. DIFFERENT LBP EOSIN and Bbox GLAND in ' num2str(i) ' and class: ' num2str(k)]);
                continue;
            else
                LBPnoMaskEosin = LBPnoMaskEosin_mask;
                VARnoMaskEosin = VARnoMaskEosin_mask;
            end
            LBPmaskEosin = double(LBPnoMaskEosin).*new_imgMask; % Mask LBP
            VARmaskEosin = double(VARnoMaskEosin).*new_imgMask; % Mask lbp VAR
            new_imgMask = logical(new_imgMask);
        [hist_LBPeosin, hist_LBP_VAReosin] = histogramOfLBP(LBPnoMaskEosin, VARnoMaskEosin, new_imgMask);

    % FINAL Texture FEATURES
    textureFeatures = [Homogeneity_c0, Contrast_c0, Energy_c0,...
        Mean_c0, Desvest_c0, Entropy_c0, Correlation_c0,...
        Homogeneity_c45, Contrast_c45, Energy_c45,...
        Mean_c45, Desvest_c45, Entropy_c45, Correlation_c45,...
        Homogeneity_h0, Contrast_h0, Energy_h0,...
        Mean_h0, Desvest_h0, Entropy_h0, Correlation_h0,...
        Homogeneity_h45, Contrast_h45, Energy_h45,...
        Mean_h45, Desvest_h45, Entropy_h45, Correlation_h45,...
        Homogeneity_e0, Contrast_e0, Energy_e0,...
        Mean_e0, Desvest_e0, Entropy_e0, Correlation_e0,...
        Homogeneity_e45, Contrast_e45, Energy_e45,...
        Mean_e45, Desvest_e45, Entropy_e45, Correlation_e45, ...
        hist_LBPcyan, hist_LBP_VARcyan, ...
        hist_LBPhemato, hist_LBP_VARhemato,...
        hist_LBPeosin, hist_LBP_VAReosin];

    textureNames = {'Homogeneity_c0', 'Contrast_c0', 'Energy_c0',...
        'Mean_c0_1', 'Mean_c0_2', 'Mean_c0_3', 'Mean_c0_4', 'Mean_c0_5', 'Mean_c0_6', 'Mean_c0_7', 'Mean_c0_8',... 
        'Desvest_c0_1', 'Desvest_c0_2', 'Desvest_c0_3', 'Desvest_c0_4', 'Desvest_c0_5', 'Desvest_c0_6', 'Desvest_c0_7', 'Desvest_c0_8', ...
        'Entropy_c0', 'Correlation_c0', ...     
        'Homogeneity_c45', 'Contrast_c45', 'Energy_c45',...
        'Mean_c45_1', 'Mean_c45_2', 'Mean_c45_3', 'Mean_c45_4', 'Mean_c45_5', 'Mean_c45_6', 'Mean_c45_7', 'Mean_c45_8', ...
        'Desvest_c45_1', 'Desvest_c45_2', 'Desvest_c45_3', 'Desvest_c45_4', 'Desvest_c45_5', 'Desvest_c45_6', 'Desvest_c45_7', 'Desvest_c45_8', ...
        'Entropy_c45', 'Correlation_c45', ...          
        'Homogeneity_h0', 'Contrast_h0', 'Energy_h0',...
        'Mean_h0_1', 'Mean_h0_2', 'Mean_h0_3', 'Mean_h0_4', 'Mean_h0_5', 'Mean_h0_6', 'Mean_h0_7', 'Mean_h0_8', ...
        'Desvest_h0_1', 'Desvest_h0_2', 'Desvest_h0_3', 'Desvest_h0_4', 'Desvest_h0_5', 'Desvest_h0_6', 'Desvest_h0_7', 'Desvest_h0_8', ...
        'Entropy_h0', 'Correlation_h0',...                
        'Homogeneity_h45', 'Contrast_h45', 'Energy_h45',...
        'Mean_h45_1', 'Mean_h45_2', 'Mean_h45_3', 'Mean_h45_4', 'Mean_h45_5', 'Mean_h45_6', 'Mean_h45_7', 'Mean_h45_8', ...
        'Desvest_h45_1', 'Desvest_h45_2', 'Desvest_h45_3', 'Desvest_h45_4', 'Desvest_h45_5', 'Desvest_h45_6', 'Desvest_h45_7', 'Desvest_h45_8', ...
        'Entropy_h45', 'Correlation_h45', ...
        'Homogeneity_e0', 'Contrast_e0', 'Energy_e0',...
        'Mean_e0_1', 'Mean_e0_2', 'Mean_e0_3', 'Mean_e0_4', 'Mean_e0_5', 'Mean_e0_6', 'Mean_e0_7', 'Mean_e0_8', ...
        'Desvest_e0_1', 'Desvest_e0_2', 'Desvest_e0_3', 'Desvest_e0_4', 'Desvest_e0_5', 'Desvest_e0_6', 'Desvest_e0_7', 'Desvest_e0_8', ...
        'Entropy_e0', 'Correlation_e0', ...
        'Homogeneity_e45', 'Contrast_e45', 'Energy_e45',...
        'Mean_e45_1', 'Mean_e45_2', 'Mean_e45_3', 'Mean_e45_4', 'Mean_e45_5', 'Mean_e45_6', 'Mean_e45_7', 'Mean_e45_8', ...
        'Desvest_e45_1', 'Desvest_e45_2', 'Desvest_e45_3', 'Desvest_e45_4','Desvest_e45_5','Desvest_e45_6', 'Desvest_e45_7', 'Desvest_e45_8', ... 
        'Entropy_e45', 'Correlation_e45', ...                
        'hist_LBPcyan1', 'hist_LBPcyan2', 'hist_LBPcyan3', 'hist_LBPcyan4', 'hist_LBPcyan5', 'hist_LBPcyan6', 'hist_LBPcyan7', 'hist_LBPcyan8', 'hist_LBPcyan9', 'hist_LBPcyan10', ...
        'hist_LBP_VARcyan1', 'hist_LBP_VARcyan2', 'hist_LBP_VARcyan3', 'hist_LBP_VARcyan4', 'hist_LBP_VARcyan5', 'hist_LBP_VARcyan6', 'hist_LBP_VARcyan7', 'hist_LBP_VARcyan8', 'hist_LBP_VARcyan9', 'hist_LBP_VARcyan10', ...
        'hist_LBPhemato1', 'hist_LBPhemato2', 'hist_LBPhemato3', 'hist_LBPhemato4', 'hist_LBPhemato5', 'hist_LBPhemato6', 'hist_LBPhemato7', 'hist_LBPhemato8', 'hist_LBPhemato9', 'hist_LBPhemato10', ...
        'hist_LBP_VARhemato1', 'hist_LBP_VARhemato2', 'hist_LBP_VARhemato3', 'hist_LBP_VARhemato4', 'hist_LBP_VARhemato5', 'hist_LBP_VARhemato6', 'hist_LBP_VARhemato7', 'hist_LBP_VARhemato8', 'hist_LBP_VARhemato9', 'hist_LBP_VARhemato10', ...
        'hist_LBPeosin1', 'hist_LBPeosin2', 'hist_LBPeosin3', 'hist_LBPeosin4', 'hist_LBPeosin5', 'hist_LBPeosin6', 'hist_LBPeosin7', 'hist_LBPeosin8', 'hist_LBPeosin9', 'hist_LBPeosin10', ...
        'hist_LBP_VAReosin1', 'hist_LBP_VAReosin2', 'hist_LBP_VAReosin3', 'hist_LBP_VAReosin4', 'hist_LBP_VAReosin5', 'hist_LBP_VAReosin6', 'hist_LBP_VAReosin7', 'hist_LBP_VAReosin8', 'hist_LBP_VAReosin9', 'hist_LBP_VAReosin10'};


%% %%%%%%%%%%%%%%%%%%%%%
%%% FRACTAL ANALYSIS %%% (15)
%%%%%%%%%%%%%%%%%%%%%%%%
    % Hurst (15)
        directions = [0,30,45,60,90];
        for jj = 1:5
            HurstCyan(jj) = fracdim2o(cyanMask,'GSE', directions(jj));
        end
        for jj = 1:5
            HurstHemato(jj) = fracdim2o(hematoMask,'GSE', directions(jj));
        end
        for jj = 1:5
            HurstEosin(jj) = fracdim2o(eosinMask,'GSE', directions(jj));
        end

    % FINAL TEXTURAL AND FRACTAL FEATURES
    fractalFeatures = [HurstCyan, HurstHemato, HurstEosin];
    fractalNames = {'HurstCyan1', 'HurstCyan2', 'HurstCyan3', 'HurstCyan4', 'HurstCyan5', ...
        'HurstHemato1', 'HurstHemato2', 'HurstHemato3', 'HurstHemato4', 'HurstHemato5', ...
        'HurstEosin1', 'HurstEosin2', 'HurstEosin3', 'HurstEosin4', 'HurstEosin5'};


        %% ALL FEATURES
        ALL_FEATURES = [shapeFeatures, contextualFeatures, textureFeatures, fractalFeatures];
        ALL_NAMES = [shapeNames, contextualNames, textureNames, fractalNames];
        if k == 1
            response = 0;
        elseif k == 2
            response = 1;
        else
            response = 2;
        end
        featMatrix(i,:) = [ALL_FEATURES, response];

        waitbar(step/leng,h,leng-step); % Show how many images remain to be saved
        step = step+1;
    end
    if k == 1
        save featMatrixArtifacts.mat featMatrix;
    elseif k == 2
        save featMatrixGrade3.mat featMatrix;
    else
        save featMatrixHealthy.mat featMatrix;
    end
    delete(h);
end
% Creation of .mat
load featMatrixArtifacts.mat;
artifacts = featMatrix;
load featMatrixGrade3.mat;
grade3 = featMatrix;
load featMatrixHealthy.mat;
healthy = featMatrix;

featMatrix = [artifacts; grade3; healthy];

save featMatrixExtracted_v1.mat featMatrix;

delete featMatrixArtifacts.mat
delete featMatrixGrade3.mat
delete featMatrixHealthy.mat

names = ALL_NAMES';
save featNames_v1.mat names;   
    

