function [featTest_mc] = Extract_Features_Test_mc(maskImg, Bbox, img, mask_black, nuclei_mask, cytoplasm_mask0, cyan)

completImage = img;
Neigh = 8;
Radius = 1;

for i = 1:length(maskImg)
    imgMask = maskImg{i};
    bb = Bbox(i,:); 
    
    %% %%%%%%%%%%%%%%%%%%%
    %%% SHAPE FEATURES %%%
    %%%%%%%%%%%%%%%%%%%%%%
    ShapeFeatures = [];
    %% GLAND
    propGland = regionprops(imgMask, 'Area', 'BoundingBox', 'Centroid', 'ConvexArea', ...
        'Eccentricity', 'EquivDiameter', 'Extent', 'Orientation', 'Perimeter', 'Solidity');
    % Features (6)
    areaGland = cat(1,propGland.Area);
    convexAreaGland = cat(1,propGland.ConvexArea);
    solidityGland = convexAreaGland/areaGland;
    eccentGland = cat(1,propGland.Eccentricity);
    equivDiamGland = cat(1,propGland.EquivDiameter);
    extentGland = cat(1,propGland.Extent);
    perimeterGland = sum(cat(1, propGland.Perimeter));
        rG=sqrt(areaGland/pi);
    if areaGland ~= 0
        roundnessGland = rG*perimeterGland/areaGland;
    else
        roundnessGland = 0;
    end
    cent = cat(1,propGland.Centroid);
    
    %% LUMEN
    lumens = mask_black;
    if length(unique(lumens))>2
        lumens = lumens>100;
    end
    lumens = imfill(lumens,'holes');
    if sum(bb) == 0
        bb(1) = 1; bb(2) = 1; bb(3) = size(img,2); bb(4) = size(img,1);
    end
    lumeni_mask = lumens(bb(2):bb(4),bb(1):bb(3));
    lumeni_mask = bwareaopen(lumeni_mask, 20);
    if max(unique(lumeni_mask)) == 0 % If some error exists, next!
        continue;
    end
    if ~isequal(size(lumeni_mask),size(imgMask))
        difRow = abs(size(lumeni_mask,1)-size(imgMask,1));
        difCol = abs(size(lumeni_mask,2)-size(imgMask,2));
        if size(lumeni_mask,1)<size(imgMask,1) & size(lumeni_mask,2)<size(imgMask,2)
            lumenix = [zeros(difRow,size(imgMask,2)-abs(size(imgMask,2)-size(lumeni_mask,2))); lumeni_mask];
            lumeni = [zeros(size(imgMask,1),difCol), lumenix];
        elseif size(lumeni_mask,1)<=size(imgMask,1) & size(lumeni_mask,2)>=size(imgMask,2)
            lumeni = [zeros(difRow,max(size(imgMask,2),size(lumeni_mask,2))); lumeni_mask];
            imgM = [zeros(max(size(lumeni_mask,1),size(imgMask,1)),difCol), imgMask]; 
            imgMask = imgM;
        elseif size(lumeni_mask,1)>=size(imgMask,1) & size(lumeni_mask,2)<=size(imgMask,2)
            lumeni = [zeros(max(size(imgMask,1),size(lumeni_mask,1)),difCol), lumeni_mask];
            imgM = [zeros(difRow,max(size(lumeni_mask,2),size(imgMask,2))); imgMask]; 
            imgMask = imgM;
        else
            lumeni = lumeni_mask;
            imgM = [zeros(difRow,size(lumeni_mask,2)-abs(size(lumeni_mask,2)-size(imgMask,2))); imgMask];
            imgMask = [zeros(size(lumeni_mask,1),difCol), imgM];
        end
    else
        lumeni = lumeni_mask;
    end
    lumenMask = lumeni.*imgMask;

    [lum_img, num_lum] = bwlabel(lumenMask);
    if num_lum>1 % if lumenMask has more than one object
        p1 = round(cent(2));
        p2 = round(cent(1));
        dist = [];
        for ll = 1:num_lum
            elem = lum_img == ll;
            pr = regionprops(elem, 'Centroid');
            centr = cat(1,pr.Centroid);
            dist(ll) = pdist2(centr,[p2,p1], 'Euclidean');
        end
        id_dist = find(dist == min(dist));
        lumenMask = lum_img == id_dist;
    end

    propLumen = regionprops(lumenMask, 'Area', 'BoundingBox', 'Centroid', 'ConvexArea',...
        'Eccentricity', 'EquivDiameter', 'Extent', 'Orientation', 'Perimeter', 'Solidity');
    % Features (7)
    areaLumen = cat(1,propLumen.Area); 
    convexAreaLumen = cat(1,propLumen.ConvexArea);
    solidityLumen = convexAreaLumen/areaLumen;
    eccentLumen = cat(1,propLumen.Eccentricity);
    equivDiamLumen = cat(1,propLumen.EquivDiameter);
    extentLumen = cat(1,propLumen.Extent);
    perimeterLumen = cat(1, propLumen.Perimeter);
        rG=sqrt(areaLumen/pi);
    if ~isempty(areaLumen)
        roundnessLumen = rG*perimeterLumen/areaLumen;
    else
        roundnessLumen = 0;
    end
    LGrelation = areaLumen/areaGland;
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
    meanDistLumen= mean(dist);
    varDistLumen = var(dist);
    
    %% NUCLEI
    nuc_mask = nuclei_mask(bb(2):bb(4),bb(1):bb(3));
    if ~isequal(size(nuc_mask),size(imgMask))
        difRow = abs(size(nuc_mask,1)-size(imgMask,1));
        difCol = abs(size(nuc_mask,2)-size(imgMask,2));
        if size(nuc_mask,1)<size(imgMask,1) & size(nuc_mask,2)<size(imgMask,2)
            nucMaski = [zeros(difRow,size(imgMask,2)-abs(size(nuc_mask,2)-size(imgMask,2))); nuc_mask];
            nucMask = [zeros(size(imgMask,1),difCol), nucMaski];
        elseif size(nuc_mask,1)<=size(imgMask,1) & size(nuc_mask,2)>=size(imgMask,2)
            nucMask = [zeros(difRow,max(size(imgMask,2),size(nuc_mask,2))); nuc_mask];
            imgM = [zeros(max(size(nuc_mask,1),size(imgMask,1)),difCol), imgMask]; 
            imgMask = imgM;
        elseif size(nuc_mask,1)>=size(imgMask,1) & size(nuc_mask,2)<=size(imgMask,2)
            nucMask = [zeros(max(size(imgMask,1),size(nuc_mask,1)),difCol), nuc_mask];
            imgM = [zeros(difRow,max(size(nuc_mask,2),size(imgMask,2))); imgMask]; 
            imgMask = imgM;
        else
            nucMask = nuc_mask;
            imgM = [zeros(difRow,size(nuc_mask,2)-abs(size(nuc_mask,2)-size(imgMask,2))); imgMask];
            imgMask = [zeros(size(nuc_mask,1),difCol), imgM];
        end
    else
        nucMask = nuc_mask;
    end

    newMask = zeros(size(img,1),size(img,2));
    [xx,yy] = find(imgMask == 1);
    for ii = 1:length(xx)
        newMask(xx(ii)+bb(2),yy(ii)+bb(1)) = 1;
    end
    newMask = imdilate(newMask,strel('disk',2));
    try
        nucleiMask = nuclei_mask(bb(2)-10:bb(4)+10,bb(1)-10:bb(3)+10);
        newImgMask = newMask(bb(2)-10:bb(4)+10,bb(1)-10:bb(3)+10);
    catch
        nucleiMask = nuclei_mask(bb(2):bb(4),bb(1):bb(3));
        newImgMask = newMask(bb(2):bb(4),bb(1):bb(3));
    end
     
    propNuclei = regionprops(nucleiMask, 'Area');
    areaNuc = cat(1,propNuclei.Area); % area of nuclei in the BB
    nucGland = nucleiMask.*newImgMask;
    propNuclei2 = regionprops(nucGland, 'Area');
    areaNuc2 = cat(1,propNuclei2.Area); % area of nuclei in the mask
    sizeBB = [size(imgMask,1),size(imgMask,2)];

  % Features (9)
    [~, numNucBB] = bwlabel(nucleiMask);
    numNucBBprop = numNucBB/(sizeBB(1)*sizeBB(2));
    [~, numNucGland] = bwlabel(nucGland);
    pixNucBB = sum(areaNuc);
    pixNucBBprop = pixNucBB/(sizeBB(1)*sizeBB(2));
    pixNucGland = sum(areaNuc2);
    if ~isempty(areaNuc2)
        numNucGlandProp = numNucGland./areaNuc2;
        pixNucGlandProp = pixNucGland./areaNuc2;
        relNumNucBBgland = numNucGland/numNucBB;
        relPixNucBBgland = pixNucGland/pixNucBB;
    else
        numNucGlandProp = 0;
        pixNucGlandProp = 0;
        relNumNucBBgland = 0;
        relPixNucBBgland = 0;
    end
    
    %% CYTOPLASM
    cytoplasm_mask = cytoplasm_mask0(bb(2):bb(4), bb(1):bb(3));
    if ~isequal(size(cytoplasm_mask),size(imgMask))
        difRow = abs(size(cytoplasm_mask,1)-size(imgMask,1));
        difCol = abs(size(cytoplasm_mask,2)-size(imgMask,2));
        if size(cytoplasm_mask,1)<size(imgMask,1) & size(cytoplasm_mask,2)<size(imgMask,2)
            cytoplasmMaski = [zeros(difRow,size(imgMask,2)-abs(size(cytoplasm_mask,2)-size(imgMask,2))); cytoplasm_mask];
            cytoplasmMask = [zeros(size(imgMask,1),difCol), cytoplasmMaski];
        elseif size(cytoplasm_mask,1)<=size(imgMask,1) & size(cytoplasm_mask,2)>=size(imgMask,2)
            cytoplasmMask = [zeros(difRow,max(size(imgMask,2),size(cytoplasm_mask,2))); cytoplasm_mask];
            imgM = [zeros(max(size(cytoplasm_mask,1),size(imgMask,1)),difCol), imgMask]; 
            imgMask = imgM;
        elseif size(cytoplasm_mask,1)>=size(imgMask,1) & size(cytoplasm_mask,2)<=size(imgMask,2)
            cytoplasmMask = [zeros(max(size(imgMask,1),size(cytoplasm_mask,1)),difCol), cytoplasm_mask];
            imgM = [zeros(difRow,max(size(cytoplasm_mask,2),size(imgMask,2))); imgMask]; 
            imgMask = imgM;
        else
            cytoplasmMask = cytoplasm_mask;
            imgM = [zeros(difRow,size(cytoplasm_mask,2)-abs(size(cytoplasm_mask,2)-size(imgMask,2))); imgMask];
            imgMask = [zeros(size(cytoplasm_mask,1),difCol), imgM];
        end
    else
        cytoplasmMask = cytoplasm_mask;
    end
    
    propCyto = regionprops(cytoplasmMask, 'Area');
    areaCyto = cat(1,propCyto.Area); % area of cytoplasm in the BB
    cytoGland = cytoplasmMask.*imgMask;
    propCyto2 = regionprops(cytoGland, 'Area');
    areaCyto2 = cat(1,propCyto2.Area); % area of cytoplasm in the mask
    % Features (5)
    pixCytoBB = sum(areaCyto);
    pixCytoBBprop = pixCytoBB/(sizeBB(1)*sizeBB(2));
    pixCytoGland = sum(areaCyto2);
    if ~isempty(areaCyto2)
        pixCytoGlandprop = pixCytoGland./areaCyto2;
        relPixCytoBBgland = pixCytoGland/pixCytoBB;
    else
        pixCytoGlandprop = 0;
        relPixCytoBBgland = 0;
    end
    
    %% NUCLEI + CYTOPLASM + LUMEN  (2)
    ring = logical(imgMask-lumenMask);
    nucRing = logical(ring.*nucMask);
    cytoRing = logical(ring.*cytoplasmMask);
    lumRing = logical(ring.*lumenMask);
    white = logical(nucRing+cytoRing+lumRing);
    pixNCL = length(find(white==1));
    if max(unique(ring))>0
        pixNCLprop = length(find(white==1))./length(find(ring == 1));
    else
        pixNCLprop = 1;
    end
    
    %% FINAL Shape + Contextual FEATURES
    ShapeFeatures = [areaGland, solidityGland, eccentGland, equivDiamGland, extentGland, roundnessGland, ...
        areaLumen,  solidityLumen, eccentLumen, equivDiamLumen, extentLumen, perimeterLumen, roundnessLumen, LGrelation, meanDistLumen, varDistLumen, ...
        numNucBB, numNucBBprop, numNucGland, pixNucBB, pixNucBBprop, numNucGlandProp, pixNucGlandProp, relNumNucBBgland, relPixNucBBgland,...
        pixCytoBB, pixCytoBBprop, pixCytoGland, pixCytoGlandprop, relPixCytoBBgland, pixNCLprop];

               
    %% %%%%%%%%%%%%%%%%%%%%%
    %%% TEXTURE FEATURES %%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    TextureFeatures = [];
    %% GLCM
    % CYAN
    % Features(12)
    cyanMask = cyan(bb(2):bb(4),bb(1):bb(3));
    glcmCyanF0 = graycomatrix(cyanMask,'Offset',[2 0]);
    statsCyanF0 = graycoprops(glcmCyanF0,{'contrast','correlation','energy','homogeneity'});
    ContrastCyanF0 = statsCyanF0.Contrast;
    CorrelationCyanF0 = statsCyanF0.Correlation;
    EnergyCyanF0 = statsCyanF0.Energy;
    HomogeneityCyanF0 = statsCyanF0.Homogeneity;
    meanGLCMcyanF0 = mean(glcmCyanF0);

    glcmCyanF45=graycomatrix(cyanMask,'Offset',[-2 2]);
    statsCyanF45 = graycoprops(glcmCyanF45,{'contrast','correlation','energy','homogeneity'});
    ContrastCyanF45 = statsCyanF45.Contrast;
    CorrelationCyanF45 = statsCyanF45.Correlation;
    HomogeneityCyanF45 = statsCyanF45.Homogeneity;

    % HEMATOXYLIN
    % Features(7)
    [H,E] = colour_deconvolution(completImage, 'H&E');
    hematoMask = H(bb(2):bb(4),bb(1):bb(3));
    glcmHematoF0 = graycomatrix(hematoMask, 'Offset', [2,0]);
    statsHematoF0 = graycoprops(glcmHematoF0, {'contrast','correlation','energy','homogeneity'});
    CorrelationHematoF0 = statsHematoF0.Correlation;
    EnergyHematoF0 = statsHematoF0.Energy;
    HomogeneityHematoF0 = statsHematoF0.Homogeneity;
    meanGLCMhematoF0 = mean(glcmHematoF0);       

    glcmHematoF45 = graycomatrix(hematoMask, 'Offset', [-2,2]);
    statsHematoF45 = graycoprops(glcmHematoF45, {'contrast','correlation','energy','homogeneity'});
    CorrelationHematoF45 = statsHematoF45.Correlation; 

    % EOSIN
    % Features (15)
    eosinMask = E(bb(2):bb(4),bb(1):bb(3));
    glcmEosinF0 = graycomatrix(eosinMask, 'Offset', [2,0]);
    statsEosinF0 = graycoprops(glcmEosinF0, {'contrast','correlation','energy','homogeneity'});
    ContrastEosinF0 = statsEosinF0.Contrast;
    CorrelationEosinF0 = statsEosinF0.Correlation;
    HomogeneityEosinF0 = statsEosinF0.Homogeneity;
    meanGLCMeosinF0 = mean(glcmEosinF0);      

    glcmEosinF45 = graycomatrix(eosinMask, 'Offset', [-2,2]);
    statsEosinF45 = graycoprops(glcmEosinF45, {'contrast','correlation','energy','homogeneity'});
    ContrastEosinF45 = statsEosinF45.Contrast;
    CorrelationEosinF45 = statsEosinF45.Correlation;
    EnergyEosinF45 = statsEosinF45.Energy;
    HomogeneityEosinF45 = statsEosinF45.Homogeneity;

    %% LBPs 
    % CYAN (20)
    maping = getmapping(Neigh,'riu2');
    new_imgMask = imgMask(Radius+1:end-Radius,Radius+1:end-Radius);
    LBPnoMaskCyan_mask = lbp(cyanMask,Radius,Neigh,maping,'i');
    VARnoMaskCyan_mask = cont(cyanMask,Radius,Neigh);
    if ~isequal(size(LBPnoMaskCyan_mask),size(new_imgMask))
        difRow = abs(size(LBPnoMaskCyan_mask,1)-size(new_imgMask,1));
        difCol = abs(size(LBPnoMaskCyan_mask,2)-size(new_imgMask,2));
        if size(LBPnoMaskCyan_mask,1)<size(new_imgMask,1) & size(LBPnoMaskCyan_mask,2)<size(new_imgMask,2)
            LBPnoMaskCyanx = [zeros(difRow,size(new_imgMask,2)-1); LBPnoMaskCyan_mask];
            LBPnoMaskCyan = [zeros(size(new_imgMask,1),difCol), LBPnoMaskCyanx];
            VARnoMaskCyanx = [zeros(difRow,size(new_imgMask,2)-1); VARnoMaskCyan_mask];
            VARnoMaskCyan = [zeros(size(new_imgMask,1),difCol), VARnoMaskCyanx];                        
        elseif size(LBPnoMaskCyan_mask,1)<=size(new_imgMask,1) & size(LBPnoMaskCyan_mask,2)>=size(new_imgMask,2)
            LBPnoMaskCyan = [zeros(difRow,max(size(new_imgMask,2),size(LBPnoMaskCyan_mask,2))); LBPnoMaskCyan_mask];
            VARnoMaskCyan = [zeros(difRow,max(size(new_imgMask,2),size(VARnoMaskCyan_mask,2))); VARnoMaskCyan_mask];
            imgM = [zeros(max(size(LBPnoMaskCyan_mask,1),size(new_imgMask,1)),difCol), new_imgMask]; 
            new_imgMask = imgM;
        elseif size(LBPnoMaskCyan_mask,1)>=size(new_imgMask,1) & size(LBPnoMaskCyan_mask,2)<=size(new_imgMask,2)
            LBPnoMaskCyan = [zeros(max(size(new_imgMask,1),size(LBPnoMaskCyan_mask,1)),difCol), LBPnoMaskCyan_mask];
            VARnoMaskCyan = [zeros(max(size(new_imgMask,1),size(VARnoMaskCyan_mask,1)),difCol), VARnoMaskCyan_mask];
            imgM = [zeros(difRow,max(size(LBPnoMaskCyan_mask,2),size(new_imgMask,2))); new_imgMask]; 
            new_imgMask = imgM;
        else
            LBPnoMaskCyan = LBPnoMaskCyan_mask;
            VARnoMaskCyan = VARnoMaskCyan_mask;
            imgM = [zeros(difRow,size(LBPnoMaskCyan_mask,2)-1); new_imgMask];
            new_imgMask = [zeros(size(LBPnoMaskCyan_mask,1),difCol), imgM];
        end
    else
        LBPnoMaskCyan = LBPnoMaskCyan_mask;
        VARnoMaskCyan = VARnoMaskCyan_mask;
    end
    LBPmaskCyan = double(LBPnoMaskCyan).*new_imgMask; % Mask LBP
    VARmaskCyan = double(VARnoMaskCyan).*new_imgMask; % Mask lbp VAR
    new_imgMask = logical(new_imgMask);
    [hist_LBPcyan, hist_LBP_VARcyan] = histogramOfLBP(LBPnoMaskCyan, VARnoMaskCyan, new_imgMask);

    % HEMATOXYLIN (20)
    LBPnoMaskHemato_mask = lbp(hematoMask,Radius,Neigh,maping,'i');
    VARnoMaskHemato_mask = cont(hematoMask,Radius,Neigh);
    if ~isequal(size(LBPnoMaskHemato_mask),size(new_imgMask))
        difRow = abs(size(LBPnoMaskHemato_mask,1)-size(new_imgMask,1));
        difCol = abs(size(LBPnoMaskHemato_mask,2)-size(new_imgMask,2));
        if size(LBPnoMaskHemato_mask,1)<size(new_imgMask,1) & size(LBPnoMaskHemato_mask,2)<size(new_imgMask,2)
            LBPnoMaskHematox = [zeros(difRow,size(new_imgMask,2)-1); LBPnoMaskHemato_mask];
            LBPnoMaskHemato = [zeros(size(new_imgMask,1),difCol), LBPnoMaskHematox];
            VARnoMaskHematox = [zeros(difRow,size(new_imgMask,2)-1); LBPnoMaskHemato_mask];
            VARnoMaskHemato = [zeros(size(new_imgMask,1),difCol), VARnoMaskHematox];                        
        elseif size(LBPnoMaskHemato_mask,1)<=size(new_imgMask,1) & size(LBPnoMaskHemato_mask,2)>=size(new_imgMask,2)
            LBPnoMaskHemato = [zeros(difRow,max(size(new_imgMask,2),size(LBPnoMaskHemato_mask,2))); LBPnoMaskHemato_mask];
            VARnoMaskHemato = [zeros(difRow,max(size(new_imgMask,2),size(LBPnoMaskHemato_mask,2))); LBPnoMaskHemato_mask];
            imgM = [zeros(max(size(LBPnoMaskHemato_mask,1),size(new_imgMask,1)),difCol), new_imgMask]; 
            new_imgMask = imgM;
        elseif size(LBPnoMaskHemato_mask,1)>=size(new_imgMask,1) & size(LBPnoMaskHemato_mask,2)<=size(new_imgMask,2)
            LBPnoMaskHemato = [zeros(max(size(new_imgMask,1),size(LBPnoMaskHemato_mask,1)),difCol), LBPnoMaskHemato_mask];
            VARnoMaskHemato = [zeros(max(size(new_imgMask,1),size(LBPnoMaskHemato_mask,1)),difCol), LBPnoMaskHemato_mask];
            imgM = [zeros(difRow,max(size(LBPnoMaskHemato_mask,2),size(new_imgMask,2))); new_imgMask]; 
            new_imgMask = imgM;
        else
            LBPnoMaskHemato = LBPnoMaskHemato_mask;
            VARnoMaskHemato = LBPnoMaskHemato_mask;
            imgM = [zeros(difRow,size(LBPnoMaskCyan_mask,2)-1); new_imgMask];
            new_imgMask = [zeros(size(LBPnoMaskCyan_mask,1),difCol), imgM];
        end
    else
        LBPnoMaskHemato = LBPnoMaskHemato_mask;
        VARnoMaskHemato = VARnoMaskHemato_mask;
    end
    LBPmaskHemato = double(LBPnoMaskHemato).*new_imgMask; % Mask LBP
    VARmaskHemato = double(VARnoMaskHemato).*new_imgMask; % Mask lbp VAR
    new_imgMask = logical(new_imgMask);
    [hist_LBPhemato, hist_LBP_VARhemato] = histogramOfLBP(LBPnoMaskHemato, VARnoMaskHemato, new_imgMask);

    % EOSIN (20)
    LBPnoMaskEosin_mask = lbp(eosinMask,Radius,Neigh,maping,'i');
    VARnoMaskEosin_mask = cont(eosinMask,Radius,Neigh);
    if ~isequal(size(LBPnoMaskEosin_mask),size(new_imgMask))
        difRow = abs(size(LBPnoMaskEosin_mask,1)-size(new_imgMask,1));
        difCol = abs(size(LBPnoMaskEosin_mask,2)-size(new_imgMask,2));
        if size(LBPnoMaskEosin_mask,1)<size(new_imgMask,1) & size(LBPnoMaskEosin_mask,2)<size(new_imgMask,2)
            LBPnoMaskEosinx = [zeros(difRow,size(new_imgMask,2)-1); LBPnoMaskEosin_mask];
            LBPnoMaskEosin = [zeros(size(new_imgMask,1),difCol), LBPnoMaskEosinx];
            VARnoMaskEosinx = [zeros(difRow,size(new_imgMask,2)-1); LBPnoMaskEosin_mask];
            VARnoMaskEosin = [zeros(size(new_imgMask,1),difCol), VARnoMaskEosinx];                        
        elseif size(LBPnoMaskEosin_mask,1)<=size(new_imgMask,1) & size(LBPnoMaskEosin_mask,2)>=size(new_imgMask,2)
            LBPnoMaskEosin = [zeros(difRow,max(size(new_imgMask,2),size(LBPnoMaskEosin_mask,2))); LBPnoMaskEosin_mask];
            VARnoMaskEosin = [zeros(difRow,max(size(new_imgMask,2),size(LBPnoMaskEosin_mask,2))); LBPnoMaskEosin_mask];
            imgM = [zeros(max(size(LBPnoMaskEosin_mask,1),size(new_imgMask,1)),difCol), new_imgMask]; 
            new_imgMask = imgM;
        elseif size(LBPnoMaskEosin_mask,1)>=size(new_imgMask,1) & size(LBPnoMaskEosin_mask,2)<=size(new_imgMask,2)
            LBPnoMaskEosin = [zeros(max(size(new_imgMask,1),size(LBPnoMaskEosin_mask,1)),difCol), LBPnoMaskEosin_mask];
            VARnoMaskEosin = [zeros(max(size(new_imgMask,1),size(LBPnoMaskEosin_mask,1)),difCol), LBPnoMaskEosin_mask];
            imgM = [zeros(difRow,max(size(LBPnoMaskEosin_mask,2),size(new_imgMask,2))); new_imgMask]; 
            new_imgMask = imgM;
        else
            LBPnoMaskEosin = LBPnoMaskEosin_mask;
            VARnoMaskEosin = LBPnoMaskEosin_mask;
            imgM = [zeros(difRow,size(LBPnoMaskCyan_mask,2)-1); new_imgMask];
            new_imgMask = [zeros(size(LBPnoMaskCyan_mask,1),difCol), imgM];
        end
    else
        LBPnoMaskEosin = LBPnoMaskEosin_mask;
        VARnoMaskEosin = VARnoMaskEosin_mask;
    end
    LBPmaskEosin = double(LBPnoMaskEosin).*new_imgMask; % Mask LBP
    VARmaskEosin = double(VARnoMaskEosin).*new_imgMask; % Mask lbp VAR
    new_imgMask = logical(new_imgMask);
    [hist_LBPeosin, hist_LBP_VAReosin] = histogramOfLBP(LBPnoMaskEosin, VARnoMaskEosin, new_imgMask);

    
    %% %%%%%%%%%%%%%%%%%%%%%
    %%% FRACTAL ANALYSIS %%%
    %%%%%%%%%%%%%%%%%%%%%%%%    
    %% Hurst (11)
    directions = [0,30,45,60,90];
    for jj = 1:3
        HurstCyan(jj) = fracdim2o(eosinMask,'GSE', directions(jj));
    end
    for jj = 1:5
        HurstHemato(jj) = fracdim2o(eosinMask,'GSE', directions(jj));
    end
    for jj = 1:3
        HurstEosin(jj) = fracdim2o(eosinMask,'GSE', directions(jj));
    end
    
   
    %% FINAL Texture and Fractal FEATURES
    TextureFeatures = [ContrastCyanF0, CorrelationCyanF0, EnergyCyanF0, HomogeneityCyanF0, ...
        meanGLCMcyanF0(1),meanGLCMcyanF0(2),meanGLCMcyanF0(3),meanGLCMcyanF0(7),meanGLCMcyanF0(8),...
        ContrastCyanF45, CorrelationCyanF45, HomogeneityCyanF45, ...
        CorrelationHematoF0, EnergyHematoF0, HomogeneityHematoF0, meanGLCMhematoF0(1), meanGLCMhematoF0(2), meanGLCMhematoF0(7), ...
        CorrelationHematoF45,...
        ContrastEosinF0, CorrelationEosinF0, HomogeneityEosinF0, meanGLCMeosinF0, ...
        ContrastEosinF45, CorrelationEosinF45, EnergyEosinF45, HomogeneityEosinF45, ...
        hist_LBPcyan, hist_LBP_VARcyan, hist_LBPhemato, hist_LBP_VARhemato,...
        hist_LBPeosin, hist_LBP_VAReosin, ...
        HurstCyan, HurstHemato, HurstEosin];

    featTest_mc(i,:) = [ShapeFeatures, TextureFeatures];
end
