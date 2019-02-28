%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% FEATURE SELECTION STAGE %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

verbose = 0;
numObser = 10;
version = 'Masks';
rng('default');

addpath('..\4_FeatureExtraction');

%% Parameters
trust = 1-10^-5; % Trust of the 99.99999% 
alpha = 1-trust;
th_cor = 0.95; % Threshold to eliminate variables which are very correlated

%% Load data
load featMatrixExtracted.mat ; % features VALUES
load featNames.mat   % features names
load fileMatrix.mat; % names of BBoxes
namesFeat = names;

%% Features, names and response
features = featMatrix;
files = fileMatrix;
% Remove errors
for i = 1:length(features)
    idx(i) = sum(features(i,1:end-1)) == 0;
end
indx = find(idx == 1);
features(indx,:) = [];
files(indx) = [];
% Change NaNs by mean
[r,c] = find(isnan(features));
for i = 1:length(r)
    gt = features(r(i),end);
    [x,~] = find(features(:,end) == gt);
    fea = features(x,c(i));
    fea(isnan(fea)) = [];
    m = mean(fea);
    features(r(i),c(i)) = m;
end
featMatrix = features;

%% Response and features
response = features(:,end);
features = features(:,1:end-1);
allFeatures = features;
            
%% Normalization of values of variables
[features,mediaZ,sigmaZ] = zscore(features);

%% If we want visually show the results:
if verbose == 0
    numObser = size(features,2);
end

numFig = ceil(size(features,2)/numObser);
infLim = 1;
featu = [];
nam = [];
mu = [];
sigma = [];
noNormFeatMatrix = [];

for k = 1:numFig
    supLim = k*numObser;
    if supLim > size(features,2)
        supLim = size(features,2);
    end
    feat = features(:,infLim:supLim);
    ff = allFeatures(:,infLim:supLim);
    name = names(infLim:supLim);
    muF = mediaZ(infLim:supLim);
    sigmaF = sigmaZ(infLim:supLim);

    %% Discriminatory ability
    % Clase
    respArtifacts = find(response == 0);
    respGrade3 = find(response == 1);
    respHealthy = find(response == 2);
    % Variables
    varArtifacts = feat(respArtifacts,:);
    varGrade3 = feat(respGrade3,:);
    varHealthy = feat(respHealthy,:);

    %% Mean and Std
    ArtifactMean = []; ArtifactStd = [];
    Grade3Mean = []; Grade3Std = [];
    HealthyMean = []; HealthyStd = [];
    for i = 1:size(varArtifacts,2)
        ArtifactMean(i) = mean(varArtifacts(:,i));
        ArtifactStd(i) = std(varArtifacts(:,i));
    end
    for i = 1:size(varGrade3,2)
        Grade3Mean(i) = mean(varGrade3(:,i));
        Grade3Std(i) = std(varGrade3(:,i));
    end
    for i = 1:size(varHealthy,2)
        HealthyMean(i) = mean(varHealthy(:,i));
        HealthyStd(i) = mean(varHealthy(:,i));
    end

    if verbose == 1
        h = figure
        errorbar(ArtifactMean,ArtifactStd, 'Color', 'b');
        hold on, errorbar(Grade3Mean,Grade3Std, 'Color', 'r');
        errorbar(HealthyMean, HealthyStd, 'Color', 'g');
        title(['Mean and Std of features ' num2str(infLim) ' - ' num2str(supLim)]);
        grid on; legend('ARTIFACTS', 'PATHOLOGICAL GLANDS', 'HEALTHY GLANDS');
        waitfor(h);
    end

    if verbose == 1
        %% Boxplot
        pos1 = (5:5:5*size(varArtifacts,2))';
        pos2 = (6.5:5:5*size(varArtifacts,2)+1.5)';
        pos3 = (8:5:5*size(varArtifacts,2)+3)';
        h = figure, boxplot(varArtifacts, 'colors', 'b', 'Positions',pos1, 'Width',0.8);
        outlier1 = findobj(gca, 'tag', 'Outliers');
        delete(outlier1)
        hold on, boxplot(varGrade3, 'colors', 'r', 'Positions',pos2, 'Width',0.8);
        outlier2 = findobj(gca, 'tag', 'Outliers');
        delete(outlier2)
        boxplot(varHealthy, 'colors', 'g', 'Positions',pos3, 'Width',0.8);
        outlier3 = findobj(gca, 'tag', 'Outliers');
        delete(outlier3)

        boxes = findobj(gca, 'Tag', 'Box');
        grid on; legend(boxes([end 11 1]), 'ARTIFACTS', 'PATHOLOGICAL GLANDS', 'HEALTHY GLANDS')
        title(['Boxplot of features ' num2str(infLim) ' - ' num2str(supLim)]);

        % Percentiles
        q3_art = prctile(varArtifacts,75);          q3_grade3 = prctile(varGrade3,75);          q3_healthy = prctile(varHealthy,75);
        q1_art = prctile(varArtifacts, 25);         q1_grade3 = prctile(varGrade3, 25);         q1_healthy = prctile(varHealthy, 25);
        RIC_art = iqr(varArtifacts);                RIC_g3 = iqr(varGrade3);                    RIC_healthy = iqr(varHealthy);
        b_inf_art = q1_art-1.5.*RIC_art;            b_inf_g3 = q1_grade3-1.5.*RIC_g3;           b_inf_healthy = q1_healthy-1.5.*RIC_healthy;
        b_sup_art = q3_art+1.5.*RIC_art;            b_sup_g3 = q3_grade3+1.5.*RIC_g3;           b_sup_healthy = q3_healthy+1.5.*RIC_healthy;

        minim = min(min([b_inf_art,b_inf_g3,b_inf_healthy])) - 0.5;
        maxim = max(max([b_sup_art,b_sup_g3,b_sup_healthy])) + 0.5;
        ylim([minim,maxim]);
        waitfor(h);
    end

    %% Study of normality of variables
    p_norm = [];
    h_norm = [];
    if verbose == 1
       h = figure
    end
    for j = 1:size(feat,2)
        [~,p_norm(j)] = kstest(feat(:,j)); 
        if p_norm(j)<=alpha % if p < 0.00001 Rejecth H0. 
            h_norm(j) = 1; % No normal variables;
        else
            h_norm(j) = 0; % Normal variables;
        end
        if verbose == 1
            %subplot(6,5,j); qqplot(feat(:,j)); title(name{j});
            sq = ceil(sqrt(numObser))*floor(sqrt(numObser));
            if sq >= numObser
                subplot(ceil(sqrt(numObser)),floor(sqrt(numObser)),j); histfit(feat(:,j)); title(name{j}); 
            else
                subplot(ceil(sqrt(numObser)),floor(sqrt(numObser))+1,j); histfit(feat(:,j)); title(name{j}); 
            end
        end
    end
    infLim = supLim+1;
    if verbose == 1
        waitfor(h)
    end
    
    %% Discriminatority depending on the variables are normal or not
    p = [];
    discrim = [];
    for j = 1:size(h_norm,2)
        if h_norm(j) == 1 % Comparison of medians (kruskalwallis)
            [p(j), v, stats] = kruskalwallis(feat(:,j),response,'off');
            if p(j)<=alpha
                discrim(j) = 1;
            else
                discrim(j) = 0;
            end
        else % Comparison of means (anova)
            [p(j), v, stats] = anova1(feat(:,j),response,'off');
            if p(j)<=alpha
                discrim(j) = 1;
            else
                discrim(j) = 0;
            end
        end
    end
    p(discrim==0) = [];
    feat(:,discrim==0) = [];
    featu = [featu feat];
    ff(:,discrim==0) = [];
    noNormFeatMatrix = [noNormFeatMatrix ff];
    name(discrim==0) = [];
    nam = [nam; name];
    muF(discrim==0) = [];
    mu = [mu muF];
    sigmaF(discrim==0) = [];
    sigma = [sigma sigmaF];

    featMatrix(:,discrim == 0) = [];
end

%% Correlation
[R,p_cor] = corrcoef(featu);
if verbose == 1
    figure, imagesc(abs(R)); colormap(jet); colorbar;
    title('Correlation of Discriminatory features');
end

%% Discard the correlated features
idx = (abs(R)>th_cor & p_cor<alpha);
mat_tri_sup = triu(idx,1);
[row,col] = find(mat_tri_sup);       
col = unique(col);
featu(:,col) = [];
noNormFeatMatrix(:,col) = [];
mu(col) = [];
sigma(col) = [];
nam(col) = [];
p(col) = [];

for i = 1:size(noNormFeatMatrix,2)
    featM = noNormFeatMatrix(:,i);
    for j = 1:size(allFeatures,2)
        featG = allFeatures(:,j);
        resp = isequal(featG,featM);
        if resp == 1
            idFeat(i) = j;
            break;
        end
    end
end

namesFeat = nam; % Names of selected variables
finalMatrix = [featu, response]; % Data matrix with features and response

        
%% Save
% files = nombres de las boundingbox
% namesFeat = nombres de las características seleccionadas
% finalMatrix = matriz de características + respuesta
% noNormFeatmatrix = matriz de características sin normalizar (y sin respuesta)
% idFeat = índice de características utilizadas
% mu = media
% sigma = desviación estándar
save featuresSelected_v1.mat files namesFeat...
    finalMatrix noNormFeatMatrix idFeat...
    mu sigma;
