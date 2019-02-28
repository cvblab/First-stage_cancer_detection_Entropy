function imout = watershedConstrained(imin,fg,bg,tamFg,tamBg)


%% Inputs
% Specify constriction (0 = Constrained; 1 = Unconstrained)
unconstrainedFg = 0; % Foreground
unconstrainedBg = 0; % Background

% Radius of structuring element (disk)
% tamFg = 2; % Foreground
% tamBg = 2; % Background 

% Scale for image resizing;
% scale = .25;

% Function for computation of priority
opcPriority = 3; % 1=mean; 2=median; 3=max(recommended)


%% Code
[a,b,c] = size(imin);

%% Imposición de mínimo
controlImage = imimposemin(imin,(fg | bg));

%% Cambio a formato double
controlImage = double(controlImage);

%% Calculo de minimos regionales y etiquetado inicial
minReg = imregionalmin(controlImage); 
% figure, imshow(minReg)
conn = 8; % Conectividad de mínimos regionales (OJO!!!!)
centerMap = bwlabel(minReg,conn); % Los píxeles NO etiquetados valen "0"
coverMapOrig = centerMap;
% Eliminación de los píxeles interiores de los mínimos regionales
borMinReg = bwmorph(minReg,'remove');
% figure, imshow(borMinReg), title('borMinReg')

%% Identificación de etiquetas de foreground y background
% Foreground
labFg = centerMap(fg);
labFg= unique(labFg);
% Background
labBg = centerMap(bg);
labBg= unique(labBg);

% return
%% OFFSETS
% Offsets de los índices de los vecinos en expansión (Igual que Hill Climbing)
se = strel('square',3);
offsetsVect = getneighbors(se);
offsetsVect(and(offsetsVect(:,1)==0,offsetsVect(:,2)==0),:)=[]; % Eliminación de pixel central
offsets = offsetsVect(:,1)+offsetsVect(:,2)*a;

% Offsets de los vecinos de la zona "Constrained" (MODIFICACIÓN RESPECTO A HILL CLIMBING)
% Foreground
if unconstrainedFg==1 % Caso unconstrained
    offsetsFg = [];
else    
    seFg = strel('disk',tamFg,0);
    offsetsFgVect = getneighbors(seFg);
    offsetsFg = offsetsFgVect(:,1)+offsetsFgVect(:,2)*a;
end

% Background
if unconstrainedBg==1 % Caso unconstrained
    offsetsBg = [];
else    
    seBg = strel('disk',tamBg,0);
    offsetsBgVect = getneighbors(seBg);
    offsetsBg = offsetsBgVect(:,1)+offsetsBgVect(:,2)*a;
end


%% centerMap
[rowIni, colIni] = find(borMinReg>0); % Bordes de los mínimos regionales
posIni = rowIni+(colIni-1)*a;
centerMap(posIni) = -centerMap(posIni);
coverMap = centerMap; % Inicialización de coverMap

% Cambios de bordes de centerMap
if unconstrainedFg == 0 && unconstrainedBg == 0 
    tamBor = max(tamFg,tamBg);
elseif unconstrainedFg == 1 && unconstrainedBg == 0 
    tamBor = tamBg;
elseif unconstrainedFg == 0 && unconstrainedBg == 1 
    tamBor = tamFg;
elseif unconstrainedFg == 1 && unconstrainedBg == 1 
    tamBor = 1;
else
    display('Error: Valores erróneos de "unconstrained"')
    return
end

   
%% coverMap
% En coverMap, las zonas de dilatación de los píxeles iniciales (si la zona 
% NO es unconstrained) se etiquetan con la ID del píxel
for i = 1:length(rowIni)
    ID = coverMapOrig(posIni(i));
    if ismember(ID,labFg) && unconstrainedFg == 0 
        offs = offsetsFgVect;
    elseif ismember(ID,labBg) && unconstrainedBg == 0
        offs = offsetsBgVect;
    else
        continue
    end
    posRow = rowIni(i)+offs(:,1);
    posCol = colIni(i)+offs(:,2);
    
    indCorr = (posRow>=1 & posRow<=a & posCol>=1 & posCol<=b);
    posRow = posRow(indCorr);
    posCol = posCol(indCorr);
    
    posLin = posRow+(posCol-1)*a;
    
    coverMap(posLin) = ID;
%     figure, imagesc(coverMap)
end

labBor = max(max(centerMap(:)),max(coverMap(:)))+1;

% Cambios de bordes de centerMap
centerMap = padarray(centerMap,[tamBor tamBor],labBor); %%%%% PADARRAY

% Cambios de bordes de coverMap
coverMap = padarray(coverMap,[tamBor tamBor],labBor); %%%%% PADARRAY

%% Bordes para controlImage
controlImage = padarray(controlImage,[tamBor tamBor],'replicate'); %%%%% PADARRAY

%% OFFSETS
% Estudio de dimensiones
[a,b,c] = size(controlImage);

% Offsets de los índices de los vecinos en expansión (Igual que Hill Climbing)
se = strel('square',3);
offsetsVect = getneighbors(se);
offsetsVect(and(offsetsVect(:,1)==0,offsetsVect(:,2)==0),:)=[]; % Eliminación de pixel central
offsets = offsetsVect(:,1)+offsetsVect(:,2)*a;

% Offsets de los vecinos de la zona "Constrained" (MODIFICACIÓN RESPECTO A HILL CLIMBING)
% Foreground
if unconstrainedFg==1 % Caso unconstrained
    offsetsFg = [];
else    
    seFg = strel('disk',tamFg);
    offsetsFgVect = getneighbors(seFg);
    offsetsFg = offsetsFgVect(:,1)+offsetsFgVect(:,2)*a;
end

% Background
if unconstrainedBg==1 % Caso unconstrained
    offsetsBg = [];
else    
    seBg = strel('disk',tamBg);
    offsetsBgVect = getneighbors(seBg);
    offsetsBg = offsetsBgVect(:,1)+offsetsBgVect(:,2)*a;
end

%% MEX
mexStart = tic;

% Imlementacion en C
[centerMap, coverMap]=cws(controlImage, centerMap, coverMap, labFg, labBg, offsetsFg, offsetsBg, opcPriority);

mextime = toc(mexStart);
mextime = ceil(mextime);


%% Eliminación de bordes, SOLO CASO CON PADARRAY
% CenterMap
centerMap(1:tamBor,:) = []; % Borde izquierdo
centerMap(end-tamBor+1:end,:) = []; % Borde derecho
centerMap(:,1:tamBor) = []; % Borde superior
centerMap(:,end-tamBor+1:end) = []; % Borde inferior

% CoverMap
coverMap(1:tamBor,:) = []; % Borde izquierdo
coverMap(end-tamBor+1:end,:) = []; % Borde derecho
coverMap(:,1:tamBor) = []; % Borde superior
coverMap(:,end-tamBor+1:end) = []; % Borde inferior

% controlImage
controlImage(1:tamBor,:) = []; % Borde izquierdo
controlImage(end-tamBor+1:end,:) = []; % Borde derecho
controlImage(:,1:tamBor) = []; % Borde superior
controlImage(:,end-tamBor+1:end) = []; % Borde inferior


%% Visualización de resultados finales
imout = coverMap;
imout(imout<0) = max(max(coverMap))+1;
RGB = label2rgb(imout);
% figure, imshow(RGB), title('coverMap final')


imout2 = centerMap;
imout2(imout2<0) = max(max(centerMap))+1;
RGB2 = label2rgb(imout2);
% figure, imshow(RGB2), title('centerMap final')

%% Imagen final
imout =  cat(3, centerMap, coverMap);
imout = max(imout,[],3);
RGB = label2rgb(imout);
% figure, imshow(RGB), title({'Final segmented image'; '(different coloured regions)'})


%% Cálculo de fronteras
% Método 1, con gradiente
% grad = imgradient(centerMap,'IntermediateDifference');
% grad = grad>0;
% grad = bwmorph(grad>0,'skel',Inf);
% % Puesta a 0 de líneas watershed en imagen de salida
% centerMap(grad) = 0;

% Método 2, residuo de dilatación
res = imdilate(imout,strel('disk',1))-imout;
imout(res>0) = 0;

% figure
% imshow(uint8(I))
% % imshow(controlImage)
% hold on
% [x, y] = find(imout==0);
% plot(y,x,'*g');

return


%% Ensanchamiento de fronteras
fronteras = imout==0;
% fronteras = imdilate(fronteras,strel('disk',2));
[x, y] = find(fronteras==1);

%% Imagen con fronteras superpuestas en verde
if size(I,3)==1
    imFin = repmat(I,1,1,3);
elseif size(I,3)==3
    imFin = I;
end

for i = 1:length(x)
    imFin(x(i),y(i),1) = 0; imFin(x(i),y(i),2) = 255; imFin(x(i),y(i),3) = 0;
end

figure, imshow(imFin), 
title({'Final segmented image'; '(watershed lines)'; ['Tiempo total: ' num2str(eltime) ' segundos'];...
      ['tamFg=' num2str(tamFg) ', tamBg=' num2str(tamBg)]})

