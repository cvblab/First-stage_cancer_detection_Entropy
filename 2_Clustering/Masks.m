function [cytoplasm_mask, stroma_mask, over2] = Masks(img, mask_black, nuclei_mask, cyan, se_open, se_dil)

rng('default');

if nargin<5
    se_open = 20;
    se_dil = 1;
end

%% Preprocess
c = cyan;
c(nuclei_mask == 1) = 0;

%% Kmeans to obtain stroma mask
[~, C] = kmeans(c(:), 3,'Distance','sqeuclidean', 'Replicates', 5);
[~, ind] = sort(sum(C,2));
C = C(ind,:); % Centroid coordinates by order
D = pdist2(c(:), C); 
[~, id] = min(D, [], 2); 
cluster = reshape(id,size(img,1),size(img,2));

%% Stroma postprocess
stroma_mask = cluster == 2;
mask_black = imfill(mask_black,'holes');
strom = stroma_mask-mask_black;
strom = strom == 1;
st = bwareaopen(strom,10);
st = imopen(st,strel('disk',se_dil));
stroma_mask = logical(st);

%% Cytoplasm postprocess
cytoplasm_mask = cluster == 3;
cytoplasm_mask = imdilate(cytoplasm_mask, strel('disk',se_dil));
cytoplasm_mask = cytoplasm_mask - mask_black;
cyto = cytoplasm_mask == 1;
cyto = cyto-nuclei_mask;
cyto = cyto-stroma_mask;
cyto = cyto == 1;
cyto = bwareaopen(cyto, se_open);
cytoplasm_mask = logical(cyto);

%% Overlay
over = imoverlay(img, cytoplasm_mask, 'green');
over2 = imoverlay(over, nuclei_mask, 'black');
% over3 = imoverlay(over2, stroma_mask, 'magenta');
% over4 = imoverlay(over3, mask_black, 'white');
