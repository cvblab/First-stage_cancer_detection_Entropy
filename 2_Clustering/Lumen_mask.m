 function [mask_lumen, mask_black] = Lumen_mask(img, s, LUMEN_SIZE, SE_DIL)    

rng('default');

if nargin<3
    LUMEN_SIZE = 20; % To remove lumen objects that are considered noise
    SE_DIL = 1;
end

%% Kmeans for lumen
[~, C] = kmeans(s(:), 3,'Distance','sqeuclidean', 'Replicates', 5);
[~, ind] = sort(sum(C,2));
C = C(ind,:); % Centroid coordinates by order
D = pdist2(s(:), C); 
[~, id] = min(D, [], 2); 
cluster = reshape(id,size(img,1),size(img,2)); 
lumen = cluster == 1; 

%% Postprocess
lumen_img = double(lumen); 
lumen_img = bwareaopen(lumen_img,LUMEN_SIZE);
se = strel('disk',SE_DIL);
mask_lumen = imdilate(lumen_img,se);

%% Background
mask_black = mask_lumen;
[lumen_img, lumen_num] = bwlabel(mask_black);
for i = 1:lumen_num
    obj = lumen_img == i;
    unic = sum(unique(obj(1,:)) + unique(obj(:,1)) + ...
        unique(obj(end,:)) + unique(obj(:,end)));
    if unic>0
        mask_black(obj) = 0; % Lumen mask processed with the background black
    end
end

%% Display image results
% figure, title(file);
% subplot(1,2,1); imshow(mat2gray(mask_lumen)); title('Mask with background');
% subplot(1,2,2); imshow(mat2gray(mask_black)); title('Mask without background');
% figure, imshow(mat2gray(img)); title(file);
% figure, imshow(mat2gray(mask_black)); title(file);