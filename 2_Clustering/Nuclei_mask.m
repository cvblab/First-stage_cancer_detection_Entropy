function [nucleous, nuclei_post] = Nuclei_mask(img, se_open, se_dil)
    
rng('default');

if nargin<3
    se_open = 20;
    se_dil = 1;
end

%% Gaussian Filter
gau_filter = fspecial('gaussian',200);
img = imfilter(img,gau_filter,'replicate'); 
img = imsharpen(img);
    
%% Parameters definition
img = double(img);
[s1,s2,s3] = size(img);

%% Obtention data
pixels = reshape(img(:),s1*s2,s3);
pixels_ind = randi(length(pixels), round(length(pixels)*0.05), 1)';
data = pixels(pixels_ind, :);

%% Kmeans 
[~, C] = kmeans(data, 4,'Distance','sqeuclidean', 'Replicates', 5);
[~, ind] = sort(sum(C,2));
C = C(ind,:); % Centroid coordinates by order
D = pdist2(pixels, C); 
[~, id] = min(D, [], 2); 
cluster_img = reshape(id,s1,s2);

%% Postprocess
nucleous = cluster_img == 1;
nuclei_mask = imdilate(nucleous, strel('disk',se_dil));
nuclei_mask = bwareaopen(nuclei_mask,se_open);
%nuclei_mask = logical(nuclei_mask);
nuclei_post = logical(nuclei_mask);
