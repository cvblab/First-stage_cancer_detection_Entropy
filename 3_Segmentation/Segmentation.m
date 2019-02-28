function [img_seg] = Segmentation(img,over2,nuclei_mask,mask_lumen,mask_black,stroma_mask)

%whole = tic;

%% Parameters watershed
tamFg = 1;
tamBg = 5;

% Format
img = double(img);                                   %img = imresize(img,0.5);
over2 = double(over2);                               %over2 = imresize(over2, 0.5);

%% Entry image (restrictions)       
imin = uint8(mat2gray(nuclei_mask)*255);             %imin = imresize(mask_nuclei,0.5);

%% External marker
mark_ext = mask_lumen-mask_black;
mark_ext(mark_ext==-1) = 0;
mark_ext = imerode(mark_ext,strel('disk',1));
mark_ext = mark_ext+stroma_mask;
mark_ext = logical(mark_ext);                        %mark_ext = imresize(mark_ext, 0.5); 
bg = false(size(imin));
pos = find(mark_ext > 0);
bg(pos) = true; 

%% Internal marker
mark_int = mask_black;
fg = false(size(imin));
pos = find(mark_int > 0);
fg(pos) = true; 

%% %%%%%%%%%%  WatershedConstraint  %%%%%%%%%%%%%%%%%
img_seg = watershedConstrained(imin,fg,bg,tamFg,tamBg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time = toc(whole)
%time = ceil(time);