function [bb_img,uni_img, fileNormal] = obtainBBoxCandidates(file, img_seg, img, mask_black)

[lumen_img,lumen_num] = bwlabel(mask_black); % Mask of candidates to lumen

for i = 1:lumen_num
    if lumen_num == 0
        continue; % If there are not any lumen... next!
    end
    bw = lumen_img == i;
    fg = logical(bw);  % Mask of each individual candidate                               

    prop = regionprops(fg,'Centroid', 'Orientation');
    c = cat(1,prop.Centroid);
    if ~isempty(c)
        for k = 1:size(c,1)
            label(k) = img_seg(round(c(k,2)),round(c(k,1)));
        end
        mask = ismember(img_seg,label); 
        mask = imfill(mask, 'holes');
        mask = imopen(mask,strel('disk',2));
        mask = bwareaopen(mask,25); % mask that contains the region of watershed in that centroid
        
        % if mask is out of the watershed region, we use the components around the candidate
        if mask(1,1) == 1 | mask(1,end) == 1 | mask(end,1) == 1 | mask(end,end) == 1
            fground = imdilate(fg,strel('disk',8));
            mask = fground.*double(img);
            prop = regionprops(fground, 'BoundingBox');
            bb = cat(1,prop.BoundingBox);
            bb = round(bb);
            bb(3) = bb(1)+bb(3);
            bb(4) = bb(2)+bb(4);
            bb(bb<1) = 1;
            bb(bb>1024) = 1024;
            bb_img{i} = img(bb(2):bb(4),bb(1):bb(3),:);
            imgUni = mask(bb(2):bb(4),bb(1):bb(3),:);
            uni_img{i} = uint8(imgUni);
            fileNormal{i} = [file(1:end-4) '_x' num2str(round(bb(1))) '_y' num2str(round(bb(2))) '.jpg'];
        
        else % we use the region watershed around the candidate
            prop = regionprops(mask,'BoundingBox');
            bb = cat(1,prop.BoundingBox);
            bb = round(bb);
            bb(3) = bb(1)+bb(3);
            bb(4) = bb(2)+bb(4);
            bb(bb<1) = 1;
            bb(bb>1024) = 1024;
            bb_img{i} = img(bb(2):bb(4),bb(1):bb(3),:);
            bb_mask = mask(bb(2):bb(4),bb(1):bb(3)); % Glándula sin enmascarar
            imgUni = double(bb_img{i}).*bb_mask; 
            uni_img{i} = uint8(imgUni); % Glandula enmascarada
            fileNormal{i} = [file(1:end-4) '_x' num2str(round(bb(1))) '_y' num2str(round(bb(2))) '.jpg'];
        end
    end
end