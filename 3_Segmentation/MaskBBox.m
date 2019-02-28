function [maskImg, BBOX] = MaskBBox(file, candidatesMask, img) 

for i = 1:length(candidatesMask)
    fileNormal = file{i};
    imgRGBmask = candidatesMask{i};
    imgMask = zeros(size(imgRGBmask,1),size(imgRGBmask,2));
    id = find(imgRGBmask(:,:,1)>15 & imgRGBmask(:,:,2) >15 & imgRGBmask(:,:,3) >15);
    imgMask(id) = 1;
    imgMask = imfill(imgMask, 'holes');
    imgMask = bwareaopen(imgMask, 20, 4);
    imgMask = logical(imgMask);
    w = ones(size(imgMask));
    pr = regionprops(w, 'Centroid');
    cent = cat(1,pr.Centroid);
    [maskimg, numimg] = bwlabel(imgMask);
    if numimg>1 % if the number of elements is greater than one we choose the center object
        dist = [];
        for t = 1:numimg
            elem = maskimg == t;
            bordesElem = bwboundaries(elem);
            bordesElem = bordesElem{1};
            distBor = [];
            for tt = 1:length(bordesElem)
                distBor(tt) = pdist2(bordesElem(tt,:),cent, 'Euclidean');
            end
            idD = find(distBor == min(distBor),1,'first');
            coord = bordesElem(idD,:);
            dist(t) = pdist2(coord,cent, 'Euclidean');
        end
        id_dist = find(dist == min(dist));
        imgMask = maskimg == id_dist;
    end
    if max(unique(imgMask)) == 0 % if some error exists, next!
        continue;
    end
    maskImg{i} = imgMask;
    
    propGland = regionprops(imgMask, 'Area', 'BoundingBox');
    ind = find(fileNormal == 'x', 1, 'last');
    indx = find(fileNormal == 'y', 1, 'last');
    bb = cat(1,propGland.BoundingBox);
    bb(1) = str2num(fileNormal(ind+1:indx-2)); 
    bb(2) = str2num(fileNormal(indx+1:end-4));
    bb(3) = bb(1)+bb(3);
    bb(4) = bb(2)+bb(4);
    bb(bb<1) = 1;
%     if bb(1)>size(img,2)
%         bb(1) = size(img,2);
%     end
%     if bb(2)>size(img,1) 
%         bb(2) = size(img,1);
%     end
    if bb(3)>size(img,2)
        bb(3) = size(img,2);
    end
    if bb(4)>size(img,1)
        bb(4) = size(img,1);
    end
    BBOX(i,:) = bb;
end