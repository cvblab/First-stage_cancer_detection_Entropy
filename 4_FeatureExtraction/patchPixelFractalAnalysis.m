function [H, alpha, P, f] = patchPixelFractalAnalysis(I, mask, mand, optH)
% @author Adrian COLOMER <adcogra@i3b.upv.es>
% @date 2017-02-01

%% Mask composed by external mask and optic disk mask
totalMask = mask.extMask.*not(mask.diskMask);

%% Grid in rows and columns
gridC = 1 : mand.stepC : size(I,2);
gridR = 1 : mand.stepR : size(I,1);

%% Avoiding warnings of parfor loop
bb1 = mand.bb1;
bb2 = mand.bb2;

%% H vector inicialization
H = zeros(length(gridC)*length(gridR),2)-1;
alpha = zeros(1,length(gridC)*length(gridR))-1;
P = [];
f = [];

ind = 1;

for c = 1 : length(gridC)
    for r = 1 : length(gridR)
        
        %Taking into account the bounds
        if (gridR(r) - floor(bb1/2) >= 1 && gridC(c) - floor(bb2/2) >= 1 &&...
                gridR(r) + floor(bb1/2) <= size(I,1) &&...
                gridC(c) + floor(bb2/2) <= size(I,2))
            winR = gridR(r) - floor(bb1/2) : gridR(r) + floor(bb1/2);
            winC = gridC(c) - floor(bb2/2) : gridC(c) + floor(bb2/2);
        else
            % H = -1 when I am out of bounds (inicialization)
            ind = ind + 1;
            continue;
        end
        
        %Image and Mask Patch
        imPatch = I(winR,winC);
        imPatchTM = totalMask(winR,winC);
        if(sum(imPatchTM(:)) == size(imPatchTM,1)*size(imPatchTM,2))
            %Computing fractal dimension for each ROI totally located
            %inside the FOV and W/O optic disk pixels
            [H(ind,1), ~] = fracdim2o(imPatch, optH.method, optH.direction(1));
            [H(ind,2), ~] = fracdim2o(imPatch, optH.method, optH.direction(2));
        else %
            ind = ind + 1;
            continue;
        end
        ind = ind + 1;
    end
end