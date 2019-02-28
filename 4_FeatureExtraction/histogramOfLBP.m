function [hist_LBP, hist_LBP_VAR] = histogramOfLBP(LBP, VAR, Mask_Adjusted)

 %% LBP histogram
 
 LBPtissue = LBP(Mask_Adjusted);
 my_hist = hist(LBPtissue,0:9);

if(sum(my_hist)>0)
    hist_LBP = my_hist/sum(my_hist);
else
    hist_LBP = my_hist;
end

 %% LBP/VAR histogram

 VARtissue = VAR(Mask_Adjusted);
 h_LBP_VAR = zeros(1,length(hist_LBP));
 for i=1:length(hist_LBP)
    idx = find(LBPtissue==(i-1));
    h_LBP_VAR(i) = sum(VARtissue(idx));
 end

if(sum(h_LBP_VAR)>0)
    hist_LBP_VAR = h_LBP_VAR/sum(h_LBP_VAR);
else
    hist_LBP_VAR = h_LBP_VAR;
end 
    
        
end