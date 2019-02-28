function [label_predictions] = obtainLabelPredictionsNeuralNetwork(score)

label_predictions = [];
sc = [];
for j = 1:size(score,2)
    v = score(:,j);
    v(v==max(v)) = 1;
    v(v~=1) = 0;
    sc(:,j) = v;
end
for k = 1:size(sc,2)
    v = sc(:,k);
    [x,y] = find(v==1);
    if x == 1
        v(1) = 0;
    elseif x==2
        v(1) = 1;
    elseif x==3
        v(1) = 2;
    end
    label_predictions(k) = v(1);
end

end