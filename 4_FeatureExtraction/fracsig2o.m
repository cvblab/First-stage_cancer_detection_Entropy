function [ H,r ] = fracsig2o( signal,varargin )
% @author Thomas JANVIER <thomas.janvier@univ-orleans.fr>
% @date 2014-02-03

%% DESCRIPTION
% FRACSIG2O Compute the fractal signature of the 2D input signal
% along a specified direction
%   ref : K. Harrar, L. Hamami, E. Lespessailles, and R. Jennane, "Piecewise Whittle estimator for trabecular bone radiograph characterization," Biomedical Signal Processing and Control, vol. 8, no. 6, pp. 657-666, Nov. 2013.
%
% INPUTS :
%
%   SIGNAL :
%       n-by-m numerical matrix
%       input signal
%
% OPTIONS :
%
%   METHOD :
%       {'ILE', VAR, BLK, ACF}
%       specify the method to compute the fractal dimension
%
%           - ILE : based on the generalization of the quadratic
%                   variations of any gaussian process
%                   copy of the matlab implementation (cf wfbmesti.m)
%                 ref : J. Istas and G. Lang, "Quadratic variations and estimation of the local Hölder index of a Gaussian process," in Annales de l'Institut Henri Poincare (B) Probability and Statistics, 1997, vol. 33, pp. 407?436.
%
%           - BLK : based on the mathematical morphology
%                   compute the slope of log(A)-vs-log(k)
%                 ref : P. Maragos and F.-K. Sun, "Measuring the Fractal Dimension of Signals: Morphological Covers and Iterative Optimization," IEEE Transactions on Signal Processing, vol. 41, no. 1, pp. 108-121, 1993.
%
%           - VAR : based on the quadratic variations of gaussian noise
%                 ref : the Scientific Committee of the GRIO (Groupe de Recherche et d?Information sur les Ostéoporoses), V. Bousson, C. Bergot, B. Sutter, P. Levitz, and B. Cortet, "Trabecular bone score (TBS): available knowledge, clinical relevance, and future prospects," Osteoporosis International, vol. 23, no. 5, pp. 1489-1501, May 2012.
%
%   DIRECTION :
%       a scalar [0,360[
%       the direction considered (the angle to the horizontal line)
%
%   SCALES :
%       a vector of integers
%       the scale factor to estimate the FS
%
% OUTPUTS :
%
%   H :
%       a vector [0,1]
%       the signal Hurst exponent for each scale
%
%   r :
%       a vector
%       the scale factor (same size as H)

%% INPUT PARSING
p = inputParser;
addRequired(p,'signal',@(x) validateattributes(x,{'numeric'},{'2d','real'}));
addOptional(p,'method','VAR',@(x) any(validatestring(x,{'VARorg','ILEorg','BLKorg','VARgrad','ILEgrad','BLKgrad','WhE','ACF'})));
addOptional(p,'direction',0,@(x) validateattributes(x,{'numeric'},{'real','scalar','>=',0}));
addOptional(p,'scales',0,@(x) validateattributes(x,{'numeric'},{'vector','integer','nonnegative'}));
parse(p, signal, varargin{:});
method = p.Results.method;

%% Rotate the image
alpha = p.Results.direction;
if alpha == 90
    signal = fliplr(signal');
elseif alpha == 180
    signal = fliplr(signal);
elseif alpha == 270
    signal = flipud(signal');
elseif alpha ~= 0
    signal = imrotation(signal,alpha,'nearest');
end

if p.Results.scales==0
    r = (2:(size(signal,2)-1)/2)';
else
    r = p.Results.scales(:);
end

%% Mathematical Morphology (blanket method)
if strcmpi(method,'BLKorg')
    se = strel('line',3,0);
    Ao = mean(mean(imdilate(signal,se)-imerode(signal,se)));
    % for each scales
    A = zeros(size(r));
    parfor i=1:length(r)
        se = strel('line',2*r(i)+1,0);
        % dilate the signal
        ub = imdilate(signal,se);
        % erode the signal
        lb = imerode(signal,se);
        % compute the area between the erosion and the dilataion
        A(i) = mean(ub(:)-lb(:));
    end
    H = (log(A)-log(Ao))./log(r);
    
    %% Mathematical Morphology (blanket method)
elseif strcmpi(method,'BLKgrad')
    se = strel('line',3,0);
    Ao = mean(mean(imdilate(signal,se)-imerode(signal,se)));
    % for each scales
    A = zeros(size(r));
    parfor i=1:length(r)
        se = strel('line',2*r(i)+1,0);
        % dilate the signal
        ub = imdilate(signal,se);
        % erode the signal
        lb = imerode(signal,se);
        % compute the area between the erosion and the dilataion
        A(i) = mean(ub(:)-lb(:));
    end
    H = (log(A(1:end-1))-log(A(2:end)))./(log(r(1:end-1))-log(r(2:end)));
    
    %% Istas & Lang method (generalized quadratic variations)
elseif strcmpi(method,'ILEorg')
    V = zeros(size(r));
    % compute Var for r=1
    Vo = mean(mean((2.*signal(:,2:end-1)-signal(:,1:end-2)-signal(:,3:end)).^2,2));
    % compute Var for all r
    parfor i=1:length(r)
        V(i) = mean(mean((2.*signal(:,r(i)+1:end-r(i))-signal(:,1:end-2*r(i))-signal(:,2*r(i)+1:end)).^2,2));
    end
    H = 0.5.*(log(V)-log(Vo))./log(r);
    
elseif strcmpi(method,'ILEgrad')
    V = zeros(size(r));
    % compute Var for all r
    parfor i=1:length(r)
        V(i) = mean(mean((2.*signal(:,r(i)+1:end-r(i))-signal(:,1:end-2*r(i))-signal(:,2*r(i)+1:end)).^2,2));
    end
    H = 0.5.*((log(V(1:end-1))-log(V(2:end)))./(log(r(1:end-1))-log(r(2:end))));
    
    %% VAR (quadratic variations)
elseif strcmpi(method,'VARorg')
    V = zeros(size(r));
    % compute Var for r=1
    Vo = mean(mean((signal(:,1:end-1)-signal(:,2:end)).^2,2));
    % compute Var for all r
    parfor i=1:length(r)
        V(i) = mean(mean((signal(:,1:end-r(i))-signal(:,1+r(i):end)).^2,2));
    end
    H = 0.5.*(log(V)-log(Vo))./log(r);
    
elseif strcmpi(method,'VARgrad')
    V = zeros(size(r));
    % compute Var for r=1
    Vo = mean(mean((signal(:,1:end-1)-signal(:,2:end)).^2,2));
    % compute Var for all r
    parfor i=1:length(r)
        V(i) = mean(mean((signal(:,1:end-r(i))-signal(:,1+r(i):end)).^2,2));
    end
    H = 0.5.*((log(V(1:end-1))-log(V(2:end)))./(log(r(1:end-1))-log(r(2:end))));
    % safeguard
    %     H = min(max(H,0),1);
    
    %% Auto Correlation Function using least squares
elseif strcmpi(method,'ACF')
    parfor i=1:length(r)
        k = r(i);
        % compute the increments = fractionnal Gaussian noise
        fGn = signal(:,k+1:end)-signal(:,1:end-k);
        [ny,nx] = size(fGn);
        % estimate the PSD of the fGn
        [P,f] = periodogram1(fGn(1,:));
        for y=2:ny
            P = P + periodogram1(fGn(y,:));
        end
        P = P./ny;
        % remove the 0 frequency (avoid computational errors)
        P = P(2:end);
        % compute variables
        nh = fix(nx/2);
        iseven = -nx+2*nh+1;
        % sample frequencies
        fe = [0:nh,nh-iseven:-1:1]';
        % autocorrelation function
        R = @(h) ((fe+k).^(2*h)-2.*fe.^(2*h)+abs(fe-k).^(2*h));
        % associated psd
        F = @(h) gamma(1-2*h).*cos(pi*h)./(2*pi*h).*R(h);
        % real part of psd
        T = @(H) imcrop(abs(fft(F(H))),[1,2,1,length(P)-1]);
        % normalization parameter
        c = @(H) mean(mean(P))/mean(mean(T(H)));
        % rescale the theoretical psd
        cT = @(H) T(H).*c(H);
        % define the Likehod Function (lf)
        lf = @(H) nansum((log(cT(H)./P)).^2); % Least Squares
        % maximize the WLF (here we minimize -WLF)
        H(i) = fminbnd(lf,0,1);
    end
    
    %% Auto Correlation Function using Whittle
elseif strcmpi(method,'WhE')
    parfor i=1:length(r)
        k = r(i);
        % compute the increments = fractionnal Gaussian noise
        fGn = signal(:,k+1:end)-signal(:,1:end-k);
        [ny,nx] = size(fGn);
        % estimate the PSD of the fGn
        [P,f] = periodogram1(fGn(1,:));
        for y=2:ny
            P = P + periodogram1(fGn(y,:));
        end
        P = P./ny;
        % remove the 0 frequency (avoid computational errors)
        P = P(2:end);
        % compute variables
        nh = fix(nx/2);
        iseven = -nx+2*nh+1;
        % sample frequencies
        fe = [0:nh,nh-iseven:-1:1]';
        % autocorrelation function
        R = @(h) ((fe+k).^(2*h)-2.*fe.^(2*h)+abs(fe-k).^(2*h));
        % associated psd
        F = @(h) gamma(1-2*h).*cos(pi*h)./(2*pi*h).*R(h);
        % real part of psd
        T = @(H) imcrop(abs(fft(F(H))),[1,2,1,length(P)-1]);
        % normalization parameter
        c = @(H) mean(mean(P))/mean(mean(T(H)));
        % rescale the theoretical psd
        cT = @(H) T(H).*c(H);
        % define the Likehod Function (lf)
        wlf = @(H) -sum(-log(cT(H))-(P./cT(H)));
        % maximize the WLF (here we minimize -WLF)
        H(i) = fminbnd(wlf,0,1);
    end
end
% safeguard
%     H = min(max(H,0),1);
H = reshape(H,1,[]);
r = reshape(r,1,[]);
end