function [R,G,B,cyan,s,H,magenta] = Chanel_color(ima)

%% RGB
R = double(ima(:,:,1)); R = R./max(unique(R));
G = double(ima(:,:,2)); G = G./max(unique(G));
B = double(ima(:,:,3)); B = B./max(unique(B));

%% CMYK
CMYK = makecform('srgb2cmyk');
a = applycform(ima,CMYK);
cyan = double(a(:,:,1));
cyan = cyan./max(unique(cyan)); %cyan
magenta = double(a(:,:,2));
magenta = magenta./max(unique(magenta)); %magenta

%% HSV
hsv = rgb2hsv(ima);
s = hsv(:,:,2);

%% Colour Deconvolution
HEtype = 'H&E';
[H,E,~] = colour_deconvolution(uint8(ima), HEtype); % Colour Deconvolution


