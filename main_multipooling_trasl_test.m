clear all;
close all;




%% Random patches

%% Parameters


n_ori = 8;
res_ori = 2*pi / n_ori;
ori = (0:(n_ori-1))*res_ori;

n_scales = 3;
taps = [11 23 47];

n_transformations = n_ori * n_scales;

n_templates = 10;

n_filters = n_templates * n_transformations;

%% Source images

load('compute_templates/pascal_filters.mat');

%% Gabors 

%load('compute_templates/gabor_filters.mat');

%% Parameters


%% Import, normalize and zero-center the input image
inputImg1 = double(rgb2gray(imread('lena.jpg','jpg')));
inputImg1  = inputImg1 - mean(mean(inputImg1));
inputImg1 = inputImg1 ./ norm(inputImg1, 1);
[inSizeXini inSizeYini] = size(inputImg1);

signatures = [];

numTraslX = 4; % Number of translations along X of the input image 
rangeTraslX = 100; % Range of translations along X of the input image
numTraslY = 4; % Number of translations along Y of the input image 
rangeTraslY = 100; % Range of translations along Y of the input image

for i = 1:numTrasl
    inputImg = imtranslate(I,[15, 25]);
    
    inputImg = imcrop( inputImg ,[ ceil( (1/2 - (1/(2*sqrt(2))))* inSizeXini) ceil( (1/2 - (1/(2*sqrt(2))))* inSizeYini) floor((1/sqrt(2))* inSizeXini) floor((1/sqrt(2))* inSizeYini) ]);
    [inSizeX inSizeY] = size(inputImg);    
    
    figure(1)
    subplot(4, 4, i);
    imshow(inputImg, []);
    % Simple layer

    filteredImg = filtering( inputImg , templates );

    % figure;
    % imshow(cell2mat(filteredImg(1,1,1)) , []);

    % Complex layer

    % Image pooling
    poolingSplitNum = 1;    % Number of splits along axes for determining pooling regions
    numBars = 20;   % Number of bars of the histogram

    L1hist = pooling( filteredImg , n_templates , poolingSplitNum , numBars , 'histogram');

    %bar(L1hist{2,1,2},L1hist{2,1,1});

    % Generate pooling area 1 signature by concatenating histograms generated
    % by the filtering with all templates
    signatures = [signatures ; horzcat(L1hist{:,1,1})];
    
end

figure(2)
plot(signatures');

%% Response visualization and quantitative comaprison