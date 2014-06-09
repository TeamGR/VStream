% This script is aimed at computing the layer 2 templates

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

n_templates = 2;
n_templates2 = 10;

n_filters = n_templates * n_transformations;

%% Source images

load('compute_templates/pascal_filters.mat');
templatesL2 = templates;
clear tamplates;

%% Gabors 

load('compute_templates/gabor_filters.mat');

%% Parameters


%% 

signatures = [];
filteredTemplatesL2 = cell(n_scales, n_templates, n_ori);
L1hist = cell(size(templatesL2,1),size(templatesL2{1,1},2));

for i = 1:n_templates2

%     inputImg = imcrop( inputImg ,[ ceil(rangeTraslX/2) ceil(rangeTraslX/2) floor(inSizeXini - rangeTraslX) floor(inSizeYini - rangeTraslY) ]);
%     [inSizeX inSizeY] = size(inputImg);    
    
    % Scale loop
    for j = 1:size(templatesL2,1)
        % Orientation loop
        for k = 1:size(templatesL2{1,1},2)
            inImg = squeeze(templatesL2{j,1}(i , k , : , :));
            filteredTemplatesL2{j, i, k} = filtering( inImg , templates );
        end
    end

    % Complex layer

    % Image pooling
    poolingSplitNum = 1;    % Number of splits along axes for determining pooling regions
    numBars = 20;   % Number of bars of the histogram
    % Scale loop
    for j = 1:size(templatesL2,1)
        % Orientation loop
        for k = 1:size(templatesL2{1,1},2)
            L1hist{j,k} = pooling( filteredTemplatesL2 , n_templates , poolingSplitNum , numBars , 'histogram');
        end
    end

    %bar(L1hist{2,1,2},L1hist{2,1,1});

    % Generate pooling area 1 signature by concatenating histograms generated
    % by the filtering with all templates
    signatures = [signatures ; horzcat(L1hist{:,1,1})];
    
end

%% Response visualization and quantitative comaprison

figure(2)
m = mean(signatures);
sd = std(signatures);
f = [ m+2*sd , flipdim(m-2*sd,2)]; 
fill([1:size(signatures,2) , size(signatures,2):-1:1] , f, [7 7 7]/8)
hold on;
plot(1:size(signatures,2) , m , 'b' , 'LineWidth',1);

max(sd)

