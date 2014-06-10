%% Experiment description
%
% Number of layers of the net: 1
% Employed dataset: SLBC
% Transformations:
%   3D rotations on a uniform background 

close all;
clear all;
addpath(genpath('.'));

%% Load templates of layer 1

%T = load('gabor_filters.mat');
% gabors = cell(1);
% gabors{1} = T.templates{1}(1,1,:,:);
% gabors = T.templates;
T = load('pascal_filters.mat');
templates=T.templates;

% set parameters for histogram computation at C1 layer
n_splits = 1; % the image is not divided in a single layer architecture
n_bins = 20; % bars of the histogram
n_scales = length(templates);
n_templates = size(templates{1},1);
n_ori = size(templates{1},2);

% %% Load templates of layer 2
% 
% T = load('templatesL2.mat');
% templatesL2_hist = T.templatesL2_hist;
% %templatesL2_moms = T.templatesL2_moms;

% Load images
imgType = '*.jpg'; % change based on image type
n_cars = 35;
cars = cell(n_cars,1);
imgPathCars = 'D:\IIT\CODICI\pipeline_overfeat\SBLC\images\car';
car_dir  = dir([imgPathCars '\' imgType]);
for idx = 1:n_cars
    cars{idx} = double(rgb2gray(imread(fullfile(imgPathCars, car_dir(idx).name))));
end

% Init responses
S1cars = cell(n_cars,1);
C1cars = cell(n_cars,1);

% Load images
imgType = '*.jpg'; % change based on image type
n_airplanes = 35;
airplanes = cell(n_airplanes,1);
imgPathAirplanes = 'D:\IIT\CODICI\pipeline_overfeat\SBLC\images\airplane';
airplane_dir  = dir([imgPathAirplanes '\' imgType]);
for idx = 1:n_airplanes
    airplanes{idx} = double(rgb2gray(imread(fullfile(imgPathAirplanes, airplane_dir(idx).name))));
end

% Init responses
S1airplanes = cell(n_airplanes,1);
C1airplanes = cell(n_airplanes,1);

% S1 responses 
% output of dotproductL1_giulia is of format
% cell(n_scales, 1)
% responses{idx_scale} = zeros(n_templates, n_ori, sizeY-(taps(idx_scale)-1), sizeX-(taps(idx_scale)-1))
for idx_car=1:n_cars
    S1cars{idx_car} = dotproductL1_giulia(cars{idx_car}, templates, n_splits);
end

% S1 responses 
% output of dotproductL1_giulia is of format
% cell(n_scales, 1)
% responses{idx_scale} = zeros(n_templates, n_ori, sizeY-(taps(idx_scale)-1), sizeX-(taps(idx_scale)-1))
for idx_airplane=1:n_airplanes
    S1airplanes{idx_airplane} = dotproductL1_giulia(airplanes{idx_airplane}, gabors, n_splits);
end

% Hist range
minima_cars = zeros(n_cars, n_scales, n_templates, n_ori);
maxima_cars = zeros(n_cars, n_scales, n_templates, n_ori);
minima_airplanes = zeros(n_cars, n_scales, n_templates, n_ori);
maxima_airplanes = zeros(n_cars, n_scales, n_templates, n_ori);
for idx_car=1:n_cars
    for idx_scale=1:n_scales
        for idx_template=1:n_templates
            for idx_ori=1:n_ori
                minima_cars(idx_car, idx_scale, idx_template, idx_ori) = min(min(S1cars{idx_car}{idx_scale}(idx_template, idx_ori, :, :)));
                maxima_cars(idx_car, idx_scale, idx_template, idx_ori) = max(max(S1cars{idx_car}{idx_scale}(idx_template, idx_ori, :, :)));
            end
        end
    end
end
for idx_image=1:n_cars
    for idx_scale=1:n_scales
        for idx_template=1:n_templates
            for idx_ori=1:n_ori
                minima_airplanes(idx_image, idx_scale, idx_template, idx_ori) = min(min(S1airplanes{idx_image}{idx_scale}(idx_template, idx_ori, :, :)));
                maxima_airplanes(idx_image, idx_scale, idx_template, idx_ori) = max(max(S1airplanes{idx_image}{idx_scale}(idx_template, idx_ori, :, :)));
            end
        end
    end
end
minima_cars = minima_cars(:);
maxima_cars = maxima_cars(:);
minima_airplanes = minima_airpanes(:);
maxima_airplanes = maxima_airplanes(:);

range = [min([minima_cars minima_airplanes]) max([maxima_cars maxima_airplanes])]; % range of the histogram

% C1 responses
% output of poolingL1_giulia is of format
% zeros(n_binsL1, n_reg, n_templatesL1, 2)
for idx_car=1:n_cars
    histograms = poolingL1_giulia(S1cars{idx_car}, n_splits, n_bins, range,  'histogram');
    signature = histograms(:, :, :, 1);
    signature = signature(:);    
    C1cars{idx_car} = signature;
end

% Output
allcars_signatures = zeros(n_cars, length(C1cars{idx_car}));
for idx_car=1:n_cars
    allcars_signatures(idx_car, :) = C1cars{idx_car};
end

% C1 responses
% output of poolingL1_giulia is of format
% zeros(n_binsL1, n_reg, n_templatesL1, 2)
for idx_airplane=1:n_airplanes
    histograms = poolingL1_giulia(S1airplanes{idx_airplane}, n_splits, n_bins, range,  'histogram');
    signature = histograms(:, :, :, 1);
    signature = signature(:);    
    C1airplanes{idx_airplane} = signature;
end

% Output
allairplanes_signatures = zeros(n_airplanes, length(C1airplanes{idx_airplane}));
for idx_airplane=1:n_airplanes
    allairplanes_signatures(idx_airplane, :) = C1airplanes{idx_airplane};
end

%% Quantitative comparison

figure

subplot(1,2,1)
mean_cars = mean(allcars_signatures);
sd_cars = std(allcars_signatures);
fcars = [ mean_cars+3*sd_cars , flipdim(mean_cars-3*sd_cars,2)]; 
fill([1:size(allcars_signatures,2) , size(allcars_signatures,2):-1:1] , fcars, [7 7 7]/8)
hold on;
plot(1:size(allcars_signatures,2) , mean_cars , 'b' , 'LineWidth',1);
grid on

subplot(1,2,2)
mean_airplanes = mean(allairplanes_signatures);
sd_airplanes = std(allairplanes_signatures);
fairplanes = [ mean_airplanes+3*sd_airplanes , flipdim(mean_airplanes-3*sd_airplanes,2)]; 
fill([1:size(allairplanes_signatures,2) , size(allairplanes_signatures,2):-1:1] , fairplanes, [7 7 7]/8)
stem(1:size(allairplanes_signatures,2) , mean_airplanes , 'r' , 'LineWidth',1);
grid on

max(sd_cars)
max(sd_airplanes)

figure
imagesc(squeeze(allcars_signatures))