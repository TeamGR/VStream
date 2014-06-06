function [ out ] = pooling( filteredImg , n_templates , poolingSplitNum , numBars , policy)
%POOLING Performs pooling on a given region according to the specified
%policy
%   Author: Raffaello Camoriano
%
%   filteredImg : Cell array of filtered images
%   numBars : Number of bars of the histogram
%   poolingSplitNum : Number of splits along axes for determining pooling regions
%   n_templates : Number of templates
%   policy : Pooling policy
%               'histogram': provides the concatenated histogram signature
%               'moments': provides the first 3 moments of the distribution
%               (mean, variance and max value)

poolRegsNum = poolingSplitNum^2;

[inSizeX inSizeY] = size(filteredImg(1,1,1));

% Side sizes of the pooling regions
regSideLenX = floor(inSizeX/poolingSplitNum);
regSideLenY = floor(inSizeY/poolingSplitNum);
regArea = regSideLenX * regSideLenY;

ma = max(max(cell2mat(filteredImg(1,1,1))));
mi = min(min(cell2mat(filteredImg(1,1,1))));

switch policy
    case 'histogram'
        
        step = (ma-mi)/numBars;
        x_hist = mi:step:ma;
        out = cell(n_templates, poolRegsNum, 2);

        % Templates loop
        for j = 1:n_templates

            % Reshape filtered image in column array form
            reshFiltI = reshape(cell2mat(filteredImg(1,j,:)),numel(cell2mat(filteredImg(1,j,:))),1);

            % Pooling regions loop
            for k = 1:poolRegsNum
                [out{j,k,1} out{j,k,2}] = hist(reshFiltI((k-1)*regArea + 1 : k*regArea), x_hist);
            end
        end
        
        
    case 'moments'
        % Mean
        out{1} = mean2(filteredImg);
        % Variance
        out{2} = mean2(filteredImg.^2) - m1^2;
        % Maximum
        out{3} = max(max(filteredImg));
end


end

