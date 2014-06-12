function visualizeExample(x)
% visualizeExample(x)
% where x is a vectorized square image. The function visualize the image.
%
% The following code visualizes the first example of the dataset MNIST_3_5
% Example:
% load('./MNIST_3_5.m');
% visualizeExample(X(1,:));
%
    colormap('gray')
    l = sqrt(numel(x));
    I = reshape(x, l, l)';
    J = I(l:-1:1,:);
    contourf(J); 
    % imagesc(J);
    axis equal; axis off;
end
