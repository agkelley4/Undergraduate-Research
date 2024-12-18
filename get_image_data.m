
function [x,z,phi] = get_image_data(file)
    im = imread(file);
    x = im2double(im(:,:,1)); % greyscale image
    
    % Set threshold
    flag2 = x < 0.9;
    x(flag2) = 0;
    
    % Edge detection using Sobel operators
    K_sobelx = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
    x_sobelx = conv2(x, K_sobelx); 
    K_sobely = [1, 2, 1; 0, 0, 0; -1, -2, -1];
    x_sobely = conv2(x, K_sobely);
    
    edges = sqrt(x_sobelx.^2 + x_sobely.^2);
    
    % Brighten edges above threshold
    flag = edges > 0.95;
    edges_bright = zeros(size(edges));
    edges_bright(flag) = edges(flag);
    
    % Set data size and scaling parameters
    data_size = [346, 394]; %size(edges_bright)
    px_2_dim = 1 / data_size(1);
    
    
    % Find edge points and calculate their positions
    
    edges_bright_inside = edges_bright(3:(end - 2), 3:(end - 2));
    [row, col] = find(edges_bright_inside >= 0.9);
    centroid_px = [mean(row), mean(col)];
    
    x = (row - centroid_px(1)) * px_2_dim;
    z = (col - centroid_px(2)) * px_2_dim;
    
    phi = atan2(x, z);
end

