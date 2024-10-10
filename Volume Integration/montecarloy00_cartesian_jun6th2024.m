clear; clc;

% Values of num_points to iterate over
num_points_values = 10000000;

% Initialize lists to store values
approx_volume_values = [];
variance_values = [];
num_points_list = [];

% Variables to store the last scatter plot data
last_points_inside_shape = [];
last_points_outside_shape = [];
last_num_points = 0;
last_approx_volume = 0;
last_variance = 0;
%volume of overall cube
cube_volume = 2^3;

% Iterate over each num_points value
for num_points = num_points_values

    % Generate num_points points within the cube [-1,1] x [-1,1] x [-1,1]
   points = generate_points(num_points);
    
   
    %test point
   % points = [0,0.3,0];
 %   points = [0.01,0.01,-.05];
  %  points = [0,0,-sqrt(3/pi)/4];
  %  points = [-sqrt(3/pi)/4, 0, sqrt(3/pi)/4];

    % Separate points inside and outside the shape
    inside_shape = is_inside_shape(points(:, 1), points(:, 2), points(:, 3));
    points_inside_shape = points(inside_shape, :);
    points_outside_shape = points(~inside_shape, :);
    
    % Count the number of points inside the shape
    num_points_inside_shape = size(points_inside_shape, 1);
    
    % Calculate the ratio of points inside the shape to the total number of points
    approx_volume = cube_volume * num_points_inside_shape / num_points;
    approx_volume_values = [approx_volume_values; approx_volume];
    
    % Calculate the variance
    variance = cube_volume * sqrt(approx_volume * (2 - approx_volume) / num_points);
    variance_values = [variance_values; variance];
    
    num_points_list = [num_points_list; num_points];
    
    % Store data for the last scatter plot
    if num_points == num_points_values(end)
        last_points_inside_shape = points_inside_shape;
        last_points_outside_shape = points_outside_shape;
        last_num_points = num_points;
        last_approx_volume = approx_volume;
        last_variance = variance;
    end
end

% Plot the last scatter plot
figure;
%scatter3(last_points_outside_shape(:, 1), last_points_outside_shape(:, 2), last_points_outside_shape(:, 3), 5, 'blue', 'filled'); hold on;
scatter3(last_points_inside_shape(:, 1), last_points_inside_shape(:, 2), last_points_inside_shape(:, 3), 5, 'red', 'filled');

% Plot the shape boundary
[x_sphere, y_sphere, z_sphere] = sphere(100);
[x_sphere, y_sphere, z_sphere] = sphere2cart(ones(size(x_sphere)), atan2(y_sphere, x_sphere), atan2(sqrt(x_sphere.^2 + y_sphere.^2), z_sphere));
%plot3(x_sphere(:), y_sphere(:), z_sphere(:), 'r', 'LineWidth', 1.5);

title(sprintf('%d Random Points and Shape in [-1, 1] x [-1, 1] x [-1, 1]\n%d Points Inside Shape\n%.5f Approximation of Volume\n%.5f Variance', ...
    last_num_points, size(last_points_inside_shape, 1), last_approx_volume, last_variance));
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
grid on;
axis equal;
hold off;

% Function to generate random points within the cube
function points = generate_points(num_points)
    points = -1 + 2 * rand(num_points, 3);
end

% Function to check if points are inside the shape
function inside = is_inside_shape(x, y, z)
    % alpha = 1/(2^(1/3));
    % beta =  1/(sqrt(3)*(2^(1/3))); %1.1547;

    alpha = 0.79;
    beta = 0.46;
    %gamma = 1;
    
    num_points = length(x);
    inside = false(num_points, 1);
    
    for i = 1:num_points
        [~, theta0, phi0] = xyz(x(i), y(i), z(i));
     %   display("theta0 phi0 = " + theta0 + " " + phi0)
        rho = alpha .* Y00(theta0,phi0) + beta .* Y10(theta0, phi0);
     %    rho = gamma .* Y1neg1(theta0,phi0);
       % rho = beta .* Y10(theta0, phi0);
     %   rho = Y22(theta0,phi0);
     
        [x0, y0, z0] = sphere2cart(rho, theta0, phi0);
      %  display("x0 y0 z0 = " + x0 + " " + y0 + " " + z0)
        t_values = [x(i)/x0, y(i)/y0, z(i)/z0];
      
        [M,I] = max(t_values);
        t = t_values(I);
        
        if t >= 0 && t <= 1 % t >= -1 if Yx2 called. create new code for those
            inside(i) = true;
        end
    end
end

% Function to convert Cartesian coordinates to spherical
function [rho0, theta0, phi0] = xyz(x, y, z)
    rho0 = sqrt(x.^2 + y.^2 + z.^2);
    theta0 = atan2(y,x);
    phi0 = acos(z./rho0);
    
end

% Function to convert spherical coordinates to Cartesian
function [x, y, z] = sphere2cart(rho, theta0, phi0)
    x = rho .* cos(theta0) .* sin(phi0);
    y = rho .* sin(theta0) .* sin(phi0);
    z = rho .* cos(phi0);
end

% Functions for Ylm
function rho = Y00(~, ~) % ~ is theta0 and phi0
    rho = 1 / (2 * sqrt(pi));
end

function rho = Y10(~, phi0) % ~ is theta0
    rho = sqrt(3 / pi) * cos(phi0) / 2;
end

function rho = Y11(theta0,phi0)
   rho = sqrt(3/pi) * sin(phi0) * cos(theta0) / 2;
end

function rho = Y1neg1(theta0,phi0)
    rho = sqrt(3/pi) * sin(phi0) * sin(theta0) / 2;
end

function rho = Y2neg2(theta0,phi0)
    rho = sqrt(15/pi) * sin(phi0)^2 * cos(theta0) * sin(theta0) / 2;
end

function rho = Y2neg1(theta0,phi0)
    rho = sqrt(15/pi) * sin(phi0) * cos(phi0) * sin(theta0) / 2;
end

function rho = Y20(theta0,phi0)
    rho = 3 * sqrt(5/pi) * cos(phi0)^2 / 4;
end

function rho = Y21(theta0,phi0)
    rho = sqrt(15/pi) * sin(phi0) * cos(phi0) * cos(theta0) / 2;
end

function rho = Y22(theta0,phi0)
    rho = sqrt(15/pi) * sin(phi0)^2 * cos(2*theta0) / 4;
end