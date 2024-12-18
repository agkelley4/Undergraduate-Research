file = 'deformed_cell.png';
[x_image,z_image,phi_image] = get_image_data(file);
% inital conditions
% u0 = beta,gamma,delta

close all

u0 = [0,0];
%u0 = [.2,.2];
%u0 = [.5,.5];
%u0 = [.3,.3647];
% u0 = [.183,.183]
% u0 = [0.0003, 1.7707];

[phi_image_sort, indices] = sort(phi_image);
x_image_sort = x_image(indices);
z_image_sort = z_image(indices);
r_image_sort = sqrt(x_image_sort.^2 + z_image_sort.^2);

% initialize count and total pts
count = 0;
N = 100000;

% Generate N points
points = generate_points(N);

%Find area_data
for k=1:N
     x_mc = points(k,1);
     z_mc = points(k,2);
     phi_mc = atan2(x_mc,z_mc);
     r_mc = sqrt((x_mc)^2 + (z_mc)^2);
     r_im_bdry = interp1(phi_image_sort,r_image_sort,phi_mc); 
    if 0 <= r_mc && r_mc <= r_im_bdry
        count = count + 1;
    end

end
area_data = 2^2 * count / N;

% optimization starts here, change monte carlo to make the rest 1 function

centroid_data = [mean(x_image_sort), mean(z_image_sort)];

data_vars = struct('phi', phi_image_sort, 'r', r_image_sort, 'points', points, 'x', x_image_sort, ...
                'z', z_image_sort, 'area_data', area_data);


% TODO
options = optimset('Display', 'final', 'PlotFcn', []);
%[u, error, exitflag, output] = fminsearch(@(u0) data_sph_error(abs(u0),data_vars), abs(u0),options)
[u, error, exitflag, output] = fmincon(@(u0) data_sph_error(u0,data_vars), u0', -1 .* eye(2), zeros(2,1), [], [], [], [],[],options);% @nonlinear_constraints, options);
% [u, error, exitflag, output] = fmincon(@(u0) data_sph_error(u0, data_vars), ...
%                                        u0', [], [], [], [], [], []), ...
%                                        @nonlinear_constraints, options);


alpha = fsolve(@(alpha) equation_to_solve(alpha, u(1), u(2)), 1, options)

u 
error = data_sph_plot(u, data_vars)

% error = data_sph_plot(u0, data_vars)

% u = [0.2079, 0.2416], error = 33.8517 with u0 = [.183,.183]
% u = [0.2177, 0.2400], error = 33.7597 with u0 = [0,0]


% Optimization function
function error = data_sph_error(u0, data_vars)
    % define vars
    beta = u0(1); 
    gamma = u0(2);

    % get scale
    [scale, x_sph_harm_centroid, z_sph_harm_centroid] = montecarlo_plane(data_vars.points, u0, data_vars.area_data);

    % get alpha
    options = optimoptions('fsolve', 'Display', 'none', 'PlotFcn', []);

    alpha = fsolve(@(alpha) equation_to_solve(alpha, beta, gamma), 1, options);

    % calculate spherical harmonic values
    theta = 0;
    r_sph_harm = alpha / (2 * sqrt(pi)) + (beta / 2) * sqrt(3 / pi) * cos(data_vars.phi) + ...
          (gamma / 2) * sqrt(3 / pi) * sin(data_vars.phi) * cos(theta); %+ ...
          %(delta / 2) * sqrt(3 / pi) * sin(data_vars.phi) * sin(theta);
    x_sph_harm = r_sph_harm .* sin(data_vars.phi);
    z_sph_harm = r_sph_harm .* cos(data_vars.phi);
    
    
    % calculate final sph harm values
    x_sph_harm_scaled = (x_sph_harm - x_sph_harm_centroid) .* scale;
    z_sph_harm_scaled = (z_sph_harm - z_sph_harm_centroid) .* scale;

    theta_sph_harm = atan2(z_sph_harm_scaled, x_sph_harm_scaled);
    [B,I] = sort(theta_sph_harm);
    x_sph_harm_rot = x_sph_harm_scaled(I);
    z_sph_harm_rot = z_sph_harm_scaled(I);

    theta_data = atan2(data_vars.z, data_vars.x);
    [B_data, I_data] = sort(theta_data);
    x_data = data_vars.x(I_data);
    z_data = data_vars.z(I_data);

    x_sph_harm_new = zeros(size(x_sph_harm_scaled));
    z_sph_harm_new = zeros(size(z_sph_harm_scaled));
    for k = 1:length(theta_data)
        x_sph_harm_new(k) = interp1(B, x_sph_harm_rot, B_data(k));
        z_sph_harm_new(k) = interp1(B, z_sph_harm_rot, B_data(k));
    end

    % calculate error from sph_harm to image
    error_list = sqrt((x_sph_harm_new - x_data).^2 + (z_sph_harm_new - z_data).^2);
    flag = isnan(error_list);
    error_list(flag) = 1;
    error = sum(error_list);
    
end

% Function to generate random points within the cube
function points = generate_points(num_points)
    % inside [-1,1]x[-1,1]
    points = -1 + 2 * rand(num_points, 2);
end


function F = equation_to_solve(alpha, beta, gamma)
        F = alpha^3 + 3 * alpha * (beta^2 + gamma^2) - 1;
end


% Final function to plot and check optimization
function error = data_sph_plot(u0, data_vars)
    % define vars
    beta = u0(1);
    gamma = u0(2);

 
    % get scale
    [scale, x_sph_harm_centroid, z_sph_harm_centroid] = montecarlo_plane(data_vars.points, u0, data_vars.area_data);
    

    % get alpha
    options = optimoptions('fsolve', 'Display', 'none', 'PlotFcn', []);

    alpha = fsolve(@(alpha) equation_to_solve(alpha, beta, gamma), 1, options);

    % calculate spherical harmonic values
    theta = 0;
    r_sph_harm = alpha / (2 * sqrt(pi)) + (beta / 2) * sqrt(3 / pi) * cos(data_vars.phi) + ...
          (gamma / 2) * sqrt(3 / pi) * sin(data_vars.phi) * cos(theta); %+ ...
          %(delta / 2) * sqrt(3 / pi) * sin(data_vars.phi) * sin(theta);


    x_sph_harm = r_sph_harm .* sin(data_vars.phi);
    z_sph_harm = r_sph_harm .* cos(data_vars.phi);

    
    x_sph_harm_scaled = (x_sph_harm - x_sph_harm_centroid) .* scale;
    z_sph_harm_scaled = (z_sph_harm - z_sph_harm_centroid) .* scale;
    % [M,I] = min(z_sph_harm_scaled);
    % x_sph_harm_rot = circshift(x_sph_harm_scaled, -I);
    % z_sph_harm_rot = circshift(z_sph_harm_scaled, -I);
    
    theta_sph_harm = atan2(z_sph_harm_scaled, x_sph_harm_scaled);
    [B,I] = sort(theta_sph_harm);
    x_sph_harm_rot = x_sph_harm_scaled(I);
    z_sph_harm_rot = z_sph_harm_scaled(I);

    theta_data = atan2(data_vars.z, data_vars.x);
    [B_data, I_data] = sort(theta_data);
    x_data = data_vars.x(I_data);
    z_data = data_vars.z(I_data);

    x_sph_harm_new = zeros(size(x_sph_harm_scaled));
    z_sph_harm_new = zeros(size(z_sph_harm_scaled));
    for k = 1:length(theta_data)
        x_sph_harm_new(k) = interp1(B, x_sph_harm_rot, B_data(k));
        z_sph_harm_new(k) = interp1(B, z_sph_harm_rot, B_data(k));
    end


    figure;
    plot(x_sph_harm_scaled, z_sph_harm_scaled, 'LineWidth', 1)
    hold on;
    % plot(data_vars.x,data_vars.z,'.')
    plot(x_data, z_data,'.')
    axis equal;

    figure
    num = 1:length(x_sph_harm_scaled);
    subplot(1,2,1)
    plot(num, x_data, '.')
    hold on
    % plot(num, x_sph_harm_scaled, '.')
    plot(num,x_sph_harm_new, '.')
    subplot(1,2,2)
    plot(num, z_data, '.')
    hold on
    % plot(num, z_sph_harm_scaled, '.')
    plot(num,z_sph_harm_new, '.')


    error_list = sqrt((x_sph_harm_new - x_data).^2 + (z_sph_harm_new - z_data).^2);
    figure
    plot(num, error_list)
    nans_list = isnan(error_list);
    hold on
    plot(num, nans_list, 'o')
    sum(error_list)
    sum(nans_list)

    figure
    plot(num, B)
    hold on
    plot(num, B_data)
    title('phi vs index')
    

    % calculate error from sph_harm to image
    flag = isnan(error_list);
    error_list(flag) = 1;
    error = sum(error_list);
end

%non linear constraint
function [c, ceq] = nonlinear_constraints(u0)
  
    % c = u0(1)^2 + u0(2)^2 - (sqrt(2)*0.4)^2;
    c = u0(1)^2 + u0(2)^2 - (0.2)^2;
    
    ceq = [];
end


