function [scale, x_centroid, z_centroid] = montecarlo_plane(points, u0, area_data)
    
    beta = abs(u0(1));
    gamma = abs(u0(2));
    count = 0; 
    x_count = 0;
    z_count = 0;

    options = optimoptions('fsolve', 'Display', 'none', 'PlotFcn', []);

    alpha = fsolve(@(alpha) equation_to_solve(alpha, beta, gamma), 1, options);

    num_points = size(points,1);
    
    % figure;
    % hold on;

        
    for i = 1:num_points
        x_mc = points(i, 1);
        z_mc = points(i, 2);
        [~, theta0, phi0] = xyz(x_mc, 0, z_mc); % send y = 0
     
        rho = alpha * Y00(theta0,phi0) + beta * Y10(theta0, phi0) + gamma * Y11(theta0, phi0);
     
    
        [x0, y0, z0] = sphere2cart(rho, theta0, phi0);
      
        t_values = [x_mc/x0, 0, z_mc/z0];
        [M, I] = max(t_values);
        t = t_values(I);
        
        if t >= 0 && t <= 1
            count = count + 1;
            x_count = x_count + x_mc;
            z_count = z_count + z_mc;
            % plot(x_mc, z_mc, 'bo');
            
        end
    
    end
    area_sph_harm = 2^2 * count / num_points;
    scale = sqrt(area_data/area_sph_harm);

    x_centroid = x_count / count;
    z_centroid = z_count / count;

    axis equal;
    % Function to convert Cartesian coordinates to spherical
    function [rho0, theta0, phi0] = xyz(x, y, z)
        rho0 = sqrt(x^2 + y^2 + z^2);
        theta0 = atan2(y,x);
        phi0 = acos(z/rho0);
        
    end
    
    % Function to convert spherical coordinates to Cartesian
    function [x, y, z] = sphere2cart(rho, theta0, phi0)
        x = rho * cos(theta0) * sin(phi0);
        y = rho * sin(theta0) * sin(phi0); % equals 0
        z = rho * cos(phi0);
    end
    
    % Functions for Ylm
    function rho = Y00(~, ~) % ~ is theta0 and phi0
        rho = 1 / (2 * sqrt(pi));
    end
    
    function rho = Y10(~, phi0) % ~ is theta0
        rho = sqrt(3 / pi) * cos(phi0) / 2;
    end
    
    function rho = Y11(theta0,phi0)
       rho = sqrt(3 / pi) * sin(phi0) * cos(theta0) / 2;
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

    % Function to solve for alpha
    function F = equation_to_solve(alpha, beta, gamma, delta)
        F = alpha^3 + 3 * alpha * (beta^2 + gamma^2) - 1;
    end

end

