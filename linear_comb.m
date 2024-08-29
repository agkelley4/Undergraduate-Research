% Constants
compare = 1 / (6 * sqrt(pi));

% Range with higher resolution
a = linspace(0, 3, 1000);
b = linspace(0, 3, 1000);

[A, B] = meshgrid(a, b);

% Functions
F = A .* (A.^2 + 3 * B.^2) / (6 * sqrt(pi));

% F2 = (A.^4 + 18*A.^2 .* B.^2 + 9*B.^4) ./ (24*B*sqrt(3*pi));
F2 = ((A + sqrt(3).*B).^4)./(48.*B.*sqrt(3*pi));
%F2 = B .* (2.*A.^2 + B.^2) .* sqrt(3/pi) / 8;


% Piecewise function
combinedF = F;
condition = A < sqrt(3) * B;
combinedF(condition) = F2(condition);

% Create figure
figure('Position', [100, 100, 600, 600]); % Changes resolution of figure
contour(A, B, combinedF, "ShowText", "on");
shading flat;

title('Plot of {\alpha}Y_{0}^{0} + {\beta}Y_{1}^{0}', 'FontSize', 21);
xlabel('{\alpha}', 'FontSize', 19);
ylabel('{\beta}', 'FontSize', 19, 'rotation', 0);

axis equal;
xlim([0 1.3]);
ylim([0, 1.3]);

% Scale image font size
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 21); % Adjusting font size

hold on;

% Find contour levels
C = contourc(a, b, combinedF, [compare compare]);
contourData = getContourMatrixData(C);
plot(contourData(:, 1), contourData(:, 2), 'r-', 'LineWidth', 3);

% Add the line alpha = sqrt(3) * beta
beta_line = linspace(0, 1.2, 1000);
alpha_line = sqrt(3) * beta_line;
plot(alpha_line, beta_line, 'b--', 'LineWidth', 3);

hold off;


% Extract points on the contour line for the given 'compare' value
contourPointsMatrix = getContourMatrixData(C);


% Write the matrix to a file
fileID = fopen('linear_comb_points.txt', 'w');
fprintf(fileID, 'alpha, beta\n');
fprintf(fileID, '%.6f, %.6f\n', contourPointsMatrix');
fclose(fileID);

function contourData = getContourMatrixData(C)
    contourData = [];
    k = 1;
    while k < size(C, 2)
        numPoints = C(2, k);
        contourData = [contourData; C(1, k + 1:k + numPoints)', C(2, k + 1:k + numPoints)'];
        k = k + numPoints + 1;
    end
end
