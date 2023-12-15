% constants
Y0 = 1 / (2 * sqrt(pi));

% range
a = linspace(0, 2, 100);
b = linspace(0, 3, 100);

[A, B] = meshgrid(a, b);

% find f
F = (3 * A .* B.^2 + A.^3)/(2*sqrt(pi));

% Create a pseudocolor plot
figure;
pcolor(A, B, F);
shading flat;
title('Plot of {\alpha}Y_{0}^{0} + {\beta}Y_{1}^{0}');
xlabel('{\alpha}');
ylabel('{\beta}');
colorbar;

% Add 20 ontour lines
hold on;
contour(A, B, F, 20, 'k'); 

% find contour levels
C = contourc(a, b, F, [Y0 Y0]);


contourData = getContourMatrixData(C);

plot(contourData(:, 1), contourData(:, 2), 'ro', 'MarkerSize', 10);

hold off;

function contourData = getContourMatrixData(C)
    contourData = [];
    k = 1;
    while k < size(C, 2)
        numPoints = C(2, k);
        contourData = [contourData; C(1, k + 1:k + numPoints)', C(2, k + 1:k + numPoints)'];
        k = k + numPoints + 1;
    end
end
