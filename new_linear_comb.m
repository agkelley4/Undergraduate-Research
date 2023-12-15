syms a b c

% expression from mathematica
expression = (a^3 + 3*a*(b^2 + c^2)) / (2*sqrt(pi)) - 1/(2*sqrt(pi));

% plot
fig = figure;
fimplicit3(expression, [0, 3, 0, 3, 0, 3], 'MeshDensity', 50);

% labeling
xlabel('{\alpha}');
ylabel('{\beta}');
zlabel('{\gamma}');
title('Plot of {\alpha}Y_0^0 + {\beta}Y_1^0 + {\gamma}Y_1^1');

% set axis limits
xlim([0, 3]);
ylim([0, 3]);
zlim([0, 3]);

grid on;
view(3);


