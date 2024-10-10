% Load model from COMSOL file 'flow_ellipsoid_axis.mph'
model = mphload('flow_ellipsoid_axis.mph');

% change mesh size
mesh = 6;
model.component('comp1').mesh('mesh1').autoMeshSize(mesh);

% conversions
rad2deg = 180/pi;
deg2rad = pi/180;

%initial conditions in form [thetap,phip,Ompy,Up]
% x0 = [45; 90; 1.15; 1.5];
% x0 = [45*deg2rad, 90*deg2rad, 1.15, 1.5];
% x0 = [44.9*deg2rad, 89.75*deg2rad, 1.08, 1.52];

%initial conditions in form [thetap,phip,Ompy,Up,Ompx, Ompz]
x0 = [45*deg2rad, 90*deg2rad, 1.15, 1.5, 0, 0];

% open a file to write results to a file
% filename = sprintf('x0_%.2f_%.2f_%.2f_%.2f_results.txt', x0);
date = string(datetime("today"));
folder = 'C:\Users\hoodflowlab\Documents\comsol examples\';
filename = strcat(folder, sprintf('x0_%.2f_%.2f_%.2f_%.2f', x0),'_mesh',...
    num2str(mesh),'_date',date,'_results.txt');
fid = fopen(filename, 'a'); % Open file for appending

tic

% fsolve
diffminchange = 1e-5;
function_tol = 1e-4;
options = optimoptions('fsolve','DiffMinChange',diffminchange,...%'StepTolerance',1e-8,...
    'FunctionTolerance',function_tol,'Display','iter');
[x0_solution,fval,exitflag,output] = fsolve(@solveModel, x0, options, model, fid);

% calculuate time to run and print to file
time = toc;
fprintf(fid, 'time: %f \n', time);
fprintf(fid, 'exit flag: %f \n', exitflag);

% close file
fclose(fid);

% Display the solution
disp('Optimal x0:');
disp(x0_solution);
disp('Fz and T:');
disp(fval);
disp('exit flag:');
disp(exitflag)

% Define a function to set parameters, run the solver, and extract results
function res = solveModel(x0, model, fid)
    % conversions
    rad2deg = 180/pi;
    deg2rad = pi/180;

    % x0
    % Set parameters
    model.param.set('thetap', num2str(x0(1)*rad2deg));
    model.param.set('phip', num2str(x0(2)*rad2deg));
    model.param.set('Ompy', num2str(x0(3)));
    model.param.set('Up', num2str(x0(4)));
    model.param.set('Ompx', num2str(x0(5)));
    model.param.set('Ompz', num2str(x0(6)));

    % Run solver
    model.sol('sol1').runAll;

    % Evaluate results
    model.result.evaluationGroup('eg1').run;
    temp = model.result.evaluationGroup('eg1').getReal();

    % Extract force and torque
    Fz = temp(3);  % Force
    T = temp(4:6);  % Torque
    % x0

    % write results to a file
    % fprintf(fid, 'x0: %f %f %f %f\n', x0);
    fprintf(fid, 'x0: %f %f %f %f %f %f\n', x0);
    fprintf(fid, 'Results: Fz = %f, Tx = %f, Ty = %f, Tz = %f\n', Fz, T);
    

    res = [Fz, T];
    %res = [Fz, T(2)];
end
