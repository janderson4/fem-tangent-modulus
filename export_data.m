% Export FEM data for CUDA testing
% Please run script.mlx before running this script

% Display all variables in workspace
disp('Variables in workspace:');
whos

% Check if required variables exist
required_vars = {'nel', 'n_int', 'E', 'nu', 'H', 'beta', 'sigma_y0', 'SIG', 'ALPHA', 'R'};
missing_vars = {};
for i = 1:length(required_vars)
    if ~exist(required_vars{i}, 'var')
        missing_vars{end+1} = required_vars{i};
    end
end

if ~isempty(missing_vars)
    error('The following variables are not defined, please run script.mlx first:\n%s', strjoin(missing_vars, '\n'));
end

% Check variable dimensions
disp('Variable dimensions:');
disp(['SIG size: ', num2str(size(SIG))]);
disp(['ALPHA size: ', num2str(size(ALPHA))]);
disp(['R size: ', num2str(size(R))]);

% Open file for writing
fid = fopen('matlab_data.txt', 'w');
if fid == -1
    error('Cannot create matlab_data.txt file');
end

% Write basic information
fprintf(fid, 'Number of elements: %d\n', nel);
fprintf(fid, 'Number of integration points: %d\n', n_int);

% Save material properties
fprintf(fid, 'Material Properties:\n');
fprintf(fid, 'E: %f\n', E);
fprintf(fid, 'nu: %f\n', nu);
fprintf(fid, 'H: %f\n', H);
fprintf(fid, 'beta: %f\n', beta);
fprintf(fid, 'sigma_y0: %f\n', sigma_y0);
fprintf(fid, '\n');

% Save data for each integration point
disp('Writing to matlab_data.txt...');
for elem = 1:nel
    if mod(elem, 50) == 0
        fprintf('Processing element %d/%d\n', elem, nel);
    end
    
    for ip = 1:n_int
        fprintf(fid, 'Element %d, IP %d:\n', elem, ip);
        
        % Write stress
        fprintf(fid, 'Stress: ');
        for i = 1:6
            fprintf(fid, '%f ', SIG(i,ip));
        end
        fprintf(fid, '\n');
        
        % Write back stress
        fprintf(fid, 'Alpha: ');
        for i = 1:6
            fprintf(fid, '%f ', ALPHA(i,ip));
        end
        fprintf(fid, '\n');
        
        % Write yield radius
        fprintf(fid, 'R: %f\n', R(ip));
        fprintf(fid, '\n');
    end
end

fclose(fid);
disp('matlab_data.txt writing completed');

% Create CUDA input file
fid = fopen('cuda_input.txt', 'w');
if fid == -1
    error('Cannot create cuda_input.txt file');
end

% Write input data in CUDA format
disp('Writing to cuda_input.txt...');
for elem = 1:nel
    if mod(elem, 50) == 0
        fprintf('Processing element %d/%d\n', elem, nel);
    end
    
    for ip = 1:n_int
        % Write stress
        fprintf(fid, 'Stress %d: ', (elem-1)*n_int + ip);
        for i = 1:6
            fprintf(fid, '%f ', SIG(i,ip));
        end
        fprintf(fid, '\n');
        
        % Write back stress
        fprintf(fid, 'Alpha %d: ', (elem-1)*n_int + ip);
        for i = 1:6
            fprintf(fid, '%f ', ALPHA(i,ip));
        end
        fprintf(fid, '\n');
        
        % Write yield radius
        fprintf(fid, 'R %d: %f\n', (elem-1)*n_int + ip, R(ip));
        fprintf(fid, '\n');
    end
end

fclose(fid);
disp('cuda_input.txt writing completed');

% Check if files were created successfully
if exist('matlab_data.txt', 'file')
    info = dir('matlab_data.txt');
    fprintf('matlab_data.txt size: %.2f MB\n', info.bytes/1024/1024);
else
    warning('matlab_data.txt was not created');
end

if exist('cuda_input.txt', 'file')
    info = dir('cuda_input.txt');
    fprintf('cuda_input.txt size: %.2f MB\n', info.bytes/1024/1024);
else
    warning('cuda_input.txt was not created');
end

disp('Data export completed');

% Test variables
test_vars 