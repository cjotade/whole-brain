addpath('../../jidt/demos/octave')  
javaaddpath('../../jidt/infodynamics.jar');

% Read Data
data_folder = '../../data/muestras/';
data_dir = dir(data_folder);
data_dir = data_dir(~ismember({data_dir.name}, {'.', '..', 'index.txt'}));

% Results Folder
results_folder = '../../results/muestras/';
mkdir(results_folder)

% How many data samples get from data
num_data = [100, 1000, 10000, 100000];

for i = 8:numel(data_dir)
    gaussian_i_folder = fullfile(data_folder, string(data_dir(i).name), filesep);
    disp("=========== Calculating metrics on " + gaussian_i_folder + " ============")
    gaussian_i_dir = dir(fullfile(gaussian_i_folder, '*.txt'));
    for n_data = num_data
        disp("-----> Working with " + n_data + " samples")
        num_samples = numel(gaussian_i_dir);
        data_result = zeros(num_samples, 3);
        for j = 1:num_samples
            filepath = fullfile(gaussian_i_dir(j).folder, filesep, gaussian_i_dir(j).name);
            data = load(filepath);
            data = data(1:n_data,:);
            TC = totalCorrelation(data);
            DTC = dualTotalCorrelation(data);
            O = TC - DTC;
            data_result(j, 1) = TC;
            data_result(j, 2) = DTC;
            data_result(j, 3) = O;
        end
        result_filepath = fullfile(results_folder, string(data_dir(i).name), string(n_data), filesep);
        mkdir(result_filepath)
        result_filename = fullfile(result_filepath, "results_samples.txt");
        dlmwrite(result_filename, data_result)
    end
end
    

