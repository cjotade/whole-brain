addpath('../../jidt/demos/octave')  
javaaddpath('../../jidt/infodynamics.jar');

% Read Data
data_folder = '../../data/firing_rates/';
data_dir = dir(data_folder);
data_dir = data_dir(~ismember({data_dir.name}, {'.', '..', 'index.json'}));

% Results Folder
results_folder = '../../results/firing_rates/';
mkdir(results_folder)

% How many data samples get from data
num_data = [100, 1000, 10000, 100000];

% Number of variables or regions
num_variables = 10;

for i = 1:numel(data_dir)
    firing_i_folder = fullfile(data_folder, string(data_dir(i).name), filesep);
    disp("=========== Calculating metrics on " + firing_i_folder + " ============")
    firing_i_dir = dir(fullfile(firing_i_folder, '*.txt'));
    for n_data = num_data
        disp("-----> Working with " + n_data + " samples")
        num_samples = numel(firing_i_dir);
        data_result = zeros(num_samples, 3);
        TC_result = zeros(num_samples, num_variables-1);
        DTC_result = zeros(num_samples, num_variables-1);
        for j = 1:num_samples
            filepath = fullfile(firing_i_dir(j).folder, filesep, firing_i_dir(j).name);
            data = load(filepath);
            data = data(1:n_data,:);
            [TC, TC_list] = totalCorrelation(data);
            [DTC, DTC_list] = dualTotalCorrelation(data);
            O = TC - DTC;
            data_result(j, 1) = TC;
            data_result(j, 2) = DTC;
            data_result(j, 3) = O;
            TC_result(j, :) = TC_list;
            DTC_result(j, :) = DTC_list;
        end
        result_filepath = fullfile(results_folder, string(data_dir(i).name), string(n_data), filesep);
        mkdir(result_filepath)
        result_filename = fullfile(result_filepath, "results_samples.txt");
        tc_filename = fullfile(result_filepath, "tc_samples.txt");
        dtc_filename = fullfile(result_filepath, "dtc_samples.txt");
        dlmwrite(result_filename, data_result)
        dlmwrite(tc_filename, TC_result)
        dlmwrite(dtc_filename, DTC_result)
    end
end