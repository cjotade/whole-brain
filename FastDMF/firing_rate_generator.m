%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example usage of DMF fMRI simulator, with parameters extracted from the Deco
% et al. Current Biology 2018 paper (see references in README).
%
% Pedro Mediano, Apr 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = [];
load data/sc_90.mat sc90
load data/mean5HT2A_bindingaal.mat mean5HT2A_aalsymm
load AAL/dmf_labels_ala_deco.mat dmf_labels

% DMF parameters
params.C         = sc90/max(sc90(:))*0.2; % structural connectivity
% Intercala regiones del cerebro (no simetrico)
params.receptors = mean5HT2A_aalsymm(:,1)/max(mean5HT2A_aalsymm(:,1)); %promedio receptores 5HT2A, aal simetrico del cerebro.
params.dt        = 0.1; % ms
params.tmax      = 250000; % time points
params.taon      = 100; % NMDA tau ms
params.taog      = 10; % GABA tau ms
params.gamma     = 0.641; % Kinetic Parameter of Excitation
params.sigma     = 0.01; % Noise SD nA
params.JN        = 0.15; % excitatory synaptic coupling nA
params.I0        = 0.382; % effective external input nA
params.Jexte     = 1.; % external->E coupling
params.Jexti     = 0.7;% external->I coupling
params.w         = 1.4; % local excitatory recurrence
params.g_e       = 0.16; % excitatory conductance
params.Ie        = 125.; % excitatory threshold for nonlineariy
params.ce        = 310.; % excitatory non linear shape parameter
params.g_i       = 0.087;
params.Ii        = 177.;
params.ci        = 615.;
params.wgaine    = 0.0; % neuromodulatory gain a 0.2 (excitatorio)
params.wgaini    = 0; % neuromodulatory gain (inhibitorio)
params.G         = 2; %  Global Coupling Parameter
params.alphas    = ones(size(params.C,1),1).*1.5; % parameter of the feedback inhibitory control
params.stren     = sum(params.C)'./2; % node strength

% Balloon-Windkessel parameters (from firing rates to BOLD signal)
params.decimate_val = 1;
params.subsamp      = 2; % number of seconds to sample bold signal
params.dtt          = 0.001; % seconds
params.tmin         = 20; % min time to take bold signal
params.n_min        = round(params.tmin/params.dtt); % min time in points
params.inic         = 1;

% Parallel computation parameters
params.batch_size = 5000;

% Run
% 2s resolucion de BOLD
nb_steps = 120000;

%% Generate data using params
data_folder = "../data/firing_rates/";
mkdir(data_folder)
regions_name = "4L_dmn"; % "dmn", "no_dmn", "visual"
save_folder = fullfile(data_folder, regions_name, filesep);
mkdir(save_folder)


%dmn_ids = [12 79 16 75 18 73 33 58 34 57];
%random_ids = [5 13 85 43 22 76 21 60 32 18];
%visual_ids = [19 72 20 71 25 66 26 65 27 64];
%4L_dmn_ids = [12 16 18 33];

regions = [12 16 18 33];

for i = 1:100
    tic;
    rates = DMF(params, nb_steps, 'rate');
    toc

    % Select regions
    data = rates(regions, 20000:end)';

    % Store in folder
    filename = fullfile(save_folder, "samples_" + string(i) + ".txt");
    dlmwrite(filename, data, " ")
end
