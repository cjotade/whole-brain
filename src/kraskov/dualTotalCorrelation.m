function [DTC, DTC_list] = dualTotalCorrelation(data)
    %% Data
    [data_len, n] = size(data);
    
    DTC_list = zeros(1, n-1);
    %% First Term
    % Estimator
    implementingClass = 'infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov2';

    % jointVariables Columns
    jointVariable1Columns = n; % array indices start from 1 in octave/matlab
    jointVariable2Columns = 1:n-1;

    firstTerm = MI_Joint(data, jointVariable1Columns, jointVariable2Columns, implementingClass);
    DTC_list(1) = firstTerm;
    %% Sum Terms
    % Estimator
    implementingClass = 'infodynamics.measures.continuous.kraskov.ConditionalMutualInfoCalculatorMultiVariateKraskov2';
    
    for j=2:n-1
        % jointVariables Columns
        jointVariable1Columns = j; % array indices start from 1 in octave/matlab
        jointVariable2Columns = 1:j-1;
        condVariable3Columns = j+1:n;

        miCondValue = MI_Cond(data, jointVariable1Columns, jointVariable2Columns, condVariable3Columns, implementingClass);

        DTC_list(j) = miCondValue;
    end

    %% DTC
    DTC = sum(DTC_list);
end