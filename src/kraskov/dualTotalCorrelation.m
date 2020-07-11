function [DTC] = dualTotalCorrelation(data)
    %% Data
    [data_len, n] = size(data);
    %% First Term
    % Estimator
    implementingClass = 'infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1';

    % jointVariables Columns
    jointVariable1Columns = n; % array indices start from 1 in octave/matlab
    jointVariable2Columns = 1:n-1;

    firstTerm = MI_Joint(data, jointVariable1Columns, jointVariable2Columns, implementingClass);

    %% Second Term
    % Estimator
    implementingClass = 'infodynamics.measures.continuous.kraskov.ConditionalMutualInfoCalculatorMultiVariateKraskov1';

    DTC_list = zeros(1, n-1);
    for j=2:n-1
        % jointVariables Columns
        jointVariable1Columns = j; % array indices start from 1 in octave/matlab
        jointVariable2Columns = 1:j-1;
        condVariable3Columns = j+1:n;

        miCondValue = MI_Cond(data, jointVariable1Columns, jointVariable2Columns, condVariable3Columns, implementingClass);

        DTC_list(j-1) = miCondValue;

    end

    %% DTC
    DTC = firstTerm + sum(DTC_list);
end