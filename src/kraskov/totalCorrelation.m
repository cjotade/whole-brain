function [TC, TC_list] = totalCorrelation(data)
    %% Data
    [data_len, n] = size(data);
    %% Terms
    implementingClass = 'infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov2';
    TC_list = zeros(1, n-1);
    for i=2:n
        % jointVariables Columns
        jointVariable1Columns = i; % array indices start from 1 in octave/matlab
        jointVariable2Columns = 1:i-1;

        miJointValue = MI_Joint(data, jointVariable1Columns, jointVariable2Columns, implementingClass);

        TC_list(i-1) = miJointValue;
    end
    %% TC
    TC = sum(TC_list);
end
