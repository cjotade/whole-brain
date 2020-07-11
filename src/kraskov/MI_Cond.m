function [miCondValue] = MI_Cond(data, jointVariable1Columns, jointVariable2Columns, condVariable3Columns, implementingClass)
    % jointVariables
    jointVariable1 = data(:, jointVariable1Columns);
    jointVariable2 = data(:, jointVariable2Columns);
    condVariable3 = data(:, condVariable3Columns);

    % MI Conditional Estimator
    miCalc = javaObject(implementingClass);

    miCalc.initialise(length(jointVariable1Columns),length(jointVariable2Columns),length(condVariable3Columns));
    miCalc.setObservations(octaveToJavaDoubleMatrix(jointVariable1), ...
        octaveToJavaDoubleMatrix(jointVariable2), ...
        octaveToJavaDoubleMatrix(condVariable3));
    miCondValue = (miCalc.computeAverageLocalOfObservations())/log(2);
end