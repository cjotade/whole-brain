function [miJointValue] = MI_Joint(data, jointVariable1Columns, jointVariable2Columns, implementingClass)
    % jointVariables
    jointVariable1 = data(:, jointVariable1Columns);
    jointVariable2 = data(:, jointVariable2Columns);

    % MI Estimator
    miCalc = javaObject(implementingClass);

    miCalc.initialise(length(jointVariable1Columns), length(jointVariable2Columns));
    miCalc.setObservations(octaveToJavaDoubleMatrix(jointVariable1), octaveToJavaDoubleMatrix(jointVariable2));
    miJointValue = (miCalc.computeAverageLocalOfObservations())/log(2);
end