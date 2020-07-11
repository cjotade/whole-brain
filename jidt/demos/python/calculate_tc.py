from jpype import *
import numpy
import sys
# Our python data file readers are a bit of a hack, python users will do better on this:
sys.path.append("../jidt/demos/python")
import readFloatsFile

# Add JIDT jar library to the path
jarLocation = "../jidt/infodynamics.jar"
# Start the JVM (add the "-Xmx" option with say 1024M if you get crashes due to not enough memory space)
startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=" + jarLocation)

# 0. Load/prepare the data:
dataRaw = readFloatsFile.readFloatsFile("../data/muestras/samples2.txt")
# As numpy array:
data = numpy.array(dataRaw)
# 1. Construct the calculator:
calcClass = JPackage("infodynamics.measures.continuous.kraskov").MutualInfoCalculatorMultiVariateKraskov1
calc = calcClass()
# 2. Set any properties to non-default values:
# No properties were set to non-default values

# Compute for all pairs:
for s in range(3):
    for d in range(3):
        # For each source-dest pair:
        if (s == d):
            continue
        source = data[:, s]
        destination = data[:, d]

        # 3. Initialise the calculator for (re-)use:
        calc.initialise()
        # 4. Supply the sample data:
        calc.setObservations(source, destination)
        # 5. Compute the estimate:
        result = calc.computeAverageLocalOfObservations()

        print("MI_Kraskov (KSG) alg. 1(col_%d -> col_%d) = %.4f nats" %
            (s, d, result))