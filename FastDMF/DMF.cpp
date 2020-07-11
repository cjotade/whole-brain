/**
 * C mex interface to fast DMF simulator to generate BOLD signals. See
 * README.md and DMF.hpp for details.
 *
 * Pedro Mediano, Apr 2020
 */

#include "mex.h"
#include "DMF.hpp"

/**
 * Cross-platform way of getting the pointer to the raw data in a mxArray. Uses
 * preprocessor directives to switch between pre- and post-2017 mex functions.
 *
 * NOTE: assumes that the array has doubles, but does not explicitly check.
 *
 */
double* safeGet(const mxArray* a) {
    double *p;
    #if MX_HAS_INTERLEAVED_COMPLEX
    p = mxGetDoubles(a);
    #else
    p = mxGetPr(a);
    #endif
    return p;
}


/**
 * Converts a Matlab struct with double fields into a std::map, keeping
 * pointers to all arrays without making any deep copies.
 *
 * @param structure_array_ptr pointer to Matlab struct
 *
 */
ParamStruct struct2map(const mxArray *structure_array_ptr) {
    mwSize total_num_of_elements;
    int number_of_fields, field_index;
    const char  *field_name;
    const mxArray *field_array_ptr;
    ParamStruct p;
    double *x;

    total_num_of_elements = mxGetNumberOfElements(structure_array_ptr);
    number_of_fields = mxGetNumberOfFields(structure_array_ptr);

    for (field_index=0; field_index<number_of_fields; field_index++)  {
        field_name = mxGetFieldNameByNumber(structure_array_ptr, field_index);
        field_array_ptr = mxGetFieldByNumber(structure_array_ptr,
                                             0, field_index);
        x = safeGet(field_array_ptr);
        p[field_name] = x;
    }

    return p;
}


/**
 * Check that input and output arguments to MEX function are valid.
 *
 * Makes sure that 1) number and type of inputs is correct, 2) struct has
 * all necessary fields, and 3) number of outputs is correct. If arguments
 * are not valid, throws a Matlab error.
 *
 * @throws Matlab error if input or output arguments are invalid
 */
void checkArguments( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]) {
    // Check that number and type of inputs is correct
    if(nrhs < 2 || nrhs > 3) {
        mexErrMsgIdAndTxt("DMF:nrhs","Two or three inputs required.");
    }

    /* make sure the first input argument is type struct */
    if (!mxIsStruct(prhs[0])) {
        mexErrMsgIdAndTxt("DMF:notStruct","First input must be a struct.");
    }
    
    /* make sure the second input argument is type double */
    if( !mxIsDouble(prhs[1]) ) {
        mexErrMsgIdAndTxt("DMF:notDouble","Second input must be type double.");
    }
    
    /* make sure the third input argument, if provided, is type char */
    if( (nrhs == 3) && !mxIsChar(prhs[2]) ) {
        mexErrMsgIdAndTxt("DMF:notChar","Third input (if provided) must be type char. This means using single quotes (\') instead of double quotes (\") where appropriate.");
    }
    

    // Check that parameter struct has the necessary fields
    std::vector<std::string> required_fields = {"C", "receptors", "dt",
            "taon", "taog", "gamma", "sigma", "JN", "I0", "Jexte", "Jexti",
            "w", "g_e", "Ie", "ce", "g_i", "Ii", "ci", "wgaine", "wgaini",
            "G", "alphas", "stren", "subsamp", "dtt", "batch_size"};
    for (const auto& field : required_fields) {
        mxArray *P = mxGetField(prhs[0], 0, field.c_str());
        if (!P) {
            std::cout << "Parameter struct does not have an entry for field: " << field;
            std::cout << "\nThe required fields for this struct are:\n";
            for (const auto& f : required_fields) {
                std::cout << f << ", ";
            }
            std::cout << "\b\b  \n\n";
            mexErrMsgIdAndTxt("DMF:struct","Missing fields in struct.");
        }
    }

    // Check that number of outputs matches input config
    size_t target_outs;
    if (nrhs == 2) {
        target_outs = 1;
    } else {
        char *buf = (char*) mxArrayToString(prhs[2]);
        std::string desired_out = std::string(buf);
        if (desired_out == "both") {
            target_outs = 2;
        } else {
            target_outs = 1;
        }
    }
    if (nlhs != target_outs) {
        mexErrMsgIdAndTxt("DMF:badOutputs","Wrong number of output arguments.");
    }

}


/**
 * Main function to call DMF from Matlab, using the C Mex interface.
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]) {
    checkArguments(nlhs, plhs, nrhs, prhs);

    // First input argument: parameter struct
    ParamStruct params =  struct2map(prhs[0]);
    size_t N = (size_t) mxGetN(mxGetField(prhs[0], 0, "C"));

    // Second input argument: number of steps
    size_t nb_steps = (size_t) safeGet(prhs[1])[0];

    // Third input argument: desired output
    bool return_rate=false, return_bold=true;
    if (nrhs > 2) {
      char *buf = (char*) mxArrayToString(prhs[2]);
      std::string desired_out = std::string(buf);
        if (desired_out == "both") {
            return_rate = true;
            return_bold = true;
        } else if (desired_out == "rate") {
            return_rate = true;
            return_bold = false;
        }
    }


    // Pre-allocate memory for results using Matlab's factory.
    // From here on these arrays are passed by reference and re-mapped,
    // avoiding any further memory allocations.
    size_t nb_steps_bold = nb_steps*params["dtt"][0]/params["subsamp"][0];
    size_t batch_size = params["batch_size"][0];

    mxArray *rate_res, *bold_res;
    rate_res = mxCreateDoubleMatrix(N, return_rate ? nb_steps : (2*batch_size), mxREAL);
    bold_res = mxCreateDoubleMatrix(N, nb_steps_bold, mxREAL);


    // Run, passing results by reference
    DMFSimulator sim(params, nb_steps, N, return_rate, return_bold);
    sim.run(safeGet(rate_res), safeGet(bold_res));


    // Copy back the results and return
    if (return_rate && return_bold) {
        plhs[0] = rate_res;
        plhs[1] = bold_res;
    } else if (return_rate) {
        plhs[0] = rate_res;
    } else {
        plhs[0] = bold_res;
    }

}

