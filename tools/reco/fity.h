//#ifndef fity_H
//#define fity_H

class fity {

  public:
//    fity(){} // Default ctor
    static int    gNVariables;        // Number of variables
    static int    gNCoefficients;     // Number of terms
    static double gDMean;             // Mean from training sample
    static double gXMean[];           // Mean from training sample
    static double gXMin[];            // Min from training sample
    static double gXMax[];            // Max from training sample
    static double gCoefficient[];     // Coefficients
    static double gCoefficientRMS[];  // Coefficients RMS
    static int    gPower[];           // Function powers
    double MDF(double* x);

};

//#endif fity_H
