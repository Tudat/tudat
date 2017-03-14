/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Anderson Jr., J.D. , Fundamentals of Aerodynamics, 3rd edition, McGraw Hill, 2001.
 *      Gentry, A., Smyth, D., and Oliver, W. . The Mark IV Supersonic-Hypersonic Arbitrary Body
 *          Program, Volume II - Program Formulation, Douglas Aircraft Company, 1973.
 *      Anderson Jr., J.D, Hypersonic and High-Temperature Gas Dynamics, 2nd edition, AIAA
 *          Education Series, 2006.
 *
 */

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/Aerodynamics/rarefiedFlowSimpleGeometryDrag.h"

namespace tudat
{
namespace aerodynamics
{

double computeDragCoefficientOfCubeInRarefiedFlow( const std::vector< double >& input )
{
    using namespace tudat;


    std::vector<long double>  numberDensities {input[0] ,input[1], input[2] ,input[3] ,input[4] ,input[6] ,input[7] ,input[8]};
    const double airspeed = input[10];
    const double freestreamTemperature = input[9];
    const double freestreamDensity = input[5];

    const double Tw = 300;
    const double Na = 6.022140857e23;
    const double kb = 1.38064852e-23;
    const double K  = 7.50e-17;

    std::vector<long double>  mass(8);
    std::vector<double> ratios(8);
    std::vector<double> m{4, 16, 28, 32, 40, 1, 14, 16};


    for (int i = 0; i < 8 ; i = i +1){
        mass[i] =  (m[i]    / 1000) / Na;
        ratios[i] = numberDensities[i] *  mass[i] / freestreamDensity;
    }

    const double num_O = numberDensities[1];
    const double P     = num_O * freestreamTemperature;
    const double a     = (K * P) / (1 + K * P);
    const double R = kb * Na;

    std::vector<double> gamma  {0.0 , 1};

    double CD_c [8][2];

    for (int i = 0 ; i < 8 ; i = i +1 ) {

        double c_mp = sqrt(2* freestreamTemperature * R / ( Na * mass[i]));
        double S    = airspeed / c_mp;
        double G    = 1 / (2 * pow(S,2));
        double Q    = 1 + G;

        double v_ratio = sqrt(0.5 * (1 + a * (( ( 4 * (R / (Na * mass[i])) * Tw) / pow(airspeed,2) ) -1 ) ) );

        for (int ii = 0 ; ii < 2 ; ii = ii +1 ) {

            double Z = 1 + erf(gamma[ii] * S);
            double P = (1 / S) * exp( - pow(gamma[ii],2) * pow(S,2) );
            CD_c[i][ii] =  (P/sqrt(mathematical_constants::PI)) + gamma[ii] * Q * Z +
                            0.5 * gamma[ii] * v_ratio * (gamma[ii] * sqrt(mathematical_constants::PI) * Z + P);
        }
    }

    double CD_front = 0;
    double CD_side = 0;

    for (int i = 0 ; i < 8 ; i = i +1){
        CD_side = ratios[i] * CD_c[i][0] + CD_side;
    }

    for (int i = 0 ; i < 8 ; i = i +1){
        CD_front = ratios[i] * CD_c[i][1] + CD_front;
    }

    double CD = 4 * CD_side + CD_front;

    return CD;
}

double computeDragCoefficientOfSphereInRarefiedFlow( const std::vector< double >& input )
{

    using namespace tudat;

    std::vector<long double>  rho {input[0] ,input[1], input[2] ,input[3] ,input[4] ,input[6] ,input[7] ,input[8]};
    const double v = input[10];
    const double T = input[9];
    const double rho_tot = input[5];

    const double Tw = 300;
    const double Na = 6.022140857e23;
    const double kb = 1.38064852e-23;
    const double K  = 7.50e-17;

    std::vector<long double>  mass(8);
    std::vector<double> ratios(8);
    std::vector<double> m{4, 16, 28, 32, 40, 1, 14, 16};

    for (int i = 0; i < 8 ; i = i +1){
        mass[i] =  (m[i]    / 1000) / Na;
        ratios[i] = rho[i] *  mass[i] / rho_tot;
    }

    const double num_O = rho[1];
    const double P     = num_O * T;
    const double a     = (K * P) / (1 + K * P);

    std::vector<double> CD_s(8);

    for (int i = 0; i < 8 ; i = i+1){
        const long double T_ki = (mass[i] * pow(v,2)) / (3 * kb);
        const long double beta = pow((mass[i] / (2 * kb * T)),0.5);
        const long double S = v * beta;
        const long double T_k = T_ki * (1 - a) + a * Tw;
        CD_s[i] = ((2*pow(S,2)+1) / (sqrt(mathematical_constants::PI) * pow(S,3))) * exp(-pow(S,2)) +
                ( (4 * pow(S,4) + 4*pow(S,2) - 1) / (2 * pow(S,4))) * erf(S) +
                ( (2*sqrt(mathematical_constants::PI)) / (3*S)) * sqrt((T_k / T_ki))  ;
    }

    double CD = 0;

    for (int i = 0 ; i < 8 ; i = i +1){
        CD = ratios[i] * CD_s[i] + CD;
    }

    return CD;
}

} // namespace aerodynamics
} // namespace tudat
