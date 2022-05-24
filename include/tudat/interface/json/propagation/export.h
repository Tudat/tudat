/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_EXPORT_H
#define TUDAT_JSONINTERFACE_EXPORT_H

#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/simulation/estimation_setup/variationalEquationsSolver.h"
#include "tudat/interface/json/propagation/variable.h"

#include "tudat/interface/json/support/valueAccess.h"
#include "tudat/interface/json/support/valueConversions.h"

namespace tudat
{

namespace json_interface
{

class ExportSettings
{
public:
    //! Constructor.
    ExportSettings( const boost::filesystem::path& outputFile,
                    const std::vector< std::shared_ptr< propagators::VariableSettings > >& variables ) :
        outputFile_( outputFile ), variables_( variables ) { }

    //! Destructor.
    virtual ~ExportSettings( ) { }

    //! Set the path of the output file the results will be written to.
    void setOutputFile( const boost::filesystem::path& outputFile )
    {
        outputFile_ = outputFile;
    }

    //! Path of the output file the results will be written to.
    boost::filesystem::path outputFile_;

    //! Variables to export.
    //! The variables will be exported to a table in which each row corresponds to an epoch,
    //! and the columns contain the values of the variables specified in this vector in the provided order.
    std::vector< std::shared_ptr< propagators::VariableSettings > > variables_;

    //! Header to be included in the first line of the output file. If empty, no header will be added.
    std::string header_ = "";

    //! Whether to include the epochs in the first column of the results table.
    bool epochsInFirstColumn_ = true;

    //! Number of significant digits for the exported results.
    unsigned int numericalPrecision_ = 15;

    //! Whether to print only the values corresponding to the initial integration step.
    bool onlyInitialStep_ = false;

    //! Whether to print only the values corresponding to the final integration step.
    bool onlyFinalStep_ = false;

    //! Whether to show, in the terminal, the indices in the output vector where variables are saved
    bool printVariableIndicesToTerminal_ = false;

};

//! Create a `json` object from a shared pointer to a `ExportSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< ExportSettings >& saveSettings );

//! Create a shared pointer to a `ExportSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< ExportSettings >& saveSettings );


//! Export results of \p dynamicsSimulator according to the settings specified in \p exportSettingsVector.
/*!
 * @copybrief exportResultsOfDynamicsSimulator
 * \param singleArcDynamicsSimulator The dynamics simulator containing the results.
 * \param exportSettingsVector The vector containing export settings (each element represents a file to be exported).
 * \param variationalEquationsAreSaved Boolean denoting whether to save the variational equations (default is false).
 * \throws std::exception If any of the requested variables is not recognized or was not stored in the results of
 * \p dynamicsSimulator.
 */
template< typename TimeType = double, typename StateScalarType = double >
void exportResultsOfDynamicsSimulator(
        const std::shared_ptr< propagators::SingleArcDynamicsSimulator< StateScalarType, TimeType > >& singleArcDynamicsSimulator,
        const std::vector< std::shared_ptr< ExportSettings > >& exportSettingsVector,
        const bool variationalEquationsAreSaved = false )
{
    using namespace propagators;
    using namespace input_output;

    std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > statesHistory =
            singleArcDynamicsSimulator->getEquationsOfMotionNumericalSolution( );
    std::map< TimeType, Eigen::VectorXd > dependentVariables =
            singleArcDynamicsSimulator->getDependentVariableHistory( );
    std::map< TimeType, double > cpuTimes =
            singleArcDynamicsSimulator->getCumulativeComputationTimeHistory( );

    for ( std::shared_ptr< ExportSettings > exportSettings : exportSettingsVector )
    {
        std::vector< std::shared_ptr< VariableSettings > > variables;
        std::vector< unsigned int > variableSizes;
        std::vector< unsigned int > variableIndices;

        // Determine number of columns (not including first column = epoch).
        unsigned int cols = 0;

        if( exportSettings->printVariableIndicesToTerminal_ )
        {
            std::cout<<"Printing data to file: "<<exportSettings->outputFile_.filename( ).string( )<<
                       "  output vectors contain:"<<std::endl;
            std::cout<<"Vector entry, Vector contents"<<std::endl;
        }

        for ( std::shared_ptr< VariableSettings > variable : exportSettings->variables_ )
        {
            unsigned int variableSize = 0;
            unsigned int variableIndex = 0;
            switch ( variable->variableType_ )
            {
            case independentVariable:
            case cpuTimeVariable:
            {
                variableSize = 1;
                break;
            }
            case stateVariable:
            {
                variableSize = statesHistory.begin( )->second.rows( );
                break;
            }
            case dependentVariable:
            {
                const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVar =
                        std::dynamic_pointer_cast< SingleDependentVariableSaveSettings >( variable );
                assertNonnullptrPointer( dependentVar );
                try
                {
                    const std::string variableID = getDependentVariableId( dependentVar );
                    try
                    {
                        variableIndex = getKeyWithValue(
                                    singleArcDynamicsSimulator->getDependentVariableIds( ), variableID );
                        try
                        {
                            variableSize = getDependentVariableSaveSize( dependentVar );
                        }
                        catch ( ... )
                        {
                            std::cerr << "Could not export the results for variable \"" << variableID << "\" "
                                      << "because its size is not known." << std::endl;
                        }
                    }
                    catch ( ... )
                    {
                        std::cerr << "Could not export the results for variable \"" << variableID << "\" "
                                  << "because the main propagator was not configured to compute this variable."
                                  << std::endl;
                    }
                }
                catch ( ... )
                {
                    std::cerr << "Could not export results for dependent variable of type "
                              << dependentVar->dependentVariableType_
                              << " because its ID is not known." << std::endl;
                }
                break;
            }
            case stateTransitionMatrix:
            {
                if( !variationalEquationsAreSaved )
                {
                    throw std::runtime_error( "Error, requested save of state transition matrix, but only equations of motion available" );
                }
                break;
            }
            case sensitivityMatrix:
            {
                if( !variationalEquationsAreSaved )
                {
                    throw std::runtime_error( "Error, requested save of sensitivity matrix, but only equations of motion available" );
                }
                break;
            }
            default:
            {
                std::cerr << "Could not export results for variable of unsupported type "
                          << variable->variableType_ << "." << std::endl;
                break;
            }
            }

            if ( variableSize > 0 )
            {
                variables.push_back( variable );
                variableSizes.push_back( variableSize );
                variableIndices.push_back( variableIndex );

                if( exportSettings->printVariableIndicesToTerminal_ )
                {
                    std::cout<<cols<<", "<<getVariableId( variable )<<std::endl;
                }

                cols += variableSize;
            }
        }

        if( exportSettings->printVariableIndicesToTerminal_ )
        {
            std::cout<<std::endl;
        }

        // Concatenate requested results
        std::map< TimeType, Eigen::VectorXd > results;
        for ( auto it = statesHistory.begin( ); it != statesHistory.end( ); ++it )
        {
            if ( ( it == statesHistory.begin( ) && ! exportSettings->onlyInitialStep_ &&
                   exportSettings->onlyFinalStep_ ) ||
                 ( it == --statesHistory.end( ) && ! exportSettings->onlyFinalStep_ &&
                   exportSettings->onlyInitialStep_ )
                 || ( ( it != statesHistory.begin( ) && it != --statesHistory.end( ) ) &&
                      ( exportSettings->onlyInitialStep_ || exportSettings->onlyFinalStep_ ) ) )
            {
                continue;
            }
            unsigned int currentIndex = 0;

            const TimeType epoch = it->first;
            Eigen::VectorXd result = Eigen::VectorXd::Zero( cols );
            for ( unsigned int i = 0; i < variables.size( ); ++i )
            {
                const std::shared_ptr< VariableSettings > variable = variables.at( i );
                const unsigned int variableSize = variableSizes.at( i );

                switch ( variable->variableType_ )
                {
                case independentVariable:
                {
                    result.segment( currentIndex, variableSize ) =
                            ( Eigen::VectorXd( 1 ) << static_cast< double >( epoch ) ).finished( );
                    break;
                }
                case cpuTimeVariable:
                {
                    result.segment( currentIndex, variableSize ) =
                            ( Eigen::VectorXd( 1 ) << cpuTimes.at( epoch ) ).finished( );
                    break;
                }
                case stateVariable:
                {
                    result.segment( currentIndex, variableSize ) = it->second.template cast< double >( );
                    break;
                }
                case dependentVariable:
                {
                    result.segment( currentIndex, variableSize ) =
                            dependentVariables.at( epoch ).segment( variableIndices.at( i ), variableSize );
                    break;
                }
                case stateTransitionMatrix:
                {
                    break;
                }
                case sensitivityMatrix:
                {
                    break;
                }
                default:
                    break;
                }

                currentIndex += variableSize;
            }
            results[ epoch ] = result;
        }

        if ( exportSettings->epochsInFirstColumn_ )
        {
            // Write results map to file.
            writeDataMapToTextFile( results,
                                    exportSettings->outputFile_,
                                    exportSettings->header_,
                                    exportSettings->numericalPrecision_ );
        }
        else
        {
            // Write results matrix to file.
            Eigen::MatrixXd resultsMatrix( results.size( ), cols );
            int currentRow = 0;
            for ( auto entry : results )
            {
                resultsMatrix.row( currentRow++ ) = entry.second.transpose( );
            }
            writeMatrixToFile( resultsMatrix,
                               exportSettings->outputFile_.filename( ).string( ),
                               exportSettings->numericalPrecision_,
                               exportSettings->outputFile_.parent_path( ),
                               " ", exportSettings->header_ );
        }
    }
}

template< typename StateScalarType = double, typename TimeType = double >
void exportResultsOfVariationalEquations(
        const std::shared_ptr< propagators::SingleArcVariationalEquationsSolver< StateScalarType, TimeType > > variationalEquationsSolver,
        const std::vector< std::shared_ptr< ExportSettings > >& exportSettingsVector )
{
    using namespace propagators;
    using namespace input_output;

    for ( std::shared_ptr< ExportSettings > exportSettings : exportSettingsVector )
    {
        for ( std::shared_ptr< VariableSettings > variable : exportSettings->variables_ )
        {
            switch ( variable->variableType_ )
            {
            case stateTransitionMatrix:
            {
                if ( exportSettings->epochsInFirstColumn_ )
                {
                    // Write results map to file.
                    writeDataMapToTextFile( variationalEquationsSolver->getNumericalVariationalEquationsSolution( )[ 0 ],
                                            exportSettings->outputFile_,
                                            exportSettings->header_,
                                            exportSettings->numericalPrecision_ );
                }
                else
                {
                    throw std::runtime_error( "Error saving state transition/sensitivity matrix without epochs not yet supported" );
                }
                break;
            }
            case sensitivityMatrix:
            {
                if ( exportSettings->epochsInFirstColumn_ )
                {
                    // Write results map to file.
                    writeDataMapToTextFile( variationalEquationsSolver->getNumericalVariationalEquationsSolution( )[ 1 ],
                                            exportSettings->outputFile_,
                                            exportSettings->header_,
                                            exportSettings->numericalPrecision_ );
                }
                else
                {
                    throw std::runtime_error( "Error saving state transition/sensitivity matrix without epochs not yet supported" );
                }
                break;
            }
            default:
                break;
            }
        }
    }
}


} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_EXPORT_H
