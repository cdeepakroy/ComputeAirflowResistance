/*=========================================================================

   Library:   TubeTK

   Copyright 2010 Kitware Inc. 28 Corporate Drive,
   Clifton Park, NY, 12065, USA.

   All rights reserved.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

=========================================================================*/

#include "tubeMessage.h"
#include "tubeCLIProgressReporter.h"
#include "itktubeSegmentTubesUsingMinimalPathFilter.h"

#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkGroupSpatialObject.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNotImageFilter.h"
#include "itkSpatialObjectReader.h"
#include "itkSpatialObjectWriter.h"
#include "itkTimeProbesCollectorBase.h"
#include "itkTubeSpatialObject.h"
#include "itkTubeSpatialObjectPoint.h"
#include "itkVesselTubeSpatialObject.h"
#include "itkVesselTubeSpatialObjectPoint.h"

#include <iostream>
#include <cmath>
#include <vector>

#include "ComputeAirflowResistanceCLP.h"

template< class TPixel, unsigned int TDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "<ModuleName>CLP.h"
//   and forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

template< class ImageType >
void WriteImage( typename ImageType::Pointer im, const std::string & name )
{
  typedef itk::ImageFileWriter< ImageType > WriterType;
  typename WriterType::Pointer writer  = WriterType::New();

  writer->SetInput( im );
  writer->SetFileName( name );
  writer->SetUseCompression( true );
  writer->Update();
}

//----------------------------------------------------------------------------
// Convert a list of point to ITK while converting them from RAS to LPS.
// Code taken from FiducialRegistration CLI
//
itk::Point< double, 3 >
ConvertPointFromRASToLPS(const itk::Point< double, 3 > & ptRAS)
{
  itk::Point< double, 3 > ptLPS;

  // convert RAS to LPS
  ptLPS[0] = -ptRAS[0];
  ptLPS[1] = -ptRAS[1];
  ptLPS[2] = ptRAS[2];

  return ptLPS;
}

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // Ensure that the input image dimension is valid
  // We only support 2D and 3D Images due to the
  // limitation of itkTubeSpatialObject
  if( VDimension != 2 && VDimension != 3 )
    {
    tube::ErrorMessage("Error: Only 2D and 3D data is currently supported.");
    return EXIT_FAILURE;
    }

  // setup progress reporting
  double progress = 0.0;

  tube::CLIProgressReporter progressReporter(
    "ComputeAirflowResistance", CLPProcessInformation );
  progressReporter.Start();
  progressReporter.Report( progress );

  // The timeCollector to perform basic profiling of algorithmic components
  itk::TimeProbesCollectorBase timeCollector;

  // Load tube segmentation mask
  tubeStandardOutputMacro( << "\n>> Loading tube segmentation mask" );

  typedef unsigned char                                MaskPixelType;
  typedef itk::Image< MaskPixelType, VDimension >      MaskImageType;
  typedef itk::ImageFileReader< MaskImageType >        MaskReaderType;

  timeCollector.Start( "Loading tube segmentation mask" );

  typename MaskReaderType::Pointer tubeMaskReader = MaskReaderType::New();

  try
    {
    tubeMaskReader->SetFileName( tubeMaskFile.c_str() );
    tubeMaskReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error loading tube segmention mask: "
      + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Loading tube segmentation mask" );
  progress = 0.1;
  progressReporter.Report( progress );

  // compute euclidean distance map of the intensity-negative of tube mask
  tubeStandardOutputMacro(
    << "\n>> Computing distance map of inverted tube mask" );

  timeCollector.Start( "Computing distance map of inverted tube mask" );

  typedef itk::NotImageFilter< MaskImageType, MaskImageType > NotFilterType;
  typename NotFilterType::Pointer notFilter = NotFilterType::New();

  notFilter->SetInput( tubeMaskReader->GetOutput() );

  typedef double                                DistanceMapPixelType;
  typedef itk::Image< DistanceMapPixelType,
    VDimension >                                DistanceMapType;
  typedef itk::DanielssonDistanceMapImageFilter<
    MaskImageType, DistanceMapType >            DistanceMapFilterType;

  typename DistanceMapFilterType::Pointer distanceMapFilter =
    DistanceMapFilterType::New();

  distanceMapFilter->SetInput( notFilter->GetOutput() );
  distanceMapFilter->SetInputIsBinary( true );
  distanceMapFilter->SetUseImageSpacing( true );
  distanceMapFilter->Update();

  // WriteImage< DistanceMapType >(distanceMapFilter->GetOutput(), "dmap.mha");

  timeCollector.Stop( "Computing distance map of inverted tube mask" );
  progress = 0.4;
  progressReporter.Report( progress );

  // Extract centerline
  tubeStandardOutputMacro( << "\n>> Extracting centerline" );

  timeCollector.Start( "Extracting centerline" );

  typedef itk::tube::SegmentTubesUsingMinimalPathFilter< VDimension,
    DistanceMapPixelType >                    CenterlineExtractionFilterType;
  typedef itk::Point< double, VDimension >    PointType;

  typename CenterlineExtractionFilterType::Pointer centerlineFilter =
    CenterlineExtractionFilterType::New();

  centerlineFilter->SetSpeedImage( distanceMapFilter->GetOutput() );

  if( startPoint.size() == 1 )
    {
    PointType startPathPoint;

    for( unsigned int i = 0; i < VDimension; i++ )
      {
      startPathPoint[i] = startPoint[0][i];
      }

    if( VDimension == 3 )
      {
      startPathPoint = ConvertPointFromRASToLPS( startPathPoint );
      }
    centerlineFilter->SetStartPoint( startPathPoint );
    }
  else
    {
    tubeErrorMacro(
      << "Error: start point must be provided and only one start point is allowed" );
    timeCollector.Stop( "Extracting centerline" );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  if( endPoint.size() == 1 )
    {
    PointType endPathPoint;

    for( unsigned int i = 0; i < VDimension; i++ )
      {
      endPathPoint[i] = endPoint[0][i];
      }

    if( VDimension == 3 )
      {
      endPathPoint = ConvertPointFromRASToLPS( endPathPoint );
      }

    centerlineFilter->SetEndPoint( endPathPoint );
    }
  else
    {
    tubeErrorMacro(
      << "Error: end point must be provided and only one end point is allowed" );
    timeCollector.Stop( "Extracting centerline" );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  if( intermediatePoints.size() > 0 )
    {
    std::vector< PointType > intermediatePathPoints;

    for( unsigned int i = 0; i < intermediatePoints.size(); i++ )
      {
      PointType curPathPoint;

      for( unsigned int j = 0; j < VDimension ; j++ )
        {
        curPathPoint[j] = intermediatePoints[i][j];
        }

      if( VDimension == 3 )
        {
        curPathPoint = ConvertPointFromRASToLPS( curPathPoint );
        }

      intermediatePathPoints.push_back( curPathPoint );
      }

    centerlineFilter->SetIntermediatePoints( intermediatePathPoints );
    }

  centerlineFilter->SetOptimizationMethod( optimizer );
  centerlineFilter->SetOptimizerTerminationValue( terminationValue );
  centerlineFilter->SetOptimizerNumberOfIterations( numberOfIterations );
  centerlineFilter->SetOptimizerStepLengthFactor( stepLengthFactor );
  centerlineFilter->SetOptimizerStepLengthRelax( stepLengthRelax );

  try
    {
    centerlineFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error extracting centerline: "
      + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Extracting centerline" );
  progress = 0.7;
  progressReporter.Report( progress );

  typedef typename CenterlineExtractionFilterType::InputSpatialObjectType
    TubeGroupType;
  typename TubeGroupType::Pointer tubeGroup = centerlineFilter->GetOutput();

  // Set tube radius using distance map
  tubeStandardOutputMacro(
    << "\n>> Setting tube radius using distance map" );

  timeCollector.Start( "Setting tube radius using distance map" );

  typedef itk::VesselTubeSpatialObject< VDimension >  TubeType;
  typedef typename TubeType::Pointer                  TubePointerType;
  typedef typename TubeGroupType::ChildrenListPointer TubeListPointerType;
  typedef typename TubeType::PointType                TubePointType;
  typedef typename TubeType::PointListType            TubePointListType;

  typedef itk::LinearInterpolateImageFunction< DistanceMapType, double >
    DistanceMapInterpolatorType;

  typename DistanceMapInterpolatorType::Pointer distanceMapInterpolator =
    DistanceMapInterpolatorType::New();
  distanceMapInterpolator->SetInputImage( distanceMapFilter->GetOutput() );

  typename MaskImageType::SpacingType spacing =
    tubeMaskReader->GetOutput()->GetSpacing();

  typename DistanceMapInterpolatorType::ContinuousIndexType prevPoint;
  double prevRadiusPhysp;

  char tubeName[] = "Tube";
  TubeListPointerType tubeList
    = tubeGroup->GetChildren( tubeGroup->GetMaximumDepth(), tubeName );

  for( typename TubeGroupType::ChildrenListType::iterator
    itTubes = tubeList->begin(); itTubes != tubeList->end(); ++itTubes )
    {
    TubePointerType curTube
      = dynamic_cast< TubeType * >( itTubes->GetPointer() );

    TubePointListType tubePointList = curTube->GetPoints();

    double airflow_resistance = 0;
    double tubeLength = 0;
    double minRadius = 0;
    double maxRadius = 0;
    double meanRadius = 0;

    double c = 8.0 * mu / M_PI;

    for( unsigned int ptId = 0; ptId < tubePointList.size(); ptId++ )
      {
      TubePointType curPosition = tubePointList[ptId].GetPosition();

      typename DistanceMapInterpolatorType::ContinuousIndexType curPoint;
      for( unsigned int i = 0; i < VDimension; i++ )
        {
        curPoint[i] = curPosition[i];
        }

      double curRadiusPhysp =
        distanceMapInterpolator->EvaluateAtContinuousIndex( curPoint );

      if( curRadiusPhysp < minRadius )
        {
        minRadius = curRadiusPhysp;
        }

      if( curRadiusPhysp > maxRadius )
        {
        maxRadius = curRadiusPhysp;
        }

      meanRadius += curRadiusPhysp;

      // Radius in TRE file is expected to be in continuous index space
      tubePointList[ptId].SetRadius( curRadiusPhysp / spacing[0] );

      // compute resistance
      if( ptId > 0 )
        {
        double dL = 0;
        for( unsigned int i = 0; i < VDimension; i++)
          {
          dL += pow( (curPoint[i] - prevPoint[i]) * spacing[i], 2 );
          }
        dL = sqrt( dL );

        tubeLength += dL;

        double r = 0.5 * (curRadiusPhysp + prevRadiusPhysp);

        airflow_resistance += c * ( dL / pow(r, 4) );
        }
      else
        {
        minRadius = curRadiusPhysp;
        maxRadius = curRadiusPhysp;
        }

      prevPoint = curPoint;
      prevRadiusPhysp = curRadiusPhysp;
      }
    curTube->SetPoints( tubePointList );

    std::cout << "TubeID = " << curTube->GetId() << std::endl
              << "\tairflow_resistance = " << airflow_resistance << std::endl
              << "\ttube_length = " << tubeLength << std::endl
              << "\tmin_radius = " << minRadius << std::endl
              << "\tmax_radius = " << maxRadius << std::endl
              << "\tmean_radius = " << meanRadius << std::endl;
    }

  timeCollector.Stop( "Setting tube radius using distance map" );
  progress = 0.9; // At about 90% done
  progressReporter.Report( progress );

  // Write tre file with radius of each point from the distance map
  tubeStandardOutputMacro( << "\n>> Writing output TRE file" );

  typedef itk::SpatialObjectWriter< VDimension > TubeWriterType;
  typename TubeWriterType::Pointer tubeWriter = TubeWriterType::New();

  timeCollector.Start( "Writing output TRE file" );

  try
    {
    tubeWriter->SetFileName( outputTREFile.c_str() );
    tubeWriter->SetInput( tubeGroup );
    tubeWriter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error writing TRE file: "
      + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Writing output TRE file" );
  progress = 1.0; // All Done
  progressReporter.Report( progress );

  // All done
  timeCollector.Report();
  return EXIT_SUCCESS;

}

// Main
int main( int argc, char * argv[] )
{
  try
    {
    PARSE_ARGS;
    }
  catch( const std::exception & err )
    {
    tube::ErrorMessage( err.what() );
    return EXIT_FAILURE;
    }

  PARSE_ARGS;

  return tube::ParseArgsAndCallDoIt( tubeMaskFile, argc, argv );

}
