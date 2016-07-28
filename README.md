TubeTK ComputeAirflowResistance Application
=============================================

#### Overview:

Computes airflow resistance of each tube

#### Usage:

```
ComputeAirflowResistance  [--returnparameterfile <std::string>]
                           [--processinformationaddress <std::string>]
                           [--xml] [--echo] [--stepLengthRelax
                           <double>] [--stepLengthFactor <double>]
                           [--numberOfIterations <int>]
                           [--terminationValue <double>] [--optimizer
                           <Gradient_Descent|Iterate_Neighborhood
                           |Regular_Step_Gradient_Descent>] [--endPoint
                           <std::vector<std::vector<float> >>] ... 
                           [--intermediatePoints
                           <std::vector<std::vector<float> >>] ... 
                           [--startPoint
                           <std::vector<std::vector<float> >>] ... 
                           [--mu <double>] [--] [--version] [-h]
                           <std::string> <std::string>


Where: 

   --returnparameterfile <std::string>
     Filename in which to write simple return parameters (int, float,
     int-vector, etc.) as opposed to bulk return parameters (image,
     geometry, transform, measurement, table).

   --processinformationaddress <std::string>
     Address of a structure to store process information (progress, abort,
     etc.). (default: 0)

   --xml
     Produce xml description of command line arguments (default: 0)

   --echo
     Echo the command line arguments (default: 0)

   --stepLengthRelax <double>
     Set Relaxation Factor. Only used with Regular Step Gradient Descent
     optimizer (default: 0.999)

   --stepLengthFactor <double>
     Optimizer Step Size factor. Only used with Iterate Neighborhood and
     Regular Step Gradient Descent optimizers (default: 0.1)

   --numberOfIterations <int>
     Maximum number of optimizer iterations. Only used with Gradient
     Descent and Regular Step Gradient Descent optimizers (default: 30000)

   --terminationValue <double>
     Minimum value to reach before Optimizer is terminated (default: 2)

   --optimizer <Gradient_Descent|Iterate_Neighborhood
      |Regular_Step_Gradient_Descent>
     Optimizer to extract path (default: Regular_Step_Gradient_Descent)

   --endPoint <std::vector<std::vector<float> >>  (accepted multiple times)
     End Point

   --intermediatePoints <std::vector<std::vector<float> >>  (accepted
      multiple times)
     Intermediate Points

   --startPoint <std::vector<std::vector<float> >>  (accepted multiple
      times)
     Start Point

   --mu <double>
     viscosity (default: 20)

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <std::string>
     (required)  Binary tube segmentation mask

   <std::string>
     (required)  Output centerline TRE file with radius of each centerline
     point set using the distance map


   Description: Compute Airflow Resistance

   Author(s): Deepak Roy Chittajallu (Kitware)

   Acknowledgements: This work is part of the TubeTK project at Kitware.
```


