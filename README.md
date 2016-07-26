TubeTK ComputeAirflowResistance Application
=============================================

#### Overview:

Computes airflow resistance of each tube

#### Usage:

```
ComputeAirflowResistance  [--returnparameterfile <std::string>]
                         [--processinformationaddress <std::string>]
                         [--xml] [--echo] [--mu <double>] [--]
                         [--version] [-h] <std::string> <std::string>
                         <std::string>


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

   --mu <double>
     viscosity (default: 20)

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <std::string>
     (required)  Input centerline TRE file

   <std::string>
     (required)  Binary tube segmentation mask

   <std::string>
     (required)  Output TRE file with radius of each centerline point set
     using the distance map


   Description: Computes Airflow Resistance of Each Tube

   Author(s): Deepak Roy Chittajallu (Kitware)
```


