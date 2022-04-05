# pyui
Simple python script to run umbrella integration analysis

## Requirements
Python packages:
- numpy
- scipy
- pandas
- argparse

## Usage
```
ui.py [-h] -c CONFIG [-m XMIN] [-M XMAX] -n NBINS -t TEMP [-u {kj,kcal}] [-v] [-d]

Arguments:
  -h, --help            show this help message and exit
  -c CONFIG, --config CONFIG
                        name of configuration file
  -m XMIN, --xmin XMIN  minimum value of interest, float
  -M XMAX, --xmax XMAX  maximum value of interest, float
  -n NBINS, --nbins NBINS
                        number of bins
  -t TEMP, --temp TEMP  simulation temperature
  -u {kj,kcal}, --units {kj,kcal}
                        force constant units of energy
  -v, --verbose         verbose output
  -d, --data_hist       print histogram for each window
```

This program uses very similar data format as [WHAM program](http://membrane.urmc.rochester.edu/?page_id=126) developed by Alan Grossfield.
In the configuraiton file every line specifies a simulation window:
<pre>
path_to_timeseries_file_1    window_center_1   spring_constant_1
path_to_timeseries_file_2    window_center_2   spring_constant_2
...
</pre>
The umbrella potential is assumed to be harmonic:

![\frac{1}{2}k(x-x_0)^2](https://latex.codecogs.com/svg.image?\frac{1}{2}k(x-x_0)^2)


Each timeseries file should have two columns of floating point numbers. The first one contains time values, second is for reaction coordinate values, e.g.:
```
1.0   1.25
2.0   1.27
3.0   1.28
...
```
Time values from the first column of timeseries file are ignored and used for compatibility with WHAM program by Alan Grossfield.
