# Unscented Kalman Filter Project
*Self-Driving Car Engineer Nanodegree Program*

This project utilizes an Unscented Kalman Filter to estimate the state of a moving object with noisy lidar and radar measurements.  It is run using the Term 2 Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases)

## Src Files
The UKF class housed in src/ukf.cpp has been reorganized for (mostly) proper member hiding.  In addition, method signatures have changed to reflect const correctness.  The call tree spreading from the ProcessMeasurement method is the result of breaking up functionality for unit testing.  With these unit tests, the code can be compared against expected answers from the lectures.

## Runtime

*Caution: If unit tests are enabled, the program will exit after the UKF constructor.*

INPUT: values provided by the simulator to the c++ program

["sensor_measurement"] => the measurment that the simulator observed (either lidar or radar)


OUTPUT: values provided by the c++ program to the simulator

["estimate_x"] <= kalman filter estimated position x
["estimate_y"] <= kalman filter estimated position y
["rmse_x"]
["rmse_y"]
["rmse_vx"]
["rmse_vy"]

---

## Dependencies
* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools](https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF` Previous versions use i/o from text files.  The current state uses i/o
from the simulator.


## Code Style

After too many years of coding, I had to make one change to [Google's C++ style guide](https://google.github.io/styleguide/cppguide.html);

```
void functionName()
{ //Opening bracket is here
}

```