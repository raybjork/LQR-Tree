# LQR-Tree

This project is a MATLAB implementation of the LQR-Tree algorithm for control of robotic systems as originally outlined in [this paper](https://groups.csail.mit.edu/robotics-center/public_papers/Tedrake09a.pdf). This algorithm seeks a series of controllers with regions of attraction that cover the controllable state space of a system.  This repository implements this algorithms on two test systems, (the underactuated torque-limited pendulum, and the cartpole) but is built in such a way that other systems could be tested using the same framework.  

We invite you to read our paper [An Exploration of Global Planning in LQR-Trees](https://github.com/raybjork/LQR-Tree/blob/master/An%20Exploration%20of%20Global%20Planning%20in%20LQR-Trees.pdf) to learn about the theory underlying this project and its implications.  

## Prerequisites

In order to run this program you will need to have the following installed on your machine:

* [MATLAB](https://www.mathworks.com/downloads/)
* [MOSEK](https://www.mosek.com/downloads/) - software used to solve optimization programs
* [Spotless](https://github.com/spot-toolbox/spotless) - toolbox to formulate programs in MOSEK 

## Running Tests

Explain how to run the automated tests for this system

## Future Work

## Authors

* **Ray Bjorkman**
* **David DePauw**

## Aknowledgments

This software was developed as a part of the class MEAM 517 (Nonlinear Control and Optimization for Underactuated Robotics) at the University of Pennsylvania.  Special thanks go to our professor, Michael Posa, for guidance in this project.  
