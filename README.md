# EEE588-Course-Project
This project deals with the design of LQR, H_infinity and LQG Controllers for the velocity and orientation control of a four wheeled robot, but with differential drive dynamics (Only rear wheel control). 

Open the pdf file that contains the report for the project. 

In the main code, the user will be asked to choose the controller of his choice, and the value of d (distance from center of gravity of robot to the midpoint of the line joining the rear wheels). If d > 0, the nominal system is stable. If d < 0, it is unstable. 

The code generates the transfer functions, step responses, and the bode plots for the sensitivity and complementary sensitivity for each of the controllers. 

The helper function contains some miscallaneous material that will assist the main code. 
