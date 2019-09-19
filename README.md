# Double Gaussian Fitter
This is a convenient script for fitting a double Gaussian signal with a 
second order polynomial background. Such fit is very common in particle physics.
However, once in a while, the fit does not converge and the user has to 
change the initial values for the parameters and hope for better results. 

Now this is no longer the case. This script will estimate the number of signal
and background events based on sideband subtraction and use these numbers as the 
initial values. This is extremely convenient for particle selection stuties since 
the number of backgrounds and signals will vary a lot based on different cuts.

The fitting is done by RooFit so it is a Maximum likelihood fit and it is a *binned fit*. 

# Author
Dazhi W. @ 2019
