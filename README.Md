# Starting point exploration for [robust glm methods

Robust [GLM](http://en.wikipedia.org/wiki/Generalized_linear_model) (see e.g. [Cantoni & Ronchetti, 2001](http://www.jstor.org/stable/2670248)) often optimize non concave quasi-likelihoods. This can lead to convergence issues or to the existence of several optimum during the fitting process. Hence the starting point of the fitting algorithm must be chosen carefully.

We suggest that a Monte-Carlo exploration of starting point, with a creterion to select the best fit. 

## Files 
 - `illustrative_example_poisson.R` : Illustration of the problem using a robust poisson GLM. 
 - `poissonglm_function.R` : Functions to compute the score function of a robust poisson GLM. 