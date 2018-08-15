# newtonnumdiffmatlab
Multivariate Newton-Raphson Method This script finds the root of a function by taking the function and its derivative, one seed value to start the iterations and one relative tolerance between successive values

newton(seed, func, dfuncdx, tol, functol, moreoutput, damping, maxiter)
% Finds the roots of the functions that func calculates using Newton's method
%    By default, the derivative is found using 5-point numerical differentiation
%    but user can input other functions in argument dfuncdx to use it for
%    evaluating the derivatives using exact derivatives or other techniques
%
%    Input arguments are:
%        - seed: The starting root guesses for evaluation using Newton's method
%        - func: An anonymous function that outputs the functions to find roots
%                    from
%    Optional input arguments are:
%        - dfuncdx=@numericdiff5. This is the method it will use for evaluating
%                    the partial derivatives. It could be provided an anonymous
%                    function so that the jacobian is analytical
%        - tol=1e-6. This is the relative error between two successive x1
%                    calculations
%        - functol=1e-8. This is the step size of the numerical difference.
%        - damping=1. This factor can be used to slow the change rate of the
%                    method, making shorter steps between two successive trials
%        - maxiter=1000. Number of total iterations before failing and exiting.
%        - moreoutput=0.
%    
%    
%    Output arguments are:
%        - x1: the roots found
%        - devi: the relative errors of the roots between last 2 iterations
%        - itera: number of iterations
%        - backupx0: for graphing purposes it saves root trials each iteration
%      Only if moreoutput is 1 (defaults to 0) it shows extra information
%        - backupdevi: errors in every iteration
%        - backupdx: correction factors every iteration
%        - backupf0: Function evaluation in every iteration
%        
%    Created by: Luis Jesus Diaz Manzo
%    2015/10/03 16:32
