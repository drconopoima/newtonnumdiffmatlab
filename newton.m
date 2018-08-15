function [x1, devi, itera, backupx0, varargout] = newton(self, seed, func, dfuncdx, tol, functol, moreoutput, damping, maxiter)
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

if nargin < 3
    error('This functions requires at least 2 input arguments');
end
if nargin < 4
    dfuncdx = self.numericdiff3;
end
if nargin < 5
    tol = 1e-6;
end
if nargin < 6
    functol = 1e-8;
end
if nargin < 7
    moreoutput = 0;
end
if nargin < 8
    damping = 1;
end
if nargin < 9
    maxiter = 1000;
end

itera = 0;
x0 = seed;
n = length(seed);
m = nargout(func);
checkn = nargin(func);
if checkn ~= n
    error('The number of seeds must be equal to the number of functions input arguments')
end
devi = realmax*ones(1,n); % This is used to ensure the first iteration will enter the loop
backupx0 = zeros(1,n);
if moreoutput == 1
    backupdevi = zeros(1,n);
    backupdx = zeros(1,n);
    backupf0 = zeros(1,m);
    backupdf0 = zeros(m,n);
end
while any(devi > tol) && itera <= maxiter
    itera = itera+1;
    x0 = num2cell(x0);
    funvectorial = zeros(1,m);
    funvectorial = num2cell(funvectorial(:));
    [funvectorial{:}] = func(x0{:});
    funvectorial = -cell2mat(funvectorial); 
    x0 = cell2mat(x0);
    jacob = dfuncdx(x0, func, functol);
    linearsolution = jacob\funvectorial;
    dx = damping*(linearsolution);   
    x1 = x0 + dx';
    devi = abs(x0 - x1)./abs(x0);
    backupx0(itera, :) = x0;
    if moreoutput == 1
        backupdevi(itera, :) = devi;
        backupdx(itera, :) = dx;
        backupf0(itera, :) = funvectorial;
        backupdf0(1+(itera-1)*(m):m+(itera-1)*(m), :) = jacob;
    end
    x0 = x1;
end
backupx0(itera+1, :) = x0;

if moreoutput == 1
    for nout = 1:nargout
        if nout==1
            varargout{nout} = {backupdevi};
        elseif nout == 2
            varargout{nout} = {backupdx};
        elseif nout == 3
            varargout{nout} = {backupf0};
        elseif nout == 4
            varargout{nout} = {backupdf0};
        end
    end
end

end