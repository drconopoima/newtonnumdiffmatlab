function outp = numericdiff3(self, seed, func, functol)
% This function outputs the derivative of all functions in func by the method of
% using a 3 point numerical differentiation
%
% Created by: Luis Jesus Diaz Manzo
% 2015/10/03 16:25

n = nargin(func);
m = nargout(func);
if nargin < 3
    error('Arguments must be at least 2')
end
if nargin < 4
    functol = 1e-8; % The size of step 
end

outp = zeros(m,n);

for argin = 1:n
    x0 = num2cell(seed);
    x2 = num2cell(seed);
    x0{argin} = x0{argin} - functol;
    x2{argin} = x2{argin} + functol;
    allargouts = zeros(1,m);
    allargouts = num2cell(allargouts);
    [allargouts{:}] = func(x0{:});
    f0 = cell2mat(allargouts);
    allargouts = num2cell(allargouts);
    [allargouts{:}] = func(x2{:});
    f2 = cell2mat(allargouts);
    outp(:, argin) = (f2' - f0') / (2 * functol);
end

end