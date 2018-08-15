function outp = numericdiff5(self, seed, func, functol)
% This function outputs the derivative of all functions in func by the method of
% using a 5 point numerical differentiation
%
% Created by: Luis Jesus Diaz Manzo
% 2015/10/03 15:24

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
    x1 = num2cell(seed);
    x3 = num2cell(seed);
    x4 = num2cell(seed);
    x0{argin} = x0{argin} - 2*functol;
    x1{argin} = x1{argin} - functol;
    x3{argin} = x3{argin} + functol;
    x4{argin} = x4{argin} + 2*functol;
    allargouts = zeros(1,m);
    allargouts = num2cell(allargouts);
    [allargouts{:}] = func(x0{:});
    f0 = cell2mat(allargouts);
    allargouts = num2cell(allargouts);
    [allargouts{:}] = func(x1{:});
    f1 = cell2mat(allargouts);
    allargouts = num2cell(allargouts);
    [allargouts{:}] = func(x3{:});
    f3 = cell2mat(allargouts);
    allargouts = num2cell(allargouts);
    [allargouts{:}] = func(x4{:});
    f4 = cell2mat(allargouts);

    outp(:, argin) = (f0' - 8*f1' + 8*f3' - f4') / (12 * functol);
end

end