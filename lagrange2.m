function obj = lagrange2(self, x, Ex, Ey, m, b)
% lagrange2(x, Ex, Ey, m, b) makes a generic heterogenous-step quadratic 
% interpolation by the Lagrange's polinomial of order two, and finds the root in
% which the curve described by pairs of points [Ex, Ey] intercepts with the 
% line of formula y = m*x + b.   
%
% This is an objective function, its aim is to be use as the function argument
% of a numerical methods to find root, it could be Newton's Method or maybe a
% wrapper for advanced methods like Matlab's fsolve()

Eyobj = m.*x + b;
max = find(Ex > x, 1, 'first');

if max >= 3 && max < length(Ex)
    obj = Eyobj - (Ey(max - 2)*(((x - Ex(max - 1))*(x - Ex(max+1)))/((Ex(max - 2) - Ex(max - 1))*(Ex(max - 2) - Ex(max+1)))) + Ey(max - 1)*(((x - Ex(max - 2))*(x - Ex(max+1)))/((Ex(max - 1) - Ex(max - 2))*(Ex(max - 1) - Ex(max+1)))) + Ey(max+1)*(((x - Ex(max - 2))*(x - Ex(max - 1)))/((Ex(max+1) - Ex(max-2))*(Ex(max+1) - Ex(max - 1))))); 
elseif max < 3
    min = max;
    obj = Eyobj - (Ey(min)*(((x - Ex(min + 1))*(x - Ex(min + 3)))/((Ex(min) - Ex(min + 1))*(Ex(min) - Ex(min + 3)))) + Ey(min + 1)*(((x - Ex(min))*(x - Ex(min + 3)))/((Ex(min + 1) - Ex(min))*(Ex(min + 1) - Ex(min + 3)))) + Ey(min + 3)*(((x - Ex(min))*(x - Ex(min + 1)))/((Ex(min + 3) - Ex(min))*(Ex(min + 3) - Ex(min + 1)))));
else
    obj = Eyobj - (Ey(max - 3)*(((x - Ex(max - 2))*(x - Ex(max)))/((Ex(max - 3) - Ex(max - 2))*(Ex(max - 3) - Ex(max)))) + Ey(max - 2)*(((x - Ex(max - 3))*(x - Ex(max)))/((Ex(max) - Ex(max - 3))*(Ex(max - 2) - Ex(max)))) + Ey(max)*(((x - Ex(max - 3))*(x - Ex(max - 2)))/((Ex(max) - Ex(max-3))*(Ex(max) - Ex(max)))));
end

end