## Copyright (C) 2025 DELL
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {} {@var{retval} =} numjacobi (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: DELL <DELL@DELL-CASA>
## Created: 2025-04-03

function [jac,f0] = numjacobi(func,x)
    % Returns the Jacobian matrix and f(x).
    % Adaptada de Numerical Methods in Engineering with MATLAB, de
    % Jaan Kiusalaas, Cambridge University Press, 2005.
    % Alejandro Casares,                         Agosto 20 de 2015
    h = 1.0e-4;
    n = length(x);
    jac = zeros(n);
    f0 = func(x(1),x(2:end));
    for i =1:n
        temp = x(i);
        x(i) = temp + h;
        f1 = func(x(1),x(2:end));
        x(i) = temp;
        jac(:,i) = (f1 - f0)/h;
    end
endfunction
