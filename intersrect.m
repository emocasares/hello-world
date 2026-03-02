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
## @deftypefn {} {@var{retval} =} intersrect (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: DELL <DELL@DELL-CASA>
## Created: 2025-04-03

% intersrect
% Función para hallar la intersección de dos rectas con ecuación y = ax + b
% function P = intersrect(ab1,ab2)
% Cada ab es un vector con los dos elementos, a y b, de una ecuación:
function P = intersrect(ab1,ab2)
  a1 = ab1(1); b1 = ab1(2);
  a2 = ab2(1); b2 = ab2(2);
  % Arma la matriz del sistema:
  A = [a1 -1; a2 -1]; V=-[ab1(2); ab2(2)];
  P = A\V;
  return
endfunction
