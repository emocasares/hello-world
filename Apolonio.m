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
## @deftypefn {} {@var{retval} =} Apolonio (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: DELL <DELL@DELL-CASA>
## Created: 2025-04-03

function [a,b,c] = Apolonio(Ma,Mb,Mc)
  % Según el teorema de Apolonio:
    a = sqrt((2*Mb^2+2*Mc^2-Ma^2)/2);
    b = sqrt((2*Ma^2+2*Mc^2-Mb^2)/2);
    c = sqrt((2*Ma^2+2*Mb^2-Mc^2)/2);
    return
endfunction
