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
## @deftypefn {} {@var{retval} =} the_bisectors_system (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: DELL <DELL@DELL-CASA>
## Created: 2025-04-03

  function [a,b,c] = the_bisectors_system(alfa,beta,gamma,a,b,c)
      %[a,b,c] = the_bisectors_system(alfa,beta,gamma,a,b,c)
      % Argumentos de entrada:
      %alfa,beta,gamma : longitudes de bisectrices de los ángulos A, B, C.
      %a, b, c : Aproximaciones iniciales a lados opuestos a ángulos A, B, C
      % Argumentos de salida:
      %a, b, c : Lados opuestos a ángulos A, B, C

      % Sistema no lineal a resolver:
      % [2bc/(b+c) ] * √[ s(s-a) / (bc) ] = α
      % [2ac/(a+c) ] * √[ s(s-b) / (ac) ] = β
      % [2ab/(a+b) ] * √[ s(s-c) / (ab) ] = γ
      % Donde s = (a+b+c)/2.

      % Las aproximaciones iniciales a a,b,c fueron obtenidas analogando
      % las bisectrices a medianas en las fórmulas de Apolonio.
      % ----------------------------------------------------------------------------------------
      t = a; y(1) = b; y(2) = c;

      % Sistema para esta función:
      f = @(t,y)...
      [(2*y(1)*y(2)/(y(1)+y(2))) * sqrt((t+y(1)+y(2))/2*((t+y(1)+y(2))/2-t)/(y(1)*y(2)))-alfa;...
       (2*t*y(2)/(t+y(2))) * sqrt((t+y(1)+y(2))/2*((t+y(1)+y(2))/2-y(1))/(t*y(2)))-beta;...
       (2*t*y(1)/(t+y(1))) * sqrt((t+y(1)+y(2))/2*((t+y(1)+y(2))/2-y(2))/(t*y(1)))-gamma];

      %% Solución numérica a partir de lo anterior

      % Método Newton: Iteraciones 3D
      fprintf('Valores iniciales (resultado de usar α,β y γ como medianas en el T. de Apolonio):\n')
      fprintf('t = %.4f, y(1) = %.4f, y(2) = %.4f\n',t,y)
      fprintf('\nSolución iterativa por Newton 3D:\n')
      fprintf('=================================\n')
      fprintf('it         t              y(1)            y(2)          f(t,y)       g(t,y))       h(t,y)\n')
      fprintf('------------------------------------------------------------------------------------------\n')

      % Solución del sistema de 3 ecuaciones no lineales.
      solved = false;
      nmax = 51; % Conveniente que sea impar porque las iteraciones tienden
                 %a hacerse cíclicas con período 2.
      tolx = 5e-10; toly = 5e-10;
      tic

      for n=1:nmax
          [J,~] = numjacobi(f,[t,y]);
          if cond(J)> 500 && n > 1
              fprintf('Condición de Jacobiano =%.3e\n',cond(J))
              fprintf('Termina el proceso iterativo sin convergencia.\n')
              return
          end

          F = f(t,y);
          ds = - J\F;
          t = t + ds(1); y = y + ds(2:end).';
          fprintf('%2d %15.8e %16.8e  %15.8e ',n,t,y)
          fprintf('%10.2e   %10.2e   %10.2e\n',F)

          if all(abs(F) < toly) || all(abs(ds) < tolx)
              solved = true;
              fprintf('\n')
              break
          end
      end
      if ~solved
          fprintf('No convergió en %d iteraciones\n',n)
          a = real(t); b = real(y(1)); c = real(y(2));
          return
      else
          a = t; b = y(1); c = y(2);
          fprintf('Convergió en %d iteraciones a:\n',n)
          fprintf('a = %19.9e\n',a)
          fprintf('b = %19.9e\n',b)
          fprintf('c = %19.9e\n',c)
      end
  endfunction
