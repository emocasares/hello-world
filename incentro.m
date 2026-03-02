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
## @deftypefn {} {@var{retval} =} incentro (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: DELL <DELL@DELL-CASA>
## Created: 2025-04-03

function [a,b,c,ap,bp,cp,mab,mbc,mca,taap,tbbp,tccp,inc,tipo,alfa,beta,gamma] = incentro(U,V,W)
% function [a,b,c,ap,bp,cp,mab,mbc,mca,taap,tbbp,tccp,inc,tipo,alfa,beta,gamma] = incentro(U,V,W)
% Función para calcular las intersecciones de las bisectrices de los
% ángulos de un triángulo ABC con los lados opuestos, dados sus
% vértices A(u1,u2),B(u3,u4),C(u5,u6) o A(u1,u2),B(v1,v2),C(w1,w2), según
% se haya dado un solo vector U o los tres puntos U,V,W como argumentos.
% Es un alias de la función bisint, pero incenter retorna además el incentro,
% el tipo de triángulo y las longitudes de las bisectrices.

    if nargin == 1
        x1 = U(1); y1 = U(2); x2 = U(3); y2 = U(4); x3 = U(5); y3 = U(6);
        V = U;
    elseif nargin == 3
        x1 = U(1); y1 = U(2);  x2 = V(1); y2 = V(2); x3 = W(1); y3 = W(2);
        V = [U(1),U(2),V(1),V(2),W(1),W(2)];
    end
    A = [x1,y1]; B = [x2,y2]; C = [x3,y3]; % Vértices leídos
    % Ángulos entre los lados:
    fconv = 180/pi;
    tol = 5e-11;
    % Pendientes de las rectas de los lados:
    mab = (y2-y1)/(x2-x1);
    mbc = (y3-y2)/(x3-x2);
    mca = (y1-y3)/(x1-x3);

    % Solución con función triangle:
    tipo = 1; % Escaleno
    [a,b,c,angA,angB,angC,area,r,R] = triangle(V);
    %   fprintf('Radio del círculo inscrito: %7.4f\n',r)
    if (abs(a-b) < tol && abs(b-c) < tol)
        %         fprintf('El triángulo dado es equilátero.\n')
        tipo = 3;
    elseif (abs(a-b) < tol || abs(b-c) < tol || abs(a-c) < tol)
        %         fprintf('El triángulo dado es isósceles.\n')
        tipo = 2;
    end
    % Angulo A (entre AB y AC):
    anga = angA*fconv;
    %      fprintf('Angulo A: %7.4f\n',anga);
    % Angulo B (entre BA y BC):
    angb = angB*fconv;
    %     fprintf('Angulo B: %7.4f\n',angb)
    % Angulo C (entre CA y CB):
    angc = angC*fconv;
    %    fprintf('Angulo C: %7.4f\n',angc)
    % Chequeo de la suma:
    if abs(anga+angb+angc - 180) >= 5e-8
        fprintf('Los ángulos no suman 180.')
        return
    end

    % Incentro:
    inc = [a*x1+b*x2+c*x3,a*y1+b*y2+c*y3]/(a+b+c);

    % Pendientes de bisectrices:
    taap = (inc(2)-y1)/(inc(1)-x1);
    tbbp = (inc(2)-y2)/(inc(1)-x2);
    tccp = (inc(2)-y3)/(inc(1)-x3);

    % El radio del círculo inscrito está calculado antes, en la
    % función triangle. (r)

    % Intersección de bisectrices con lados opuestos:
    % Se toman en cuenta los casos de lados verticales, que no tienen
    % ecuaciones de la forma y = mx + p

    % Bisectriz a angA:
    %fprintf('Pendiente de AA'': %7.4f\n',taap)
    [~,ordora] = ecuarecta2p(A,inc);
    [~,ordobc] = ecuarecta2p(B,C);
    if isnan(ordobc)
        bp(1) = B(1); bp(2) = taap*bp(1)+ordora;
    else
        ap = intersrect([taap,ordora],[mbc,ordobc]);
    end
    %fprintf('A'': (%7.4f,%7.4f)\n',ap(1),ap(2))
    alfa = fdist(A(1),A(2),ap(1),ap(2));

    % Bisectriz a angB:
    %fprintf('Pendiente de BB'': %7.4f\n',tbbp)
    [~,ordorb] = ecuarecta2p(B,inc);
    [~,ordoac] = ecuarecta2p(A,C);
    if isnan(ordoac)
        bp(1) = A(1); bp(2) = tbbp*bp(1)+ordorb;
    else
        bp = intersrect([tbbp,ordorb],[mca,ordoac]);
    end
    %fprintf('B'': (%7.4f,%7.4f)\n',bp(1),bp(2))
    beta = fdist(B(1),B(2),bp(1),bp(2));

    % Bisectriz a angC:
    %fprintf('Pendiente de CC'': %7.4f\n',tccp)
    [~,ordorc] = ecuarecta2p(C,inc);
    [~,ordoab] = ecuarecta2p(A,B);
    if isnan(ordoab)
        cp(1) = A(1); cp(2) = tccp*cp(1)+ordorc;
    else
        cp = intersrect([tccp,ordorc],[mab,ordoab]);
    end
    %fprintf('C'': (%7.4f,%7.4f)\n',cp(1),cp(2))
    gamma = fdist(C(1),C(2),cp(1),cp(2));
    return
endfunction
