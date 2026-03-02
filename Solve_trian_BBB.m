% Solve_trian_BBB                     ACM                      29 Mar 2025
% Integrador de todos los procesos para resolver un triángulo dadas las
% longitudes de sus bisectrices, alfa, beta y gama.
% 1).- Entrada de las tres longitudes, alfa, beta, gama.
% 2).- Aplicación del teorema de Apolonio como si fueran medianas, lo cual
%      da una aproximación inicial de los lados a, b y c.
% 3).- Usando esa aproximación, hallar los lados resolviendo el sistema no
%      lineal dado por Ricardo Ramírez, mediante el método de Newton para
%      sistemas no lineales.
% 4).- Si Newton tiene éxito, entrar a graficar el triángulo resuelto, di-
%      bujando el incentro, la bisectrices y el círculo inscrito.

% La solución de un triángulo conocidas sus bisectrices, es un problema muy
% complicado por métodos algebraicos o geométricos, a diferencia de lo que
% sucede cuando se conocen las medianas, las alturas o las mediatrices del
% triángulo y se desea hallar el resto de características (resolverlo).
% El sistema que se plantea, gracias a Ricardo Ramírez, en el caso de las
% bisectrices α,β y γ, aparece más abajo, en la función the_bisectors_system:

run carga

close all
clc
%warning('OFF','all')
tol = 5e-11; % Para el método de Newton.
tolyi = 5e-14; % Para parte imaginaria en los lados calculados.
tols = 5e-11; % Para las bisectrices recalculadas.
spc(1:30)=' '; sbr(1:75)='=';
NS = 'abc'; NV = 'ABC'; NB = 'αβγ'; BB = zeros(1,3);
fprintf('%sSolución de un triángulo dadas las tres bisectrices:\n',spc(1:10))
fprintf('%s\n',sbr)
% Pedir alfa, beta y gama como datos.
fprintf('\n%sIngrese las longitudes de las bisectrices:\n',spc(1:14))
ask = {' primera ','siguiente','siguiente'};
for ib = 1:3
    intext = sprintf('%sLongitud de la %s bisectriz:',spc(1:15),ask{ib});
    long = input(intext);
    BB(ib) = long;
end
BB = sort(BB,'descend'); % La bisectrices se ordenan en orden decreciente.
for ib = 1:3
    fprintf('%s%s = %7.4f\n',spc,NB(ib),BB(ib))
end
[ordB,idxB] = sort(BB,'descend'); % idxB: índice para mantener correspondencias luego
fprintf('\n')
tic
% Aplica el teorema de Apolonio:
[a,b,c] = Apolonio(BB(1),BB(2),BB(3));
if (a ~= real(a) || b ~= real(b)|| c ~= real(c))
    fprintf('%s¡No puede haber un triángulo con esas bisectrices!\n',spc(1:7))
    fprintf('%s Proceso terminado.\n',spc(1:23))
    return
end
% Resuelve el sistema no lineal por Newton
[a,b,c] = the_bisectors_system(BB(1),BB(2),BB(3),a,b,c);
toc
% Chequeo de que ningún lado puede ser mayor que la suma de los otros dos:
if (a > b+c || b > a+c || c > a+b)
    fprintf('Los lados generados no pueden formar un triángulo:\n')
    fprintf('a = %10.4f, b = %10.4f, c = %10.4f\n',a,b,c)
    fprintf('Se abandona el proceso en este punto.\n')
else
    % Grafica el resultado obtenido con incentro, bisectrices y círculo inscrito:
    S = sort([a,b,c]); % Los pongo en orden creciente
    a = S(1);b = S(2);c = S(3);
    fprintf('\n')
    % Calculo angA, angB y angC en radianes:
    angA = law_cos(a,b,c);angB = law_cos(b,c,a);angC = pi - (angA+angB);
    % Los paso a grados:
    fconv = 180/pi;
    anga = fconv*angA;angb = fconv*angB;angc = fconv*angC;
    % Angulo A (entre AB y CA):
    fprintf('Angulo A: %7.4f°\n',anga);
    % Angulo B (entre BC y BA):
    fprintf('Angulo B: %7.4f°\n',angb)
    % Angulo C (entre CA y BC):
    fprintf('Angulo C: %7.4f°\n',angc)
    % Adjudicar coordenadas a los vértices.
    % Criterio: poner el lado más largo (c) sobre el eje x, desde el origen.
    B = [0,0]; A = [c,0];
    % Calcula las coordenadas de A:
    C(1) = a*cos(angB); C(2) = a*sin(angB);
    V = real([A(1),A(2),B(1),B(2),C(1),C(2)]);
    for iv = 1:3
        jv1 = (iv-1)*2 + 1;
        fprintf('Vértice %s:(%7.4f,%7.4f)\n',NV(iv),V(jv1),V(jv1+1))
    end
    % Traza el triángulo ACB encontrado: (Antihorario)
    figure; hold on
    XE = [V(1),V(3),V(5),V(1)];YE = [V(2),V(4),V(6),V(2)];
    margen = 0.1515*(max(XE)-min(XE));
    xmin = min(XE)-margen; xmax = max(XE)+margen;xlim([xmin xmax]);
    ymin = min(YE)-margen; ymax = max(YE)+margen;ylim([ymin ymax]);
    desp = 0.0121*(xmax-xmin);
    line(XE,YE);
    % Nombres de los vértices:
    text(A(1)+0.5*desp, A(2)-desp,'A')
    text(B(1)-desp, B(2)-desp,'B')
    text(C(1), C(2)+desp,'C')
    % Nombres de los lados:
    text(0.5*(B(1)+C(1))-desp,0.5*(B(2)+C(2))+desp,'a')
    text(0.5*(A(1)+C(1)),0.5*(A(2)+C(2))+desp,'b')
    text(0.5*(A(1)+B(1)),0.5*(A(2)+B(2))-desp,'c')
    xlabel('Abscisas');ylabel('Ordenadas')
    title({'Círculo inscrito en un triángulo generado a partir de sus bisectrices.',...
        'Bisectrices e incentro.'})
    % Escribir en el gráfico las longitudes obtenidas:
    % Bisectrices: BB: datos. bb: recalculadas
    [a,b,c,ap,bp,cp,mab,mbc,mca,taap,tbbp,tccp,inc,tipo,bb(3),bb(1),bb(2)] = incentro(V);
    % Halla los ángulos:
    [~,~,~,angA,angB,angC,~,r,~] = triangle(V);
    fprintf('Lados: a = %7.4f, b = %7.4f, c = %7.4f\n',a,b,c)
    texto = sprintf(' a = %7.4f\n b = %7.4f\n c = %7.4f',a,b,c);
    text((min(XE)+max(XE)-3*margen)/2,ymax-0.4*margen,texto,'FontSize',8)
    % Dibujar el incentro y las bisectrices.
    fconv = 180/pi;
    anga = fconv*angA;angb = fconv*angB;angc = fconv*angC;
    % Diámetro del círculo inscrito:
    d_ins = 2*r;
    % Reordeno para asignar las correspondencias correctas:
    bb = sort(bb,'descend');
    % Longitudes de las bisectrices:
    texto = sprintf(' α = %7.4f\n β = %7.4f\n γ = %7.4f',bb(1),bb(2),bb(3));%'αβγ'
        text((min(XE)+max(XE)+margen)/2,ymax-0.4*margen,texto,'FontSize',8)
    for ib = 1:3
        epsilon(ib) = (bb(ib)-ordB(ib))/bb(ib);
    end
    epsilon = epsilon(idxB);
    bb = bb(idxB);
    % Se escriben las bisectrices recalculadas
    fprintf('Longitud de alfa  = AA'': %7.4f (Error relativo: %12.5e)\n',bb(1),epsilon(1))
    fprintf('Longitud de beta  = BB'': %7.4f (Error relativo: %12.5e)\n',bb(2),epsilon(2))
    fprintf('Longitud de gamma = CC'': %7.4f (Error relativo: %12.5e)\n',bb(3),epsilon(3))
    fprintf('Diámetro del círculo inscrito: %7.4f\n',d_ins)
    fprintf('Incentro: [%7.4f,%7.4f]\n',inc(1),inc(2))
    % Dibuja las bisectrices y el círculo inscrito.
    % De A:
    x = [A(1) inc(1)]; y = [A(2) inc(2)]; line(x,y,'color','red','LineWidth',1)
    % De B:
    x = [B(1) inc(1)]; y = [B(2) inc(2)]; line(x,y,'color','red','LineWidth',1)
    % De C:
    x = [C(1) inc(1)]; y = [C(2) inc(2)]; line(x,y,'color','red','LineWidth',1)
    p = plot(inc(1),inc(2),'k.');
    %p.MarkerSize = 3;
    %p.Marker = 'o';
    %p.MarkerFaceColor = 'g';%Para el incentro
    viscircles([inc(1),inc(2)],r,'color','red','LineWidth',1,'LineStyle','--');
    % Coloca la identificación negra sobre las bisectrices:
    text(0.5*(inc(1)+A(1))-0.5*desp,0.5*(inc(2)+A(2))+desp,'α','color','black')
    text(0.5*(inc(1)+B(1))-desp,0.5*(inc(2)+B(2))+desp,'β','color','black')
    text(0.5*(inc(1)+C(1))-1.5*desp,0.5*(inc(2)+C(2))+desp,'γ','color','black')
    if abs(epsilon(1)) > tols || abs(epsilon(2)) > tols || abs(epsilon(2)) > tols
        texto = sprintf('La imprecisión de los resultados excede la tolerancia de error.');
        text((min(XE)+max(XE)/2-2.1*margen),ymin,texto,'FontSize',8)
    end
    fprintf('Proceso concluido\n')
    axis square; axis  equal;
    grid on; grid minor
    %dock
end

