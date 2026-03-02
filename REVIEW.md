# Revisión del programa Solve_trian_BBB (Resolver triángulo dadas sus bisectrices)

## Resumen

El programa resuelve un triángulo a partir de las longitudes de sus tres bisectrices
usando el método de Newton para sistemas no lineales. La aproximación inicial se obtiene
tratando las bisectrices como medianas mediante el Teorema de Apolonio.

Se encontraron **9 problemas**, clasificados por severidad.

---

## BUGS CRÍTICOS (el programa falla o da resultados incorrectos)

### BUG 1: Funciones `ecuarecta2p` y `fdist` no incluidas en el repositorio

**Archivo:** `incentro.m` (líneas 93, 94, 101, 105, 106, 113, 117, 118, 125)

Las funciones `ecuarecta2p` (ecuación de recta por dos puntos) y `fdist` (distancia entre
dos puntos) se usan en `incentro.m` pero no están incluidas en el repositorio. Sin ellas,
la fase de graficación y verificación falla con:

```
'ecuarecta2p' undefined
```

**Funciones faltantes necesarias:**

```matlab
function [m, b] = ecuarecta2p(P1, P2)
  % Ecuación de la recta por dos puntos: y = mx + b
  if P2(1) == P1(1)
    m = Inf; b = NaN;
  else
    m = (P2(2) - P1(2)) / (P2(1) - P1(1));
    b = P1(2) - m * P1(1);
  end
endfunction

function d = fdist(x1, y1, x2, y2)
  d = sqrt((x2-x1)^2 + (y2-y1)^2);
endfunction
```

---

### BUG 2: Variable `bp` en vez de `ap` en `incentro.m` (línea 96)

**Archivo:** `incentro.m`, líneas 95-99

Cuando el lado BC es vertical (`isnan(ordobc)`), el punto de intersección de la bisectriz
desde A con el lado BC se escribe incorrectamente en la variable `bp` en vez de `ap`:

```matlab
% CÓDIGO ACTUAL (INCORRECTO):
if isnan(ordobc)
    bp(1) = B(1); bp(2) = taap*bp(1)+ordora;  % <<< Debería ser ap, no bp
else
    ap = intersrect([taap,ordora],[mbc,ordobc]);
end
```

**Efecto:** Cuando el lado BC es vertical, `ap` queda sin definir y la línea 101
(`alfa = fdist(A(1),A(2),ap(1),ap(2))`) causa un error:
```
'ap' undefined
```

**Corrección:**
```matlab
if isnan(ordobc)
    ap(1) = B(1); ap(2) = taap*ap(1)+ordora;
else
    ap = intersrect([taap,ordora],[mbc,ordobc]);
end
```

**Verificado:** Se reprodujo el error con triángulo A=(4,0), B=(0,0), C=(0,3).

---

### BUG 3: Bisectrices verticales no manejadas en `incentro.m`

**Archivo:** `incentro.m`, líneas 91-125

El código solo maneja el caso cuando un **lado** del triángulo es vertical (pendiente
infinita), pero NO maneja cuando una **bisectriz** es vertical (es decir, cuando el
vértice y el incentro tienen la misma coordenada x).

**Ejemplo reproductor:** Triángulo isósceles simétrico respecto al eje y:
- A=(4,0), B=(0,0), C=(2,3): la bisectriz desde B es vertical → `gamma = NaN`
- A=(3,0), B=(0,0), C=(1.5, 2.598): la bisectriz desde C es vertical → `gamma = NaN`

Cuando `tccp = Inf` (o `taap`/`tbbp`), `ecuarecta2p` devuelve `ordorc = NaN`, y
`intersrect([Inf, NaN], ...)` produce resultados `NaN`.

**Corrección necesaria:** Agregar manejo para bisectrices verticales. Por ejemplo:
```matlab
if isinf(tccp)
    cp(1) = C(1);
    cp(2) = mab * cp(1) + ordoab;
else
    % ... código existente ...
end
```

---

## BUGS MENORES (resultados incorrectos en ciertos casos)

### BUG 4: `epsilon(2)` revisado dos veces, `epsilon(3)` nunca revisado

**Archivo:** `Solve_trian_BBB.m`, línea 152

```matlab
% CÓDIGO ACTUAL (INCORRECTO):
if abs(epsilon(1)) > tols || abs(epsilon(2)) > tols || abs(epsilon(2)) > tols

% CORRECCIÓN:
if abs(epsilon(1)) > tols || abs(epsilon(2)) > tols || abs(epsilon(3)) > tols
```

Error de copy-paste. La tercera bisectriz nunca se verifica contra la tolerancia.

---

### BUG 5: Variable `tolyi` declarada pero nunca usada

**Archivo:** `Solve_trian_BBB.m`, línea 26

```matlab
tolyi = 5e-14; % Para parte imaginaria en los lados calculados.
```

Esta variable se declara con un comentario que indica su propósito (verificar la parte
imaginaria de los lados), pero nunca se usa en el código. Después de Newton, no hay
verificación de que los lados sean reales.

**Efecto demostrado:** Con bisectrices de un triángulo casi degenerado (1, 1, 1.999),
Newton converge pero devuelve lados con parte imaginaria pequeña:
```
imag(a) = -1.34e-14, imag(b) = -1.34e-14, imag(c) = 1.90e-14
```

**Corrección sugerida:** Agregar después de la convergencia de Newton:
```matlab
if abs(imag(a)) > tolyi || abs(imag(b)) > tolyi || abs(imag(c)) > tolyi
    fprintf('Los lados tienen componente imaginaria significativa.\n')
    return
end
a = real(a); b = real(b); c = real(c);
```

---

### BUG 6: No se verifica que los lados sean positivos

**Archivo:** `Solve_trian_BBB.m`, líneas 57-60

Se verifica la desigualdad del triángulo (`a > b+c`, etc.) pero no se verifica que los
lados sean positivos (`a > 0, b > 0, c > 0`). Newton podría converger a valores
negativos que satisfagan las ecuaciones algebraicamente.

---

## PROBLEMAS DE CONVERGENCIA

### CASO 7: Oscilación cíclica con bisectrices de triángulos muy obtuso

Cuando una bisectriz es mucho menor que las otras (triángulo muy obtuso), Newton no
converge sino que oscila entre dos valores con período 2. Esto se observa con:

- Bisectrices (5, 5, 0.1): converge en 7 iteraciones (límite)
- Bisectrices (5, 5, 0.01): **NO converge en 51 iteraciones** (oscilación clara)
- Bisectrices (5, 5, 0.001): **NO converge en 51 iteraciones**

El autor ya notó este comportamiento (comentario en línea 63: "las iteraciones tienden a
hacerse cíclicas con período 2"). Posibles mejoras:
- Usar **amortiguamiento** (damped Newton): `x = x + λ·ds` con `0 < λ < 1`
- Usar el **promedio de las dos últimas iteraciones** cuando se detecta oscilación
- Cambiar a un método más robusto como Levenberg-Marquardt

---

## PROBLEMAS DE DISEÑO / ROBUSTEZ

### ISSUE 8: Bisectrices dibujadas solo hasta el incentro, no hasta el lado opuesto

**Archivo:** `Solve_trian_BBB.m`, líneas 138-142

Las bisectrices se dibujan desde cada vértice hasta el incentro:
```matlab
x = [A(1) inc(1)]; y = [A(2) inc(2)]; line(x,y,...)
```

Pero una bisectriz va desde el vértice hasta el punto donde interseca el lado opuesto
(pasando por el incentro). Los puntos `ap`, `bp`, `cp` (calculados en `incentro`)
representan esas intersecciones y deberían usarse:
```matlab
x = [A(1) ap(1)]; y = [A(2) ap(2)]; line(x,y,...)
```

---

### ISSUE 9: `intersrect` no detecta rectas paralelas

**Archivo:** `intersrect.m`

Cuando las dos rectas son paralelas (misma pendiente), la matriz del sistema es singular.
Octave no lanza error sino que devuelve una pseudo-solución (vía pseudoinversa) que es
matemáticamente incorrecta. Debería verificarse:

```matlab
if abs(a1 - a2) < eps
    error('Las rectas son paralelas, no hay intersección');
end
```

---

## RESUMEN DE CASOS DE FALLA ENCONTRADOS

| Caso | Bisectrices (α, β, γ) | Resultado |
|------|----------------------|-----------|
| Triángulo normal (3-4-5) | 4.2164, 3.3541, 2.4244 | Converge perfectamente en 4 iteraciones |
| Equilátero | 0.866, 0.866, 0.866 | Converge en 3 iteraciones |
| Isósceles | 4.879, 4.879, 4.000 | Converge en 4 iteraciones |
| Casi degenerado (1,1,1.999) | 1.333, 1.333, 0.032 | Converge en 9 iter pero con parte imaginaria residual |
| Muy obtuso (γ muy pequeña) | 5, 5, 0.01 | **NO converge** (oscilación cíclica) |
| Extremo | 5, 5, 0.001 | **NO converge** (oscilación cíclica) |
| Lado BC vertical | (recalculadas) | **ERROR: 'ap' undefined** (Bug 2) |
| Bisectriz vertical | (recalculadas) | **NaN en bisectrices** (Bug 3) |
| Imposible | 10, 1, 1 | Detectado correctamente por Apolonio |
| Escala grande | 1000, 800, 600 | Converge en 4 iteraciones |
| Escala pequeña | 0.001, 0.0008, 0.0006 | Converge en 6 iteraciones |

---

## CONCLUSIÓN

El núcleo matemático (Newton + Apolonio) funciona bien para triángulos "normales" y
converge rápidamente (3-4 iteraciones). Los problemas principales son:

1. **Funciones faltantes** (`ecuarecta2p`, `fdist`) que impiden la ejecución completa
2. **Bug de nombre de variable** (`bp` vs `ap`) que crashea con lados verticales
3. **Bisectrices verticales no manejadas** en triángulos simétricos
4. **No convergencia** para triángulos muy obtusos (bisectriz γ muy pequeña)
5. **Errores de copy-paste** menores (`epsilon(2)` duplicado, `tolyi` sin usar)
