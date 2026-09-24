# Propuestas de funciones constitutivas $\Phi(\xi)$ y $M(\lambda)$ para el modelo Mt-ETV

**Marco:** energía libre de Helmholtz conformacional con dos tensores de conformación ($\mathbf c^C$, $\mathbf c^R$) y una variable estructural escalar $\lambda$, bajo las restricciones de la Segunda Ley derivadas en el borrador.

---

## 0. Errata y verificación de las derivadas del borrador

Antes de proponer nuevas funciones se verificaron las derivadas de la Sección "Thermodynamic Determination of the Coupling Function $\Phi(\xi)$".

### 0.1 Factor espurio $m$ en $\partial\hat f_c/\partial\lambda$

El borrador escribe

$$\frac{\partial\hat f_c}{\partial\lambda} = -\frac{G_R m}{2\rho}\lambda^{m-1} + \frac{\partial\Phi}{\partial\xi}e^{-I_R}.$$

Sin embargo, el término correspondiente en la energía libre es $-\dfrac{G_R}{2m\rho}\lambda^m$, cuya derivada es

$$\frac{\partial}{\partial\lambda}\left(-\frac{G_R}{2m\rho}\lambda^{m}\right) = -\frac{G_R}{2m\rho}\,m\,\lambda^{m-1} = -\frac{G_R}{2\rho}\lambda^{m-1},$$

**sin** el factor $m$. La nota de Andrés en la sección anterior es la correcta. La consecuencia es que la condición de punto crítico en $(\lambda=1, I_R=-3)$ da

$$\boxed{\;a = \frac{e^{-3m}}{3m+1}\;}\qquad\text{y no}\qquad a=\frac{m\,e^{-3m}}{3m+1}.$$

Para $m=3/2$: $a = \dfrac{2e^{-9/2}}{11} \approx 2.0198\times 10^{-3}$, y no $\approx 3.03\times10^{-3}$. **Este valor aparece en la implementación simbólica y debe corregirse**, pues escala toda la contribución $\Phi$ a $\partial\hat f_c/\partial\lambda$ y a $\dot\lambda_{\rm relax}$.

### 0.2 Segunda derivada

Con $a$ corregida, un cálculo directo (o el criterio general (A4) de la §1) da

$$\frac{\partial^2\hat f_c}{\partial\lambda^2}\bigg|_{\lambda=1,\,I_R=-3} = \frac{G_R}{2\rho}\,\frac{m}{3m+1}\;>\;0,$$

y **no** $\dfrac{G_R m^2}{2\rho(3m+1)}$. Para $m=3/2$ el valor es $\dfrac{3}{11}\dfrac{G_R}{2\rho}$, no $\dfrac{9}{22}$.

*Verificación cerrada.* Poniendo $G_R/2\rho=1$, $m=3/2$, $I_R=-3$ fijo, la parte $\lambda$-dependiente de $\hat f_c$ es
$F(\lambda)=\lambda^{3/2}\left[-\tfrac{4}{33}+\tfrac{2}{11}\ln\lambda\right]$, de donde $F'(\lambda)=\tfrac{3}{11}\lambda^{1/2}\ln\lambda$ y $F''(1)=\tfrac{3}{11}$. Nótese además que $F'(\lambda)<0$ para todo $\lambda\in(0,1)$: el mínimo en $\lambda=1$ es **global** en $(0,1]$, no sólo local. Este es un resultado fuerte que conviene reportar explícitamente en el artículo.

### 0.3 Errata tipográfica en $\dot\lambda_{\rm relax}$

En la expresión de Andrés

$$\dot\lambda_\text{relax} = \frac{\lambda^{m-1}e^{-(I_R+3)}}{\tau_\lambda(3m+1)}\left[(3m+1)e^{m(I_R+3)}+mI_R-m\ln\lambda-1\right]$$

el prefactor debe ser $e^{-m(I_R+3)}$ (no $e^{-(I_R+3)}$) para que sea consistente con el corchete. La forma correcta y simplificada es

$$\dot\lambda_{\rm relax} = \frac{\lambda^{m-1}}{\tau_\lambda(3m+1)}\Big[(3m+1)-e^{-m(I_R+3)}\big(m\ln\lambda - mI_R+1\big)\Big].$$

### 0.4 Terminología: mínimo, no máximo

La nota "one wants $\hat f_c$ to have a **maximum** with respect to $\lambda$" es un lapsus. Termodinámicamente lo que se requiere es un **mínimo** de la energía libre de Helmholtz en el estado en reposo, que es exactamente lo que garantiza $\partial^2\hat f_c/\partial\lambda^2>0$. La Sección 5.4 ya lo dice correctamente; conviene homogeneizar el texto.

---

## 1. Marco unificado: normalización y condiciones de admisibilidad

### 1.1 Normalización

Se propone escribir la función libre como

$$\Phi(\xi) \;=\; \frac{G_R}{2\rho}\,\varphi(\xi),\qquad \xi = \lambda e^{-I_R},\qquad I_R = \log|\mathbf c^R| - \operatorname{tr}(\mathbf c^R),$$

de modo que $\varphi$ es **adimensional**. La energía libre conformacional queda

$$\hat f_c = \frac{G_C}{2\rho}\Big[\operatorname{tr}(\mathbf c^C)-\log|\mathbf c^C|\Big] - \frac{G_R}{2m\rho}\lambda^m + \frac{G_R}{2\rho}\varphi(\xi).$$

Es útil introducir el **invariante de deformación de los rouleaux**

$$D_R \;\equiv\; -(I_R+3) \;=\; \operatorname{tr}(\mathbf c^R)-\log|\mathbf c^R|-3 \;\geq\; 0,$$

que se anula si y sólo si $\mathbf c^R=\boldsymbol\delta$ (desigualdad AM–GM para tensores SPD, la misma que se usa en el borrador). Entonces $\xi = \lambda\,e^{3}e^{D_R}$ y el estado de equilibrio ($\lambda=1$, $\mathbf c^R=\boldsymbol\delta$) corresponde a $\xi_{\rm eq}=e^3\approx 20.086$.

**Dominio admisible:** $\xi\in(0,\infty)$, con $\xi\to 0$ cuando $\lambda\to 0$. Este punto es crucial: cualquier condición sobre $\varphi'$ debe verificarse **hasta $\xi=0^+$**, no sólo en un entorno de $\xi_{\rm eq}$.

### 1.2 Derivadas fundamentales

$$\frac{\partial\hat f_c}{\partial\lambda} = \frac{G_R}{2\rho}\Big[\varphi'(\xi)\,e^{-I_R} - \lambda^{m-1}\Big],\qquad
\frac{\partial\hat f_c}{\partial I_R} = -\frac{G_R}{2\rho}\,\xi\,\varphi'(\xi),$$

$$\frac{\partial\hat f_c}{\partial \mathbf c^R} = \frac{\partial\hat f_c}{\partial I_R}\Big[(\mathbf c^R)^{-1}-\boldsymbol\delta\Big].$$

### 1.3 Las cuatro condiciones de admisibilidad

Toda propuesta debe satisfacer:

| | Condición | Origen |
|---|---|---|
| **(A1)** | $\varphi'(\xi)\geq 0\quad \forall\,\xi>0$ | Segunda Ley: $\dfrac{\partial\hat f_c}{\partial\mathbf c^R}:\dot{\mathbf c}^R_{\rm relax}\leq 0$ |
| **(A2)** | $M \geq 0$ | Segunda Ley (Onsager): $\dfrac{\partial\hat f_c}{\partial\lambda}\dot\lambda_{\rm relax} = -M\left(\dfrac{\partial\hat f_c}{\partial\lambda}\right)^2\leq 0$ |
| **(A3)** | $e^{3}\varphi'(e^{3}) = 1$ | Punto crítico en $(\lambda=1,\ \mathbf c^R=\boldsymbol\delta)$ |
| **(A4)** | $\displaystyle \sigma \equiv \frac{d\ln\varphi'}{d\ln\xi}\bigg|_{\xi=e^3} > m-1$ | Estabilidad: mínimo de $\hat f_c$ en $\lambda=1$ |

**Demostración de (A1).** Con $\dot{\mathbf c}^R_{\rm relax}=\tau_R^{-1}(\boldsymbol\delta-\mathbf c^R)$,

$$\frac{\partial\hat f_c}{\partial\mathbf c^R}:\dot{\mathbf c}^R_{\rm relax}
= -\frac{G_R}{2\rho}\,\frac{\xi\,\varphi'(\xi)}{\tau_R}\Big[\operatorname{tr}\big((\mathbf c^R)^{-1}\big)+\operatorname{tr}(\mathbf c^R)-6\Big].$$

Como $\xi\geq0$, $\tau_R>0$ y el corchete es $\geq 0$, la desigualdad $\leq 0$ equivale exactamente a $\varphi'(\xi)\geq0$. $\square$

**Demostración de (A3) y (A4).** De (1.2), $\partial\hat f_c/\partial\lambda=0$ en $(\lambda=1,I_R=-3)$ exige $\varphi'(e^3)e^3=1$. Derivando otra vez,

$$\frac{\partial^2\hat f_c}{\partial\lambda^2}\bigg|_{\rm eq} = \frac{G_R}{2\rho}\Big[e^{6}\varphi''(e^3)-(m-1)\Big]
= \frac{G_R}{2\rho}\Big[\underbrace{e^{3}\frac{\varphi''(e^3)}{\varphi'(e^3)}}_{=\;\sigma}\cdot\underbrace{e^{3}\varphi'(e^3)}_{=\;1}-(m-1)\Big],$$

de modo que se obtiene la identidad compacta

$$\boxed{\;\frac{\partial^2\hat f_c}{\partial\lambda^2}\bigg|_{\rm eq} = \frac{G_R}{2\rho}\Big[\sigma-(m-1)\Big].\;}\tag{$\star$}$$

La condición de mínimo es entonces $\sigma>m-1$. Para $m=3/2$: $\sigma>1/2$. $\square$

**Verificación del ($\star$) con el ansatz del borrador.** Para $\varphi_S(\xi)=a\,\xi^m\ln\xi$ se tiene $\varphi_S'=a\xi^{m-1}(m\ln\xi+1)$, luego $\sigma = (m-1)+\dfrac{m}{3m+1}$, y ($\star$) da $\dfrac{G_R}{2\rho}\dfrac{m}{3m+1}$ — coincidiendo con §0.2 y confirmando ambos cálculos por vía independiente.

### 1.4 Condición adicional (deseable, no obligatoria): estabilidad global

Para que en reposo el sistema evolucione monótonamente hacia $\lambda=1$ conviene exigir

$$\textbf{(A5)}\qquad e^{3}\varphi'(\lambda e^{3}) < \lambda^{m-1}\qquad \forall\,\lambda\in(0,1),$$

que garantiza $\partial\hat f_c/\partial\lambda<0$ en todo el intervalo y elimina barreras de nucleación espurias.

### 1.5 Producción de entropía

Con las elecciones anteriores, las tres contribuciones son:

$$\dot\sigma_C = \frac{G_C^2}{2T\rho\,\eta_C(\dot\gamma)}\Big[\operatorname{tr}(\mathbf c^C)+\operatorname{tr}((\mathbf c^C)^{-1})-6\Big]\geq0,$$

$$\dot\sigma_R = \frac{G_R}{2T\,\tau_R}\,\xi\,\varphi'(\xi)\Big[\operatorname{tr}(\mathbf c^R)+\operatorname{tr}((\mathbf c^R)^{-1})-6\Big]\geq0,$$

$$\dot\sigma_\lambda = \frac{\rho}{T}M\left(\frac{\partial\hat f_c}{\partial\lambda}\right)^2\geq0.$$

Obsérvese que $\xi\varphi'(\xi)$ es precisamente la cantidad normalizada a $1$ en el equilibrio por (A3): la condición de punto crítico admite la lectura física de **normalizar la intensidad del acoplamiento elástico estructura–deformación en el estado de reposo**.

### 1.6 Estructura general de $\dot\lambda_{\rm relax}$

Escribiendo $M = \dfrac{2\rho}{G_R\tau_\lambda}\,\mu$, con $\mu\geq0$ adimensional:

$$\dot\lambda_{\rm relax} = -M\frac{\partial\hat f_c}{\partial\lambda} = \frac{\mu}{\tau_\lambda}\Big[\lambda^{m-1}-\varphi'(\xi)e^{-I_R}\Big].$$

---

## 2. Las cinco propuestas

---

## Propuesta 1 — Red elástica de agregados con ley de potencias

### 2.1.1 Formulación matemática

$$\boxed{\;\varphi(\xi)=\frac{e^{-3n}}{n}\,\xi^{\,n},\qquad n>m\;}
\qquad\Longleftrightarrow\qquad
\Phi(\xi)=\frac{G_R\,e^{-3n}}{2\rho\,n}\,\xi^{\,n}.$$

$$\boxed{\;M(\lambda)=\frac{2\rho}{G_R\,\tau_\lambda}\,\lambda^{1-m}(1-\lambda)\;}$$

**Verificación de (A1).** $\varphi'(\xi)=e^{-3n}\xi^{\,n-1}>0$ para todo $\xi>0$, $n>0$. Se cumple **estrictamente en todo el dominio**, sin restricción de subdominios.

**Verificación de (A3).** $e^{3}\varphi'(e^{3}) = e^{3}\cdot e^{-3n}\cdot e^{3(n-1)} = e^{3}e^{-3}=1$. La condición se satisface **automáticamente por construcción**: no hay constante de ajuste que determinar. Ésta es la ventaja estructural de la propuesta.

**Verificación de (A4).** $\sigma = \dfrac{d\ln\varphi'}{d\ln\xi}=n-1$, luego por ($\star$)

$$\frac{\partial^2\hat f_c}{\partial\lambda^2}\bigg|_{\rm eq} = \frac{G_R}{2\rho}(n-m)>0 \iff n>m.$$

**Verificación de (A5).** $e^3\varphi'(\lambda e^3)=\lambda^{n-1}<\lambda^{m-1}$ para $\lambda<1$ si y sólo si $n>m$. La misma condición garantiza estabilidad local **y** global.

**Verificación de (A2) y cinética resultante.** De (1.2),

$$\frac{\partial\hat f_c}{\partial\lambda}=\frac{G_R}{2\rho}\lambda^{m-1}\Big[\lambda^{\,n-m}e^{\,nD_R}-1\Big],$$

donde se usó $e^{-nI_R}=e^{3n}e^{nD_R}$. Sustituyendo en (1.6):

$$\boxed{\;\dot\lambda_{\rm relax}=\frac{1-\lambda}{\tau_\lambda}\Big[\,1-\lambda^{\,n-m}e^{\,nD_R}\Big].\;}$$

Como $M\geq0$ para $\lambda\in[0,1]$, la Segunda Ley se satisface idénticamente vía Onsager:
$\dfrac{\partial\hat f_c}{\partial\lambda}\dot\lambda_{\rm relax}=-M\left(\dfrac{\partial\hat f_c}{\partial\lambda}\right)^2\leq0$.

### 2.1.2 Justificación e intuición física

La cinética resultante es notablemente transparente:

- **Término de agregación:** $\dfrac{1-\lambda}{\tau_\lambda}$, idéntico al término browniano/de agregación en reposo del t-ETV original.
- **Término de ruptura:** el factor $\lambda^{n-m}e^{nD_R}$. Como $D_R\geq0$ siempre, con igualdad si y sólo si los rouleaux están **indeformados**, la ruptura estructural está controlada exponencialmente por la deformación *interna* de los agregados, no por el corte impuesto. Esto materializa la cadena causal $\dot\gamma\to\mathbf c^R\to\lambda$ que el borrador ya identifica como la diferencia conceptual clave con el t-ETV.

Físicamente, $\varphi\propto\xi^n$ es la energía elástica de una red de agregados cuya rigidez escala como una potencia de la "extensión estructural" $\xi=\lambda e^{-I_R}$ (producto del grado de agregación por la deformación acumulada). Es la forma estándar de la elasticidad de redes de agregados fractales/coloidales, donde el exponente $n$ está relacionado con la dimensión fractal del agregado y con la conectividad de la red. La condición $n>m$ tiene lectura directa: **la energía elástica almacenada por el acoplamiento estructura–deformación debe crecer más rápido con la estructura que la energía de agregación pura** $-\lambda^m$; de otro modo el estado agregado no es estable.

**Parsimonia:** un único parámetro nuevo $n$ reemplaza a los tres parámetros empíricos $\{\tau_{\rm aggr},\tau_{\rm break},d\}$ del t-ETV.

### 2.1.3 Respaldo bibliográfico

- Beris, A. N., & Edwards, B. J. (1994). *Thermodynamics of Flowing Systems: With Internal Microstructure*. Oxford University Press. (Formalismo de bracket generalizado; energías libres de tipo potencia en tensores de conformación.)
- Larson, R. G., & Wei, Y. (2019). A review of thixotropy and its rheological modeling. *Journal of Rheology, 63*(3), 477–501. (Revisión de cinéticas estructurales de tipo potencia.)
- Mewis, J., & Wagner, N. J. (2009). Thixotropy. *Advances in Colloid and Interface Science, 147–148*, 214–227.
- Stephanou, P. S., & Georgiou, G. C. (2018). A nonequilibrium thermodynamics perspective of thixotropy. *The Journal of Chemical Physics, 149*(24), 244902. https://doi.org/10.1063/1.5049397

---

## Propuesta 2 — Extensibilidad finita del rouleau (tipo FENE): ruptura por deformación crítica

### 2.2.1 Formulación matemática

$$\boxed{\;\varphi(\xi) = -\,b\,\xi_c\ln\!\left(1-\frac{\xi}{\xi_c}\right),\qquad 0\le\xi<\xi_c,\qquad b=(1-u)e^{-3},\quad u\equiv\frac{e^3}{\xi_c}\;}$$

$$\boxed{\;M(\lambda,\mathbf c^R)=\frac{2\rho}{G_R\tau_\lambda}\,\lambda^{1-m}(1-\lambda)\left(1-\frac{\xi}{\xi_c}\right)\;}$$

**Verificación de (A1).** $\varphi'(\xi)=\dfrac{b}{1-\xi/\xi_c}>0$ para $b>0$ y $\xi<\xi_c$. $\;\checkmark$

**Verificación de (A3).** $e^{3}\varphi'(e^3)=\dfrac{b\,e^3}{1-u}=1 \Rightarrow b=(1-u)e^{-3}>0$ para $u<1$. $\;\checkmark$

**Verificación de (A4).**
$$\sigma = \frac{\xi\varphi''}{\varphi'}\bigg|_{e^3} = \frac{\xi/\xi_c}{1-\xi/\xi_c}\bigg|_{e^3}=\frac{u}{1-u}>m-1 \iff u>\frac{m-1}{m}.$$
Para $m=3/2$: $u>1/3$, es decir $\xi_c<3e^3\approx 60.3$. Combinando con $u<1$:
$$\boxed{\;e^3<\xi_c<\frac{m}{m-1}e^3\;}\qquad (m=3/2:\ \ 20.09<\xi_c<60.26).$$

**Verificación de (A2) y cinética resultante.** La elección de $M$ **cancela la singularidad**:

$$\dot\lambda_{\rm relax}=\frac{(1-\lambda)}{\tau_\lambda}\Big[\Big(1-\tfrac{\xi}{\xi_c}\Big)-b\,e^{-I_R}\lambda^{1-m}\Big]
=\frac{(1-\lambda)}{\tau_\lambda}\Big[1-u\lambda e^{D_R}-(1-u)\lambda^{1-m}e^{D_R}\Big].$$

En reposo ($D_R=0$) y $\lambda=1$: $1-u-(1-u)=0$. $\;\checkmark$ $M\geq0$ mientras $\xi<\xi_c$. $\;\checkmark$

### 2.2.2 Justificación e intuición física

Esta es la única propuesta que incorpora un **umbral de deformación crítica** para el agregado. La divergencia $\varphi\to\infty$ cuando $\xi\to\xi_c^-$ representa que un rouleau no puede sostener deformaciones arbitrarias: existe una elongación máxima más allá de la cual la energía necesaria es infinita y la estructura debe romperse. La restricción $\xi<\xi_c$ se traduce en

$$\lambda\, e^{D_R} < \frac{1}{u},$$

es decir, **el producto grado de agregación × deformación de los agregados está acotado**. Es la traducción continua de la observación experimental de que los doblete/rouleaux se disgregan por encima de un esfuerzo de corte crítico.

**Propiedad emergente: biestabilidad estructural.** Un análisis de (A5) muestra que esta forma **no** satisface la estabilidad global. Definiendo $H(\lambda)=\dfrac{(1-u)}{(1-u\lambda)\lambda^{m-1}}$, se tiene $H(1)=1$ y $H(0^+)=\infty$, de modo que existe un $\lambda_b\in(0,1)$ con $H(\lambda_b)=1$: un **punto crítico inestable** que separa el estado agregado ($\lambda=1$, mínimo) del estado disperso. Ejemplo para $m=3/2$:

| $u$ | $\xi_c$ | $\lambda_b$ (barrera) | $D_R^{\max}$ en $\lambda=1$ |
|---|---|---|---|
| $0.50$ | $40.2$ | $\approx 0.35$ | $0.693$ |
| $0.75$ | $26.8$ | $\approx 0.09$ | $0.288$ |
| $0.90$ | $22.3$ | $\approx 0.010$ | $0.105$ |

Esto **no es necesariamente un defecto**: la biestabilidad es el mecanismo canónico de la bifurcación de viscosidad e histéresis en fluidos tixotrópicos con esfuerzo de fluencia, fenomenología que la sangre exhibe. Pero implica que el modelo predecirá histéresis y sensibilidad a condiciones iniciales, y que $u$ debe elegirse cerca del extremo superior del rango admisible para reducir la barrera. **Advertencia práctica:** valores de $u$ próximos a $1$ restringen severamente la deformación admisible ($D_R^{\max}=-\ln u$ en $\lambda=1$), por lo que debe verificarse numéricamente que la trayectoria no alcance la singularidad. Es la propuesta de mayor riqueza física pero también la más delicada de integrar.

### 2.2.3 Respaldo bibliográfico

- Warner, H. R., Jr. (1972). Kinetic theory and rheology of dilute suspensions of finitely extendible dumbbells. *Industrial & Engineering Chemistry Fundamentals, 11*(3), 379–387.
- Bird, R. B., Curtiss, C. F., Armstrong, R. C., & Hassager, O. (1987). *Dynamics of Polymeric Liquids, Vol. 2: Kinetic Theory* (2nd ed.). Wiley.
- Owens, R. G. (2006). A new microstructure-based constitutive model for human blood. *Journal of Non-Newtonian Fluid Mechanics, 140*(1–3), 57–70.
- Moyers-González, M., Owens, R. G., & Fang, J. (2008). A non-homogeneous constitutive model for human blood. Part 1. Model derivation and steady flow. *Journal of Fluid Mechanics, 617*, 327–354.
- Coussot, P., Nguyen, Q. D., Huynh, H. T., & Bonn, D. (2002). Viscosity bifurcation in thixotropic, yielding fluids. *Journal of Rheology, 46*(3), 573–589.

---

## Propuesta 3 — Entropía combinatoria de agregación (ansatz de Stephanou regularizado)

### 2.3.1 Motivación: un defecto de signo en el ansatz actual

El ansatz del borrador, en la normalización de §1.1, es $\varphi_S(\xi)=a\,\xi^m\ln\xi$, de donde

$$\varphi_S'(\xi)=a\,\xi^{m-1}\big(m\ln\xi+1\big)\;<\;0\quad\text{para}\quad \xi<e^{-1/m}.$$

Como $\xi=\lambda e^{3}e^{D_R}\geq \lambda e^3$, la violación de (A1) ocurre para

$$\lambda < e^{-3-1/m}\qquad (m=3/2:\ \ \lambda<2.56\times10^{-2}),$$

y para agregados deformados el umbral es aún menor. **Es una región estrecha pero no vacía**, y justamente la que se visita en protocolos de inversión de flujo y corte alto, donde $\lambda\to0$. Cabe notar que la Segunda Ley sólo exige que la **suma** de las tres contribuciones sea no negativa, por lo que la violación podría estar compensada; pero el borrador adopta explícitamente la forma fuerte (por mecanismo), y con ella el ansatz falla en ese subdominio. La corrección natural es regularizar el logaritmo.

### 2.3.2 Formulación matemática

$$\boxed{\;\varphi(\xi)=a\,\xi^{m}\ln(1+\xi),\qquad a=\frac{e^{-3m}}{m\,L+r},\quad L\equiv\ln(1+e^3),\ \ r\equiv\frac{e^3}{1+e^3}\;}$$

$$\boxed{\;M(\lambda)=\frac{2\rho}{G_R\tau_\lambda}\,\lambda^{1-m}(1-\lambda)\;}$$

**Verificación de (A1).**
$$\varphi'(\xi)=a\left[m\,\xi^{m-1}\ln(1+\xi)+\frac{\xi^{m}}{1+\xi}\right]\;>\;0\quad\forall\,\xi>0,\ a>0.$$
Ambos sumandos son estrictamente positivos en $(0,\infty)$; la regularización $\ln\xi\to\ln(1+\xi)$ **elimina por completo el defecto de signo**. $\;\checkmark$

**Verificación de (A3).** Evaluando en $\xi=e^3$ (con $L\approx3.04858$, $r\approx0.95258$):
$$e^{3}\varphi'(e^3)=a\,e^{3m}\big[mL+r\big]=1\;\Longrightarrow\; a=\frac{e^{-3m}}{mL+r}.$$
Para $m=3/2$: $mL+r=5.52545$ y $a=e^{-4.5}/5.52545\approx 2.0105\times10^{-3}$. (Notablemente próximo al valor corregido del ansatz original, $2.0198\times10^{-3}$: la regularización casi no perturba el equilibrio, como debe ser.) $\;\checkmark$

**Verificación de (A4).** Con $\varphi''$ calculado por regla del producto,
$$\sigma=\frac{m(m-1)L+2mr-r^{2}}{mL+r}.$$
Para $m=3/2$:
$$\sigma=\frac{2.28644+2.85774-0.90741}{5.52545}=\frac{4.23677}{5.52545}\approx0.7668\;>\;\tfrac12=m-1.\;\checkmark$$
Por ($\star$): $\dfrac{\partial^2\hat f_c}{\partial\lambda^2}\Big|_{\rm eq}=\dfrac{G_R}{2\rho}(0.7668-0.5)=0.2668\,\dfrac{G_R}{2\rho}>0$. Para otros valores de $m$ la condición $\sigma(m)>m-1$ debe verificarse numéricamente (se cumple holgadamente en el rango $m\in[1,2]$).

**Cinética resultante.**
$$\dot\lambda_{\rm relax}=\frac{1-\lambda}{\tau_\lambda}\Big[1-a\lambda^{1-m}e^{-I_R}\Big(m\xi^{m-1}\ln(1+\xi)+\tfrac{\xi^{m}}{1+\xi}\Big)\Big],\qquad \xi=\lambda e^{3}e^{D_R}.$$

### 2.3.3 Justificación e intuición física

El factor $\ln(1+\xi)$ es la forma canónica de una **entropía combinatoria de mezcla** en el régimen diluido regularizado: en la teoría de Flory–Huggins y en los tratamientos de tipo Smoluchowski de la distribución de tamaños de agregados, el término $c\ln c$ se sustituye por $c\ln(1+c)$ para mantener finita la contribución cuando la población de agregados tiende a cero. El prefactor $\xi^m$ pondera esa entropía por el número de contactos célula–célula dentro de un rouleau de "tamaño" $\lambda^m$ — exactamente el escalado que el t-ETV ya asume en $\boldsymbol\tau^R\propto\lambda^m$.

Interpretación: $\Phi$ mide el costo entrópico de organizar los eritrocitos en pilas orientadas. Al deformarse los agregados ($D_R\uparrow$), $\xi$ crece, la entropía conformacional disponible cae, y el sistema es empujado hacia la desagregación. La ventaja frente a la Propuesta 1 es que preserva íntegramente la motivación estadística del ansatz de Stephanou que ya se está usando, con la mínima modificación necesaria para satisfacer (A1) globalmente.

**Parsimonia:** cero parámetros libres nuevos ($a$ queda fijada por (A3)). Es la propuesta **más parsimoniosa** de las cinco.

### 2.3.4 Respaldo bibliográfico

- Tsimouri, I. Ch., Stephanou, P. S., & Mavrantzas, V. G. (2018). A constitutive rheological model for agglomerating blood derived from nonequilibrium thermodynamics. *Physics of Fluids, 30*(3), 030710. https://doi.org/10.1063/1.5016913
- Stephanou, P. S., & Georgiou, G. C. (2018). A nonequilibrium thermodynamics perspective of thixotropy. *The Journal of Chemical Physics, 149*(24), 244902.
- Germann, N., Cook, L. P., & Beris, A. N. (2013). Nonequilibrium thermodynamic modeling of the structure and rheology of concentrated wormlike micellar solutions. *Journal of Non-Newtonian Fluid Mechanics, 196*, 51–57.
- Flory, P. J. (1942). Thermodynamics of high polymer solutions. *The Journal of Chemical Physics, 10*(1), 51–61.

---

## Propuesta 4 — Energía de adhesión por depleción con saturación (Neu–Meiselman)

### 2.4.1 Formulación matemática

$$\boxed{\;\varphi(\xi)=w\left[\xi+\xi_0\big(e^{-\xi/\xi_0}-1\big)\right],\qquad w=\frac{e^{-3}}{1-e^{-u}},\quad u\equiv\frac{e^3}{\xi_0}\;}$$

$$\boxed{\;M(\lambda)=\frac{2\rho}{G_R\tau_\lambda}\,\lambda^{2-m}(1-\lambda)\;}$$

**Verificación de (A1).**
$$\varphi'(\xi)=w\big[1-e^{-\xi/\xi_0}\big]\;\geq\;0\quad\forall\,\xi\geq0,\ w>0.$$
Además $\varphi'(0)=0$: **sin estructura no hay energía de adhesión**, lo cual es la condición de consistencia física correcta en el límite $\lambda\to0$. $\;\checkmark$

**Verificación de (A3).** $e^{3}\varphi'(e^3)=e^3w(1-e^{-u})=1\Rightarrow w=\dfrac{e^{-3}}{1-e^{-u}}>0$. $\;\checkmark$

**Verificación de (A4).**
$$\sigma=\frac{\xi\varphi''}{\varphi'}\bigg|_{e^3}=\frac{(\xi/\xi_0)e^{-\xi/\xi_0}}{1-e^{-\xi/\xi_0}}\bigg|_{e^3}=\frac{u}{e^{u}-1}\;>\;m-1.$$
La función $u\mapsto u/(e^u-1)$ es monótona decreciente desde $1$ (en $u\to0$) hasta $0$. Para $m=3/2$ la condición es $u/(e^u-1)>1/2$, cuya raíz es $u^*\approx1.256$; por tanto

$$\boxed{\;\xi_0>\frac{e^3}{u^*}\approx16.0\;}\qquad (m=3/2).$$

Nótese que la restricción es interpretable: $\sigma\in(0,1)$ acota esta familia a $m<2$, condición satisfecha por el valor $m=3/2$ del t-ETV.

**Cinética resultante.**
$$\dot\lambda_{\rm relax}=\frac{\lambda(1-\lambda)}{\tau_\lambda}\Big[1-w\,\lambda^{1-m}e^{-I_R}\big(1-e^{-\xi/\xi_0}\big)\Big].$$

### 2.4.2 Justificación e intuición física

La forma **saturante** $\varphi'(\xi)=w(1-e^{-\xi/\xi_0})$ es la firma de una energía de adhesión por depleción. En el modelo de Neu y Meiselman, la atracción entre eritrocitos surge del gradiente de presión osmótica generado cuando se solapan las capas de depleción de las macromoléculas plasmáticas (fibrinógeno, dextranos) adyacentes al glicocálix. Esa energía de interacción **satura** una vez que el solapamiento es completo: agregar más área de contacto ya no incrementa la afinidad. La escala $\xi_0$ es precisamente el análogo continuo del espesor de la capa de depleción reducido al espacio de la variable estructural.

La movilidad $M\propto\lambda^{2-m}(1-\lambda)$, que da el factor $\lambda(1-\lambda)$ en la cinética, es la **movilidad degenerada** estándar de los modelos de campo de fase (Cahn–Hilliard/Allen–Cahn) para variables confinadas a $[0,1]$. Garantiza estructuralmente que $\lambda=0$ y $\lambda=1$ sean estados invariantes de la parte relajacional, eliminando de raíz los desbordamientos numéricos del parámetro estructural.

**Advertencia honesta:** con esta movilidad $\lambda=0$ es absorbente para la parte relajacional (no hay nucleación desde el cero exacto). Físicamente $\lambda=0$ nunca se alcanza en tiempo finito, pero conviene inicializar con $\lambda_0>0$ y verificar la robustez numérica.

### 2.4.3 Respaldo bibliográfico

- Neu, B., & Meiselman, H. J. (2002). Depletion-mediated red blood cell aggregation in polymer solutions. *Biophysical Journal, 83*(5), 2482–2490.
- Baskurt, O. K., Neu, B., & Meiselman, H. J. (2011). *Red Blood Cell Aggregation*. CRC Press.
- Chien, S. (1970). Shear dependence of effective cell volume as a determinant of blood viscosity. *Science, 168*(3934), 977–979.
- Cahn, J. W., & Hilliard, J. E. (1958). Free energy of a nonuniform system. I. Interfacial free energy. *The Journal of Chemical Physics, 28*(2), 258–267.
- Elliott, C. M., & Garcke, H. (1996). On the Cahn–Hilliard equation with degenerate mobility. *SIAM Journal on Mathematical Analysis, 27*(2), 404–423.

---

## Propuesta 5 — Adsorción cooperativa tipo Langmuir + movilidad activada por deformación

> Ésta es la propuesta que **recupera la fenomenología del t-ETV desde primeros principios**, reemplazando el factor empírico $|\tau_{\rm aggr}\dot\gamma|^d$ por una función de variables de estado.

### 2.5.1 Formulación matemática

$$\boxed{\;\varphi(\xi)=A\left[\xi-\xi_c\ln\!\left(1+\frac{\xi}{\xi_c}\right)\right],\qquad A=(\xi_c+e^3)\,e^{-6},\quad \xi_c>e^3\;}$$

$$\boxed{\;M(\lambda,\mathbf c^R)=\frac{2\rho}{G_R\tau_\lambda}\,\lambda^{1-m}(1-\lambda)\Big[1+(\beta D_R)^{d/2}\Big],\qquad \beta>0,\ d>0\;}$$

**Verificación de (A1).** $\varphi'(\xi)=\dfrac{A\,\xi}{\xi_c+\xi}\geq0$ para todo $\xi\geq0$, $A>0$, con $\varphi'(0)=0$ y saturación $\varphi'\to A$ para $\xi\gg\xi_c$. $\;\checkmark$

**Verificación de (A3).** $e^{3}\varphi'(e^3)=\dfrac{A\,e^{6}}{\xi_c+e^3}=1\Rightarrow A=(\xi_c+e^3)e^{-6}>0$. $\;\checkmark$

**Verificación de (A4).** $\sigma=\dfrac{\xi_c}{\xi_c+\xi}\bigg|_{e^3}=\dfrac{\xi_c}{\xi_c+e^3}>m-1$. Para $m=3/2$: $\xi_c>e^3\approx20.09$. $\;\checkmark$

**Verificación de (A5) — estabilidad global.** Con $G(\lambda)\equiv e^3\varphi'(\lambda e^3)/\lambda^{m-1}=\dfrac{(\xi_c+e^3)\lambda^{2-m}}{\xi_c+\lambda e^3}$, se comprueba que $G'(\lambda)>0$ para $\xi_c>e^3$ y $m<2$, con $G(0)=0$, $G(1)=1$. Por tanto $G<1$ en $(0,1)$ y $\partial\hat f_c/\partial\lambda<0$ en todo el intervalo: **el mínimo en $\lambda=1$ es global y no hay barrera de nucleación**. $\;\checkmark$

**Verificación de (A2).** $M\geq0$ trivialmente, pues $D_R\geq0$, $\beta>0$, $\lambda\in[0,1]$. $\;\checkmark$

**Cinética resultante.**
$$\dot\lambda_{\rm relax}=\frac{(1-\lambda)}{\tau_\lambda}\Big[1+(\beta D_R)^{d/2}\Big]\left[1-\frac{(\xi_c+e^3)\lambda^{2-m}e^{D_R}}{\xi_c+\lambda e^{3}e^{D_R}}\right].$$

### 2.5.2 Justificación e intuición física

**(a) La función $\Phi$.** $\varphi'(\xi)=A\xi/(\xi_c+\xi)$ es exactamente una **isoterma de Langmuir**: la fracción de sitios de unión ocupados en la superficie eritrocitaria por las macromoléculas puenteantes. Para $\xi\ll\xi_c$ la respuesta es lineal (régimen diluido); para $\xi\gg\xi_c$ satura (superficie cubierta). La constante $\xi_c$ es la constante de disociación efectiva del proceso de unión célula–célula.

**(b) La movilidad y la recuperación del t-ETV.** Éste es el punto central. El t-ETV original modula la agregación con el factor empírico $\big[1+|\tau_{\rm aggr}\dot\gamma|^{d}\big]$. Esto es termodinámicamente incómodo, porque $\dot{\boldsymbol\gamma}$ **no es una variable de estado**: aparecer dentro de un coeficiente de relajación (disipativo) mezcla la cinemática impuesta con la respuesta interna del material. La sustitución propuesta es usar el invariante de estado $D_R$.

*Demostración de que el reemplazo es exacto al orden dominante.* Para un tensor de conformación upper-convected en corte simple estacionario con tiempo de relajación $\tau_R$ y $\mathrm{Wi}=\tau_R\dot\gamma$:

$$c^R_{12}=\mathrm{Wi},\quad c^R_{11}=1+2\mathrm{Wi}^2,\quad c^R_{22}=c^R_{33}=1
\;\Longrightarrow\; \operatorname{tr}\mathbf c^R=3+2\mathrm{Wi}^2,\quad |\mathbf c^R|=1+\mathrm{Wi}^2.$$

Luego, a orden dominante en $\mathrm{Wi}$,
$$D_R=\operatorname{tr}\mathbf c^R-\ln|\mathbf c^R|-3 = 2\mathrm{Wi}^2-\mathrm{Wi}^2+O(\mathrm{Wi}^4)=(\tau_R\dot\gamma)^2+O(\mathrm{Wi}^4),$$

de modo que
$$\big(\beta D_R\big)^{d/2}=\beta^{d/2}(\tau_R\dot\gamma)^{d}=\big|\tau_{\rm aggr}\dot\gamma\big|^{d}\quad\text{con}\quad \boxed{\;\tau_{\rm aggr}=\sqrt{\beta}\;\tau_R\;}.$$

**Se recupera exactamente el término empírico del t-ETV, con $\tau_{\rm aggr}$ ya no como parámetro libre sino ligado al tiempo de relajación de los rouleaux por un único factor adimensional $\beta$.** En flujos transitorios o con inversión, $D_R$ y $\dot\gamma^2$ dejan de ser proporcionales, y la formulación en $D_R$ es la físicamente correcta (respeta el retardo microestructural que el propio borrador identifica como $\dot\gamma\to\mathbf c^R\to\lambda$).

### 2.5.3 Respaldo bibliográfico

- Armstrong, M., Pincot, A., Jariwala, S., Horner, J., Wagner, N., & Beris, A. (2022). Tensorial formulations for improved thixotropic viscoelastic modeling of human blood. *Journal of Rheology, 66*(2), 327–347. https://doi.org/10.1122/8.0000346
- Öttinger, H. C. (2005). *Beyond Equilibrium Thermodynamics*. Wiley. (Dependencia admisible de la matriz disipativa respecto de las variables de estado.)
- Grmela, M., & Öttinger, H. C. (1997). Dynamics and thermodynamics of complex fluids. I. Development of a general formalism. *Physical Review E, 56*(6), 6620–6632.
- Dullaert, K., & Mewis, J. (2006). A structural kinetics model for thixotropy. *Journal of Non-Newtonian Fluid Mechanics, 139*(1–2), 21–30.
- Stephanou, P. S. (2020). A constitutive hemorheological model addressing both the deformability and aggregation of red blood cells. *Physics of Fluids, 32*(10), 103103. (Erratum: *Physics of Fluids, 33*, 039901, 2021.)

---

## 3. Tabla comparativa

| # | $\varphi(\xi)$ | $\varphi'(\xi)$ | $\sigma$ | Condición de estabilidad | Par. nuevos | (A5) global | Riesgo numérico |
|---|---|---|---|---|---|---|---|
| **1** | $\frac{e^{-3n}}{n}\xi^n$ | $e^{-3n}\xi^{n-1}$ | $n-1$ | $n>m$ | 1 ($n$) | Sí | Bajo |
| **2** | $-b\xi_c\ln(1-\xi/\xi_c)$ | $\frac{b}{1-\xi/\xi_c}$ | $\frac{u}{1-u}$ | $\frac{m-1}{m}<u<1$ | 1 ($\xi_c$) | No (biestable) | Alto (singularidad) |
| **3** | $a\,\xi^m\ln(1+\xi)$ | $a\big[m\xi^{m-1}L_\xi+\frac{\xi^m}{1+\xi}\big]$ | $\frac{m(m-1)L+2mr-r^2}{mL+r}$ | $\approx0.767>0.5$ ($m{=}3/2$) | **0** | Sí (verificar num.) | Bajo |
| **4** | $w[\xi+\xi_0(e^{-\xi/\xi_0}-1)]$ | $w(1-e^{-\xi/\xi_0})$ | $\frac{u}{e^u-1}$ | $\xi_0\gtrsim16.0$ ($m{=}3/2$) | 1 ($\xi_0$) | Sí (verificar num.) | Bajo ($\lambda{=}0$ absorbente) |
| **5** | $A[\xi-\xi_c\ln(1+\xi/\xi_c)]$ | $\frac{A\xi}{\xi_c+\xi}$ | $\frac{\xi_c}{\xi_c+e^3}$ | $\xi_c>e^3$ | 2 ($\xi_c,\beta$) | **Sí (demostrado)** | Bajo |

Recuento de parámetros frente al t-ETV (que usa $\tau_{\rm aggr}$, $\tau_{\rm break}$, $d$): las propuestas 1, 2 y 4 reducen de $11$ a $9$; la propuesta 3 reduce a $8$; la propuesta 5 a $10$ pero con la ventaja de la equivalencia asintótica demostrada.

---

## 4. Recomendación de implementación

**Orden sugerido de exploración numérica:**

1. **Propuesta 1** ($\varphi=\tfrac{e^{-3n}}{n}\xi^n$, $n>m$). Es la más robusta: (A1) y (A3) se cumplen idénticamente sin ajuste, (A4) y (A5) colapsan a la misma desigualdad $n>m$, y la cinética resultante es analíticamente transparente. Barrido sugerido: $n\in\{1.75,\,2,\,2.5,\,3\}$ para $m=3/2$. Debería ser directo de implementar en el marco simbólico existente, sustituyendo únicamente la expresión de $\Phi$ y $M$.

2. **Propuesta 3** (Stephanou regularizado). Preserva la motivación estadística ya adoptada, corrige el defecto de signo en $\lambda\lesssim0.026$, y no introduce parámetros nuevos. Es la comparación natural con los resultados que ya se tienen.

3. **Propuesta 5** para el artículo. Es la que da el argumento más fuerte: demuestra que el término empírico $|\tau_{\rm aggr}\dot\gamma|^d$ del t-ETV **emerge** como el límite de pequeña deformación de una movilidad de Onsager dependiente del estado, con $\tau_{\rm aggr}=\sqrt\beta\,\tau_R$. Es un resultado publicable por sí mismo.

4. **Propuestas 2 y 4** como exploración de fenomenología adicional (histéresis/bifurcación de viscosidad y saturación por depleción, respectivamente).

### 4.1 Observación sobre la acotación de $\lambda$

Independientemente de la elección de $\Phi$ y $M$, la ecuación de evolución

$$\dot\lambda = \mathbf g:\dot{\boldsymbol\gamma}+\dot\lambda_{\rm relax},\qquad \mathbf g:\dot{\boldsymbol\gamma}=-2\lambda\dot\gamma_0 c^R_{12}$$

**no garantiza por sí sola $\lambda\in[0,1]$**. En inversión de flujo, $\dot\gamma_0$ y $c^R_{12}$ pueden tener signos opuestos transitoriamente, con lo que el término cinemático se vuelve positivo y puede empujar $\lambda>1$. Ninguna elección de $M$ lo evita, porque $M$ sólo actúa sobre $\dot\lambda_{\rm relax}$.

El borrador ya deja abierta la puerta a la solución, al indicar que $\mathbf g$ puede ser "*a tensor proportional to this expression*". Basta tomar

$$\mathbf g = -h(\lambda)\,(\mathbf c^R-\boldsymbol\delta),\qquad h(\lambda)=\lambda(1-\lambda),$$

con lo cual el término cinemático se anula idénticamente en $\lambda=0$ y $\lambda=1$, y **la acotación queda garantizada estructuralmente**, sin recortes ni `Piecewise` en la integración numérica.

El costo es que cambia la EDP característica. Repitiendo la derivación de la §"Exploration" con $h$ general:

$$\frac{\partial \hat f_c}{\partial I_R} + h(\lambda)\frac{\partial \hat f_c}{\partial\lambda} = -\frac{G_R}{2\rho}\lambda^m,$$

cuyas características son $\dfrac{dI_R}{ds}=1$, $\dfrac{d\lambda}{ds}=h(\lambda)$, con invariante $I_R-H(\lambda)$ donde $H'=1/h$. Para $h=\lambda$ se recupera $H=\ln\lambda$ e $\xi=\lambda e^{-I_R}$, el caso del borrador. Para $h=\lambda(1-\lambda)$ se obtiene $H=\ln\!\big(\tfrac{\lambda}{1-\lambda}\big)$ y el invariante generalizado

$$\tilde\xi = \frac{\lambda}{1-\lambda}\,e^{-I_R},$$

que es la **variable logit** de la fracción estructural — una elección con excelente pedigrí en modelos de transición de fase y ocupación de sitios. Toda la maquinaria de las cinco propuestas se traslada intacta sustituyendo $\xi\to\tilde\xi$; sólo hay que recalcular el punto de equilibrio, que ahora es $\tilde\xi\to\infty$ cuando $\lambda\to1$, lo que exige reformular la condición (A3) como un límite (por ejemplo, exigiendo el punto crítico en $\lambda=1-\epsilon$, o eligiendo $\varphi$ con la asintótica adecuada). Merece un análisis aparte, pero resuelve simultáneamente el problema de acotación y da un invariante estructural más natural.

---

## 5. Referencias

### 5.1 Formato APA

Armstrong, M., Pincot, A., Jariwala, S., Horner, J., Wagner, N., & Beris, A. (2022). Tensorial formulations for improved thixotropic viscoelastic modeling of human blood. *Journal of Rheology, 66*(2), 327–347. https://doi.org/10.1122/8.0000346

Barnes, H. A. (1997). Thixotropy—A review. *Journal of Non-Newtonian Fluid Mechanics, 70*(1–2), 1–33.

Baskurt, O. K., Neu, B., & Meiselman, H. J. (2011). *Red blood cell aggregation*. CRC Press.

Beris, A. N., & Edwards, B. J. (1994). *Thermodynamics of flowing systems: With internal microstructure*. Oxford University Press.

Bird, R. B., Curtiss, C. F., Armstrong, R. C., & Hassager, O. (1987). *Dynamics of polymeric liquids, Vol. 2: Kinetic theory* (2nd ed.). Wiley.

Cahn, J. W., & Hilliard, J. E. (1958). Free energy of a nonuniform system. I. Interfacial free energy. *The Journal of Chemical Physics, 28*(2), 258–267.

Chien, S. (1970). Shear dependence of effective cell volume as a determinant of blood viscosity. *Science, 168*(3934), 977–979.

Coussot, P., Nguyen, Q. D., Huynh, H. T., & Bonn, D. (2002). Viscosity bifurcation in thixotropic, yielding fluids. *Journal of Rheology, 46*(3), 573–589.

Dullaert, K., & Mewis, J. (2006). A structural kinetics model for thixotropy. *Journal of Non-Newtonian Fluid Mechanics, 139*(1–2), 21–30.

Elliott, C. M., & Garcke, H. (1996). On the Cahn–Hilliard equation with degenerate mobility. *SIAM Journal on Mathematical Analysis, 27*(2), 404–423.

Flory, P. J. (1942). Thermodynamics of high polymer solutions. *The Journal of Chemical Physics, 10*(1), 51–61.

Germann, N., Cook, L. P., & Beris, A. N. (2013). Nonequilibrium thermodynamic modeling of the structure and rheology of concentrated wormlike micellar solutions. *Journal of Non-Newtonian Fluid Mechanics, 196*, 51–57.

Grmela, M., & Öttinger, H. C. (1997). Dynamics and thermodynamics of complex fluids. I. Development of a general formalism. *Physical Review E, 56*(6), 6620–6632.

Larson, R. G., & Wei, Y. (2019). A review of thixotropy and its rheological modeling. *Journal of Rheology, 63*(3), 477–501.

Mewis, J., & Wagner, N. J. (2009). Thixotropy. *Advances in Colloid and Interface Science, 147–148*, 214–227.

Mewis, J., & Wagner, N. J. (2012). *Colloidal suspension rheology*. Cambridge University Press.

Moyers-González, M., Owens, R. G., & Fang, J. (2008). A non-homogeneous constitutive model for human blood. Part 1. Model derivation and steady flow. *Journal of Fluid Mechanics, 617*, 327–354.

Neu, B., & Meiselman, H. J. (2002). Depletion-mediated red blood cell aggregation in polymer solutions. *Biophysical Journal, 83*(5), 2482–2490.

Öttinger, H. C. (2005). *Beyond equilibrium thermodynamics*. Wiley.

Owens, R. G. (2006). A new microstructure-based constitutive model for human blood. *Journal of Non-Newtonian Fluid Mechanics, 140*(1–3), 57–70.

Stephanou, P. S. (2020). A constitutive hemorheological model addressing both the deformability and aggregation of red blood cells. *Physics of Fluids, 32*(10), 103103. (Erratum: *Physics of Fluids, 33*, 039901, 2021.)

Stephanou, P. S., & Georgiou, G. C. (2018). A nonequilibrium thermodynamics perspective of thixotropy. *The Journal of Chemical Physics, 149*(24), 244902. https://doi.org/10.1063/1.5049397

Tsimouri, I. Ch., Stephanou, P. S., & Mavrantzas, V. G. (2018). A constitutive rheological model for agglomerating blood derived from nonequilibrium thermodynamics. *Physics of Fluids, 30*(3), 030710. https://doi.org/10.1063/1.5016913

Warner, H. R., Jr. (1972). Kinetic theory and rheology of dilute suspensions of finitely extendible dumbbells. *Industrial & Engineering Chemistry Fundamentals, 11*(3), 379–387.

### 5.2 BibTeX

```bibtex
@article{armstrong_tensorial_2022,
  author  = {Armstrong, Matthew and Pincot, Andre and Jariwala, Soham and Horner, Jeff and Wagner, Norman and Beris, Antony},
  title   = {Tensorial formulations for improved thixotropic viscoelastic modeling of human blood},
  journal = {Journal of Rheology}, volume = {66}, number = {2}, pages = {327--347}, year = {2022},
  doi     = {10.1122/8.0000346}
}

@article{tsimouri_constitutive_2018,
  author  = {Tsimouri, Ioanna Ch. and Stephanou, Pavlos S. and Mavrantzas, Vlasis G.},
  title   = {A constitutive rheological model for agglomerating blood derived from nonequilibrium thermodynamics},
  journal = {Physics of Fluids}, volume = {30}, number = {3}, pages = {030710}, year = {2018},
  doi     = {10.1063/1.5016913}
}

@article{stephanou_nonequilibrium_2018,
  author  = {Stephanou, Pavlos S. and Georgiou, Georgios C.},
  title   = {A nonequilibrium thermodynamics perspective of thixotropy},
  journal = {The Journal of Chemical Physics}, volume = {149}, number = {24}, pages = {244902}, year = {2018},
  doi     = {10.1063/1.5049397}
}

@article{stephanou_hemorheological_2020,
  author  = {Stephanou, Pavlos S.},
  title   = {A constitutive hemorheological model addressing both the deformability and aggregation of red blood cells},
  journal = {Physics of Fluids}, volume = {32}, number = {10}, pages = {103103}, year = {2020}
}

@article{germann_nonequilibrium_2013,
  author  = {Germann, Natalie and Cook, L. Pamela and Beris, Antony N.},
  title   = {Nonequilibrium thermodynamic modeling of the structure and rheology of concentrated wormlike micellar solutions},
  journal = {Journal of Non-Newtonian Fluid Mechanics}, volume = {196}, pages = {51--57}, year = {2013}
}

@article{owens_microstructure_2006,
  author  = {Owens, Robert G.},
  title   = {A new microstructure-based constitutive model for human blood},
  journal = {Journal of Non-Newtonian Fluid Mechanics}, volume = {140}, number = {1--3}, pages = {57--70}, year = {2006}
}

@article{moyersgonzalez_nonhomogeneous_2008,
  author  = {Moyers-Gonzalez, Miguel and Owens, Robert G. and Fang, Jiannong},
  title   = {A non-homogeneous constitutive model for human blood. {Part} 1. {Model} derivation and steady flow},
  journal = {Journal of Fluid Mechanics}, volume = {617}, pages = {327--354}, year = {2008}
}

@article{neu_depletion_2002,
  author  = {Neu, Bj{\"o}rn and Meiselman, Herbert J.},
  title   = {Depletion-mediated red blood cell aggregation in polymer solutions},
  journal = {Biophysical Journal}, volume = {83}, number = {5}, pages = {2482--2490}, year = {2002}
}

@book{baskurt_rbc_aggregation_2011,
  author    = {Baskurt, Oguz K. and Neu, Bj{\"o}rn and Meiselman, Herbert J.},
  title     = {Red Blood Cell Aggregation}, publisher = {CRC Press}, year = {2011}
}

@article{chien_shear_1970,
  author  = {Chien, Shu},
  title   = {Shear dependence of effective cell volume as a determinant of blood viscosity},
  journal = {Science}, volume = {168}, number = {3934}, pages = {977--979}, year = {1970}
}

@article{mewis_thixotropy_2009,
  author  = {Mewis, Jan and Wagner, Norman J.},
  title   = {Thixotropy},
  journal = {Advances in Colloid and Interface Science}, volume = {147--148}, pages = {214--227}, year = {2009}
}

@article{larson_review_2019,
  author  = {Larson, Ronald G. and Wei, Yufei},
  title   = {A review of thixotropy and its rheological modeling},
  journal = {Journal of Rheology}, volume = {63}, number = {3}, pages = {477--501}, year = {2019}
}

@article{dullaert_structural_2006,
  author  = {Dullaert, Konraad and Mewis, Jan},
  title   = {A structural kinetics model for thixotropy},
  journal = {Journal of Non-Newtonian Fluid Mechanics}, volume = {139}, number = {1--2}, pages = {21--30}, year = {2006}
}

@article{barnes_thixotropy_1997,
  author  = {Barnes, Howard A.},
  title   = {Thixotropy---a review},
  journal = {Journal of Non-Newtonian Fluid Mechanics}, volume = {70}, number = {1--2}, pages = {1--33}, year = {1997}
}

@book{beris_thermodynamics_1994,
  author    = {Beris, Antony N. and Edwards, Brian J.},
  title     = {Thermodynamics of Flowing Systems: With Internal Microstructure},
  publisher = {Oxford University Press}, year = {1994}
}

@book{ottinger_beyond_2005,
  author    = {{\"O}ttinger, Hans Christian},
  title     = {Beyond Equilibrium Thermodynamics}, publisher = {Wiley}, year = {2005}
}

@article{grmela_dynamics_1997,
  author  = {Grmela, Miroslav and {\"O}ttinger, Hans Christian},
  title   = {Dynamics and thermodynamics of complex fluids. {I}. {Development} of a general formalism},
  journal = {Physical Review E}, volume = {56}, number = {6}, pages = {6620--6632}, year = {1997}
}

@article{warner_kinetic_1972,
  author  = {Warner, Harold R.},
  title   = {Kinetic theory and rheology of dilute suspensions of finitely extendible dumbbells},
  journal = {Industrial \& Engineering Chemistry Fundamentals}, volume = {11}, number = {3}, pages = {379--387}, year = {1972}
}

@article{cahn_free_1958,
  author  = {Cahn, John W. and Hilliard, John E.},
  title   = {Free energy of a nonuniform system. {I}. {Interfacial} free energy},
  journal = {The Journal of Chemical Physics}, volume = {28}, number = {2}, pages = {258--267}, year = {1958}
}

@article{elliott_cahnhilliard_1996,
  author  = {Elliott, Charles M. and Garcke, Harald},
  title   = {On the {Cahn--Hilliard} equation with degenerate mobility},
  journal = {SIAM Journal on Mathematical Analysis}, volume = {27}, number = {2}, pages = {404--423}, year = {1996}
}

@article{flory_thermodynamics_1942,
  author  = {Flory, Paul J.},
  title   = {Thermodynamics of high polymer solutions},
  journal = {The Journal of Chemical Physics}, volume = {10}, number = {1}, pages = {51--61}, year = {1942}
}

@article{coussot_viscosity_2002,
  author  = {Coussot, Philippe and Nguyen, Q. D. and Huynh, H. T. and Bonn, Daniel},
  title   = {Viscosity bifurcation in thixotropic, yielding fluids},
  journal = {Journal of Rheology}, volume = {46}, number = {3}, pages = {573--589}, year = {2002}
}
```


