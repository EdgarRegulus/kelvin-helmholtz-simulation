# kelvin-helmholtz-simulation
Este repositorio contiene el **proyecto final de simulación numérica** sobre la expansión de un **remanente de supernova (RSN)**, desarrollado mediante un **código en Fortran 95** y análisis posterior de los resultados.  

El proyecto modela la evolución temporal de la onda de choque producida por una supernova y compara la tasa de expansión con el modelo teórico de **Taylor–Sedov**, incorporando y evaluando el efecto del **enfriamiento radiativo**.

---

## Objetivo del proyecto
Desarrollar y validar una simulación numérica en Fortran de la expansión de un remanente de supernova durante 50 000 años, resolviendo las ecuaciones de Euler en dos dimensiones y comparando cuantitativamente la evolución del radio del frente de choque con el modelo de Taylor–Sedov, para evaluar el impacto del enfriamiento radiativo mediante análisis de datos y visualización de resultados.

El núcleo del proyecto es un **solver numérico** implementado en Fortran, basado en:

- Método de **MacCormack** (esquema predictor–corrector de segundo orden).
- Condición de estabilidad **Courant–Friedrichs–Lewy (CFL)**.
- Implementación de **viscosidad artificial** para controlar oscilaciones numéricas.
- Manejo explícito de condiciones de frontera y simetría.

Se resuelven de forma conservativa:
- Masa
- Momento
- Energía

---

## Análisis de datos y postprocesado

Los resultados de la simulación se almacenan en archivos de salida que posteriormente son analizados con **scripts externos**, permitiendo:

- Extracción del radio del frente de choque en función del tiempo.
- Comparación cuantitativa con modelos teóricos.
- Generación de gráficas de evolución temporal.
- Análisis comparativo entre simulaciones con y sin enfriamiento radiativo.

Este flujo de trabajo reproduce un **pipeline típico de análisis de datos científicos**: simulación → extracción → limpieza → análisis → visualización.

---

## Habilidades demostradas

Este proyecto es relevante para perfiles de **programación científica y análisis de datos**, ya que demuestra:

- Programación estructurada en **Fortran** (alto desempeño computacional).
- Implementación de **métodos numéricos** para ecuaciones diferenciales parciales.
- Uso de estructuras de datos y control de memoria.
- Automatización del análisis mediante scripts.
- Interpretación física de datos simulados.
- Comparación entre modelos teóricos y resultados numéricos.
- Documentación técnica clara de un proyecto computacional.

---

## Estructura del repositorio

```text
supernova-expansion-simulation/
.
├── output/                # Salidas de la simulación sin enfriamiento radiativo y visualizaciones
├── output_cooling/        # Salidas de la simulación con enfriamiento radiativo y visualizaciones
│
├── Euler2D.f95            # Código principal del solver hidrodinámico 2D
├── Supernova.f95          # Definición del problema físico específico
├── parameters.mod         # Parámetros físicos y numéricos globales
├── constants.mod          # Constantes físicas (GAM, AMU, KB, etc.)
├── globals.mod            # Variables globales compartidas
│
├── Euler                  # Ejecutable generado tras la compilación
└── Euler2D.f95.o          # Archivo objeto generado por el compilador
```

## `Euler2D.f95` — Secuencia de la simulación

1. **Inicialización del entorno**
   - Se importan los módulos de constantes físicas y parámetros globales.
   - Se define la malla cartesiana uniforme en dos dimensiones y se asigna memoria para las variables conservadas.

2. **Configuración física**
   - Se leen los parámetros numéricos (CFL, solver, condiciones de frontera).
   - Se establecen las opciones de salida (formato ASCII o binario) y el directorio de escritura.

3. **Aplicación de la condición inicial**
   - Se asigna densidad, momento y energía en cada celda física.
   - Se inicializan las celdas fantasma de acuerdo con las condiciones de frontera seleccionadas.

4. **Bucle temporal**
   - Se calcula el paso de tiempo usando la condición de Courant–Friedrichs–Lewy (CFL).
   - Se evalúan los flujos hidrodinámicos mediante el solver seleccionado (Lax o MacCormack).
   - Se actualizan las variables conservadas en toda la malla.

5. **Procesos físicos adicionales**
   - Si está activado, se aplica el enfriamiento radiativo del gas.
   - Se considera la geometría axisimétrica cuando corresponde.

6. **Salida de resultados**
   - Se escriben archivos de salida en intervalos regulares definidos por `DTOUT`.
   - La simulación finaliza al alcanzar el tiempo total de integración (`TFIN`).

## `Supernova.f95`

El módulo del problema define una simulación de la expansión de un remanente de supernova en un medio interestelar homogéneo. Se especifican las propiedades físicas del entorno (densidad, temperatura, composición, velocidad inicial) y de la explosión (energía total, masa eyectada, radio inicial y centro geométrico). La condición inicial (IC_custom) construye un perfil esférico de la supernova, asignando densidad, presión y un campo de velocidades radial lineal dentro del radio inicial, mientras que fuera de este se mantiene el medio interestelar en equilibrio. La energía total se reparte entre componentes cinética y térmica, garantizando conservación de masa y energía.

---
