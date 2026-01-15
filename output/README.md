## output/

Esta carpeta contiene las **salidas numéricas de la simulación hidrodinámica sin enfriamiento radiativo**. Incluye archivos generados en distintos tiempos de integración, utilizados para analizar la evolución del remanente bajo un modelo adiabático y para producir visualizaciones comparativas.


## `plot_SNR.py`

Este script se encarga de **leer y visualizar las salidas binarias** generadas por el solver hidrodinámico de Euler en 2D. Para cada instante de tiempo, reconstruye las variables físicas en la malla y genera mapas bidimensionales de **densidad**, **magnitud de la velocidad**, **presión** y **temperatura**. El programa automatiza la lectura secuencial de archivos, convierte las unidades a escalas astrofísicas interpretables y guarda las figuras resultantes, facilitando el análisis comparativo de la evolución temporal del remanente de supernova.

## `radio_SNR.py` — Evolución temporal del remanente

Este script estima y grafica el **radio de expansión del remanente de supernova (RSN)** a partir de las salidas de la simulación hidrodinámica. Para cada instante, lee los campos de velocidad, calcula su magnitud y localiza la **discontinuidad más pronunciada** asociada al frente de choque, a partir de la cual se determina el radio respecto al centro de la explosión. El radio obtenido se almacena como función del tiempo y se compara con la **solución analítica de Sedov–Taylor**, permitiendo validar numéricamente el comportamiento autosimilar del remanente y analizar la coherencia física del modelo.

Por limitaciones de almacenamiento, este repositorio **no incluye todas las salidas generadas por la simulación** (Tan solo la primera y última salida). Para reproducir completamente los resultados y obtener los archivos de salida necesarios para el análisis y las visualizaciones, es necesario **compilar y ejecutar el código principal `Euler2D.f95`**, el cual genera automáticamente los datos utilizados por los scripts de post-proceso. Ejecute el código **sin enfriamiento radiativo**, verificando que `USE_COOLING = 0` en el archivo de parámetros. En esta configuración se generarán **50 salidas**, correspondientes a una salida cada **mil años** de evolución del sistema.


