! ==============================================================================
! Euler2D.f95
!
! por Edgar Muñoz Sánchez
! 4 / dic / 2023


! Programa que resuelve las ecuaciones de Euler de la dinámica de gases en 2D
! usando los métodos de Lax y Macormack.

! Esta versión utiliza un archivo externo en el que se define el problema a
! resolver, definiendo ahí constantes globales, condiciones iniciales y
! condiciones de frontera (si necesarias).
!
! ==============================================================================

! ==============================================================================
 include 'Supernova.f95'

! ==============================================================================
! VARIABLES GLOBALES
! ============================================================================

module globals
  
  use parameters

  real*8    U(NEQ, 0:NX+1, 0:NY+1)    ! Variables conservadas actuales
  real*8   UP(NEQ, 0:NX+1, 0:NY+1)    ! Variables conservadas "avanzadas"
  real*8    F(NEQ, 0:NX+1, 0:NY+1)    ! Flujos físicos
  real*8    G(NEQ, 0:NX+1, 0:NY+1)    ! Flujos físicos
  real*8   UT(NEQ, 0:NX+1, 0:NY+1)    ! Us temporales (sólo usadas en Macormack)
  real*8 PRIM(NEQ, 0:NX+1, 0:NY+1)    ! Variables primitivas

  real*8  :: time      ! Tiempo actual
  integer :: it        ! Iteración actual
  real*8  :: dt        ! Paso de tiempo
  integer :: nout      ! Número de la siguiente salida
  real*8  :: tout      ! Tiempo de la siguiente salida

  ! Para medir tiempo de ejecución
  integer :: clock_start, clock_count, clock_rate, clock_max

end module globals

! ==============================================================================
! SUBRUTINAS
! ==============================================================================

! Devuelve la coordenada X de la celda i
subroutine xcoord(i, x)
  use parameters, only: X1, DX
  implicit none
  integer, intent(in) :: i
  real*8, intent(out) :: x
  x = X1 + i*DX
end subroutine xcoord

! Devuelve la coordenada Y de la celda j
subroutine ycoord(j, y)
  use parameters, only: Y1, DY
  implicit none
  integer, intent(in) :: j
  real*8, intent(out) :: y
  y = Y1 + j*DY
end subroutine ycoord

! ==============================================================================
! INICIALIZACIONES
! ==============================================================================

! Inicializaciones genéricas
subroutine initmain()

  use globals, only: time, it, nout, tout

  ! Inicializaciones de variables globales del código
  time = 0.0
  it = 0
  nout = 0
  tout = 0.0

end subroutine initmain

! ==============================================================================
! CONDICIONES INICIALES
! ==============================================================================

! Impone las condiciones iniciales
subroutine initflow(UU)
  
  use parameters
  implicit none
  real*8, intent(inout) :: UU(NEQ, 0:NX+1, 0:NY+1)

  ! Inicializar variables conservadas U de acuerdo al problema elegido
  ! Esta función debe estar definida en el archivo del problema
  call IC_custom(UU)

end subroutine initflow

! ==============================================================================
! CONDICIONES DE FRONTERA
! ==============================================================================

! Aplicar condiciones de frontera a celdas fantasma
! El arreglo pasado es al que aplicaremos las BCs
subroutine boundary(UU)
  
  use parameters
  use globals, only: time
  implicit none
  real*8, intent(inout) :: UU(NEQ, 0:NX+1, 0:NY+1)

  integer :: i, j, e

  ! BC a la izquierda
  if (BC_LEFT == BC_FREEFLOW) then
    do j=1,NY
      do e=1,NEQ
        UU(e,0,j) = UU(e,1,j)
      end do
    end do
  else if (BC_LEFT == BC_INFLOW) then
    ! Especificar valores aquí
    do j=1,NY
      UU(1,0,j) = 0.0
      UU(2,0,j) = 0.0
      UU(3,0,j) = 0.0
      UU(4,0,j) = 0.0
    end do
  else if (BC_LEFT == BC_REFLECTIVE) then
    do j=1,NY
      UU(1,0,j) =  UU(1,1,j)
      UU(2,0,j) = -UU(2,1,j)
      UU(3,0,j) =  UU(3,1,j)
      UU(4,0,j) =  UU(4,1,j)
    end do
  else if (BC_LEFT == BC_PERIODIC) then
    do j=1,NY
      do e=1,NEQ
        UU(e,0,j) = UU(e,NX,j)
      end do
    end do
  end if

  ! BC a la derecha
  if (BC_RIGHT == BC_FREEFLOW) then
    do j=1,NY
      do e=1,NEQ
        UU(e,NX+1,j) = UU(e,NX,j)
      end do
    end do
  else if (BC_RIGHT == BC_INFLOW) then
    do j=1,NY
      UU(1,NX+1,j) = 0.0
      UU(2,NX+1,j) = 0.0
      UU(3,NX+1,j) = 0.0
      UU(4,NX+1,j) = 0.0
    end do
  else if (BC_RIGHT == BC_REFLECTIVE) then
    do j=1,NY
      UU(1,NX+1,j) =  UU(1,NX,j)
      UU(2,NX+1,j) = -UU(2,NX,j)
      UU(3,NX+1,j) =  UU(3,NX,j)
      UU(4,NX+1,j) =  UU(4,NX,j)
    end do
  else if (BC_RIGHT == BC_PERIODIC) then
    do j=1,NY
      do e=1,NEQ
        UU(e,NX+1,j) = UU(e,1,j)
      end do
    end do
  end if

  ! BC inferior
  if (BC_BOTTOM == BC_FREEFLOW) then
    do i=1,NX
      do e=1,NEQ
        UU(e,i,0) = UU(e,i,1)
      end do
    end do
  else if (BC_BOTTOM == BC_INFLOW) then
    ! Especificar valores aquí
    do i=1,NX
      UU(1,i,0) = 0.0
      UU(2,i,0) = 0.0
      UU(3,i,0) = 0.0
      UU(4,i,0) = 0.0
    end do
  else if (BC_BOTTOM == BC_REFLECTIVE) then
    do i=1,NX
      UU(1,i,0) =  UU(1,i,1)
      UU(2,i,0) =  UU(2,i,1)
      UU(3,i,0) = -UU(3,i,1)
      UU(4,i,0) =  UU(4,i,1)
    end do
  else if (BC_BOTTOM == BC_PERIODIC) then
    do i=1,NX
      do e=1,NEQ
        UU(e,i,0) = UU(e,i,NY)
      end do
    end do
  end if

  ! BC superior
  if (BC_TOP == BC_FREEFLOW) then
    do i=1,NX
      do e=1,NEQ
        UU(e,i,NY+1) = UU(e,i,NY)
      end do
    end do
  else if (BC_TOP == BC_INFLOW) then
    do i=1,NX
      UU(1,i,NY+1) = 0.0
      UU(2,i,NY+1) = 0.0
      UU(3,i,NY+1) = 0.0
      UU(4,i,NY+1) = 0.0
    end do
  else if (BC_TOP == BC_REFLECTIVE) then
    do i=1,NX
      UU(1,i,NY+1) =  UU(1,i,NY)
      UU(2,i,NY+1) =  UU(2,i,NY)
      UU(3,i,NY+1) = -UU(3,i,NY)
      UU(4,i,NY+1) =  UU(4,i,NY)
    end do
  else if (BC_TOP == BC_PERIODIC) then
    do i=1,NX
      do e=1,NEQ
        UU(e,i,NY+1) = UU(e,i,1)
      end do
    end do
  end if

  ! Aplicar las condiciones de frontera especificados por el usuario
  ! en el archivo de parémetros importado
  call BC_custom(UU, time)

end subroutine boundary

! ==============================================================================
! FLUJOS FÍSICOS Y CONVERSIONES DE PRIMITIVAS
! ==============================================================================

! Calcula las primitivas en toda la malla a partir de las U pasadas
! Esto incluye las celdas fantasma
subroutine update_prim(UU, PP)
  
  use parameters
  implicit none
  real*8, intent(inout) :: UU(NEQ, 0:NX+1, 0:NY+1)
  real*8, intent(out) :: PP(NEQ, 0:NX+1, 0:NY+1)

  integer :: i, j

  do j=0,NY+1
    do i=0,NX+1
      PP(1,i,j) = UU(1,i,j)
      PP(2,i,j) = UU(2,i,j)/UU(1,i,j)
      PP(3,i,j) = UU(3,i,j)/UU(1,i,j)
      PP(4,i,j) = (GAM-1.0)*(UU(4,i,j) - 0.5*(UU(2,i,j)**2 + UU(3,i,j)**2)/UU(1,i,j))
      if (PP(4,i,j) <= 0) then
        PP(4,i,j) = 1D-30
        UU(4,i,j) = 0.5*PP(1,i,j)*(PP(2,i,j)**2+PP(3,i,j)**2) + PP(4,i,j)/(GAM-1.0);
      end if
    end do
  end do

end subroutine update_prim

! ==============================================================================

! Calcular los flujos físicos F -- Ecuaciones de Euler 1D
! No olvidar calcular F en las celdas fantasma!
subroutine Euler_fluxes(PP, FF, GG)
  
  use parameters
  implicit none
  real*8, intent(in) :: PP(NEQ, 0:NX+1, 0:NY+1)
  real*8, intent(out) :: FF(NEQ, 0:NX+1, 0:NY+1)
  real*8, intent(out) :: GG(NEQ, 0:NX+1, 0:NY+1)

  real*8 :: E
  integer :: i, j

  do j=0,NY+1
    do i=0,NX+1
      E = 0.5*PP(1,i,j)*(PP(2,i,j)**2 + PP(3,i,j)**2) + PP(4,i,j)/(GAM-1)
      ! Flujos en X
      FF(1,i,j) = PP(1,i,j) * PP(2,i,j)
      FF(2,i,j) = PP(1,i,j) * PP(2,i,j)**2 + PP(4,i,j)
      FF(3,i,j) = PP(1,i,j) * PP(2,i,j) * PP(3,i,j)
      FF(4,i,j) = PP(2,i,j) * (E + PP(4,i,j))
      ! Flujos en Y
      GG(1,i,j) = PP(1,i,j) * PP(3,i,j)
      GG(2,i,j) = PP(1,i,j) * PP(2,i,j) * PP(3,i,j)
      GG(3,i,j) = PP(1,i,j) * PP(3,i,j)**2 + PP(4,i,j)
      GG(4,i,j) = PP(3,i,j) * (E + PP(4,i,j))
    end do
  end do

end subroutine Euler_fluxes

! ==============================================================================
! PASO DE TIMEPO
! ==============================================================================

! Calcula el paso de tiempo resultante de la condición CFL
subroutine timestep(PP, dt)
  
  use parameters
  use globals, only: time, tout
  implicit none
  real*8, intent(in) :: PP(NEQ, 0:NX+1, 0:NY+1)
  real*8, intent(out) :: dt

  real*8 :: cs, u_plus_cs, v_plus_cs, max_speed_x, max_speed_y, dtx, dty
  integer :: i, j

  ! Determinamos el máximo valor de u + cs en ambas dimensiones
  max_speed_x = 0.0
  max_speed_y = 0.0
  do i=1,NX
    do j=1,NY
      cs = sqrt(GAM*PP(4,i,j)/PP(1,i,j))
      u_plus_cs = abs(PP(2,i,j)) + cs
      v_plus_cs = abs(PP(3,i,j)) + cs
      if (u_plus_cs > max_speed_x) then
        max_speed_x = u_plus_cs
      end if
      if (v_plus_cs > max_speed_y) then
        max_speed_y = v_plus_cs
      end if
    end do
  end do

  ! Pasos de tiempo usando el criterio CFL, en cada dimensión
  dtx = CFL * DX / max_speed_x
  dty = CFL * DY / max_speed_y

  ! Usamos el menor de los pasos de tiempo
  ! El factor de 1/sqrt(2) es necesario en 2D para que la condición de
  ! estabilidad siga siendo CFL < 1
  dt = min(dtx, dty) / sqrt(2.0)

  if (time + dt > tout) then
    dt = tout - time
  end if

end subroutine timestep

! ==============================================================================
! SOLVERS NUMÉRICOS
! ==============================================================================

! Método de Lax-Friedrichs
subroutine Lax(UU, PP, FF, GG, UUP)
  
  use parameters
  use globals, only: dt
  implicit none
  real*8, intent(in) :: UU(NEQ, 0:NX+1, 0:NY+1)
  real*8, intent(in) :: PP(NEQ, 0:NX+1, 0:NY+1)
  real*8, intent(inout) :: FF(NEQ, 0:NX+1, 0:NY+1)
  real*8, intent(inout) :: GG(NEQ, 0:NX+1, 0:NY+1)
  real*8, intent(out) :: UUP(NEQ, 0:NX+1, 0:NY+1)

  integer :: i, j, e

  ! Calcular los flujos físicos (ecs Euler)
  call Euler_fluxes(PP, FF, GG)

  ! Aplicar método de Lax
  do j=1,NY
    do i=1,NX
      do e=1,NEQ
        UUP(e,i,j) = &
          0.25*(UU(e,i+1,j) + UU(e,i-1,j) + UU(e,i,j+1) + UU(e,i,j-1)) &
          - dt/(2*DX) * (FF(e,i+1,j) - FF(e,i-1,j)) &
          - dt/(2*DY) * (GG(e,i,j+1) - GG(e,i,j-1))
      end do
    end do
  end do

end subroutine Lax

! ==============================================================================

! Método de Macormack
subroutine Macormack(UU, PP, FF, GG, UUP)
  
  use parameters
  use globals, only: dt, UT
  implicit none
  real*8, intent(in) :: UU(NEQ, 0:NX+1, 0:NY+1)
  real*8, intent(inout) :: PP(NEQ, 0:NX+1, 0:NY+1)
  real*8, intent(inout) :: FF(NEQ, 0:NX+1, 0:NY+1)
  real*8, intent(inout) :: GG(NEQ, 0:NX+1, 0:NY+1)
  real*8, intent(out) :: UUP(NEQ, 0:NX+1, 0:NY+1)

  integer :: i, j, e

  ! Calcular los flujos físicos (ecs Euler)
  call Euler_fluxes(PP, FF, GG)

  ! Paso predictor: actualizamos las UT con flujos hacia adelante
  do j=1,NY
    do i=1,NX
      do e=1,NEQ
        UT(e,i,j) = UU(e,i,j) &
          - dt/DX * (FF(e,i+1,j) - FF(e,i,j)) &
          - dt/DY * (GG(e,i,j+1) - GG(e,i,j))
      end do
    end do
  end do

  ! Aplicamos las BCs a las UT
  call boundary(UT)

  ! Actualizar las primitivas de nuevo usando las UT esta vez
  call update_prim(UT, PP)

  ! Re-calculamos los flujos F y G pero usando las primitivas actualizadas
  call Euler_fluxes(PP, FF, GG)

  ! Paso corrector: obtenemos las UP usando U, UT y F actualizados
  do j=1,NY
    do i=1,NX
      do e=1,NEQ
        UUP(e,i,j) = 0.5*(UU(e,i,j) + UT(e,i,j)) &
          - 0.5*dt/DX * (FF(e,i,j) - FF(e,i-1,j)) &
          - 0.5*dt/DY * (GG(e,i,j) - GG(e,i,j-1))
      end do
    end do
  end do

end subroutine Macormack

! ==============================================================================

! Volcar las UPs sobre las Us "finales" del paso de tiempo -- sin viscosidad
subroutine step(UU, UUP)
  
  use parameters
  implicit none
  real*8, intent(out) :: UU(NEQ, 0:NX+1, 0:NY+1)
  real*8, intent(in) :: UUP(NEQ, 0:NX+1, 0:NY+1)

  integer :: i, j, e

  do j=0,NY+1
    do i=0,NX+1
      do e=1,NEQ
        UU(e,i,j) = UUP(e,i,j)
      end do
    end do
  end do

end subroutine step

! ==============================================================================

! Volcar las UPs sobre las Us "finales" del paso de tiempo, aplicando
! viscosidad artificial donde haya máximos o mínimos locales
subroutine step_viscosity(UU, UUP)
  
  use parameters
  implicit none
  real*8, intent(out) :: UU(NEQ, 0:NX+1, 0:NY+1)
  real*8, intent(in) :: UUP(NEQ, 0:NX+1, 0:NY+1)

  integer :: i, j, e

  do j=1,NY
    do i=1,NX
      do e=1,NEQ
        UU(e,i,j) = UUP(e,i,j) &
         + ETA*(UUP(e,i+1,j) + UUP(e,i-1,j) + UUP(e,i,j+1) + UUP(e,i,j-1) &
         - 4*UUP(e,i,j))
        ! Versión selectiva
        ! Aplicamos la viscosidad sólo donde hay mínimos/máximos locales
        ! En las demás celdas simplemente copiamos las UP sobre las U
        ! ETA debe ser estrictamente menor que 1/2
        ! if (((UUP(e,i+1,j)-UUP(e,i,j))*(UUP(e,i,j)-UU(e,i-1,j)) < 0) .or. &
        ! ((UUP(e,i,j+1)-UU(e,i,j))*(UUP(e,i,j)-UU(e,i,j-1)) < 0)) then
        !   UU(e,i,j) = UUP(e,i,j) &
        !     + ETA*(UUP(e,i+1,j) + UUP(e,i-1,j) + UUP(e,i,j+1) + UUP(e,i,j-1) &
        !     - 4*UUP(e,i,j))
        ! else
        !   UU(e,i,j) = UUP(e,i,j)
        ! end if
      end do
    end do
  end do

  ! Lo de arriba no toca las celdas fantasma, pero éstas también deben
  ! ser copiadas (sin aplicar viscosidad)
  do j=1,NY
    do e=1,NEQ
      UU(e,0,j) = UUP(e,0,j)
      UU(e,NX+1,j) = UUP(e,NX+1,j)
    end do
  end do
  do i=1,NX
    do e=1,NEQ
      UU(e,i,0) = UUP(e,i,0)
      UU(e,i,NY+1) = UUP(e,i,NY+1)
    end do
  end do

end subroutine step_viscosity

! =============================================================================
! ENFRIAMIENTO RADIATIVO
! =============================================================================

! Calcula la temperatura del gas a partir de su presión y densidad
subroutine temperature(pres, dens, T)

  use parameters
  implicit none
  real*8, intent(in) :: pres
  real*8, intent(in) :: dens
  real*8, intent(out) :: T

  T = mu0*AMU/KB*pres/dens
  if (T > 1e4) then
    T = mui*AMU/KB*pres/dens
  end if

end subroutine temperature

! -----------------------------------------------------------------------------

! La función de enfriamiento Λ(T) da la tasa de emisión de radiación, como
! función de la temperatura, de un plasma a temperaturas arriba de 10^4 K.
!
! Esto resulta principalmente de la suma de líneas atómicas de emisión
! (principalmente hidrógeno, carbono, oxígeno, nitrógeno, neón y hierro) y
! emisión libre-libre (bremsstrahlung) que domina a altas temperaturas.
!
! Las unidades son erg cm^3 s^-1. Multiplicar esta cantidad por n_e*n_ion
! resulta en la tasa de pérdida de energía, en erg cm^-3 s^-1.
!
! La función implementada aquí está dada por interpolaciones simples:
! T < 1e4 K: no hay enfriamiento
! 1e4 <= T <= 1e5: rampa lineal desde Lambda = 1e-25 hasta 1e-21
! T >= 1e5: fórmula de Draine & Wood (1991): líneas + libre-libre
!
! Una función más realista usa coeficientes obtenidos de cálculos muy
! detallados de procesos atómicos de emisión en plasmas, y en general
! varía apreciablemente según de la metalicidad del gas y puede depender
! de su densidad. Ver Dalgarno & McCray (1972), Sutherland & Dopita (1993),
! y la base de datos CHIANTI (https://www.chiantidatabase.org/).
subroutine Lambda(T, L)

  use parameters
  implicit none
  real*8, intent(in) :: T
  real*8, intent(out) :: L

  real*8, parameter :: T0 = 1.0D4, T1 = 1.0D5, T2 = 1.0D7
  real*8, parameter :: L0 = 1D-25, L1 = 1D-21, L2 = 1D-23
  real*8, parameter :: m = (L1-L0)/(T1-T0)
  
  if (T < T0) then
    L = 0.0;
  else if ((T >= T0) .and. (T <= T1)) then
    L = L0 + m*(T-T0)
  else
    L = L2*((T2/T)+sqrt(T/T2))
  end if

end subroutine Lambda

! -----------------------------------------------------------------------------

! Aplica enfriamiento radiativo a todas las celdas
! Esto se hace calculando la cantidad de energía térmica perdida por
! radiación en el intervalo de tiempo dt y reduciendo la presión concorde
! Además, para evitar errores de presiones negativas, el gas no se dejará
! enfriar a menos de 10^3 K.
subroutine cooling(UU, PP)

  use parameters
  use globals, only: dt
  implicit none
  real*8, intent(inout) :: UU(NEQ, 0:NX+1, 0:NY+1)
  real*8, intent(inout) :: PP(NEQ, 0:NX+1, 0:NY+1)

  real*8, parameter :: Tmin = 1000.0
  real*8 :: T, n, tau, L, Ekin, Eth0, cool_factor
  integer :: i, j

  ! Celdas físicas
  do i=1,NX
    do j=1,NY

      ! Energía cinética y térmica de la celda
      Eth0 = PP(4,i,j)/(GAM-1)
      Ekin = UU(4,i,j) - Eth0

      ! Temperatura de la celda
      call temperature(PP(4,i,j), PP(1,i,j), T)

      ! Función de enfriamiento
      call Lambda(T, L)

      if (L > 0) then
     
        ! Densidad numérica de la celda (suponiendo ionización completa)
        n = PP(1,i,j)/(mui*AMU)

        ! Eth0 / n^2 Λ(T) es la escala de tiempo de enfriamiento
        tau = Eth0 / (n**2 * L)

        ! Factor de enfriamiento
        cool_factor = exp(-dt/tau)

        ! Limitamos el factor de enfriamiento para minimizar errores
        if (cool_factor < 0.5) then
          cool_factor = 0.5;
          !printf("Cooling limit reached!\n")
          !printf("\nCell %i, %i\n", i, j)
          !printf("rho = %e\n", PRIM[i][j][0])
          !printf("n = %e\n", n)
          !printf("P = %e\n", PRIM[i][j][3])
          !printf("T = %E\n", T)
          !printf("Pnew = %e\n", PRIM[i][j][3]*cool_factor)
          !printf("cool_factor = %e\n", cool_factor)
          !printf("Eth_old = %E\n", Eth0)
        end if
        
        ! Actualizar la presión y la energía total de la celda
        PP(4,i,j) = PP(4,i,j) * cool_factor
        UU(4,i,j) = Ekin + PP(4,i,j)/(GAM-1)

      end if
    
    end do
  end do

end subroutine cooling
  
! =============================================================================
! OTROS
! =============================================================================

! Agrega los términos geométricos cilíndricos para simulaciones axisimétricas
!
! Se considera que los ejes x y y corresponden realmente a los ejes z y r
! de un sistema de coordenadas cilíndrico 3D, y que el problema tiene
! simetría en torno al eje z. Así podemos simular un corte 2D de un problema
! que realmente es 3D (aprovechando esa simetría).
!
! Escribiendo las ecs de Euler 3D en coordenadas cilíndricas salen términos 
! adicionales (llamados "geométricos") que no contienen derivadas, así que
! son movidos al lado derecho de las ecuaciones y se consideran términos
! fuente. Se suman, multiplicados por dt, a las conservadas al final de cada
! paso.
subroutine add_cyl_geom_terms(UU)

  use parameters
  use globals, only: dt
  implicit none
  real*8, intent(inout) :: UU(NEQ, 0:NX+1, 0:NY+1)

  real*8 :: z, r, rho, u, v, P, E
  integer :: i, j

  do i=1,NX
    call xcoord(i, z)
    do j=1,NY
      call ycoord(j, r)
      rho = UU(1,i,j)
      u = UU(2,i,j)/UU(1,i,j)
      v = UU(3,i,j)/UU(1,i,j)
      P = (GAM-1)*(UU(4,i,j) - 0.5*rho*(u**2 + v**2))
      E = UU(4,i,j)
      UU(1,i,j) = UU(1,i,j) + dt * (-rho*v/r)
      UU(2,i,j) = UU(2,i,j) + dt * (-rho*u*v/r)
      UU(3,i,j) = UU(3,i,j) + dt * (-rho*v*v/r)
      UU(4,i,j) = UU(4,i,j) + dt * (-v*(E+P)/r)
    end do
  end do

end subroutine add_cyl_geom_terms

! ==============================================================================
! OUTPUT DE DATOS
! ==============================================================================
  
! Escribe a disco el estado de la simulación
subroutine output(PRIM)
  
  use parameters
  use globals, only: time, it, dt, nout, tout
  implicit none
  real*8, intent(in) :: PRIM(NEQ, 0:NX+1, 0:NY+1)

  character(len=256) :: fname
  real*8 :: x, y, tsc
  integer :: i, j, e

  if (OUTPUT_FORMAT == OUTPUT_FORMAT_ASCII) then

    ! Generar el nombre del archivo de salida
    write(fname, "(a,a,i4.4,a)") TRIM(OUT_DIR), "output_", nout, ".txt"

    ! Abrir el archivo
    open(unit=10, file=fname, status='unknown')

    ! Escribir la posición, densidad, velocidad y presión a disco,
    ! separadas por espacios
    do i=1,NX
      call xcoord(i, x)
      do j=1,NY
        call ycoord(j, y)
        write(10,*) x, " ", y, " ", PRIM(1,i,j), " ", PRIM(2,i,j), " ", &
          PRIM(3,i,j), " ", PRIM(4,i,j)
      end do
    end do

    ! Cerrar archivo
    close(10)

  else if (OUTPUT_FORMAT == OUTPUT_FORMAT_BINARY) then

    ! Generar el nombre del archivo de salida
    write(fname, "(a,a,i4.4,a)") trim(OUT_DIR), "output_", nout, ".dat"

    ! Abrir el archivo
    open(unit=10, file=fname, access="stream")

    ! Escribir tiempo actual de simulación
    write(10) time

    ! Escribir NEQ, NX y NY
    write(10) NEQ, NX, NY

    ! Escribir X1, X2, Y1, Y2
    write(10) X1, X2, Y1, Y2

    ! Escribir la posición, densidad, velocidad y presión a disco,
    ! separadas por espacios
    do e=1,NEQ
      do i=1,NX
        do j=1,NY
          write(10) prim(e,i,j)
        end do
      end do
    end do

    ! Cerrar archivo
    close(10)

  end if

  if (time > YR) then
    tsc = YR
  else
    tsc = 1.0
  end if
  write(*,'(a,i0,a,f12.6,a,i0,a,es10.3)') "Salida ", nout, ", t=", time/tsc, ", it=", it, ", dt=", dt/tsc, fname

  ! Avanzar variables para el siguiente output
  nout = nout + 1
  tout = nout * DTOUT

end subroutine output

! ==============================================================================
! PROGRAMA PRINCIPAL
! ==============================================================================

program Euler1D
  
  use parameters
  use globals
  implicit none

  ! Inicializaciones generales
  call initmain()

  ! Condición inicial e inicializaciones
  call initflow(U)

  ! Opcional: Agregar términos fuente geométricos (coords cilíndricas)
  if (AXISYMMETRIC) then
    call add_cyl_geom_terms(U)    
  end if

  ! Llenar celdas fantasma
  call boundary(U)

  ! Calcular las primitivas
  call update_prim(U, PRIM)

  ! Escribir condición inicial a disco
  call output(PRIM)

  ! Tiempo de inicio de la simulación
  call system_clock(clock_start, clock_rate, clock_max)
  
  ! Bucle principal
  do while (time <= TFIN)

    ! Calcular el paso de tiempo
    call timestep(PRIM, dt)

    ! Aplicar el método numérico para calcular las UP
    if (SOLVER == SOLVER_LAX) then
      call Lax(U, PRIM, F, G, UP)
    else if (SOLVER == SOLVER_MACORMACK) then
      call Macormack(U, PRIM, F, G, UP)
    end if

    ! Opcional: Agregar términos fuente geométricos axisimétricos
    if (AXISYMMETRIC) then
      call add_cyl_geom_terms(UP)
    end if

    ! Opcional: Aplicar enfriamiento radiativo, si aplicable
    if (USE_COOLING) then
      call cooling(UP, PRIM)
    end if

    ! Aplicar condiciones de frontera a las UP recién calculadas
    call boundary(UP)

    ! Avanzar las U, con viscosidad para Macormack
    if (SOLVER == SOLVER_LAX) then
      call step(U, UP)
    else if (SOLVER == SOLVER_MACORMACK) then
      call step_viscosity(U, UP)
    end if

    ! Actualizar las primitivas PRIM usando las nuevas U
    call update_prim(U, PRIM)

    ! Avanzar el estado de la simulación
    time = time + dt
    it = it + 1

    ! Escribir a disco
    if (time >= tout) then
      call output(PRIM)
    end if

  end do

  ! Imprimir tiempo transcurrido
  call system_clock(clock_count, clock_rate, clock_max)
  write(*,'(a,i5,a,f10.5,a)') "Se calcularon ", it, " iteraciones en ", (clock_count-clock_start)/real(clock_rate), " s"

end program Euler1D

! ==============================================================================
