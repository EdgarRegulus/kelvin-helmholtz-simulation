! ==============================================================================
! CONSTANTES GLOBALES
! ==============================================================================
module constants
  implicit none

  ! Constantes con nombre (a usar como opciones más abajo) -- No modificar

  ! Solvers numéricos
  integer, parameter :: SOLVER_LAX = 1
  integer, parameter :: SOLVER_MACORMACK = 2

  ! Tipos de condiciones de frontera
  integer, parameter ::  BC_FREEFLOW = 1
  integer, parameter ::  BC_INFLOW = 2
  integer, parameter ::  BC_REFLECTIVE = 3
  integer, parameter ::  BC_PERIODIC = 4

  ! Tipos de formato de salida
  integer, parameter :: OUTPUT_FORMAT_ASCII = 1
  integer, parameter :: OUTPUT_FORMAT_BINARY = 2

  ! ----------------------------------------------------------------------------

  ! Constantes astrofísicas y fundamentales
  real*8, parameter :: AMU  = 1.660538782e-24    ! Unidad de masa atómica (g)
  real*8, parameter :: KB   = 1.380650400e-16    ! Constante de Boltzmann (erg/K)
  real*8, parameter :: PC   = 3.085677588e+18    ! Parsec (cm)
  real*8, parameter :: AU   = 1.495978707e+13    ! Unidad astronómica (cm)
  real*8, parameter :: YR   = 3.155673600e+7     ! Año sideral terrestre (s)
  real*8, parameter :: KYR  = 3.155673600e+10    ! Mil años en segundos (s)
  real*8, parameter :: MSUN = 1.988920000e+33    ! Masa solar (g)
  real*8, parameter :: KPS  = 1.0e5              ! Un km/s en cm/s
  
  real*8, parameter :: PI = acos(-1.0)  

end module constants

! ==============================================================================
! DEFINICIÓN DEL PROBLEMA
! ==============================================================================
module parameters
  
  use constants
  implicit none

  ! Parámetros de la simulación
  integer, parameter :: NEQ = 4         ! Número de ecuaciones
  integer, parameter :: NX = 500        ! Tamaño de la malla
  integer, parameter :: NY = 250        ! Tamaño de la malla
  real*8, parameter :: X1 = -30*PC      ! Coordenada física del extremo izquierdo
  real*8, parameter :: X2 = 30*PC       ! Coordenada física del extremo derecho
  real*8, parameter :: Y1 = 0.0         ! Coordenada física del extremo izquierdo
  real*8, parameter :: Y2 = 30*PC       ! Coordenada física del extremo derecho
  real*8, parameter :: TFIN = 50e3*YR   ! Tiempo final de integración
  real*8, parameter :: DTOUT = 1000*YR  ! Intervalo para escribir a disco
  
  real*8, parameter :: CFL = 0.9        ! Parametro de Courant
  real*8, parameter :: GAM = 5.0/3.0    ! Razón de capacidades caloríficas

  ! Viscosidad artficial, sólo para Macormack
  real*8, parameter :: ETA = 0.05

  ! Solucionador numérico
  ! integer, parameter :: SOLVER = SOLVER_LAX
  integer, parameter ::  SOLVER = SOLVER_MACORMACK

  ! Condiciones de frontera, usar las constantes siguientes:
  ! BC_FREEFLOW: salida libre (gradiente cero)
  ! BC_INFLOW: entrada debe ser especificada en boundary()
  ! BC_REFLECTIVE: reflectiva
  ! BC_PERIODIC: periódica
  integer, parameter :: BC_LEFT   = BC_FREEFLOW
  integer, parameter :: BC_RIGHT  = BC_FREEFLOW
  integer, parameter :: BC_BOTTOM = BC_REFLECTIVE
  integer, parameter :: BC_TOP    = BC_FREEFLOW

  ! Directorio donde escribir las salidas (usar "./" para dir actual)
  ! Debe terminar en una diagonal '/'
  character(len=128), parameter :: OUT_DIR = "output_cooling/"

  ! Formato de salidas: ASCII o binario
  ! integer, parameter :: OUTPUT_FORMAT = OUTPUT_FORMAT_ASCII
  integer, parameter :: OUTPUT_FORMAT = OUTPUT_FORMAT_BINARY

  ! Enfriamiento radiativo del gas
  logical, parameter :: USE_COOLING = .true.

  ! Axisimetría en torno al eje X
  logical, parameter :: AXISYMMETRIC = .true.

  ! Constantes calculadas
  real*8, parameter :: DX = (X2-X1)/NX      ! Espaciamiento en X de la malla
  real*8, parameter :: DY = (Y2-Y1)/NY      ! Espaciamiento en Y de la malla

  ! ============================================================================
  ! PARÄMETROS DEL PROBLEMA: REMANENTES DE SUPERNOVA

  ! Parámetros del medio ambiente
  real*8, parameter :: n_ism = 1.0    ! Densidad numérica (cm^-3)
  real*8, parameter :: mu0 = 1.3      ! Masa por partícula, gas neutro (amu)
  real*8, parameter :: mui = 0.61     ! Masa por partícula, gas ionizado (amu)
  real*8, parameter :: T_ism = 100    ! Temperatura (K)
  real*8, parameter :: u_ism = 0.0    ! Velocidad en x (cm/s)
  real*8, parameter :: v_ism = 0.0    ! Velocidad en y (cm/s)

  ! Parámetros de la supernova
  real*8, parameter :: RSN = 1*PC        ! Radio inicial de la explosión (cm)
  real*8, parameter :: ESN = 1.0D51      ! Energía de la explosión (erg)
  real*8, parameter :: MSN = 1.0*MSUN    ! Masa eyectada (g)
  real*8, parameter :: SN_xc = 0.0     ! Coord x del centro de la explosión (cm)
  real*8, parameter :: SN_yc = 0.0     ! Coord y del centro de la explosión (cm)

contains

  ! ============================================================================
  ! CONDICIÓN INICIAL
  ! ============================================================================
  
  ! Condición inicial definida por el usuario -- A MODIFICAR
  ! La función debe fijar el valor de las variables conservadas en todas
  ! las celdas físicas (no necesario llenar celdas fantasma)
  ! Se puede usar x(i) = X1 + i*DX y y(i) = Y1 + j*DY para obtener
  ! las coords físicas de una celda
  subroutine IC_custom(UU) 

    implicit none
    real*8, intent(inout) :: UU(NEQ, 0:NX+1, 0:NY+1)

    ! Valores calculados del medio ambiente
    real*8, parameter :: rho_ism = n_ism * mu0 * AMU
    real*8, parameter :: P_ism = n_ism*KB*T_ism
    real*8, parameter :: E_ism = 0.5*rho_ism*(u_ism**2+v_ism**2) + P_ism/(GAM-1)

    ! Parámetros calculados del remanente
    real*8, parameter :: chi = 0.5              ! Fracción de energía cinética a total
    real*8, parameter :: Ekin = chi*ESN         ! Energía cinética
    real*8, parameter :: Eth = (1-chi)*ESN      ! Energía térmica
    real*8, parameter :: rho_SN = MSN/(4.0*PI/3.0*RSN**3)   ! Densidad interior
    real*8, parameter :: vmax = sqrt(10.0/3.0*Ekin/MSN)       ! Velocidad en el borde
    real*8, parameter :: pres = (GAM-1)*Eth/(4.0*PI/3.0*RSN**3)   ! Presión interior
    real*8, parameter :: rho = rho_SN + rho_ism
  
    real*8 :: dist, vel_mag, u, v, x, y
    integer :: i, j

    write(*,*) "Parámetros del medio interestelar"
    write(*,*) "rho_ism = ",  rho_ism, " g/cm^3"
    write(*,*) "P_ism = ", P_ism, " dyn/cm^2"
    write(*,*) "T_ism = ", T_ism, " K"

    write(*,*) "Parámetros de la supernova"
    write(*,*) "RSN = ", RSN/PC, "pc"
    write(*,*) "MSN = ", MSN/MSUN, " Msun"
    write(*,*) "ESN = ", ESN, " erg"
    write(*,*) "Ekin = ", Ekin, " erg"
    write(*,*) "Eth = ", Eth, " erg"
    write(*,*) "xc = ", SN_xc/PC, " pc, yc = ", SN_yc/PC, " pc"
    write(*,*) "rho_SN = ", rho_SN, " g/cm^3"
    write(*,*) "vmax = ", vmax/KPS, " km/s"
    write(*,*) "pres = ", pres, " dyn/cm^2"

    do i=1, NX
      x = X1 + i*DX
      do j=1, NY
        y = Y1 + j*DY

        ! Distancia al centro de la explosión
        dist = sqrt((x-SN_xc)**2 + (y-SN_yc)**2)

        if (dist <= RSN) then

          ! Magnitud y componentes de la velocidad
          if (dist == 0) then
            vel_mag = 0
            u = 0
            v = 0
          else
            vel_mag = (dist/RSN)*vmax
            u = vel_mag*(x-SN_xc)/dist
            v = vel_mag*(y-SN_yc)/dist
          end if

          ! Llenar conservadas
          UU(1,i,j) = rho
          UU(2,i,j) = rho * u
          UU(3,i,j) = rho * v
          UU(4,i,j) = 0.5*rho*vel_mag**2 + pres/(GAM-1)

        else

          ! Llenar conservadas
          UU(1,i,j) = rho_ism
          UU(2,i,j) = rho_ism * u_ism
          UU(3,i,j) = rho_ism * v_ism
          UU(4,i,j) = E_ism

        end if

      end do
    end do

  end subroutine IC_custom

  ! ============================================================================
  ! CONDICIÓN DE FRONTERA "CUSTOM"
  ! ============================================================================

  subroutine BC_custom(UU, time) 

    implicit none
    real*8, intent(inout) :: UU(NEQ, 0:NX+1, 0:NY+1)
    real*8, intent(in) :: time

  end subroutine BC_custom

end module parameters
