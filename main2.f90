MODULE SUBS

CONTAINS

SUBROUTINE INITIALIZE(RHO, U, V, RHO_STAR, U_STAR, V_STAR, RHO_1, U_1, V_1)
    IMPLICIT NONE
    
    DOUBLE PRECISION, DIMENSION(:,:) :: RHO, U, V, RHO_STAR, U_STAR, V_STAR, RHO_1, U_1, V_1
    INTEGER :: NX, NY
    
    NX = SIZE(U,2)
    NY = SIZE(U,1)
    
    RHO = 1000
    U = 0
    V = 0
    
    RHO_STAR = 1000
    U_STAR = 0
    V_STAR = 0
    
    RHO_1 = 1000
    U_1 = 0
    V_1 = 0

END SUBROUTINE INITIALIZE

SUBROUTINE APPLY_VELOCITY_BC(RHO, U, V, RHO_STAR, U_STAR, V_STAR, RHO_1, U_1, V_1, U0)
    IMPLICIT NONE
    
    DOUBLE PRECISION, DIMENSION(:,:) :: RHO, U, V, RHO_STAR, U_STAR, V_STAR, RHO_1, U_1, V_1
    INTEGER :: NX, NY
    DOUBLE PRECISION :: U0
    
    NX = SIZE(U,2)
    NY = SIZE(U,1)
    
    U(:,1) = 0      
    U(:,NX) = 0
    U(1,:) = 0
    U(NY,:) = U0
    
    !U_STAR(:,1) = 0      
    !U_STAR(:,NX) = 0
    !U_STAR(1,:) = 0
    !U_STAR(NY,:) = U0
    
    U_1(:,1) = 0      
    U_1(:,NX) = 0
    U_1(1,:) = 0
    U_1(NY,:) = U0
    
    V(:,1) = 0      
    V(:,NX) = 0
    V(1,:) = 0
    V(NY,:) = 0
    
    !V_STAR(:,1) = 0      
    !V_STAR(:,NX) = 0
    !V_STAR(1,:) = 0
    !V_STAR(NY,:) = 0
    
    V_1(:,1) = 0      
    V_1(:,NX) = 0
    V_1(1,:) = 0
    V_1(NY,:) = 0
    
END SUBROUTINE APPLY_VELOCITY_BC


SUBROUTINE PREDICTOR(RHO, U, V, RHO_STAR, U_STAR, V_STAR, RHO_1, U_1, V_1, DX, DY, DT, MU, C, U0)
    IMPLICIT NONE
    
    DOUBLE PRECISION, DIMENSION(:,:) :: RHO, U, V, RHO_STAR, U_STAR, V_STAR, RHO_1, U_1, V_1
    INTEGER :: NX, NY
    
    DOUBLE PRECISION :: DX
    DOUBLE PRECISION :: DY
    DOUBLE PRECISION :: DT
    DOUBLE PRECISION :: C
    DOUBLE PRECISION :: MU
    DOUBLE PRECISION :: U0
    
    INTEGER :: ROW, COL
    DOUBLE PRECISION :: C1, C2, C3, C4, C5
    
    NX = SIZE(U,2)
    NY = SIZE(U,1)


    
    C1 = DT/DX
    C2 = DT/DY
    C3 = MU*DT/DX**2
    C4 = MU*DT/DY**2
    C5 = MU*DT/(12*DX*DY)
    
    ! RHO: internal points
    DO ROW = 2,NY-1
        DO COL = 2,NX-1
            RHO_STAR(ROW,COL) = RHO(ROW,COL) - C1*(RHO(ROW,COL+1)*U(ROW,COL+1) - RHO(ROW,COL)*U(ROW,COL)) &
                                     - C2*(RHO(ROW+1,COL)*V(ROW+1,COL) - RHO(ROW,COL)*V(ROW,COL))
        END DO    
    END DO  

    ! RHO: left surface
    COL = 1
    DO ROW = 1,NY
        RHO_STAR(ROW,COL) = RHO(ROW,COL) - ( DT/(2*DX) ) * ( -RHO(ROW,COL+2)*U(ROW,COL+2) + 4*RHO(ROW,COL+1)*U(ROW,COL+1) &
                                                             - 3*RHO(ROW,COL)*U(ROW,COL))
    END DO
    
    ! RHO: right surface
    COL = NX
    DO ROW = 1,NY
        RHO_STAR(ROW,COL) = RHO(ROW,COL) + ( DT/(2*DX) ) * ( -RHO(ROW,COL-2)*U(ROW,COL-2) + 4*RHO(ROW,COL-1)*U(ROW,COL-1) &
                                                             - 3*RHO(ROW,COL)*U(ROW,COL))
    END DO
    
    ! RHO: bottom surface
    ROW = 1
    DO COL = 1,NX
        RHO_STAR(ROW,COL) = RHO(ROW,COL) - ( DT/(2*DY) ) * ( -RHO(ROW+2,COL)*V(ROW+2,COL) + 4*RHO(ROW+1,COL)*V(ROW+1,COL) &
                                                             - 3*RHO(ROW,COL)*V(ROW,COL))
    END DO
    
    ! ROW: top surface
    ROW = NY
    DO COL = 2,NX-1
        RHO_STAR(ROW,COL) = RHO(ROW,COL) - ( DT*U0/(2*DX) ) * ( RHO(ROW,COL+1) - RHO(ROW,COL-1) ) &
                                         + ( DT/(2*DY) ) * ( -RHO(ROW-2,COL)*V(ROW-2,COL) + 4*RHO(ROW-1,COL)*V(ROW-1,COL) &
                                                             - 3*RHO(ROW,COL)*V(ROW,COL))
    END DO
    
    ! U, V: internal points
    DO ROW = 2,NY-1
        DO COL = 2,NX-1      
            U_STAR(ROW,COL) = ( RHO(ROW,COL)*U(ROW,COL) &
                                - C1*(RHO(ROW,COL+1)*(U(ROW,COL+1)**2+C**2) - RHO(ROW,COL)*(U(ROW,COL)**2+C**2)) &
                                - C2*(RHO(ROW+1,COL)*U(ROW+1,COL)*V(ROW+1,COL) - RHO(ROW,COL)*U(ROW,COL)*V(ROW,COL)) &
                                + (4.0/3.0)*C3*(U(ROW,COL+1)-2*U(ROW,COL)+U(ROW,COL-1)) &
                                + C4*(U(ROW+1,COL)-2*U(ROW,COL)+U(ROW-1,COL)) &
                                + C5*(V(ROW+1,COL+1)+V(ROW-1,COL-1)-V(ROW-1,COL+1)-V(ROW+1,COL-1))  )/RHO_STAR(ROW,COL)
            
            V_STAR(ROW,COL) = ( RHO(ROW,COL)*V(ROW,COL) &
                                - C1*(RHO(ROW,COL+1)*U(ROW,COL+1)*V(ROW,COL+1) - RHO(ROW,COL)*U(ROW,COL)*V(ROW,COL)) &
                                - C2*(RHO(ROW+1,COL)*(V(ROW+1,COL)**2+C**2) - RHO(ROW,COL)*(V(ROW,COL)**2+C**2)) &
                                + C3*(V(ROW,COL+1)-2*V(ROW,COL)+V(ROW,COL-1)) &
                                + (4.0/3.0)*C4*(V(ROW+1,COL)-2*V(ROW,COL)+V(ROW-1,COL)) &
                                + C5*(U(ROW+1,COL+1)+U(ROW-1,COL-1)-U(ROW-1,COL+1)-U(ROW+1,COL-1))  )/RHO_STAR(ROW,COL)
        END DO
    END DO
    
END SUBROUTINE PREDICTOR




SUBROUTINE CORRECTOR(RHO, U, V, RHO_STAR, U_STAR, V_STAR, RHO_1, U_1, V_1, DX, DY, DT, MU, C, U0)
    IMPLICIT NONE
    
    DOUBLE PRECISION, DIMENSION(:,:) :: RHO, U, V, RHO_STAR, U_STAR, V_STAR, RHO_1, U_1, V_1
    INTEGER :: NX, NY
    
    DOUBLE PRECISION :: DX
    DOUBLE PRECISION :: DY
    DOUBLE PRECISION :: DT
    DOUBLE PRECISION :: C
    DOUBLE PRECISION :: MU
    DOUBLE PRECISION :: U0
    
    INTEGER :: ROW, COL
    DOUBLE PRECISION :: C1, C2, C3, C4, C5
    
    NX = SIZE(U,2)
    NY = SIZE(U,1)
 
    C1 = DT/DX
    C2 = DT/DY
    C3 = MU*DT/DX**2
    C4 = MU*DT/DY**2
    C5 = MU*DT/(12*DX*DY)
    
    
    ! RHO: internal points
    DO ROW = 2,NY-1
        DO COL = 2,NX-1
            RHO_1(ROW,COL) = ( RHO(ROW,COL) + RHO_STAR(ROW,COL) &
                            - C1*(RHO_STAR(ROW,COL)*U_STAR(ROW,COL) - RHO_STAR(ROW,COL-1)*U_STAR(ROW,COL-1)) &
                            - C2*(RHO_STAR(ROW,COL)*V_STAR(ROW,COL) - RHO_STAR(ROW-1,COL)*V_STAR(ROW-1,COL)) )/2
        END DO    
    END DO  

    ! RHO: left surface
    COL = 1
    DO ROW = 1,NY
        RHO_1(ROW,COL) = 0.5 * ( RHO(ROW,COL) + RHO_STAR(ROW,COL) &
                               - ( DT/(2*DX) ) * ( -RHO_STAR(ROW,COL+2)*U_STAR(ROW,COL+2)  &
                                                   + 4*RHO_STAR(ROW,COL+1)*U_STAR(ROW,COL+1) &
                                                   - 3*RHO_STAR(ROW,COL)*U_STAR(ROW,COL)) )
    END DO
    
    ! RHO: right surface
    COL = NX
    DO ROW = 1,NY
        RHO_1(ROW,COL) = 0.5 * ( RHO(ROW,COL) + RHO_STAR(ROW,COL) &
                               + ( DT/(2*DX) ) * ( -RHO_STAR(ROW,COL-2)*U_STAR(ROW,COL-2)  &
                                                   + 4*RHO_STAR(ROW,COL-1)*U_STAR(ROW,COL-1) &
                                                   - 3*RHO_STAR(ROW,COL)*U_STAR(ROW,COL)) )
    END DO
    
    ! RHO: bottom surface
    ROW = 1
    DO COL = 1,NX
        RHO_1(ROW,COL) = 0.5 * ( RHO(ROW,COL) + RHO_STAR(ROW,COL) &
                               - ( DT/(2*DY) ) * ( -RHO_STAR(ROW+2,COL)*V_STAR(ROW+2,COL)  &
                                                   + 4*RHO_STAR(ROW+1,COL)*V_STAR(ROW+1,COL) &
                                                   - 3*RHO_STAR(ROW,COL)*V_STAR(ROW,COL)) )
    END DO
    
    ! ROW: top surface
    ROW = NY
    DO COL = 2,NX-1
        RHO_1(ROW,COL) = 0.5 * ( RHO(ROW,COL) + RHO_STAR(ROW,COL) &
                               - ( DT*U0/(2*DX) ) * ( RHO_STAR(ROW,COL+1) - RHO_STAR(ROW,COL-1) ) &
                               + ( DT/(2*DY) ) * ( -RHO_STAR(ROW-2,COL)*V_STAR(ROW-2,COL)  &
                                                   + 4*RHO_STAR(ROW-1,COL)*V_STAR(ROW-1,COL) &
                                                   - 3*RHO_STAR(ROW,COL)*V_STAR(ROW,COL)) )
    END DO
    
    ! U, V: internal points
    DO ROW = 2,NY-1
        DO COL = 2,NX-1      
            U_1(ROW,COL) = ( RHO(ROW,COL)*U(ROW,COL) + RHO_STAR(ROW,COL)*U_STAR(ROW,COL) &
                                - C1*(RHO_STAR(ROW,COL)*(U_STAR(ROW,COL)**2+C**2) &
                                            - RHO_STAR(ROW,COL-1)*(U_STAR(ROW,COL-1)**2+C**2)) &
                                - C2*(RHO_STAR(ROW,COL)*U_STAR(ROW,COL)*V_STAR(ROW,COL) &
                                            - RHO_STAR(ROW-1,COL)*U_STAR(ROW-1,COL)*V_STAR(ROW-1,COL)) &
                                + (4.0/3.0)*C3*(U_STAR(ROW,COL+1)-2*U_STAR(ROW,COL)+U_STAR(ROW,COL-1)) &
                                + C4*(U_STAR(ROW+1,COL)-2*U_STAR(ROW,COL)+U_STAR(ROW-1,COL)) &
                                + C5*(V_STAR(ROW+1,COL+1)+V_STAR(ROW-1,COL-1)-V_STAR(ROW-1,COL+1)-V_STAR(ROW+1,COL-1))  &
                                )/ (2*RHO_1(ROW,COL))
            
            V_1(ROW,COL) = ( RHO(ROW,COL)*V(ROW,COL) + RHO_STAR(ROW,COL)*V_STAR(ROW,COL) &
                                - C1*(RHO_STAR(ROW,COL)*U_STAR(ROW,COL)*V_STAR(ROW,COL) &
                                            - RHO_STAR(ROW,COL-1)*U_STAR(ROW,COL-1)*V_STAR(ROW,COL-1)) &
                                - C2*(RHO_STAR(ROW,COL)*(V_STAR(ROW,COL)**2+C**2) &
                                            - RHO_STAR(ROW-1,COL)*(V_STAR(ROW-1,COL)**2+C**2)) &
                                + C3*(V_STAR(ROW,COL+1)-2*V_STAR(ROW,COL)+V_STAR(ROW,COL-1)) &
                                + (4.0/3.0)*C4*(V_STAR(ROW+1,COL)-2*V_STAR(ROW,COL)+V_STAR(ROW-1,COL)) &
                                + C5*(U_STAR(ROW+1,COL+1)+U_STAR(ROW-1,COL-1)-U_STAR(ROW-1,COL+1)-U_STAR(ROW+1,COL-1))  &
                                )/ (2*RHO_1(ROW,COL))
        END DO
    END DO
    
END SUBROUTINE CORRECTOR

END MODULE SUBS




PROGRAM NS
    USE commonUtil
    USE SUBS
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NX = 16
    INTEGER, PARAMETER :: NY = 16
    DOUBLE PRECISION, PARAMETER :: DX = 1
    DOUBLE PRECISION, PARAMETER :: DY = 1
    DOUBLE PRECISION, PARAMETER :: DT = 1E-4
    DOUBLE PRECISION, PARAMETER :: C = 1500    ! speed of sound
    DOUBLE PRECISION, PARAMETER :: MU = 1E-3   ! viscosity
    
    DOUBLE PRECISION, PARAMETER :: U0 = 1E-5
    
    DOUBLE PRECISION, DIMENSION(NY,NX) :: RHO, U, V, RHO_STAR, U_STAR, V_STAR, RHO_1, U_1, V_1
    
    
    INTEGER :: I
    INTEGER :: ROW,COL
    CHARACTER(LEN=10) :: INT_STR, STR
    CHARACTER(LEN=10) :: ARG
    INTEGER :: NT
    
    
    CALL INITIALIZE(RHO, U, V, RHO_STAR, U_STAR, V_STAR, RHO_1, U_1, V_1)
    
    CALL APPLY_VELOCITY_BC(RHO, U, V, RHO_STAR, U_STAR, V_STAR, RHO_1, U_1, V_1, U0)
    
    CALL GETARG(1,ARG)
    READ(ARG,*) NT
    
    DO I=1,NT
        CALL PREDICTOR(RHO, U, V, RHO_STAR, U_STAR, V_STAR, RHO_1, U_1, V_1, DX, DY, DT, MU, C, U0)       ! CALCULATE U_STAR
        CALL APPLY_VELOCITY_BC(RHO, U, V, RHO_STAR, U_STAR, V_STAR, RHO_1, U_1, V_1, U0)
        CALL CORRECTOR(RHO, U, V, RHO_STAR, U_STAR, V_STAR, RHO_1, U_1, V_1, DX, DY, DT, MU, C, U0)       ! CALCULATE U_1
        CALL APPLY_VELOCITY_BC(RHO, U, V, RHO_STAR, U_STAR, V_STAR, RHO_1, U_1, V_1, U0)
        
        ! REPLACE U WITH U_1
        RHO = RHO_1
        U = U_1
        V = V_1
        
        WRITE(INT_STR,'(I)') I
        STR = 'U_'//INT_STR//'.dat'
        !CALL write2DArrayToFile(U, STR)
        
    END DO
    
    
    CALL write2DArrayToFile(U, 'U.dat')
    CALL write2DArrayToFile(V, 'V.dat')
    CALL write2DArrayToFile(RHO, 'RHO.dat')
    
    
END PROGRAM



