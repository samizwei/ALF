!  Copyright (C) 2016 - 2020 The ALF project
!
!     The ALF project is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     The ALF project is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with ALF.  If not, see http://www.gnu.org/licenses/.
!
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version
!
!     MODIFIED: new model file (not part of the original ALF distribution).

!--------------------------------------------------------------------
!> @brief
!> Two-band model of Christensen et al., PRL 125, 247001 (2020), with a
!> selectable inter-orbital spin interaction:
!>
!>   H = H_0 - Ham_U * sum_i sum_{a in A} (S^a_i)^2 ,
!>   S^a_i = sum_{s,s'} c^dag_{is} sigma^a_{ss'} d_{is'} + h.c.
!>
!>   Ham_Int = "Ising"      : A = {z}        (paper, Eq. 3)
!>   Ham_Int = "XY"         : A = {x, y}
!>   Ham_Int = "Heisenberg" : A = {x, y, z}
!>
!> The same prefactor Ham_U multiplies every component.
!>
!> Representation: original basis (no gauge transformation), N_FL = 1, N_SUN = 1,
!> four orbitals per site:  1 = c_up, 2 = c_dn, 3 = d_up, 4 = d_dn.
!> Each (S^a_i)^2 is decoupled with its own discrete HS field (type 2),
!> g = sqrt(Dtau*Ham_U).
!>
!> Sign problem: the antiunitary T = (d -> -d) * (i sigma^y K) satisfies T^2 = -1
!> and commutes with H_0 and with every S^a_i (both factors flip S^a). Hence the
!> fermion determinant is real and non-negative for all three couplings (Kramers).
!--------------------------------------------------------------------

    submodule (Hamiltonian_main) ham_TwoBand_Spin_smod

      Use Operator_mod
      Use WaveFunction_mod
      Use Lattices_v3
      Use MyMats
      Use Random_Wrap
      Use Files_mod
      Use Matrix
      Use Observables
      Use Fields_mod
      Use Predefined_Hoppings
      use runtime_error_mod

      Implicit none

      type, extends(ham_base) :: ham_TwoBand_Spin
      contains
        ! Set Hamiltonian-specific procedures
        procedure, nopass :: Ham_Set
        procedure, nopass :: Alloc_obs
        procedure, nopass :: Obser
        procedure, nopass :: ObserT
#ifdef HDF5
        procedure, nopass :: write_parameters_hdf5
#endif
      end type ham_TwoBand_Spin

      !#PARAMETERS START# VAR_lattice
      Character (len=64) :: Model = ''  ! Value irrelevant
      Character (len=64) :: Lattice_type = 'Square'  ! Possible Values: 'Square'
      Integer            :: L1 = 6   ! Length in direction a_1
      Integer            :: L2 = 6   ! Length in direction a_2
      !#PARAMETERS END#

      !#PARAMETERS START# VAR_TwoBand_Spin
      Character (len=64)    :: Ham_Int   = 'Heisenberg' ! Interaction: 'Ising', 'XY' or 'Heisenberg'
      real(Kind=Kind(0.d0)) :: Ham_T     = 1.d0      ! Nearest-neighbour hopping t
      real(Kind=Kind(0.d0)) :: Ham_delta = 0.4d0     ! Hopping anisotropy delta
      real(Kind=Kind(0.d0)) :: Ham_mu    = -2.d0     ! Band offset mu: eps_c ~ -mu, eps_d ~ +mu (Eq. 1)
      real(Kind=Kind(0.d0)) :: Ham_chem  = 0.d0      ! Additional global chemical potential -Ham_chem*N
      real(Kind=Kind(0.d0)) :: Ham_U     = 1.d0      ! H_int = -U sum_i sum_a (S^a_i)^2
      real(Kind=Kind(0.d0)) :: Dtau      = 0.1d0     ! Thereby Ltrot=Beta/dtau
      real(Kind=Kind(0.d0)) :: Beta      = 5.d0      ! Inverse temperature
      !logical              :: Projector = .false.   ! Whether the projective algorithm is used
      real(Kind=Kind(0.d0)) :: Theta     = 10.d0     ! Projection parameter
      !logical              :: Symm      = .true.    ! Whether symmetrization takes place
      Integer               :: N_part    = -1        ! Total particles in trial wave function. If < 0 -> 2*L1*L2
      !#PARAMETERS END#

      Integer, parameter :: Norb_site = 4            ! c_up, c_dn, d_up, d_dn

      Type (Lattice),       target :: Latt
      Type (Unit_cell),     target :: Latt_unit      ! Four orbitals per unit cell
      Type (Unit_cell),     target :: Latt_unit_S    ! One "orbital" per unit cell, for S^a and pair correlations
      Integer, allocatable :: List(:,:), Invlist(:,:)
      Real (Kind=Kind(0.d0)), allocatable :: H0_mat(:,:)     ! Single-particle Hamiltonian (Ndim x Ndim)

      Integer :: N_comp                               ! Number of spin components in H_int
      Integer :: Comp(3)                              ! Which components (1=x, 2=y, 3=z)
      Complex (Kind=Kind(0.d0)) :: M_S(Norb_site,Norb_site,3)   ! On-site matrices of S^x, S^y, S^z
      Real    (Kind=Kind(0.d0)) :: S_pair(2)          ! Pair form factor: P = 2 (c_up c_dn - d_up d_dn)

    contains

      module Subroutine Ham_Alloc_TwoBand_Spin
        allocate(ham_TwoBand_Spin::ham)
      end Subroutine Ham_Alloc_TwoBand_Spin

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_TwoBand_Spin_read_write_parameters.F90"

!--------------------------------------------------------------------
      Subroutine Ham_Set

#if defined (MPI) || defined(TEMPERING)
          Use mpi
#endif
          Implicit none

          integer                :: nf, unit_info
          Character (len=64)     :: file_info

#ifdef MPI
          Integer        :: ierr, Isize, Irank, irank_g, isize_g, igroup
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
#endif

          call read_parameters()

          If (L1 < 2) then
             Write(error_unit,*) 'TwoBand_Spin: L1 has to be >= 2. For one-dimensional lattices set L2=1'
             Call Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          endif
          If (Ham_U < 0.d0) then
             Write(error_unit,*) 'TwoBand_Spin: Warning, Ham_U < 0 leads to a complex HS coupling (sign problem).'
          endif

          select case (trim(Ham_Int))
          case ('Ising', 'ising')
             N_comp = 1; Comp(1:1) = [3]
          case ('XY', 'xy')
             N_comp = 2; Comp(1:2) = [1, 2]
          case ('Heisenberg', 'heisenberg')
             N_comp = 3; Comp(1:3) = [1, 2, 3]
          case default
             Write(error_unit,*) 'TwoBand_Spin: Ham_Int must be Ising, XY or Heisenberg, got: ', trim(Ham_Int)
             Call Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          end select

          Call Set_spin_matrices

          Ltrot  = nint(beta/dtau)
          Thtrot = 0
          if (Projector) Thtrot = nint(theta/dtau)
          Ltrot  = Ltrot+2*Thtrot
          N_SUN  = 1
          N_FL   = 1

          Call Ham_Latt
          if (N_part < 0) N_part = 2*Latt%N     ! Total density n_c + n_d = 2 per site
          Call Ham_Hop
          Call Ham_V
          if (Projector) Call Ham_Trial()

#ifdef MPI
          If (Irank_g == 0) then
#endif
             File_info = "info"
#if defined(TEMPERING)
             write(File_info,'(A,I0,A)') "Temp_",igroup,"/info"
#endif
             Open(newunit=unit_info, file=file_info, status="unknown", position="append")
             Write(unit_info,*) '====================================='
             Write(unit_info,*) 'Model is      : TwoBand_Spin  (H_0 - U sum_i sum_a S^a_i S^a_i)'
             Write(unit_info,*) 'Interaction   : ', trim(Ham_Int)
             Write(unit_info,*) 'Lattice is    : ', Lattice_type
             Write(unit_info,*) 'L1            : ', L1
             Write(unit_info,*) 'L2            : ', L2
             Write(unit_info,*) '# unit cells  : ', Latt%N
             Write(unit_info,*) '# of orbitals : ', Latt_unit%Norb
             Write(unit_info,*) 'Ndim          : ', Ndim
             Write(unit_info,*) 'Symm. decomp  : ', Symm
             if (Projector) then
                Write(unit_info,*) 'Projective version'
                Write(unit_info,*) 'Theta         : ', Theta
                Write(unit_info,*) 'Tau_max       : ', beta
                Write(unit_info,*) '# of particles: ', N_part
             else
                Write(unit_info,*) 'Finite temperture version'
                Write(unit_info,*) 'Beta          : ', Beta
             endif
             Write(unit_info,*) 'dtau,Ltrot_eff: ', dtau,Ltrot
             Write(unit_info,*) 'N_SUN         : ', N_SUN
             Write(unit_info,*) 'N_FL          : ', N_FL
             Write(unit_info,*) 't             : ', Ham_T
             Write(unit_info,*) 'delta         : ', Ham_delta
             Write(unit_info,*) 'mu            : ', Ham_mu
             Write(unit_info,*) 'Ham_chem      : ', Ham_chem
             Write(unit_info,*) 'Ham_U         : ', Ham_U
             if (Projector) then
                Do nf = 1,N_FL
                   Write(unit_info,*) 'Degen of right trial wave function: ', WF_R(nf)%Degen
                   Write(unit_info,*) 'Degen of left  trial wave function: ', WF_L(nf)%Degen
                enddo
             endif
             Close(unit_info)
#ifdef MPI
          Endif
#endif

        end Subroutine Ham_Set

!--------------------------------------------------------------------
!> @brief
!> On-site matrices M^a with S^a_i = sum_{kl} x^dag_{ik} M^a_{kl} x_{il},
!> orbital order (c_up, c_dn, d_up, d_dn):  M^a = [[0, sigma^a], [sigma^a, 0]].
!--------------------------------------------------------------------
        Subroutine Set_spin_matrices

          Implicit none
          Complex (Kind=Kind(0.d0)) :: sig(2,2,3), Z0, Z1, ZI
          Integer :: a, s, sp

          Z0 = cmplx(0.d0, 0.d0, kind(0.d0))
          Z1 = cmplx(1.d0, 0.d0, kind(0.d0))
          ZI = cmplx(0.d0, 1.d0, kind(0.d0))
          sig = Z0
          sig(1,2,1) =  Z1 ; sig(2,1,1) =  Z1          ! sigma^x
          sig(1,2,2) = -ZI ; sig(2,1,2) =  ZI          ! sigma^y
          sig(1,1,3) =  Z1 ; sig(2,2,3) = -Z1          ! sigma^z

          M_S = Z0
          Do a = 1,3
             Do s = 1,2
                Do sp = 1,2
                   M_S(s,   2+sp, a) = sig(s,sp,a)     ! c^dag_s sigma d_s'
                   M_S(2+s, sp,   a) = sig(s,sp,a)     ! d^dag_s sigma c_s'  (h.c., sigma Hermitian)
                Enddo
             Enddo
          Enddo

          S_pair(1) =  1.d0                             ! c band
          S_pair(2) = -1.d0                             ! d band

        end Subroutine Set_spin_matrices

!--------------------------------------------------------------------
        Subroutine Ham_Latt

          Implicit none
          Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)
          Integer :: I, no, nc

          If (Lattice_type /= "Square") then
             Write(error_unit,*) 'The TwoBand_Spin model is only defined for the square lattice'
             Call Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          Endif

          a1_p(1) =  1.0d0  ; a1_p(2) =  0.d0
          a2_p(1) =  0.0d0  ; a2_p(2) =  1.d0
          L1_p    =  dble(L1)*a1_p
          L2_p    =  dble(L2)*a2_p
          Call Make_Lattice( L1_p, L2_p, a1_p, a2_p, Latt )

          Latt_unit%Norb    = Norb_site
          Latt_unit%N_coord = 2
          allocate(Latt_unit%Orb_pos_p(Latt_unit%Norb,2))
          Do no = 1,Norb_site
             Latt_unit%Orb_pos_p(no,:) = [0.d0, 0.d0]
          Enddo

          Latt_unit_S%Norb    = 1
          Latt_unit_S%N_coord = 2
          allocate(Latt_unit_S%Orb_pos_p(1,2))
          Latt_unit_S%Orb_pos_p(1,:) = [0.d0, 0.d0]

          Ndim = Latt%N*Latt_unit%Norb
          allocate(List(Ndim,2), Invlist(Latt%N,Latt_unit%Norb))
          nc = 0
          Do no = 1,Latt_unit%Norb
             Do I = 1,Latt%N
                nc = nc + 1
                List(nc,1)    = I
                List(nc,2)    = no
                Invlist(I,no) = nc
             Enddo
          Enddo

        end Subroutine Ham_Latt

!--------------------------------------------------------------------
!> @brief
!> Real-space hopping of Eq. (1), spin independent.
!> c-band (orbitals 1,2): t_x = t+delta, t_y = t-delta, on-site -mu
!> d-band (orbitals 3,4): t_x = t-delta, t_y = t+delta, on-site +mu
!--------------------------------------------------------------------
        Subroutine Set_H0_mat(H, Delta_mod)

          Implicit none
          Real (Kind=Kind(0.d0)), Intent(Out) :: H(:,:)
          Real (Kind=Kind(0.d0)), Intent(In)  :: Delta_mod   ! Small bond modulation (only for the trial wave function)

          Integer :: I, Ix, Iy, no, nb, I1, J1
          Real (Kind=Kind(0.d0)) :: tx(2), ty(2), eps(2), X

          tx (1) = Ham_T + Ham_delta ; ty (1) = Ham_T - Ham_delta ; eps(1) = -Ham_mu - Ham_chem
          tx (2) = Ham_T - Ham_delta ; ty (2) = Ham_T + Ham_delta ; eps(2) =  Ham_mu - Ham_chem

          H = 0.d0
          Do I = 1,Latt%N
             X  = 1.d0 + Delta_mod*cos(acos(-1.d0)*real(Latt%list(I,1) + Latt%list(I,2),Kind(0.d0)))
             Ix = Latt%nnlist(I,1,0)
             Iy = Latt%nnlist(I,0,1)
             Do no = 1,Norb_site
                nb = (no+1)/2                    ! band: 1 = c, 2 = d
                I1 = Invlist(I,no)
                ! Bonds are added (+=) such that L=2 correctly yields -2 t cos(k)
                J1 = Invlist(Ix,no)
                H(I1,J1) = H(I1,J1) - tx(nb)*X
                H(J1,I1) = H(J1,I1) - tx(nb)*X
                If (L2 > 1) then
                   J1 = Invlist(Iy,no)
                   H(I1,J1) = H(I1,J1) - ty(nb)
                   H(J1,I1) = H(J1,I1) - ty(nb)
                Endif
                H(I1,I1) = H(I1,I1) + eps(nb)
             Enddo
          Enddo

        end Subroutine Set_H0_mat

!--------------------------------------------------------------------
        Subroutine Ham_Hop

          Implicit none
          Integer :: nf, I, J

          allocate(H0_mat(Ndim,Ndim))
          Call Set_H0_mat(H0_mat, 0.d0)

          allocate(Op_T(1,N_FL))
          do nf = 1,N_FL
             Call Op_make(Op_T(1,nf),Ndim)
             Do I = 1,Ndim
                Do J = 1,Ndim
                   Op_T(1,nf)%O(I,J) = cmplx(H0_mat(I,J), 0.d0, kind(0.D0))
                Enddo
                Op_T(1,nf)%P(I) = I
             Enddo
             Op_T(1,nf)%g      = -Dtau
             Op_T(1,nf)%alpha  =  cmplx(0.d0,0.d0, kind(0.D0))
             Call Op_set(Op_T(1,nf))
          enddo

        end Subroutine Ham_Hop

!--------------------------------------------------------------------
!> @brief
!> Slater determinant ground state of H_0 (small bond modulation to lift
!> degeneracies), N_part particles in total.
!--------------------------------------------------------------------
        Subroutine Ham_Trial()

          Implicit none
          Integer                              :: nf, I, n
          Real (Kind=Kind(0.d0)), allocatable  :: H(:,:), U0(:,:), E0(:)
          Real (Kind=Kind(0.d0)), parameter    :: Delta_mod = 0.01d0

          If (N_part < 1 .or. N_part >= Ndim) then
             Write(error_unit,*) 'TwoBand_Spin: N_part has to satisfy 0 < N_part < Ndim'
             Call Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          Endif

          Allocate(WF_L(N_FL),WF_R(N_FL))
          do nf = 1,N_FL
             Call WF_alloc(WF_L(nf),Ndim,N_part)
             Call WF_alloc(WF_R(nf),Ndim,N_part)
          enddo

          Allocate(H(Ndim,Ndim), U0(Ndim,Ndim), E0(Ndim))
          Call Set_H0_mat(H, Delta_mod)
          Call Diag(H,U0,E0)

          Do nf = 1,N_FL
             do n = 1,N_part
                do I = 1,Ndim
                   WF_L(nf)%P(I,n) = U0(I,n)
                   WF_R(nf)%P(I,n) = U0(I,n)
                enddo
             enddo
             WF_L(nf)%Degen = E0(N_part+1) - E0(N_part)
             WF_R(nf)%Degen = E0(N_part+1) - E0(N_part)
          enddo

          Deallocate(H, U0, E0)

        end Subroutine Ham_Trial

!--------------------------------------------------------------------
!> @brief
!> -U (S^a_i)^2  ->  one type-2 operator per site and component,
!> exp(Dtau U (S^a)^2) ~ sum_l gamma_l exp( sqrt(Dtau U) eta_l S^a ).
!--------------------------------------------------------------------
        Subroutine Ham_V

          Implicit none
          Integer :: I, nc, n, k, l

          Allocate(Op_V(Latt%N*N_comp,1))
          n = 0
          Do nc = 1,N_comp
             Do I = 1,Latt%N
                n = n + 1
                Call Op_make(Op_V(n,1), Norb_site)
                Do k = 1,Norb_site
                   Op_V(n,1)%P(k) = Invlist(I,k)
                   Do l = 1,Norb_site
                      Op_V(n,1)%O(k,l) = M_S(k,l,Comp(nc))
                   Enddo
                Enddo
                Op_V(n,1)%g     = SQRT(CMPLX(Dtau*Ham_U, 0.d0, kind(0.D0)))
                Op_V(n,1)%alpha = cmplx(0.d0, 0.d0, kind(0.D0))
                Op_V(n,1)%type  = 2
                Call Op_set( Op_V(n,1) )
             Enddo
          Enddo

        end Subroutine Ham_V

!--------------------------------------------------------------------
        Subroutine Alloc_obs(Ltau)

          Implicit none
          Integer, Intent(In) :: Ltau
          Integer    ::  i, N, Nt
          Character (len=64) ::  Filename
          Character (len=:), allocatable ::  Channel

          ! Scalar observables
          Allocate ( Obs_scal(5) )
          Do I = 1,Size(Obs_scal,1)
             select case (I)
             case (1)
                N = 1;   Filename = "Kin"
             case (2)
                N = 1;   Filename = "Pot"
             case (3)
                N = 1;   Filename = "Part"
             case (4)
                N = 1;   Filename = "Ener"
             case (5)
                N = 2;   Filename = "Part_cd"
             case default
                Write(error_unit,*) ' Error in Alloc_obs '
             end select
             Call Obser_Vec_make(Obs_scal(I),N,Filename)
          enddo

          ! Equal time correlators
          Allocate ( Obs_eq(6) )
          Nt = 1
          Channel = '--'
          Do I = 1,Size(Obs_eq,1)
             select case (I)
             case (1)
                Filename = "Green"
                Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit,   Channel, dtau)
             case (2)
                Filename = "Den"
                Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit,   Channel, dtau)
             case (3)
                Filename = "SpinX"
                Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_S, Channel, dtau)
             case (4)
                Filename = "SpinY"
                Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_S, Channel, dtau)
             case (5)
                Filename = "SpinZ"
                Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_S, Channel, dtau)
             case (6)
                Filename = "Pair"
                Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_S, Channel, dtau)
             case default
                Write(error_unit,*) ' Error in Alloc_obs '
             end select
          enddo

          If (Ltau == 1) then
             Allocate ( Obs_tau(6) )
             Nt = Ltrot+1-2*Thtrot
             Do I = 1,Size(Obs_tau,1)
                select case (I)
                case (1)
                   Channel = 'P' ; Filename = "Green"
                   If (Projector) Channel = 'T0'
                   Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit,   Channel, dtau)
                case (2)
                   Channel = 'PH'; Filename = "Den"
                   If (Projector) Channel = 'T0'
                   Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit,   Channel, dtau)
                case (3)
                   Channel = 'PH'; Filename = "SpinX"
                   If (Projector) Channel = 'T0'
                   Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit_S, Channel, dtau)
                case (4)
                   Channel = 'PH'; Filename = "SpinY"
                   If (Projector) Channel = 'T0'
                   Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit_S, Channel, dtau)
                case (5)
                   Channel = 'PH'; Filename = "SpinZ"
                   If (Projector) Channel = 'T0'
                   Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit_S, Channel, dtau)
                case (6)
                   Channel = 'PP'; Filename = "Pair"
                   If (Projector) Channel = 'T0'
                   Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit_S, Channel, dtau)
                case default
                   Write(error_unit,*) ' Error in Alloc_obs '
                end select
             enddo
          endif

        End Subroutine Alloc_obs

!--------------------------------------------------------------------
!> @brief
!> Equal-time observables.  GR(i,j) = <x_i x^dag_j>,  GRC(i,j) = <x^dag_i x_j>.
!> Wick (single flavor):
!>  <A_I B_J> = sum M_ab N_cd [ GRC(a,b) GRC(c,d) + GRC(a,d) GR(b,c) ]
!>  <x^dag_1 x^dag_2 x_3 x_4> = GRC(1,4) GRC(2,3) - GRC(1,3) GRC(2,4)
!--------------------------------------------------------------------
        subroutine Obser(GR,Phase,Ntau, Mc_step_weight)

          Use Predefined_Obs

          Implicit none

          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
          Integer,                   INTENT(IN) :: Ntau
          Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight

          Complex (Kind=Kind(0.d0)), allocatable :: GRC(:,:,:)
          Complex (Kind=Kind(0.d0)) :: ZP, ZS, Zkin, Zpot, Zrho, Zc, Zd, ZSS, ZPair
          Complex (Kind=Kind(0.d0)) :: ZA(3)
          Integer :: I, J, I1, imj, nf, a, nc, no
          Integer :: IS(Norb_site), JS(Norb_site)

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          ZS = ZS*Mc_step_weight

          allocate(GRC(Ndim,Ndim,N_FL))
          Do nf = 1,N_FL
             Do I = 1,Ndim
                Do J = 1,Ndim
                   GRC(I, J, nf) = -GR(J, I, nf)
                Enddo
                GRC(I, I, nf) = 1.D0 + GRC(I, I, nf)
             Enddo
          Enddo

          Do I = 1,Size(Obs_scal,1)
             Obs_scal(I)%N         =  Obs_scal(I)%N + 1
             Obs_scal(I)%Ave_sign  =  Obs_scal(I)%Ave_sign + Real(ZS,kind(0.d0))
          Enddo

          ! <H_0>
          Zkin = cmplx(0.d0, 0.d0, kind(0.D0))
          Do I = 1,Ndim
             Do J = 1,Ndim
                If (H0_mat(I,J) /= 0.d0) Zkin = Zkin + H0_mat(I,J)*GRC(I,J,1)
             Enddo
          Enddo
          Obs_scal(1)%Obs_vec(1) = Obs_scal(1)%Obs_vec(1) + Zkin*ZP*ZS

          ! -U sum_i sum_a <(S^a_i)^2>
          Zpot = cmplx(0.d0, 0.d0, kind(0.D0))
          Do I = 1,Latt%N
             IS = Invlist(I,:)
             Do nc = 1,N_comp
                Zpot = Zpot + SS_eq(IS, IS, Comp(nc))
             Enddo
          Enddo
          Zpot = -Ham_U*Zpot
          Obs_scal(2)%Obs_vec(1) = Obs_scal(2)%Obs_vec(1) + Zpot*ZP*ZS

          ! Particle numbers
          Zc = cmplx(0.d0, 0.d0, kind(0.D0))
          Zd = cmplx(0.d0, 0.d0, kind(0.D0))
          Do I = 1,Latt%N
             Do no = 1,2
                Zc = Zc + GRC(Invlist(I,no),  Invlist(I,no),  1)
                Zd = Zd + GRC(Invlist(I,no+2),Invlist(I,no+2),1)
             Enddo
          Enddo
          Zrho = Zc + Zd
          Obs_scal(3)%Obs_vec(1) = Obs_scal(3)%Obs_vec(1) + Zrho*ZP*ZS
          Obs_scal(4)%Obs_vec(1) = Obs_scal(4)%Obs_vec(1) + (Zkin + Zpot)*ZP*ZS
          Obs_scal(5)%Obs_vec(1) = Obs_scal(5)%Obs_vec(1) + Zc*ZP*ZS
          Obs_scal(5)%Obs_vec(2) = Obs_scal(5)%Obs_vec(2) + Zd*ZP*ZS

          ! Equal-time correlations
          Call Predefined_Obs_eq_Green_measure( Latt, Latt_unit, List, GR, GRC, N_SUN, ZS, ZP, Obs_eq(1) )
          Call Predefined_Obs_eq_Den_measure  ( Latt, Latt_unit, List, GR, GRC, N_SUN, ZS, ZP, Obs_eq(2) )

          Do I1 = 3,6
             Obs_eq(I1)%N        = Obs_eq(I1)%N + 1
             Obs_eq(I1)%Ave_sign = Obs_eq(I1)%Ave_sign + Real(ZS,kind(0.d0))
          Enddo
          Do I = 1,Latt%N
             IS = Invlist(I,:)
             Do a = 1,3
                ZA(a) = Bil(IS, a)
                Obs_eq(2+a)%Obs_Latt0(1) = Obs_eq(2+a)%Obs_Latt0(1) + ZA(a)*ZP*ZS
             Enddo
             Do J = 1,Latt%N
                JS  = Invlist(J,:)
                imj = Latt%imj(I,J)
                Do a = 1,3
                   ZSS = SS_eq(IS, JS, a)
                   Obs_eq(2+a)%Obs_Latt(imj,1,1,1) = Obs_eq(2+a)%Obs_Latt(imj,1,1,1) + ZSS*ZP*ZS
                Enddo
                ZPair = Pair_eq(IS, JS)
                Obs_eq(6)%Obs_Latt(imj,1,1,1) = Obs_eq(6)%Obs_Latt(imj,1,1,1) + ZPair*ZP*ZS
             Enddo
          Enddo

          deallocate(GRC)

        contains

          !> <S^a_I> = sum_kl M_kl <x^dag_k x_l>
          Complex (Kind=Kind(0.d0)) function Bil(IS, a)
            Integer, Intent(In) :: IS(Norb_site), a
            Integer :: k, l
            Bil = cmplx(0.d0, 0.d0, kind(0.D0))
            Do k = 1,Norb_site
               Do l = 1,Norb_site
                  If (M_S(k,l,a) /= cmplx(0.d0,0.d0,kind(0.d0))) &
                       & Bil = Bil + M_S(k,l,a)*GRC(IS(k),IS(l),1)
               Enddo
            Enddo
          end function Bil

          !> <S^a_I S^a_J>
          Complex (Kind=Kind(0.d0)) function SS_eq(IS, JS, a)
            Integer, Intent(In) :: IS(Norb_site), JS(Norb_site), a
            Integer :: k, l, m, n
            Complex (Kind=Kind(0.d0)) :: Zcon
            Zcon = cmplx(0.d0, 0.d0, kind(0.D0))
            Do k = 1,Norb_site
               Do l = 1,Norb_site
                  If (M_S(k,l,a) == cmplx(0.d0,0.d0,kind(0.d0))) cycle
                  Do m = 1,Norb_site
                     Do n = 1,Norb_site
                        If (M_S(m,n,a) == cmplx(0.d0,0.d0,kind(0.d0))) cycle
                        Zcon = Zcon + M_S(k,l,a)*M_S(m,n,a)*GRC(IS(k),JS(n),1)*GR(IS(l),JS(m),1)
                     Enddo
                  Enddo
               Enddo
            Enddo
            SS_eq = Bil(IS,a)*Bil(JS,a) + Zcon
          end function SS_eq

          !> <P^dag_I P_J>,  P_I = 2 sum_o s_o x_{I o up} x_{I o dn},  s_c = +1, s_d = -1
          Complex (Kind=Kind(0.d0)) function Pair_eq(IS, JS)
            Integer, Intent(In) :: IS(Norb_site), JS(Norb_site)
            Integer :: o, op, au, ad, bu, bd
            Pair_eq = cmplx(0.d0, 0.d0, kind(0.D0))
            Do o = 1,2
               au = IS(2*o-1); ad = IS(2*o)
               Do op = 1,2
                  bu = JS(2*op-1); bd = JS(2*op)
                  ! <x^dag_{a dn} x^dag_{a up} x_{b up} x_{b dn}>
                  Pair_eq = Pair_eq + S_pair(o)*S_pair(op) * &
                       & ( GRC(ad,bd,1)*GRC(au,bu,1) - GRC(ad,bu,1)*GRC(au,bd,1) )
               Enddo
            Enddo
            Pair_eq = 4.d0*Pair_eq
          end function Pair_eq

        end Subroutine Obser

!--------------------------------------------------------------------
!> @brief
!> Time-displaced observables.
!>  <x^dag_a(tau) x_b(tau)> = delta_ab - GTT(b,a),  <x^dag_a(0) x_b(0)> = delta_ab - G00(b,a)
!>  <x^dag_a(tau) x_d(0)>   = -G0T(d,a),            <x_b(tau) x^dag_c(0)> = GT0(b,c)
!--------------------------------------------------------------------
        Subroutine ObserT(NT,  GT0,G0T,G00,GTT, PHASE,  Mc_step_weight)

          Use Predefined_Obs

          Implicit none

          Integer         , INTENT(IN) :: NT
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase
          Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight

          Complex (Kind=Kind(0.d0)) :: ZP, ZS, ZSS, ZPair, ZAI, ZAJ, Zcon
          Integer :: I, J, I1, imj, a, k, l, m, n, o, op, au, ad, bu, bd
          Integer :: IS(Norb_site), JS(Norb_site)

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          ZS = ZS * Mc_step_weight

          Call Predefined_Obs_tau_Green_measure( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT, N_SUN, ZS, ZP, Obs_tau(1) )
          Call Predefined_Obs_tau_Den_measure  ( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT, N_SUN, ZS, ZP, Obs_tau(2) )

          If (NT == 0) then
             Do I1 = 3,6
                Obs_tau(I1)%N        = Obs_tau(I1)%N + 1
                Obs_tau(I1)%Ave_sign = Obs_tau(I1)%Ave_sign + Real(ZS,kind(0.d0))
             Enddo
          Endif

          Do I = 1,Latt%N
             IS = Invlist(I,:)
             Do J = 1,Latt%N
                JS  = Invlist(J,:)
                imj = Latt%imj(I,J)
                Do a = 1,3
                   ZAI  = cmplx(0.d0, 0.d0, kind(0.D0))
                   ZAJ  = cmplx(0.d0, 0.d0, kind(0.D0))
                   Zcon = cmplx(0.d0, 0.d0, kind(0.D0))
                   Do k = 1,Norb_site
                      Do l = 1,Norb_site
                         If (M_S(k,l,a) == cmplx(0.d0,0.d0,kind(0.d0))) cycle
                         ZAI = ZAI - M_S(k,l,a)*GTT(IS(l),IS(k),1)
                         ZAJ = ZAJ - M_S(k,l,a)*G00(JS(l),JS(k),1)
                         Do m = 1,Norb_site
                            Do n = 1,Norb_site
                               If (M_S(m,n,a) == cmplx(0.d0,0.d0,kind(0.d0))) cycle
                               Zcon = Zcon - M_S(k,l,a)*M_S(m,n,a)*G0T(JS(n),IS(k),1)*GT0(IS(l),JS(m),1)
                            Enddo
                         Enddo
                      Enddo
                   Enddo
                   ZSS = ZAI*ZAJ + Zcon
                   Obs_tau(2+a)%Obs_Latt(imj,NT+1,1,1) = Obs_tau(2+a)%Obs_Latt(imj,NT+1,1,1) + ZSS*ZP*ZS
                   If (J == 1) Obs_tau(2+a)%Obs_Latt0(1) = Obs_tau(2+a)%Obs_Latt0(1) + ZAI*ZP*ZS
                Enddo

                ! <P^dag_I(tau) P_J(0)>, <x^dag_1(tau) x^dag_2(tau) x_3(0) x_4(0)>
                !   = C(1,4) C(2,3) - C(1,3) C(2,4),  C(p,q) = -G0T(q,p)
                ZPair = cmplx(0.d0, 0.d0, kind(0.D0))
                Do o = 1,2
                   au = IS(2*o-1); ad = IS(2*o)
                   Do op = 1,2
                      bu = JS(2*op-1); bd = JS(2*op)
                      ZPair = ZPair + S_pair(o)*S_pair(op) * &
                           & ( G0T(bd,ad,1)*G0T(bu,au,1) - G0T(bu,ad,1)*G0T(bd,au,1) )
                   Enddo
                Enddo
                ZPair = 4.d0*ZPair
                Obs_tau(6)%Obs_Latt(imj,NT+1,1,1) = Obs_tau(6)%Obs_Latt(imj,NT+1,1,1) + ZPair*ZP*ZS
             Enddo
          Enddo

        end Subroutine OBSERT

    end submodule ham_TwoBand_Spin_smod
