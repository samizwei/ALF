!  Copyright (C) 2022 The ALF project
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
!     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
!
!     - It is our hope that this program makes a contribution to the scientific community. Being
!       part of that community we feel that it is reasonable to require you to give an attribution
!       back to the original authors if you have benefitted from this program.
!       Guidelines for a proper citation can be found on the project's homepage
!       http://alf.physik.uni-wuerzburg.de
!
!     - We require the preservation of the above copyright notice and this license in all original files.
!
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
!
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version


!--------------------------------------------------------------------
!> @author
!> ALF-project (model file written for the two-band S^z S^z model)
!>
!> @brief
!> Two-band Hubbard model with interband S^z S^z interaction, Eq. (3) of
!> M. H. Christensen et al., Phys. Rev. Lett. 125, 247001 (2020):
!>
!> \f[ H = H_0 - U \sum_i S^z_i S^z_i , \qquad
!>     S^z_i = \sum_{\alpha\beta} c^\dagger_{i\alpha} \sigma^z_{\alpha\beta} d_{i\beta} + h.c. \f]
!>
!> with the real-space form of Eq. (1)  (see Infos/Fourier_transform_H0.pdf)
!>
!> \f[ H_0 = - \sum_{i\sigma} \left[ (t+\delta) c^\dagger_{i\sigma} c_{i+x,\sigma}
!>                                 + (t-\delta) c^\dagger_{i\sigma} c_{i+y,\sigma}
!>                                 + (t-\delta) d^\dagger_{i\sigma} d_{i+x,\sigma}
!>                                 + (t+\delta) d^\dagger_{i\sigma} d_{i+y,\sigma} + h.c. \right]
!>           - \mu \sum_{i\sigma} \left( n^c_{i\sigma} - n^d_{i\sigma} \right)
!>           - \mu_0 \sum_{i\sigma} \left( n^c_{i\sigma} + n^d_{i\sigma} \right) \f]
!>
!> The last term (Ham_chem = mu_0) is an optional global chemical potential, zero in the paper.
!>
!> @details
!> Implementation details:
!>  - Square lattice, two orbitals per unit cell: orbital 1 = c, orbital 2 = d.
!>  - S^z_i conserves spin, S^z_i = A_{i,up} - A_{i,dn} with A_{i,s} = c^dag_{is} d_{is} + h.c.
!>    The gauge transformation  d_{i,dn} -> - d_{i,dn}  leaves H_0 invariant and maps
!>    A_{i,dn} -> -A_{i,dn}.  In the new basis
!>                 H = H_0 - U \sum_i [ \sum_s ( c^dag_{is} d_{is} + h.c. ) ]^2
!>    which is SU(2) symmetric in the spin ("color") index. We therefore use
!>    N_FL = 1, N_SUN = 2: the weight is det(M)^2 with real M -> no sign problem for U > 0.
!>  - The interaction is a perfect square of a single-body operator (type=2 vertex,
!>    Predefined_Int_V_SUN) with g = sqrt(Dtau*U), alpha = 0.
!>  - Observables are measured in the gauge-transformed basis. Gauge invariant quantities
!>    (densities, S^z, the pair field P_+- below) are identical to the original basis.
!>    Exception: the c-d off-diagonal elements of "Green" equal sum_s s <c^dag_s d_s>
!>    in the original basis (spin-resolved sign), the c-c and d-d elements are unchanged.
!>
!> Observables:
!>  Scalar   : Kin  = <H_0> (incl. on-site terms), Pot = -U sum_i <(S^z_i)^2>,
!>             Part = <N>, Ener = Kin + Pot, Part_cd = (<N_c>, <N_d>)
!>  Eq./tau  : Green (orbital resolved), Den (orbital resolved),
!>             SpinZ = <S^z_i S^z_j>,   Pair = <P^dag_i P_j>  with
!>             P_i = 2 ( c_{i,up} c_{i,dn} - d_{i,up} d_{i,dn} )   (s+- pairing of the paper)
!>
!--------------------------------------------------------------------

    submodule (Hamiltonian_main) ham_Hubbard_Sz_smod

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
      Use LRC_Mod
      use runtime_error_mod

      Implicit none

      type, extends(ham_base) :: ham_Hubbard_Sz
      contains
        ! Set Hamiltonian-specific procedures
        procedure, nopass :: Ham_Set
        procedure, nopass :: Alloc_obs
        procedure, nopass :: Obser
        procedure, nopass :: ObserT
#ifdef HDF5
        procedure, nopass :: write_parameters_hdf5
#endif
      end type ham_Hubbard_Sz

      !#PARAMETERS START# VAR_lattice
      Character (len=64) :: Model = ''  ! Value irrelevant
      Character (len=64) :: Lattice_type = 'Square'  ! Possible Values: 'Square'
      Integer            :: L1 = 6   ! Length in direction a_1
      Integer            :: L2 = 6   ! Length in direction a_2
      !#PARAMETERS END#

      !#PARAMETERS START# VAR_Hubbard_Sz
      real(Kind=Kind(0.d0)) :: Ham_T     = 1.d0      ! Nearest-neighbour hopping t
      real(Kind=Kind(0.d0)) :: Ham_delta = 0.4d0     ! Hopping anisotropy delta
      real(Kind=Kind(0.d0)) :: Ham_mu    = -2.d0     ! Band offset mu: eps_c ~ -mu, eps_d ~ +mu (Eq. 1)
      real(Kind=Kind(0.d0)) :: Ham_chem  = 0.d0      ! Additional global chemical potential -Ham_chem*N
      real(Kind=Kind(0.d0)) :: Ham_U     = 1.d0      ! Interaction H_int = -U sum_i (S^z_i)^2
      real(Kind=Kind(0.d0)) :: Dtau      = 0.1d0     ! Thereby Ltrot=Beta/dtau
      real(Kind=Kind(0.d0)) :: Beta      = 5.d0      ! Inverse temperature
      !logical              :: Projector = .false.   ! Whether the projective algorithm is used
      real(Kind=Kind(0.d0)) :: Theta     = 10.d0     ! Projection parameter
      !logical              :: Symm      = .true.    ! Whether symmetrization takes place
      Integer               :: N_part    = -1        ! Particles per spin in trial wave function. If < 0 -> L1*L2
      !#PARAMETERS END#

      Type (Lattice),       target :: Latt
      Type (Unit_cell),     target :: Latt_unit      ! Two orbitals (c,d) per unit cell
      Type (Unit_cell),     target :: Latt_unit_S    ! One "orbital" per unit cell, for S^z and pair correlations
      Integer, allocatable :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell
      Real (Kind=Kind(0.d0)), allocatable :: H0_mat(:,:)  ! Single-particle Hamiltonian (per spin)

    contains

      module Subroutine Ham_Alloc_Hubbard_Sz
        allocate(ham_Hubbard_Sz::ham)
      end Subroutine Ham_Alloc_Hubbard_Sz

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_Hubbard_Sz_read_write_parameters.F90"

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the Hamiltonian
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

          ! From dynamically generated file "Hamiltonian_Hubbard_Sz_read_write_parameters.F90"
          call read_parameters()

          If (L1 < 2) then
             Write(error_unit,*) 'Hubbard_Sz: L1 has to be >= 2. For one-dimensional lattices set L2=1'
             Call Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          endif
          If (Ham_U < 0.d0) then
             Write(error_unit,*) 'Hubbard_Sz: Warning, Ham_U < 0 leads to a complex HS coupling (sign problem).'
          endif

          Ltrot  = nint(beta/dtau)
          Thtrot = 0
          if (Projector) Thtrot = nint(theta/dtau)
          Ltrot  = Ltrot+2*Thtrot
          N_SUN  = 2   ! Spin, after the gauge transformation d_dn -> -d_dn
          N_FL   = 1

          ! Setup the Bravais lattice
          Call Ham_Latt

          if (N_part < 0) N_part = Latt%N   ! Total density n_c + n_d = 2 per site

          ! Setup the hopping / single-particle part
          Call Ham_Hop

          ! Setup the interaction.
          Call Ham_V

          ! Setup the trial wave function, in case of a projector approach
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
             Write(unit_info,*) 'Model is      : Hubbard_Sz  (H_0 - U sum_i S^z_i S^z_i)'
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
!> Square lattice with two orbitals (c and d) on the same site.
!--------------------------------------------------------------------
        Subroutine Ham_Latt

          Implicit none
          Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)
          Integer :: I, no, nc

          If (Lattice_type /= "Square") then
             Write(error_unit,*) 'The Hubbard_Sz model is only defined for the square lattice'
             Call Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          Endif

          a1_p(1) =  1.0d0  ; a1_p(2) =  0.d0
          a2_p(1) =  0.0d0  ; a2_p(2) =  1.d0
          L1_p    =  dble(L1)*a1_p
          L2_p    =  dble(L2)*a2_p
          Call Make_Lattice( L1_p, L2_p, a1_p, a2_p, Latt )

          ! Unit cell for orbital-resolved quantities: orbital 1 = c, orbital 2 = d
          Latt_unit%Norb    = 2
          Latt_unit%N_coord = 2
          allocate(Latt_unit%Orb_pos_p(Latt_unit%Norb,2))
          Latt_unit%Orb_pos_p(1,:) = [0.d0, 0.d0]
          Latt_unit%Orb_pos_p(2,:) = [0.d0, 0.d0]

          ! Unit cell for on-site composite operators (S^z_i, P_i)
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
!> Real-space hopping matrix of Eq. (1), identical for both spins.
!> c-band: t_x = t+delta, t_y = t-delta, on-site -mu
!> d-band: t_x = t-delta, t_y = t+delta, on-site +mu
!--------------------------------------------------------------------
        Subroutine Set_H0_mat(H, Delta_mod)

          Implicit none
          Real (Kind=Kind(0.d0)), Intent(Out) :: H(:,:)
          Real (Kind=Kind(0.d0)), Intent(In)  :: Delta_mod   ! Small bond modulation (only for the trial wave function)

          Integer :: I, Ix, Iy, no, I1, J1
          Real (Kind=Kind(0.d0)) :: tx(2), ty(2), eps(2), X

          tx (1) = Ham_T + Ham_delta ; ty (1) = Ham_T - Ham_delta ; eps(1) = -Ham_mu - Ham_chem
          tx (2) = Ham_T - Ham_delta ; ty (2) = Ham_T + Ham_delta ; eps(2) =  Ham_mu - Ham_chem

          H = 0.d0
          Do I = 1,Latt%N
             X  = 1.d0 + Delta_mod*cos(acos(-1.d0)*real(Latt%list(I,1) + Latt%list(I,2),Kind(0.d0)))
             Ix = Latt%nnlist(I,1,0)
             Iy = Latt%nnlist(I,0,1)
             Do no = 1,Latt_unit%Norb
                I1 = Invlist(I,no)
                ! Bonds are added (+=) such that L=2 correctly yields -2 t cos(k)
                J1 = Invlist(Ix,no)
                H(I1,J1) = H(I1,J1) - tx(no)*X
                H(J1,I1) = H(J1,I1) - tx(no)*X
                If (L2 > 1) then
                   J1 = Invlist(Iy,no)
                   H(I1,J1) = H(I1,J1) - ty(no)
                   H(J1,I1) = H(J1,I1) - ty(no)
                Endif
                H(I1,I1) = H(I1,I1) + eps(no)
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
!> Slater determinant ground state of H_0 (with a small bond modulation
!> to lift degeneracies), N_part particles per spin.
!--------------------------------------------------------------------
        Subroutine Ham_Trial()

          Implicit none

          Integer                              :: nf, I, n
          Real (Kind=Kind(0.d0)), allocatable  :: H(:,:), U0(:,:), E0(:)
          Real (Kind=Kind(0.d0)), parameter    :: Delta_mod = 0.01d0

          If (N_part < 1 .or. N_part >= Ndim) then
             Write(error_unit,*) 'Hubbard_Sz: N_part has to satisfy 0 < N_part < Ndim'
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
!> H_int = -U sum_i [ sum_s (c^dag_{is} d_{is} + h.c.) ]^2   (gauge-transformed basis)
!> Predefined_Int_V_SUN(V) implements -V/N_SUN [ sum_s (c^dag_i c_j + h.c.) ]^2,
!> hence V = U*N_SUN  (g = sqrt(Dtau*U), alpha = 0, type = 2).
!--------------------------------------------------------------------
        Subroutine Ham_V

          Use Predefined_Int
          Implicit none

          Integer :: I

          Allocate(Op_V(Latt%N,N_FL))
          Do I = 1,Latt%N
             Call Predefined_Int_V_SUN( Op_V(I,1), Invlist(I,1), Invlist(I,2), N_SUN, Dtau, &
                  &                     Ham_U*real(N_SUN,kind(0.d0)) )
          Enddo

        end Subroutine Ham_V

!--------------------------------------------------------------------
!> @brief
!> Specifiy the equal time and time displaced observables
!--------------------------------------------------------------------
        Subroutine Alloc_obs(Ltau)

          Implicit none
          !>  Ltau=1 if time displaced correlations are considered.
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
          Allocate ( Obs_eq(4) )
          Do I = 1,Size(Obs_eq,1)
             Nt = 1
             Channel = '--'
             select case (I)
             case (1)
                Filename = "Green"
                Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit,   Channel, dtau)
             case (2)
                Filename = "Den"
                Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit,   Channel, dtau)
             case (3)
                Filename = "SpinZ"
                Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_S, Channel, dtau)
             case (4)
                Filename = "Pair"
                Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_S, Channel, dtau)
             case default
                Write(error_unit,*) ' Error in Alloc_obs '
             end select
          enddo

          If (Ltau == 1) then
             ! Time-displaced correlators
             Allocate ( Obs_tau(4) )
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
                   Channel = 'PH'; Filename = "SpinZ"
                   If (Projector) Channel = 'T0'
                   Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit_S, Channel, dtau)
                case (4)
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
!> @author
!> ALF Collaboration
!>
!> @brief
!> Computes equal time observables
!> @details
!> @param [IN] Gr   Complex(:,:,:)
!> \verbatim
!>  Green function: Gr(I,J,nf) = <c_{I,nf } c^{dagger}_{J,nf } > on time slice ntau
!> \endverbatim
!> @param [IN] Phase   Complex
!> \verbatim
!>  Phase
!> \endverbatim
!> @param [IN] Ntau Integer
!> \verbatim
!>  Time slice
!> \endverbatim
!-------------------------------------------------------------------
        subroutine Obser(GR,Phase,Ntau, Mc_step_weight)

          Use Predefined_Obs

          Implicit none

          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
          Integer,                   INTENT(IN) :: Ntau
          Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight

          !Local
          Complex (Kind=Kind(0.d0)), allocatable :: GRC(:,:,:)
          Complex (Kind=Kind(0.d0)) :: ZP, ZS, ZN, Zkin, Zpot, Zrho, Zc, Zd, ZA_I, ZA_J, ZZ, ZPair
          Integer :: I, J, I1, J1, imj, nf, n1, n2
          Integer :: IS(2), JS(2)    ! (c,d) indices of sites I and J

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))

          ZS = ZS*Mc_step_weight
          ZN = cmplx(dble(N_SUN), 0.d0, kind(0.D0))

          allocate(GRC(Ndim,Ndim,N_FL))

          Do nf = 1,N_FL
             Do I = 1,Ndim
                Do J = 1,Ndim
                   GRC(I, J, nf) = -GR(J, I, nf)
                Enddo
                GRC(I, I, nf) = 1.D0 + GRC(I, I, nf)
             Enddo
          Enddo
          ! GRC(i,j,nf) = < c^{dagger}_{i,nf } c_{j,nf } >

          ! Compute scalar observables.
          Do I = 1,Size(Obs_scal,1)
             Obs_scal(I)%N         =  Obs_scal(I)%N + 1
             Obs_scal(I)%Ave_sign  =  Obs_scal(I)%Ave_sign + Real(ZS,kind(0.d0))
          Enddo

          ! Kinetic energy (incl. on-site band offsets): <H_0> = N_SUN sum_{ij} H0(i,j) <c^dag_i c_j>
          Zkin = cmplx(0.d0, 0.d0, kind(0.D0))
          Do I = 1,Ndim
             Do J = 1,Ndim
                If (H0_mat(I,J) /= 0.d0) Zkin = Zkin + H0_mat(I,J)*GRC(I,J,1)
             Enddo
          Enddo
          Zkin = Zkin*ZN
          Obs_scal(1)%Obs_vec(1) = Obs_scal(1)%Obs_vec(1) + Zkin*ZP*ZS

          ! Potential energy: -U sum_i <(S^z_i)^2>
          Zpot = cmplx(0.d0, 0.d0, kind(0.D0))
          Do I = 1,Latt%N
             IS = Invlist(I,:)
             Zpot = Zpot + SzSz_eq(IS, IS)
          Enddo
          Zpot = -Ham_U*Zpot
          Obs_scal(2)%Obs_vec(1) = Obs_scal(2)%Obs_vec(1) + Zpot*ZP*ZS

          ! Particle numbers
          Zc = cmplx(0.d0, 0.d0, kind(0.D0))
          Zd = cmplx(0.d0, 0.d0, kind(0.D0))
          Do I = 1,Latt%N
             Zc = Zc + GRC(Invlist(I,1),Invlist(I,1),1)
             Zd = Zd + GRC(Invlist(I,2),Invlist(I,2),1)
          Enddo
          Zc   = Zc*ZN
          Zd   = Zd*ZN
          Zrho = Zc + Zd
          Obs_scal(3)%Obs_vec(1) = Obs_scal(3)%Obs_vec(1) + Zrho*ZP*ZS
          Obs_scal(4)%Obs_vec(1) = Obs_scal(4)%Obs_vec(1) + (Zkin + Zpot)*ZP*ZS
          Obs_scal(5)%Obs_vec(1) = Obs_scal(5)%Obs_vec(1) + Zc*ZP*ZS
          Obs_scal(5)%Obs_vec(2) = Obs_scal(5)%Obs_vec(2) + Zd*ZP*ZS

          ! Equal-time correlations
          Call Predefined_Obs_eq_Green_measure( Latt, Latt_unit, List, GR, GRC, N_SUN, ZS, ZP, Obs_eq(1) )
          Call Predefined_Obs_eq_Den_measure  ( Latt, Latt_unit, List, GR, GRC, N_SUN, ZS, ZP, Obs_eq(2) )

          Do I1 = 3,4
             Obs_eq(I1)%N        = Obs_eq(I1)%N + 1
             Obs_eq(I1)%Ave_sign = Obs_eq(I1)%Ave_sign + Real(ZS,kind(0.d0))
          Enddo
          Do I = 1,Latt%N
             IS   = Invlist(I,:)
             ZA_I = ZN*( GRC(IS(1),IS(2),1) + GRC(IS(2),IS(1),1) )   ! <S^z_I>
             Do J = 1,Latt%N
                JS  = Invlist(J,:)
                imj = Latt%imj(I,J)
                ! <S^z_I S^z_J>
                ZZ  = SzSz_eq(IS, JS)
                ! <P^dag_I P_J>,  P_I = 2 sum_o x_{I o up} x_{I o dn}  (gauge basis)
                ZPair = cmplx(0.d0, 0.d0, kind(0.D0))
                Do n1 = 1,2
                   Do n2 = 1,2
                      ZPair = ZPair + GRC(IS(n1),JS(n2),1)**2
                   Enddo
                Enddo
                ZPair = 4.d0*ZPair
                Obs_eq(3)%Obs_Latt(imj,1,1,1) = Obs_eq(3)%Obs_Latt(imj,1,1,1) + ZZ   *ZP*ZS
                Obs_eq(4)%Obs_Latt(imj,1,1,1) = Obs_eq(4)%Obs_Latt(imj,1,1,1) + ZPair*ZP*ZS
             Enddo
             Obs_eq(3)%Obs_Latt0(1) = Obs_eq(3)%Obs_Latt0(1) + ZA_I*ZP*ZS
          Enddo

          deallocate(GRC)

        contains

          !> <S^z_I S^z_J> for N_SUN identical colors, S^z_I = sum_s (x^dag_{Ic s} x_{Id s} + h.c.)
          !> = N^2 <A_I><A_J> + N sum_{(a,b),(c,e)} <x^dag_a x_e> <x_b x^dag_c>
          Complex (Kind=Kind(0.d0)) function SzSz_eq(IS, JS)
            Integer, Intent(In) :: IS(2), JS(2)
            Integer :: n1, n2, a, b, c, e
            Complex (Kind=Kind(0.d0)) :: ZA_I, ZA_J, Zcon
            ZA_I = GRC(IS(1),IS(2),1) + GRC(IS(2),IS(1),1)
            ZA_J = GRC(JS(1),JS(2),1) + GRC(JS(2),JS(1),1)
            Zcon = cmplx(0.d0, 0.d0, kind(0.D0))
            Do n1 = 1,2
               a = IS(n1); b = IS(3-n1)
               Do n2 = 1,2
                  c = JS(n2); e = JS(3-n2)
                  Zcon = Zcon + GRC(a,e,1)*GR(b,c,1)
               Enddo
            Enddo
            SzSz_eq = ZN*ZN*ZA_I*ZA_J + ZN*Zcon
          end function SzSz_eq

        end Subroutine Obser

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Computes time displaced  observables
!> @details
!> @param [IN] NT, Integer
!> \verbatim
!>  Imaginary time
!> \endverbatim
!> @param [IN] GT0, GTT, G00, GTT,  Complex(:,:,:)
!> \verbatim
!>  Green functions:
!>  GT0(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(0  )>
!>  G0T(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(tau)>
!>  G00(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(0  )>
!>  GTT(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(tau)>
!> \endverbatim
!> @param [IN] Phase   Complex
!> \verbatim
!>  Phase
!> \endverbatim
!-------------------------------------------------------------------
        Subroutine ObserT(NT,  GT0,G0T,G00,GTT, PHASE,  Mc_step_weight)

          Use Predefined_Obs

          Implicit none

          Integer         , INTENT(IN) :: NT
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase
          Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight

          !Locals
          Complex (Kind=Kind(0.d0)) :: ZP, ZS, ZN, ZA_I, ZA_J, Zcon, ZZ, ZPair
          Integer :: I, J, I1, imj, n1, n2, a, b, c, e
          Integer :: IS(2), JS(2)

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          ZS = ZS * Mc_step_weight
          ZN = cmplx(dble(N_SUN), 0.d0, kind(0.D0))

          Call Predefined_Obs_tau_Green_measure( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT, N_SUN, ZS, ZP, Obs_tau(1) )
          Call Predefined_Obs_tau_Den_measure  ( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT, N_SUN, ZS, ZP, Obs_tau(2) )

          If (NT == 0) then
             Do I1 = 3,4
                Obs_tau(I1)%N        = Obs_tau(I1)%N + 1
                Obs_tau(I1)%Ave_sign = Obs_tau(I1)%Ave_sign + Real(ZS,kind(0.d0))
             Enddo
          Endif

          Do I = 1,Latt%N
             IS   = Invlist(I,:)
             ! <A_I(tau)> per color,  <x^dag_a x_b> = delta_ab - G(b,a)
             ZA_I = -GTT(IS(2),IS(1),1) - GTT(IS(1),IS(2),1)
             Do J = 1,Latt%N
                JS   = Invlist(J,:)
                imj  = Latt%imj(I,J)
                ZA_J = -G00(JS(2),JS(1),1) - G00(JS(1),JS(2),1)
                ! <S^z_I(tau) S^z_J(0)>:  connected part  <x^dag_a(tau) x_e(0)> <x_b(tau) x^dag_c(0)>
                !                                       = -G0T(e,a) * GT0(b,c)
                Zcon = cmplx(0.d0, 0.d0, kind(0.D0))
                Do n1 = 1,2
                   a = IS(n1); b = IS(3-n1)
                   Do n2 = 1,2
                      c = JS(n2); e = JS(3-n2)
                      Zcon = Zcon - G0T(e,a,1)*GT0(b,c,1)
                   Enddo
                Enddo
                ZZ = ZN*ZN*ZA_I*ZA_J + ZN*Zcon
                ! <P^dag_I(tau) P_J(0)> = 4 sum_{o,o'} <x^dag_{Io}(tau) x_{Jo'}>^2 = 4 sum G0T(Jo',Io)^2
                ZPair = cmplx(0.d0, 0.d0, kind(0.D0))
                Do n1 = 1,2
                   Do n2 = 1,2
                      ZPair = ZPair + G0T(JS(n2),IS(n1),1)**2
                   Enddo
                Enddo
                ZPair = 4.d0*ZPair
                Obs_tau(3)%Obs_Latt(imj,NT+1,1,1) = Obs_tau(3)%Obs_Latt(imj,NT+1,1,1) + ZZ   *ZP*ZS
                Obs_tau(4)%Obs_Latt(imj,NT+1,1,1) = Obs_tau(4)%Obs_Latt(imj,NT+1,1,1) + ZPair*ZP*ZS
             Enddo
             Obs_tau(3)%Obs_Latt0(1) = Obs_tau(3)%Obs_Latt0(1) + ZN*ZA_I*ZP*ZS
          Enddo

        end Subroutine OBSERT

    end submodule ham_Hubbard_Sz_smod
