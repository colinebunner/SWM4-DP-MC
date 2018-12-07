module energy
  implicit none
  private
  save
  public :: energy_nonbond, energy_drude, energy_qreal, energy_qkspace, energy_qexclude, energy_qmimage

  real*8,parameter::N_Avogadro=6.02214129E23&
   ,k_B=1.3806488E-23

  real*8,parameter::e_unit=1.602176565E-19&
   ,eps_0=8.854187817E-12&
   ,eXV_to_K=e_unit/k_B& !< e times V to Kelvin
   ,qqfact=(e_unit**2)*1E10/(3.1415926535897932*4.0d0)/eps_0/k_B !< conversion factor for coulomb interactions (including e-->C, m-->Ang, fourpi*epsilon_0, J-->K)

contains

  subroutine energy_nonbond(nmol,nunit,rcut,boxlength,lpbc,xcoords,ycoords,zcoords,sig,eps,enb)
    integer,intent(in) :: nmol, nunit
    logical,intent(in) :: lpbc
    real*8,intent(in)  :: rcut, xcoords(:,:), ycoords(:,:), zcoords(:,:),boxlength,sig,eps
    real*8,intent(out) :: enb

    integer :: imol, jmol
    real*8  :: rxi, ryi, rzi, rxij, ryij, rzij, r2, sigr2, sigr6, rc2

    enb = 0.0d0

    rc2 = rcut*rcut

    ! Loop over unique pairs of molecules
    first_mol: do imol = 1, nmol-1

      ! Store oxygen coordinate (others have no LJ contribution)
      rxi   = xcoords(imol,1)
      ryi   = ycoords(imol,1)
      rzi   = zcoords(imol,1)

      second_mol: do jmol = imol+1, nmol
       
        rxij = rxi - xcoords(jmol,1)
        ryij = ryi - ycoords(jmol,1)
        rzij = rzi - zcoords(jmol,1)
        if (lpbc) then
          rxij = rxij - boxlength*nint(rxij/boxlength)
          ryij = ryij - boxlength*nint(ryij/boxlength)
          rzij = rzij - boxlength*nint(rzij/boxlength)
        end if
        r2 = rxij*rxij + ryij*ryij + rzij*rzij

        if (r2.lt.rc2) then
          ! Non-bonded energy (no Coulomb)
          sigr2 = (sig*sig)/r2
          sigr6 = sigr2*sigr2*sigr2
          enb   = enb + 4.0d0*eps*sigr6*(sigr6-1.0d0)
        end if
      end do second_mol
    end do first_mol

  return

  end subroutine energy_nonbond

  subroutine energy_drude(nmol,nunit,lrestrain,xcoords,ycoords,zcoords,kdrude,edrude)
    integer,intent(in) :: nmol, nunit
    logical,intent(in) :: lrestrain
    real*8,intent(in)  :: xcoords(:,:), ycoords(:,:), zcoords(:,:), kdrude
    real*8,intent(out) :: edrude

    integer :: imol
    real*8  :: rxij, ryij, rzij, r2

    edrude = 0.0d0

    ! Loop over all molecules, calculate Drude oscillator energy
    molecules: do imol = 1, nmol

      rxij = xcoords(imol,5) - xcoords(imol,1)
      ryij = ycoords(imol,5) - ycoords(imol,1)
      rzij = zcoords(imol,5) - zcoords(imol,1)

      r2 = rxij*rxij + ryij*ryij + rzij*rzij

      edrude = edrude + 0.5d0*kdrude*r2

      ! Optional restraining terms. Quartic and sextic terms with force
      ! constants 10.0 and 100.0 times larger, respectively.
      if (lrestrain) then
        edrude = edrude + 500.0d0*kdrude*r2*r2 + 5000.0d0*kdrude*r2*r2*r2
      end if

    end do molecules

    return

  end subroutine energy_drude

  ! Computes real-space portion of Coulomb energy
  subroutine energy_qreal(nmol,nunit,boxlength,rcut,kalp,xcoords,ycoords,zcoords,qbeads,eqreal,lovr)
    integer,intent(in) :: nmol, nunit
    real*8,intent(in)  :: xcoords(:,:), ycoords(:,:), zcoords(:,:), qbeads(5),&
                         &boxlength,kalp,rcut
    real*8,intent(out) :: eqreal
    logical,intent(out) :: lovr

    integer :: imol, jmol, iunit, junit
    real*8  :: rxi,ryi,rzi,rxij,ryij,rzij,rij

    ! Set total real space energy to 0
    eqreal = 0.0d0

    lovr = .false.

    ! Loop over unique pairs of molecules
    first_mol: do imol = 1, nmol-1
      second_mol: do jmol = imol+1, nmol
      
        ! Subloop over interacting pairs of charges. All 5 beads charged in this
        ! model. 
        first_unit: do iunit = 1, 5

          rxi = xcoords(imol,iunit)
          ryi = ycoords(imol,iunit)
          rzi = zcoords(imol,iunit)

          second_unit: do junit = 1, 5

            rxij = rxi - xcoords(jmol,junit)
            ryij = ryi - ycoords(jmol,junit)
            rzij = rzi - zcoords(jmol,junit)

            rxij = rxij - boxlength*nint(rxij/boxlength)
            ryij = ryij - boxlength*nint(ryij/boxlength)
            rzij = rzij - boxlength*nint(rzij/boxlength)

            rij = sqrt(rxij*rxij + ryij*ryij + rzij*rzij)

            if (rij.lt.1.2e0) then
              lovr = .true.
              return
            end if

            if (rij.lt.rcut) then
              eqreal = eqreal + (qbeads(iunit)*qbeads(junit)*erfc(kalp*rij))/rij
            end if

          end do second_unit
        end do first_unit
      end do second_mol
    end do first_mol

    eqreal = eqreal*qqfact

    return

  end subroutine energy_qreal

  ! Computes k-space portion of Coulomb energy
  subroutine energy_qkspace(nmol,nunit,boxlength,kalp,nkvec,xcoords,ycoords,zcoords,qbeads,eqrecip)
    integer,intent(in) :: nmol, nunit, nkvec
    real*8,intent(in)  :: xcoords(:,:), ycoords(:,:), zcoords(:,:), qbeads(5),&
                         &boxlength,kalp
    real*8,intent(out) :: eqrecip

    integer :: imol, iunit, nx, ny, nz, nymin, nzmin
    real*8  :: kx,ky,kz,ksq,sumr,sumi,tpl,arg,prefac,vol,kalpsq

    tpl    = (2.0d0*3.1415926535897932)/boxlength

    ! Set total k-space energy to 0. Set real and imaginary contributions to 0.
    eqrecip = 0.0d0
    vol = boxlength*boxlength*boxlength

    kalpsq = kalp*kalp

    ! Loop over k vectors
    do nx = 0, nkvec

      if (nx.eq.0) then
        nymin = 0
      else
        nymin = -nkvec
      end if

      do ny = nymin, nkvec

        if (nx.eq.0 .and. ny.eq.0) then
          nzmin = 1
        else
          nzmin = -nkvec
        end if

        do nz = nzmin, nkvec

          sumr   = 0.0d0
          sumi   = 0.0d0
         
          kx = real(nx)*tpl
          ky = real(ny)*tpl
          kz = real(nz)*tpl
          ksq = kx*kx+ky*ky+kz*kz

          prefac = exp(-ksq/(4.0d0*kalpsq)) / (ksq*vol)

          first_mol: do imol = 1, nmol
            first_unit: do iunit = 1, 5
              arg  = kx*xcoords(imol,iunit)+ky*ycoords(imol,iunit)+kz*zcoords(imol,iunit)
              sumr = sumr + cos(arg)*qbeads(iunit)
              sumi = sumi + sin(arg)*qbeads(iunit)
            end do first_unit
          end do first_mol

          print *, nx, ny, nz, (sumr*sumr+sumi*sumi)*prefac*qqfact

          eqrecip = eqrecip + (sumr*sumr+sumi*sumi)*prefac

        end do
      end do
    end do

    eqrecip = eqrecip*qqfact

    return

  end subroutine energy_qkspace

  ! Removes intramolecular charge interaction between offdiagonal terms (i.e.
  ! not self-correction)
  subroutine energy_qexclude(nmol,nunit,kalp,xcoords,ycoords,zcoords,qbeads,eqexclude)
    integer,intent(in) :: nmol, nunit
    real*8,intent(in)  :: xcoords(:,:), ycoords(:,:), zcoords(:,:), qbeads(5),&
                         &kalp
    real*8,intent(out) :: eqexclude

    integer :: imol, iunit, junit
    real*8  :: bxi,byi,bzi,bx,by,bz,bdist

    eqexclude = 0.0d0

    first_mol: do imol = 1, nmol
      first_unit: do iunit = 1, 4
        bxi = xcoords(imol,iunit)
        byi = ycoords(imol,iunit)
        bzi = zcoords(imol,iunit)
        second_unit: do junit = iunit+1, 5
          bx = xcoords(imol,junit) - bxi
          by = ycoords(imol,junit) - byi
          bz = zcoords(imol,junit) - bzi
          bdist = sqrt(bx*bx+by*by+bz*bz)
          eqexclude = eqexclude + qbeads(iunit)*qbeads(junit)*(erf(kalp*bdist)/bdist)
        end do second_unit
      end do first_unit
    end do first_mol

    eqexclude = eqexclude*qqfact

    return

  end subroutine energy_qexclude
  
  ! Hedging my bets and coding a subroutine that just treats Coulomb
  ! interactions with minimum image convention
  subroutine energy_qmimage(nmol,nunit,boxlength,lpbc,xcoords,ycoords,zcoords,qbeads,eqmimage,lovr)
    integer,intent(in)  :: nmol, nunit
    logical,intent(in)  :: lpbc
    real*8,intent(in)   :: xcoords(:,:), ycoords(:,:), zcoords(:,:), qbeads(5),&
                         &boxlength
    real*8,intent(out)  :: eqmimage
    logical,intent(out) :: lovr

    integer :: imol, jmol, iunit, junit
    real*8  :: rxi,ryi,rzi,rxij,ryij,rzij,rij

    ! Set total real space energy to 0
    eqmimage = 0.0d0

    ! Since all beads have charges, we can use this subroutine to implement
    ! hard-cutoff on moves
    lovr = .false.

    ! Loop over unique pairs of molecules
    first_mol: do imol = 1, nmol-1
      second_mol: do jmol = imol+1, nmol
      
        ! Subloop over interacting pairs of charges. All 5 beads charged in this
        ! model. 
        first_unit: do iunit = 1, 5

          rxi = xcoords(imol,iunit)
          ryi = ycoords(imol,iunit)
          rzi = zcoords(imol,iunit)

          second_unit: do junit = 1, 5

            rxij = rxi - xcoords(jmol,junit)
            ryij = ryi - ycoords(jmol,junit)
            rzij = rzi - zcoords(jmol,junit)
            if (lpbc) then
              rxij = rxij - boxlength*nint(rxij/boxlength)
              ryij = ryij - boxlength*nint(ryij/boxlength)
              rzij = rzij - boxlength*nint(rzij/boxlength)
            end if
            rij = sqrt(rxij*rxij + ryij*ryij + rzij*rzij)

            if (rij.lt.1.2) then
              lovr = .true.
              return
            end if

            eqmimage = eqmimage + (qbeads(iunit)*qbeads(junit))/rij

          end do second_unit
        end do first_unit
      end do second_mol
    end do first_mol

    eqmimage = eqmimage*qqfact

    return
  end subroutine energy_qmimage

end module
