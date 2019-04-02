module tuv_photolysis
  use machine,           only: rk => kind_phys

  use module_prates_tuv, only: nj, nlambda_start, wc, wl, nwave, is_full_tuv
  use module_prates_tuv, only: o2_xs, so2_xs, o3_xs_tab, no2_xs_tab
  use module_prates_tuv, only: calc_tuv_prates, calc_tuv_init
  use module_prates_tuv, only: xs_int ! for temperature dependent (level dependent) o3_xs and no2_xs
  use tuv_subs, only: tuv_radfld
  use module_xsections,  only: o3xs, no2xs_jpl06a
  
  implicit none
  private
  public :: tuv_photolysis_init
  public :: tuv_photolysis_run
  public :: tuv_photolysis_finalize
  public :: tuv_jnames
  public :: tuv_n_phot
  
  integer :: cld_od_opt=1
  logical :: has_aer_ra_feedback = .false.
  real    :: dobsi = 0.
  integer :: nlev 
  integer :: nlyr
  
  real, parameter :: kboltz= 1.38064852e-16 ! boltzmann constant (erg/K)
  real, parameter :: R=2.8704e6       ! gas constant (erg/g/K)
  real, parameter :: g=980.616        ! grav acceleration (cm/sec2)

  real, parameter :: tuv_n_phot = 113
  character(len=16), parameter :: tuv_jnames(tuv_n_phot) = &
    (/'j_o2            ' &
    , 'j_o1d           ' &
    , 'j_o3p           ' &
    , 'j_no2           ' &
    , 'j_no3_a         ' &
    , 'j_no3_b         ' &
    , 'j_n2o5_a        ' &
    , 'j_n2o5_b        ' &
    , 'j_hno2          ' &
    , 'j_hno3          ' &
    , 'j_hno4          ' &
    , 'j_h2o2          ' &
    , 'j_chbr3         ' &
    , 'j_ch3cho_a      ' &
    , 'j_ch3cho_b      ' &
    , 'j_ch3cho_c      ' &
    , 'j_c2h5cho       ' &
    , 'j_gly_a         ' &
    , 'j_gly_b         ' &
    , 'j_gly_c         ' &
    , 'j_mgly          ' &
    , 'j_ch3coch3      ' &
    , 'j_ch3ooh        ' &
    , 'j_ch3ono2       ' &
    , 'j_pan_a         ' &
    , 'j_pan_b         ' &
    , 'j_ccl2o         ' &
    , 'j_ccl4          ' &
    , 'j_cclfo         ' &
    , 'j_cf2o          ' &
    , 'j_cf2clcfcl2    ' &
    , 'j_cf2clcf2cl    ' &
    , 'j_cf3cf2cl      ' &
    , 'j_ccl3f         ' &
    , 'j_ccl2f2        ' &
    , 'j_ch3br         ' &
    , 'j_ch3ccl3       ' &
    , 'j_ch3cl         ' &
    , 'j_cloo          ' &
    , 'j_cf3chcl2      ' &
    , 'j_cf3chfcl      ' &
    , 'j_ch3cfcl2      ' &
    , 'j_ch3cf2cl      ' &
    , 'j_cf3cf2chcl2   ' &
    , 'j_cf2clcf2chfcl ' &
    , 'j_chclf2        ' &
    , 'j_ho2           ' &
    , 'j_cf2bf2        ' &
    , 'j_cf2brcl       ' &
    , 'j_cf3br         ' &
    , 'j_cf2brcf2br    ' &
    , 'j_n2o           ' &
    , 'j_clono2_a      ' &
    , 'j_clono2_b      ' &
    , 'j_brono2_a      ' &
    , 'j_brono2_b      ' &
    , 'j_cl2           ' &
    , 'j_glyald_a      ' &
    , 'j_glyald_b      ' &
    , 'j_glyald_c      ' &
    , 'j_biacetyl      ' &
    , 'j_mvk           ' &
    , 'j_macr          ' &
    , 'j_ch3cocooh     ' &
    , 'j_ch3ch2ono2    ' &
    , 'j_ch3chono2ch3  ' &
    , 'j_ch2ohch2ono2  ' &
    , 'j_ch3coch2ono2  ' &
    , 'j_bnit1         ' &
    , 'j_cloocl        ' &
    , 'j_hyac_a        ' &
    , 'j_hyac_b        ' &
    , 'j_hobr          ' &
    , 'j_bro           ' &
    , 'j_br2           ' &
    , 'j_no3_aq_a      ' &
    , 'j_no3_aq_b      ' &
    , 'j_no3_aq_c      ' &
    , 'j_mek           ' &
    , 'j_ppn_a         ' &
    , 'j_ppn_b         ' &
    , 'j_hoch2ooh      ' &
    , 'j_acrol         ' &
    , 'j_ch3coooh      ' &
    , 'j_amine         ' &
    , 'j_clo_a         ' &
    , 'j_clo_b         ' &
    , 'j_clno2         ' &
    , 'j_brno          ' &
    , 'j_brno2         ' &
    , 'j_brono_a       ' &
    , 'j_brono_b       ' &
    , 'j_hocl          ' &
    , 'j_nocl          ' &
    , 'j_oclo          ' &
    , 'j_brcl          ' &
    , 'j_ch3oono2      ' &
    , 'j_bnit2         ' &
    , 'j_clono         ' &
    , 'j_hcl           ' &
    , 'j_ch2o_r        ' &
    , 'j_ch2o_m        ' &
    , 'j_ch3cooh       ' &
    , 'j_ch3ocl        ' &
    , 'j_chcl3         ' &
    , 'j_c2h5ono2      ' &
    , 'j_nc3h7ono2     ' &
    , 'j_1c4h9ono2     ' &
    , 'j_2c4h9ono2     ' &
    , 'j_perfluoro     ' &
    , 'j_i2            ' &
    , 'j_io            ' &
    , 'j_ioh           ' &
    /)

contains

!> \section arg_table_tuv_photolysis_init Argument Table
!! | local_name | standard_name             | long_name                 | units   | rank | type      | kind      | intent | optional |
!! |------------|---------------------------|---------------------------|---------|------|-----------|-----------|--------|----------|
!! | nlevels    | num_levels_for_photolysis | number of column layers   | count   |    0 | integer   |           | in     | F        |
!! ! input_root | tuv_inputdata_root        | root dir for TUV inputs   | none    |    0 | character | len=*     | in    | F        |
!! | errmsg     | ccpp_error_message        | CCPP error message        | none    |    0 | character | len=512   | out    | F        |
!! | errflg     | ccpp_error_flag           | CCPP error flag           | flag    |    0 | integer   |           | out    | F        |
!!
  subroutine tuv_photolysis_init(nlevels, input_root, errmsg, errflg)
    use params_mod, only: input_data_root

    integer, intent(in) :: nlevels
    character(len=*),   intent(in ) :: input_root
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    character(len=300) :: xsqy_filepath
    logical :: full_tuv

    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    input_data_root = trim(input_root)
    xsqy_filepath = trim(input_root)//'/wrf_tuv_xsqy.nc'
    full_tuv=.true.

 
    nlev = nlevels
    nlyr = nlev-1

    call calc_tuv_init(xsqy_filepath, full_tuv, tuv_jnames)

  end subroutine tuv_photolysis_init

!> \section arg_table_tuv_photolysis_run Argument Table
!! | local_name | standard_name                         | long_name                      | units     | rank | type      | kind      | intent | optional |
!! |------------|---------------------------------------|--------------------------------|-----------|------|-----------|-----------|--------|----------|
!! | zenith     | solar_zenith                          | solar zenith angle             | degrees   |    0 | real      | kind_phys | in     | F        |
!! | albedo     | surface_albedo                        | surface albedo                 | none      |    0 | real      | kind_phys | in     | F        |
!! | press_mid  | layer_pressure                        | mid-point layer pressure       | Pa        |    1 | real      | kind_phys | in     | F        |
!! | press_int  | layer_interface_pressure              | layer interface pressure       | Pa        |    1 | real      | kind_phys | in     | F        |
!! | alt        | layer_altitude                        | mid-point layer altitude       | km        |    1 | real      | kind_phys | in     | F        |
!! | temp       | layer_temperature                     | mid-point layer temperature    | K         |    1 | real      | kind_phys | in     | F        |
!! | o2vmr      | O2_vmr_col                            | O2 volume mixing ratio column  | mole/mole |    1 | real      | kind_phys | in     | F        |
!! | o3vmr      | O3_vmr_col                            | O3 volume mixing ratio column  | mole/mole |    1 | real      | kind_phys | in     | F        |
!! | so2vmr     | SO2_vmr_col                           | SO2 volume mixing ratio column | mole/mole |    1 | real      | kind_phys | in     | F        |
!! | no2vmr     | NO2_vmr_col                           | NO2 volume mixing ratio column | mole/mole |    1 | real      | kind_phys | in     | F        |
!! | prates     | photolysis_rates_col                  | photolysis rates column        | s-1       |    2 | real      | kind_phys | out    | F        |
!! | o3totcol   | ozone_column_density                  | total ozone column density     | DU        |    0 | real      | kind_phys | out    | F        |
!! | errmsg     | ccpp_error_message                    | CCPP error message             | none      |    0 | character | len=512   | out    | F        |
!! | errflg     | ccpp_error_flag                       | CCPP error flag                | flag      |    0 | integer   |           | out    | F        |
!!
  subroutine tuv_photolysis_run( zenith, albedo, press_mid, press_int, alt, temp, o2vmr, o3vmr, so2vmr, no2vmr, prates, o3totcol, errmsg, errflg )

    real(rk), intent(in) :: zenith
    real(rk), intent(in) :: albedo
    real(rk), intent(in) :: press_int(:)
    real(rk), intent(in) :: press_mid(:)
    real(rk), intent(in) :: alt(:)  ! km
    real(rk), intent(in) :: temp(:) ! K
    real(rk), intent(in) :: o2vmr(:)
    real(rk), intent(in) :: o3vmr(:)
    real(rk), intent(in) :: so2vmr(:)
    real(rk), intent(in) :: no2vmr(:)
    real(rk), intent(out) :: prates(:,:) ! /sec
    real(rk), intent(out) :: o3totcol ! total # molecules / cm2 
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    real :: zlev(nlev) ! km 
    real :: tlev(nlev)
    real :: dz_cm(nlyr)   ! thickness of each layer (cm)
    real :: airdens(nlev) ! # molecules / cm3 at each level
    real :: aircol(nlyr)  ! # molecules / cm2 in each layer
    real :: o2col(nlyr)  
    real :: o3col(nlyr) 
    real :: so2col(nlyr)
    real :: no2col(nlyr)
    real :: dpress(nlyr)
    real :: zen
    real :: alb
    real, allocatable :: tuv_prates(:,:)
    integer :: j, k
    real, parameter :: du_fac = 1./2.687e16

    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    dpress(1:nlyr) = press_mid(2:nlyr+1) - press_mid(1:nlyr)
    do k=1,nlev
       airdens(k) = 10.*press_mid(k)/(kboltz*temp(k))
    end do
    do k=1,nlyr
       aircol(k) = 10.*dpress(k)*R/(kboltz*g)
       o3col(k) = 0.5*(o3vmr(k)+o3vmr(k+1))*aircol(k)
       o2col(k) = 0.5*(o2vmr(k)+o2vmr(k+1))*aircol(k)
       so2col(k) = 0.5*(so2vmr(k)+so2vmr(k+1))*aircol(k)
       no2col(k) = 0.5*(no2vmr(k)+no2vmr(k+1))*aircol(k)
    end do

    aircol(1:nlyr)  =  aircol(nlyr:1:-1)
    airdens(1:nlyr) =  airdens(nlyr:1:-1)

    o3col(1:nlyr)  = o3col(nlyr:1:-1)
    o2col(1:nlyr)  = o2col(nlyr:1:-1)
    so2col(1:nlyr) = so2col(nlyr:1:-1)
    no2col(1:nlyr) = no2col(nlyr:1:-1)

    o3totcol = sum(o3col)*du_fac
 
    tlev(nlev:1:-1) = temp(1:nlev)
    zlev(nlev:1:-1) = alt(1:nlev) ! km

    allocate(tuv_prates(nlev, nj))

    zen = zenith
    alb = albedo
    call tuv_photolysis_calc( zen, alb, zlev, tlev, airdens, aircol, o2col, o3col, so2col, no2col, tuv_prates )
    do j=1,nj
       prates(1:nlev,j) = tuv_prates(nlev:1:-1,j)
    end do
    deallocate(tuv_prates)

  end subroutine tuv_photolysis_run
  
!> \section arg_table_tuv_photolysis_finalize Argument Table
!! | local_name | standard_name                         | long_name                      | units     | rank | type      | kind      | intent | optional |
!! |------------|---------------------------------------|--------------------------------|-----------|------|-----------|-----------|--------|----------|
!! | errmsg     | ccpp_error_message                    | CCPP error message             | none      |    0 | character | len=512   | out    | F        |
!! | errflg     | ccpp_error_flag                       | CCPP error flag                | flag      |    0 | integer   |           | out    | F        |
!!
  subroutine tuv_photolysis_finalize( errmsg, errflg )

    !--- arguments
    character(len=512), intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

  end subroutine tuv_photolysis_finalize

  subroutine tuv_photolysis_calc( zenith, albedo_in, zlev, tlev, airdens, aircol, o2col, o3col, so2col, no2col, tuv_prates )

    real, intent(in) :: zenith
    real, intent(in) :: albedo_in
    real, intent(in) :: zlev(:)
    real, intent(in) :: tlev(:)
    real, intent(in) :: airdens(:) ! # molecules / cm3 in each layer
    real, intent(in) :: aircol(:)  ! # molecules / cm2 in each layer
    real, intent(in) :: o2col(:)
    real, intent(in) :: o3col(:)
    real, intent(in) :: so2col(:)
    real, intent(in) :: no2col(:)
    real, intent(out) :: tuv_prates(:,:) ! /sec

    real :: tauaer300(nlev) ! aerosol properties
    real :: tauaer400(nlev)
    real :: tauaer600(nlev)
    real :: tauaer999(nlev)
    real :: waer300(nlev)
    real :: waer400(nlev)
    real :: waer600(nlev)
    real :: waer999(nlev)
    real :: gaer300(nlev)
    real :: gaer400(nlev)
    real :: gaer600(nlev)
    real :: gaer999(nlev)

    real :: albedo(nwave)
    
    real :: dtaer(nlyr,nwave), omaer(nlyr,nwave), gaer(nlyr,nwave)
    real :: dtcld(nlyr,nwave), omcld(nlyr,nwave), gcld(nlyr,nwave)
    real :: dt_cld(nlyr)
    
    real :: qll(nlev) ! cld water content (g/m3)
    real :: cldfrac(nlev)
    real :: o3_xs(nwave,nlev)
    real :: no2_xs(nwave,nlev)
    real :: srb_o2_xs(nwave,nlev)
    real :: radfld(nwave,nlev)
    real :: efld(nlev,nwave)
    real :: e_dir(nlev,nwave)
    real :: e_dn(nlev,nwave)
    real :: e_up(nlev,nwave)
    real :: dir_fld(nlev,nwave)
    real :: dwn_fld(nlev,nwave)
    real :: up_fld(nlev,nwave)
    real :: o3_xs_tpose(nlev,nwave)
    real :: no2_xs_tpose(nlev,nwave)
    
    if( .not. is_full_tuv ) then
       call xs_int( o3_xs, tlev, o3_xs_tab )
       call xs_int( no2_xs, tlev, no2_xs_tab )
    else
       call o3xs( nlev,tlev,nwave+1,wl,o3_xs_tpose )
       call no2xs_jpl06a( nlev,tlev,nwave+1,wl,no2_xs_tpose )
       o3_xs  = transpose( o3_xs_tpose )
       no2_xs = transpose( no2_xs_tpose )
    endif

    qll=0.0
    cldfrac=0.0
    tauaer300=0.0
    tauaer400=0.0
    tauaer600=0.0
    tauaer999=0.0
    waer300=1.0
    waer400=1.0
    waer600=1.0
    waer999=1.0
    gaer300=0.0
    gaer400=0.0
    gaer600=0.0
    gaer999=0.0

    albedo(:) = albedo_in
    
    call tuv_radfld( nlambda_start, cld_od_opt, cldfrac, nlyr, nwave, &
         zenith, zlev, albedo, &
         aircol, o2col, o3col, so2col, no2col, &
         tauaer300, tauaer400, tauaer600, tauaer999, &
         waer300, waer400, waer600, waer999, &
         gaer300, gaer400, gaer600, gaer999, &
         dtaer, omaer, gaer, dtcld, omcld, gcld, &
         has_aer_ra_feedback, &
         qll, dobsi, o3_xs, no2_xs, o2_xs, &
         so2_xs, wl(1), wc, tlev, srb_o2_xs, radfld, efld, &
         e_dir, e_dn, e_up, &
         dir_fld, dwn_fld, up_fld, dt_cld )

    call calc_tuv_prates(1,nlev,nlev, wl,wc, tlev, airdens, radfld, srb_o2_xs, nj, tuv_prates )

  end subroutine tuv_photolysis_calc

end module tuv_photolysis
