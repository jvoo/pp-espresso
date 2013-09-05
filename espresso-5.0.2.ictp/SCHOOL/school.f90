!-----------------------------------------------------------------------
program school
  !-----------------------------------------------------------------------
  !
  ! read in PWSCF data in XML format using IOTK lib
  ! then prepare matrices for GWL calculation
  ! 
  ! input:  namelist "&inputpp", with variables
  !   prefix       prefix of input files saved by program pwscf
  !   outdir       temporary directory where files resides
  !   pp_file      output file. If it is omitted, a directory 
  !                "prefix.export/" is created in outdir and
  !                some output files are put there. Anyway all the data 
  !                are accessible through the "prefix.export/index.xml" file which
  !                contains implicit pointers to all the other files in the
  !                export directory. If reading is done by the IOTK library
  !                all data appear to be in index.xml even if physically it
  !                is not. 
  !   uspp_spsi    using US PP if set .TRUE. writes S | psi > 
  !                and | psi > separately in the output file 
  !   single_file  one-file output is produced
  !   ascii        ....
  !
  !   pseudo_dir   pseudopotential directory
  !   psfile(:)    name of the pp file for each species 
  !    

  use io_files,  ONLY : prefix, tmp_dir, outdir
  use io_files,  ONLY : psfile, pseudo_dir
  use io_global, ONLY : stdout, ionode, ionode_id
  USE mp_global,     ONLY: mp_startup,mpime,kunit
  USE environment,   ONLY: environment_start
  USE mp, ONLY : mp_bcast
  use ldaU, ONLY : lda_plus_u
  use scf, only : vrs, vltot, v, kedtau
  USE fft_base,             ONLY : dfftp
  use pwcom, only : doublegrid, nspin
  use uspp, ONLY : okvan
  use realus, ONLY : qpointlist
  
  implicit none
  character(len=9) :: code = 'SCHOOL'
  integer :: ios, kunittmp
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  character(len=200) :: pp_file
  logical :: uspp_spsi, ascii, single_file, raw


  NAMELIST /inputschool/ prefix

  CALL mp_startup ( )
  CALL environment_start ( code )

  prefix='export'
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  IF ( ionode ) THEN
     CALL input_from_file ( )
     READ(5,inputschool,IOSTAT=ios)
     IF (ios /= 0) CALL errore ('SCHOOL', 'reading inputschool namelist', ABS(ios) )
  endif

  tmp_dir = trimcheck( outdir )
  CALL mp_bcast( outdir, ionode_id )
  CALL mp_bcast( tmp_dir, ionode_id )
  CALL mp_bcast( prefix, ionode_id )

  call read_file

  call openfile_school

#if defined __PARA
  kunittmp = kunit
#else
  kunittmp = 1
#endif

  pp_file= ' '
  uspp_spsi = .FALSE.
  ascii = .FALSE.
  single_file = .FALSE.
  raw = .FALSE.


  call read_export(pp_file,kunittmp,uspp_spsi, ascii, single_file, raw)


  call summary()
  
  CALL print_ks_energies()

  CALL hinit0()
!                                                                                               
  if(lda_plus_u) then
    CALL init_ns()
  endif
  CALL set_vrs(vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid )

  IF ( okvan) CALL qpointlist()

  CALL check_wfcs

  call stop_pp

  stop
end program school



