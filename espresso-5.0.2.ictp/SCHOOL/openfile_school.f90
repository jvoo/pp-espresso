!
!----------------------------------------------------------------------------
SUBROUTINE openfile_school()
  !----------------------------------------------------------------------------
  !
  ! ... This routine opens all files needed to the self consistent run,
  ! ... sets various file names, units, record lengths
  !
  USE kinds,          ONLY : DP
  USE wvfct,          ONLY : nbnd, npwx, igk
  use control_flags,  ONLY:  twfcollect
  USE io_files,       ONLY : prefix, iunwfc, nwordwfc, iunat, iunsat, nwordatwfc, diropn,iunigk,seqopn
  USE noncollin_module, ONLY : npol
  USE ldaU,             ONLY : lda_plus_u
  USE basis,            ONLY : natomwfc
  USE ions_base,        ONLY : nat, ityp
  USE noncollin_module,   ONLY : noncolin
  USE uspp_param,         ONLY : n_atom_wfc
  !
  IMPLICIT NONE
  !
  LOGICAL       :: exst
  !
  !
  twfcollect=.false.
  !
  ! ... nwordwfc is the record length for the direct-access file
  ! ... containing wavefunctions
  !
  nwordwfc = nbnd * npwx * npol
  !
  CALL diropn( iunwfc, 'wfc', 2*nwordwfc, exst )
  !
  IF ( .NOT. exst ) THEN
     call errore ('openfile_school','file '//TRIM( prefix )//'.wfc'//' not found',1)     
  END IF

  CALL seqopn( iunigk, 'igk', 'UNFORMATTED', exst )

  ! ... Needed for LDA+U
  !
  ! ... iunat  contains the (orthogonalized) atomic wfcs 
  ! ... iunsat contains the (orthogonalized) atomic wfcs * S
  ! ... iunocc contains the atomic occupations computed in new_ns
  ! ... it is opened and closed for each reading-writing operation  
  !
  natomwfc = n_atom_wfc( nat, ityp, noncolin )
  nwordatwfc = 2*npwx*natomwfc*npol
  !
  IF ( lda_plus_u ) then
     CALL diropn( iunat,  'atwfc',  nwordatwfc, exst )
     IF ( .NOT. exst ) THEN
        call errore ('openfile_school','file '//TRIM( prefix )//'.atwfc'//' not found',1)
     END IF

     CALL diropn( iunsat, 'satwfc', nwordatwfc, exst )
     IF ( .NOT. exst ) THEN
        call errore ('openfile_school','file '//TRIM( prefix )//'.satwfc'//' not found',1)
     END IF
  END IF
  !

  RETURN
  !
END SUBROUTINE openfile_school
