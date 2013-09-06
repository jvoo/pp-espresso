subroutine check_wfcs
  
  USE io_files,             ONLY : prefix, iunwfc, nwordwfc, iunigk, find_free_unit
  USE kinds,    ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : evc,psic
  USE wvfct,                ONLY : nbnd,npw,npwx,igk
  USE gvect, ONLY : gstart
  USE gvecs, ONLY : nls,nlsm
  USE mp, ONLY : mp_sum, mp_barrier
  USE io_global, ONLY : stdout,ionode
  USE fft_base,         ONLY : dffts
  USE fft_interfaces,ONLY : fwfft, invfft
  USE klist,                ONLY : nks,ngk,xk
  USE io_files, ONLY : iunigk
  USE lsda_mod,   ONLY : lsda, nspin,current_spin,isk
  USE becmod,        ONLY : bec_type, becp, calbec,allocate_bec_type, deallocate_bec_type 
  USE uspp,     ONLY : nkb, vkb, becsum, nhtol, nhtoj, indv, okvan, qq
  USE uspp_param, ONLY : upf, nh
  USE ions_base,  ONLY : nat, nsp, ityp

  implicit none

  real(kind=DP), allocatable :: omat(:,:)
  real(kind=DP), allocatable :: rwfcs(:,:)
  integer :: i,j
  integer :: ik, ios
  complex(kind=DP), allocatable :: zmat(:,:)
  complex(kind=DP), allocatable :: cwfcs(:,:)
  integer :: is
  INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0

  allocate(omat(nbnd,nbnd))
  if(gamma_only) then

!loop on spin
     do is=1,nspin
!read in wfcs
        CALL davcio(evc,2*nwordwfc,iunwfc,is,-1)
!US part 
        if(okvan) then
           CALL allocate_bec_type (nkb,nbnd,becp)
            IF ( nkb > 0 ) &
                 CALL init_us_2( npw, igk, xk(1,ik), vkb )
            CALL calbec( npw, vkb, evc, becp )
        endif
!do inner products
        call dgemm('T','N',nbnd,nbnd,2*npw,2.d0,evc,2*npwx,evc,2*npwx,0.d0,omat,nbnd)
        if(gstart==2) then
           do i=1,nbnd
              do j=1,nbnd
                 omat(i,j)=omat(i,j)-dble(evc(1,i)*conjg(evc(1,j)))
              enddo
           enddo
        endif
        call mp_sum(omat)
!US part
        if(okvan) then
           ijkb0 = 0
           DO nt = 1, nsp
              IF ( upf(nt)%tvanp ) THEN
                 DO na = 1, nat
                    IF ( ityp(na) == nt ) THEN
                       
                       DO jh = 1, nh(nt)
                          jkb = ijkb0 + jh
                          DO ih = 1, nh(nt)
                             ikb = ijkb0 + ih
                             do i=1,nbnd
                                do j=1,nbnd
                                   omat(i,j)=omat(i,j)+qq(ih,jh,nt)*becp%r(ikb,i)*becp%r(jkb,j)
                                enddo
                             enddo
                             
                          END DO
                       END DO
                       ijkb0 = ijkb0 + nh(nt)
                    END IF
                 END DO
              ELSE
                 DO na = 1, nat
                    IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
                 END DO
              END IF
           END DO
          
        endif

!
        if(ionode) then
           do i=1,nbnd
              do j=1,nbnd
                 write(stdout,*) 'CHECK_WFCS IS',is,i,j,omat(i,j)
              enddo
           enddo
           call flush_unit(stdout)
        endif
!do the same in real space
        allocate(rwfcs(dffts%nnr,nbnd))
        do i = 1, nbnd, 2
           psic(1:dffts%nnr)=0.d0
           if ( i < nbnd ) then
         !                                                                               
                ! ... two ffts at the same time                                                 
                !                                                                               
              psic(nls(1:npw))  = evc(1:npw,i) + &
                   ( 0.D0, 1.D0 ) * evc(1:npw,i+1)
              psic(nlsm(1:npw)) = CONJG( evc(1:npw,i) - &
                   ( 0.D0, 1.D0 ) * evc(1:npw,i+1) )
              !                                                                               
           else
              !                                                                                       
              psic(nls(1:npw))  = evc(1:npw,i)
              psic(nlsm(1:npw)) = CONJG( evc(1:npw,i) )
              !                                                                                    
           endif
           CALL invfft ('Wave', psic, dffts)
           rwfcs(1:dffts%nnr,i)= DBLE(psic(1:dffts%nnr))
           if(i/=nbnd)  rwfcs(1:dffts%nnr,i+1)= DIMAG(psic(1:dffts%nnr))

        enddo
      
        call dgemm('T','N',nbnd,nbnd,dffts%nnr,1.d0,rwfcs,dffts%nnr,rwfcs,dffts%nnr,0.d0,omat,nbnd)
        call mp_sum(omat)
        omat(1:nbnd,1:nbnd)=omat(1:nbnd,1:nbnd)/(dffts%nr1*dffts%nr2*dffts%nr3)
        if(ionode) then
           do i=1,nbnd
              do j=1,nbnd
                 write(stdout,*) 'CHECK_RWFCS',i,j,omat(i,j)
              enddo
           enddo
           call flush_unit(stdout)
        endif


        deallocate(rwfcs)
     
     if(okvan) CALL deallocate_bec_type ( becp )
  enddo
  else
     allocate(zmat(nbnd,nbnd))
     if (nks>1) rewind (unit = iunigk)
     if(okvan)  CALL allocate_bec_type (nkb,nbnd,becp)
     do ik=1,nks
        if (lsda) current_spin = isk (ik)
        npw = ngk (ik)
        IF ( nks > 1 ) READ( iunigk ) igk
       
        call davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)
        !US part
        if(okvan) then
           IF ( nkb > 0 ) &
                CALL init_us_2( npw, igk, xk(1,ik), vkb )
           CALL calbec( npw, vkb, evc, becp )
        endif

        call zgemm('C','N',nbnd,nbnd,npw,(1.d0,0.d0),evc,npwx,evc,npwx,(0.d0,0.d0),zmat,nbnd)
        call mp_sum(zmat)
        if(okvan) then
           ijkb0 = 0
           DO nt = 1, nsp
              IF ( upf(nt)%tvanp ) THEN
                 DO na = 1, nat
                    IF ( ityp(na) == nt ) THEN

                       DO jh = 1, nh(nt)
                          jkb = ijkb0 + jh
                          DO ih = 1, nh(nt)
                             ikb = ijkb0 + ih
                             do i=1,nbnd
                                do j=1,nbnd
                                   zmat(i,j)=zmat(i,j)+qq(ih,jh,nt)*conjg(becp%k(ikb,i))*becp%k(jkb,j)
                                enddo
                             enddo

                          END DO
                       END DO
                       ijkb0 = ijkb0 + nh(nt)
                    END IF
                 END DO
              ELSE
                 DO na = 1, nat
                    IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
                 END DO
              END IF
           END DO

        endif

        if(ionode) then
           do i=1,nbnd
              do j=1,nbnd
                 write(stdout,*) 'CHECK_WFCS is',current_spin,ik,i,j,zmat(i,j)
              enddo
           enddo
           call flush_unit(stdout)
        endif


     enddo
     allocate(cwfcs(1:dffts%nnr,nbnd))
     if (nks>1) rewind (unit = iunigk)
     do ik=1,nks
        npw = ngk (ik)
        IF ( nks > 1 ) READ( iunigk ) igk
        call davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)
        do i=1,nbnd
           psic(1:dffts%nnr)=0.d0
           psic(nls(igk(1:npw)))=evc(1:npw,i)
           CALL invfft ('Wave', psic, dffts)
           cwfcs(1:dffts%nnr,i)=psic(1:dffts%nnr)
        enddo
        call zgemm('C','N',nbnd,nbnd,dffts%nnr,(1.d0,0.d0),cwfcs,dffts%nnr,cwfcs,dffts%nnr,(0.d0,0.d0),zmat,nbnd)
        call mp_sum(zmat)
        zmat(1:nbnd,1:nbnd)=zmat(1:nbnd,1:nbnd)/(dffts%nr1*dffts%nr2*dffts%nr3)
        if(ionode) then
           do i=1,nbnd
              do j=1,nbnd
                 write(stdout,*) 'RCHECK_WFCS IS',current_spin, ik,i,j,zmat(i,j)
              enddo
           enddo
           call flush_unit(stdout)
        endif
     enddo

     deallocate(cwfcs)
     deallocate(zmat)
  endif
  deallocate(omat)
  if(okvan)CALL deallocate_bec_type (becp)
  return
end subroutine check_wfcs
