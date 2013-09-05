
subroutine check_wfcs_R
    USE fft_base, ONLY			: dffts
    USE fft_interfaces, ONLY		: fwfft, invfft


    USE io_global, ONLY 		: ionode
    USE io_files, ONLY 			: iunwfc, nwordwfc
    USE kinds, ONLY 			: DP
    USE control_flags, ONLY		: gamma_only
    USE wvfct, ONLY			: nbnd, npw, npwx
    USE gvect, ONLY			: gstart
    USE mp, ONLY			: mp_sum


    USE wavefunctions_module, ONLY	: evc,psic
    USE gvecs, ONLY			: nls, nlsm


    implicit none

    real(kind=DP), allocatable		:: omat(:,:), rwfcs(:,:)
    integer				:: i, j
    real(kind=DP)			:: offmax


    allocate(omat(nbnd,nbnd))

    if(gamma_only) then
	!read in wfcs
	CALL davcio(evc, 2*nwordwfc, iunwfc, 1, -1)
	!do inner products
	CALL dgemm('T', 'N', nbnd, nbnd, 2*npw, 2.d0, evc, 2*npwx, evc, 2*npwx, 0.d0, omat, nbnd)
	if(gstart==2) then
	    do i=1,nbnd
		do j=1,nbnd
		    omat(i,j)=omat(i,j)-dble(evc(1,i)*conjg(evc(1,j)))
		enddo
	    enddo
	endif
	call mp_sum(omat)
    endif

    if(ionode) then
	do i=1,nbnd
	    WRITE(*,*) i, omat(i,i)
	enddo
	offmax = 0d0
	do i=1,nbnd
	    do j=1,nbnd
		if((i.ne.j).and.(abs(omat(i,j))>offmax)) offmax = abs(omat(i,j))
	    enddo
	enddo
	WRITE(*,*) 'maximum off-diagonal element =', offmax
    endif

    deallocate(omat)

end subroutine check_wfcs_R
