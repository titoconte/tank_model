C23456789
	subroutine atualize(nth,nr,nz,w2,w1,uT1,uT2,uR1,uR2,eta1,eta2)
		implicit none
		integer nth,nr,nz
		real, dimension(nth,nr,nz) :: w1,w2,uT1,uT2
		real, dimension(nth,nr,nz) :: uR1,uR2
		real, dimension(nth,nr) :: eta1,eta2
		w1=w2
		uR1=uR2
		uT1=uT2
		eta1=eta2
	return
	end subroutine atualize

