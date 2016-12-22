csajhsakfjks
	subroutine boundaries(fco,eta2,w2,uR2,uT2,nr,kmare,nth,nz)
		implicit none
		integer k,nr,nth,nz
		real fco
		real, dimension(nth,nr,nz) :: w2,uT2,uR2
		real, dimension(nth,nr) :: eta2,kmare
		! Condicoes de contorno com extrapolacao linear
		uT2(:,nr,:)=fco
		do k=1,nth
			eta2(k,1)=(2*eta2(k,2)-eta2(k,3))*kmare(k,1)
			eta2(k,nr)=(2*eta2(k,nr-1)-eta2(k,nr-2))*kmare(k,nr)
			w2(k,1,2:nz-1)=(2*w2(k,2,2:nz-1)-w2(k,3,2:nz-1))*kmare(k,1)
			w2(k,nr,2:nz-1)=(2*w2(k,nr-1,2:nz-1)-w2(k,nr-2,2:nz-1))
			w2(k,nr,2:nz-1)=w2(k,nr,2:nz-1)*kmare(k,nr)
			uR2(k,1,2:nz-1)=(2*uR2(k,2,2:nz-1)-uR2(k,3,2:nz-1))*kmare(k,1)
			uR2(k,nr,2:nz-1)=(2*uR2(k,nr-1,2:nz-1)-uR2(k,nr-2,2:nz-1))
			uR2(k,nr,2:nz-1)=uR2(k,nr,2:nz-1)*kmare(k,nr)
		end do   
	return
	end subroutine boundaries
