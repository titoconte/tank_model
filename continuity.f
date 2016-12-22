C23456789
	subroutine cont(nth,nr,nz,w1,eta1,eta2,w2,uR1,uT1,dr,dth,dt,kmare)
		implicit none
		integer j,k,l,nth,nr,nz
		real Z,forcR,forcT,dz,dth,dt,dr,forc
		real, dimension(nth,nr,nz) :: w1,w2,uT1,uR1,rho1
		real, dimension(nth,nr) :: eta1,eta2,kmare,bate
		do j=2,nr-1
			do k=1,nth-1
				if (kmare(k,j)>0) then
					forcR=0
					forcT=0
					dz=(eta1(k,j)+bate(k,j))/nz
					do l=2,nz-1
						forcR=(uR1(k,j,l)-uR1(k,j-1,l))/dr+forcR 
						forcT=(uT1(k+1,j,l)-uT1(k,j,l))/dth+forcT
						forc =-dz*(forcR+forcT)/dr
						w2(k,j,l)=w1(k,j,l)+forc
					end do
					eta2(k,j) = eta1(k,j)-w2(k,j,1)*dt
				end if
			end do
			if (kmare(nth,j)>0) then
				forcR=0
				forcT=0
				dz=(eta1(nth,j)+bate(k,j))/nz
				do l=2,nz-1
					forcR=(uR1(nth,j,l)-uR1(nth,j-1,l))/dr+forcR 
					forcT=(uT1(1,j,l)-uT1(nth,j,l))/dth+forcT
					forc=-dz*(forcR+forcT)/dr
					w2(nth,j,l)=w1(nth,j,l)+forc
				end do
				eta2 = eta1-w2(nth,j,1)*dt
			end if
		end do
	return
	end subroutine cont
