C23456789
c	subroutine theta_movement(nth,nr,nz,w2,eta2,w2,uR1,uT1,uR2,dr,dtheta,dt,rfric,fco,difz,g,rho1,bate,batT)
c		implicit none
c		integer j,k,l,nth,nr,nz,lfim
c		real Z,forc,dz,dtheta,dt,dr,difz,difzdz,coriolis,dp,difusao,atrito_fundo,adveccao,umedR
c		real drdzinf,drdzsup,rfric
c		real, dimension(nth,nr,nz) :: w2,uT1,uT2,uR2,rho1
c		real, dimension(nth,nr) :: eta1,eta2
c		do j=2,nr-1
c			if (kmare(k,j)>0) then
c					bat3=[batT(1,j) bate(1,j+1) bate(1,j)]
c					dz=(eta2(1,j)+Z)/nz
c					lfim=min(bat3)/dz-1
c					difzdz=difz/dz
c			    		do l=2,lfim
c						umedR=(uR2(1,j,l)+uR2(1,j+1,l)+uR2(nth,j+1,l)+uR2(nth,j,l))/4
c						drdzinf=(uT1(1,j,l+1)-uT1(1,j,l)).*difzdz
c						drdzsup=(uT1(1,j,l)-uT1(1,j,l-1)).*difzdz
c						dp=g.*(eta2(1,j,1)-eta2(nth,j,1))./(dtheta*R)
c						coriolis = fco*umedR*rho1(1,j,l)
c						atrito_fundo = rfric*uR1(1,j,l)
c						difusao=(drdzinf-drdzsup)/dz(1,j)
c						adveccao = uR2(1,j,l)*((uT1(1,j,l)-uT1(nth,j,l))/dr)+(uT1(1,j,l)/Raio)*((uT1(1,j+1,l)-uT1(1,j,l))/dtheta)-(uT1(1,j,l)*uR1(1,j,l))/Raio+w2(1,j,l)*((uT1(1,j,l)-uT1(1,j,l-1))./dz)
c						forc=coriolis-dp+difusao-atrito_fundo+adveccao
c						uT2(1,j,l)=uT1(1,j,l)*rho1(1,j,l)+forc*dt
c					end do
c				end if
c			do k=2,nth
c				if (kmare(k,j)>0) then
c					bat3=[batT(k,j) bate(k,j+1) bate(k,j)]
c					dz=(eta2(k,j)+Z)/nz
c					lfim=min(bat3)/dz-1
c					difzdz=difz/dz
c			    		do l=2,lfim
c						umedR=(uR2(k,j,l)+uR2(k,j+1,l)+uR2(k-1,j+1,l)+uR2(k-1,j,l))/4
c						drdzinf=(uT1(k,j,l+1)-uT1(k,j,l)).*difzdz
c						drdzsup=(uT1(k,j,l)-uT1(k,j,l-1)).*difzdz
c						dp=g.*(eta2(k,j,1)-eta2(k-1,j,1))./(dtheta*Raio)
c						coriolis = fco*umedR*rho1(k,j,l)
c						atrito_fundo = rfric*uR1(k,j,l)
c						difusao=(drdzinf-drdzsup)/dz(k,j)
c						adveccao = uR2(k,j,l)*((uT1(k,j,l)-uT1(k-1,j,l))/dr)+(uT1(k,j,l)/Raio)*((uT1(k,j+1,l)-uT1(k,j,l))/dtheta)-(uT1(k,j,l)*uR1(k,j,l))/Raio+w2(k,j,l)*((uT1(k,j,l)-uT1(k,j,l-1))./dz)
c						forc=coriolis-dp+difusao-atrito_fundo+adveccao
c						uT2(k,j,l)=uT1(k,j,l)*rho1(k,j,l)+forc*dt
c					end do
c				end if
c			end do
c		end do
c	return
c	end subroutine theta_movement

