C
      subroutine r_movement(nth,nr,nz,w2,eta2,w2,uR1,uT1,uR2,dr,dth,dt
&,rfric,fco,difz,g,rho1,bate,batR,Raio)
		implicit none
		integer j,k,l,nth,nr,nz,lfim
		real Z,forc,dz,dth,dt,dr,difz,difzdz,coriolis,dp
		real difusao,atrito_fundo,adveccao,umedT,Raio
		real drdzinf,drdzsup,rfric
		real, dimension(nth,nr,nz) :: w2,uT1,uR1,uR2,rho1
		real, dimension(nth,nr) :: eta1,eta2,batR,bate
		real, dimension(3) :: bat3
		do j=2,nr-1
			if (kmare(1,j)>0) then
					bat3=[batR(1,j),bate(nth,j),bate(1,j)]
					dz=(eta2(1,j)+Z)/nz
					lfim=min(bat3)/dz-1
					difzdz=difz/dz
			    		do l=2,lfim
						umedT=(uT1(1,j,l)+uT1(1,j+1,l)+uT1(nth,j+1,l)+uT1(nth,j,l))/4
						drdzinf=(uR1(1,j,l+1)-uR1(1,j,l))*difzdz
						drdzsup=(uR1(1,j,l)-uR1(1,j,l-1))*difzdz
						dp=g*(eta2(1,j)-eta2(nth,j))/dr
						coriolis = fco*umedT*rho1(1,j,l)
						atrito_fundo = rfric*uR1(1,j,l)
						difusao=(drdzinf-drdzsup)/dz(1,j)
     adveccao = uR1(1,j,l)*((uR1(1,j,l)-uR1(nth,j,l))/dr)+
 &(uT1(1,j,l)/Raio)*((uR1(1,j+1,l)-uR1(1,j,l))/dth)-(uT1(1,j,l)^2)/Raio+
 &w2(1,j,l)*((uR1(1,j,l)-uR1(1,j,l-1))/dz)
						forc=coriolis-dp+difusao-atrito_fundo+adveccao
						uR2(1,j,l)=uR1(1,j,l)*rho1(1,j,l)+forc*dt
					end do
				end if
			do k=2,nth
				if (kmare(k,j)>0) then
					bat3=[batR(k,j),bate(k-1,j),bate(k,j)]
					dz=(eta2(k,j)+Z)/nz
					lfim=min(bat3)/dz-1
					difzdz=difz/dz
			    		do l=2,lfim
						umedT=(uT1(k,j,l)+uT1(k,j+1,l)+uT1(k-1,j+1,l)+uT1(k-1,j,l))/4
						drdzinf=(uR1(k,j,l+1)-uR1(k,j,l))*difzdz
						drdzsup=(uR1(k,j,l)-uR1(k,j,l-1))*difzdz
						dp=g*(eta2(k,j)-eta2(k-1,j))/dr
						coriolis = fco*umedT*rho1(k,j,l)
						atrito_fundo = rfric*uR1(k,j,l)
						difusao=(drdzinf-drdzsup)/dz(k,j)
		adveccao = uR1(k,j,l)*((uR1(k,j,l)-uR1(k-1,j,l))/dr)+(uT1(k,j,l)/Raio)*((uR1(k,j+1,l)-uR1(k,j,l))/dth)-(uT1(k,j,l)^2)/Raio+w2(k,j,l)*((uR1(k,j,l)-uR1(k,j,l-1))/dz)
						forc=coriolis-dp+difusao-atrito_fundo+adveccao
						uR2(k,j,l)=uR1(k,j,l)*rho1(k,j,l)+forc*dt
					end do
				end if
			end do
		end do
	return
	end subroutine r_movement

