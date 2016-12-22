C23456789
	program tank_model
		implicit none
		! declaracoes
		external atualize
		external boundaries
		external cont
		external r_movement
c		external theta_movement
		real g,rpm,Raio,Z,rho,difz,difT,rfric,fco,dth,dr,dt,tempo
		integer ntime,nr,nth,nz,fsave,pnumber,kplot
		integer ii,jj,n
		real, allocatable, dimension(:,:,:) :: w1,w2,uT1,uT2,uR1,uR2,rho1
		real, allocatable, dimension(:,:) :: bat,batR,batT,kmarT,kmare
		real, allocatable, dimension(:,:) :: eta1,eta2,kmarR,bate
		character(len=1024) :: fname

		! ler inputs
		open(5, file = "params.txt")
		read(5,*) g,dt,ntime,nr,nth,nz,difz,difT,rfric,fsave,rpm,rho,Raio
		print *,g
		close(5)

		! alocar matrizes
		allocate(w1(nth,nr,nz))
		allocate(uR1(nth,nr,nz))
		allocate(uT1(nth,nr,nz))
		allocate(rho1(nth,nr,nz))
		allocate(eta1(nth,nr))
		allocate(w2(nth,nr,nz))
		allocate(uR2(nth,nr,nz))
		allocate(uT2(nth,nr,nz))
		allocate(eta2(nth,nr))
		allocate(bat(nth+1,nr+1))
		allocate(bate(nth+1,nr+1))
		allocate(batT(nth,nr))
		allocate(batR(nth,nr))
		allocate(kmarT(nth,nr))
		allocate(kmarR(nth,nr))
		allocate(kmare(nth,nr))

		! rotacao
		fco = rpm/60  

		! deltas R  e theta
		dth = 360/nth
		dr = Raio/nr

		! zerar matrizes
		w1(:,:,:) = 0
		uR1(:,:,:) = 0
		uT1(:,:,:) = 0
		rho1(:,:,:)=0
		eta1(:,:)=0
		w2(:,:,:) = 0
		uR2(:,:,:) = 0
		uT2(:,:,:) = 0
		eta2(:,:)=0
		! forcante de rotacao
		uT1(:,nr,:)=fco
		uT2(:,nr,:)=fco

		!ler batimetria
		open(11,file='dep.dep')
		do  jj=1,nth+1,1
			read(11,*) bat(jj,:)
		end do
		close(11)

		do jj=1,nth
			do ii=1,nr
				bate(jj,ii) = bat(jj,ii)
				kmare(jj,ii)=abs(bate(jj,ii))/bate(jj,ii)
				batT(jj,ii)=(bat(jj,ii)+bat(jj+1,ii))/2
				kmarR(jj,ii)=abs(batT(jj,ii))/batT(jj,ii)
				batR(jj,ii)=(bat(jj,ii)+bat(jj,ii+1))/2
				kmarR(jj,ii)=abs(batR(jj,ii))/batR(jj,ii)
			end do
		end do


		! loop no tempo
		kplot=0
		do n=2,ntime
			tempo=n*dt
			kplot=kplot+1
			!Eq da Continuidade
			call cont(nth,nr,nz,w1,eta1,eta2,w2,uR1,uT1,dr,dth,dt,kmare)
			call r_movement(nth,nr,nz,w2,eta2,w2,uR1,uT1,uR2,dr,dth,dt, &
				rfric,fco,difz,g,rho1,bate,batR,Raio)
c			call theta_movement(nth,nr,nz,w2,eta2,w2,uR1,uT1,uR2,dr,dth,dt,rfric,fco,difz,g,rho1,bate,batT)
			call boundaries(fco,eta2,w2,uR2,uT2,nr,kmare,nth,nz)
			call atualize(nth,nr,nz,w2,w1,uT1,uT2,uR1,uR2,eta1,eta2)
			
			if (kplot==fsave) then
				pnumber=pnumber+1
				write(fname,"(A5,I2)") "eta_", pnumber			
				open(unit = 10, status='replace',form='unformatted',file=fname)
				write(10)  ! write the data in array x to the file
				close(10) ! close the file
				write(fname,"(A5,I2)") "v_", pnumber			
				open(unit = 11, status='replace',form='unformatted',file=fname)
				write(11)  ! write the data in array x to the file
				close(11) ! close the file
				write(fname,"(A5,I2)") "u_", pnumber			
				open(unit = 12, status='replace',form='unformatted',file=fname)
				write(12)  ! write the data in array x to the file
				close(12) ! close the file
				write(fname,"(A5,I2)") "w_", pnumber			
				open(unit = 15, status='replace',form='unformatted',file=fname)
				write(15)  ! write the data in array x to the fil
				close(15) ! close the file
				kplot=0
			end if
		end do
	end program tank_model
