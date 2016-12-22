%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modelo Numerico hidrodinamico 3D linear (decaimento explicito) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc

%Configuracao da grade (pontos tipo eta, u, v)
%NORTE (k)
%    V V V V V
%4  UEUEUEUEUE
%    V V V V V
%3  UEUEUEUEUE
%    V V V V V
%2  UEUEUEUEUE
%    V V V V V
%1  UEUEUEUEUE
%    1 2 3 4 5 LESTE (j)

% desconsiderar primeira coluna e ultima  linha da grade acima
% de modo a ter ponto tipo eta em todo o contorno

%Sistema Internacional de unidades (SI)

% Constantes iniciais
g=9.8;          % aceleracao da gravidade
dt=1;
ntime = 1000;   % passos de tempo
nr=100;        % espacamento da grade em x
ntheta=72;     % espacamento da grade em y
nz=10;          % espacamento de grade em z
difz=0.001;     % coeficiente de difusao na vertical
difT=0.001;	% difusão temperatura
rfric=0.002;    % coeficiente de friccao no fundo
freqplot=20;    % frequencia de plotagem
dens=1000;      % densidade media da agua do mar
rpm=10;         % rotação
R = 0.2215;     % raio
Z = 0.1;        % profundidade  
rho=1000;

% roação
fco = rpm/60;  

% deltas
dtheta = 360/ntheta;
dr = R/nr;


%tamanho da grade e batimetria (chaves para pontos maritimos)
nr2=nr*2;
ntheta2=ntheta*2;
bat(1:ntheta2,1:nr2)=Z;

bat(:,1)=0;bat(:,2)=0;
bat(:,nr2)=0;bat(:,nr2-1)=0;

bat(1,:)=0; bat(2,:)=0;
bat(ntheta2,:)=0;bat(ntheta2-1,:)=0;

bate=bat(1:2:end-1,2:2:end);
batR=bat(1:2:end-1,1:2:end-1);
batT=bat(2:2:end,2:2:end);
kmare=bate*0;
kmare(bate>0)=1;
kmarR=batR*0;
kmarR(batR>0)=1;
kmarT=batT*0;
kmarT(batT>0)=1;

% Condicoes iniciais de repouso (1 - valores atuais, 2 - renovados)
w1=zeros(ntheta,nr,nz);
uR1=zeros(ntheta,nr,nz);
uT1=zeros(ntheta,nr,nz);
eta1=zeros(ntheta,nr,1);
rho1=zeros(ntheta,nr,nz)*rho;
%T1=zeros(ntheta,nr,nz)*T0;
w2=zeros(ntheta,nr,nz);
uR2=zeros(ntheta,nr,nz);
uT2=zeros(ntheta,nr,nz);
eta1=zeros(ntheta,nr,1);
%T2=zeros(ntheta,nr,nz)*T0;
uT1(:,nr,:)=fco;
uT2(:,nr,:)=fco;
% %Condicoes de vento e calculo das tensoes de cisalhamento na superficie
% dens_ar=1.025;
% fric=2.6*1E-3;
% 
% uwind(1:ntheta,1:nr)=7.0711;
% vwind(1:ntheta,1:nr)=7.0711;
% wwind=sqrt(uwind.^2+vwind.^2);
% taux=fric*dens_ar.*uwind.*wwind;
% tauy=fric*dens_ar.*vwind.*wwind;
% 
% vento=sqrt(uwind.^2+vwind.^2);
% ventoma=max(vento);
% ventomax=max(ventoma);

% %plotando os ventos
% figure(2)
% quiver(X,Y,uwind,vwind,'LineWidth',2)
% title(['Vento - intensidade maxima ',...
%     num2str(ventomax),' m/s'],'fontsize',12)
% axis equal
% axis([xgrid(1) xgrid(nr) ygrid(1) ygrid(ntheta)])
% xlabel('DISTANCIA (m) EW','fontsize',12)
% ylabel('DISTANCIA (m) NS','fontsize',12)
%print -djpeg fig_vento

%Loop no tempo
kplot=1;
for n=2:ntime
   tempo=n*dt;
   kplot=kplot+1;
    %Eq da Continuidade
    for j=2:nr-1
        for k=2:ntheta-1
           if kmare(k,j)>0
               forcR=0;
               forcT=0;
	       dz=(eta1+Z)/nz;
               for l=2:nz-1
               forcR=(uR1(k,j+1,l)-uR1(k,j,l))/dr+forcR; 
               forcT=(uT1(k,j,l)-uT1(k-1,j,l))/dtheta+forcT;
               w2(k,j,l)=w1(k,j,l)-dz(k,j)*(forcR+forcT)/dr;
               eta2 = eta1-w2(k,j,1)*dt;
               end
           end
        end
    end
  
  %Eq. do movimento em r
    for j=2:nr-1
        for k=2:ntheta-1
            if (kmarR(k,j)*kmare(k,j)*kmare(k,j-1))>0
		    bat3=[batR(k,j) bate(k,j) bate(k,j-1)];
		    batmin=min(bat3);
		    lfim=batmin/dz(k,j)-1;
		    difzdz=difz/dz(k,j);
		    for l=2:lfim
		    	umedT=(uT1(k,j,l)+uT1(k,j-1,l)+uT1(k+1,j-1,l)+uT1(k+1,j,l))/4;
		    	drdzinf=(uR1(k,j,l+1)-uR1(k,j,l)).*difzdz;
		    	drdzsup=(uR1(k,j,l)-uR1(k,j,l-1)).*difzdz;
			dp=g.*(eta2(k,j,1)-eta2(k,j-1,1))./dr;
			coriolis = fco*umedT*rho1(k,j,l);
			atrito_fundo = rfric*uR1(k,j,l);
			difusao=(drdzinf-drdzsup)/dz(k,j);
			%adveccao = uR1(k,j,l)*((uR1(k,j,l)-uR1(k,j-1,l))/dr)...
			%+(uT1(k,j,l)/R)*((uR1(k,j,l)-uR1(k,j-1,l))/dtheta)-(uT1(k,j,l).^2)/R...
			%+w2(k,j,l)*((uR1(k,j,l)-uR1(k,j-1,l))./dzdz(k,j))
		    	forc=coriolis-dp+difusao-atrito_fundo;%+adveccao;
		    	uR2(k,j,l)=uR1(k,j,l)*rho1(k,j,l)+forc*dt;
		    end
            end
            
         end
    end
% Movimento em theta
    for j=2:nr-1
        for k=2:ntheta-1
			if (kmarT(k,j)*kmare(k,j)*kmare(k,j-1))>0
		    bat3=[batT(k,j) bate(k,j) bate(k,j-1)];
		    batmin=min(bat3);
		    lfim=batmin/dz(k,j)-1;
		    difzdz=difz/dz(k,j);
		    for l=2:lfim
		    	umedR=(uR2(k,j,l)+uR2(k,j-1,l)+uR2(k+1,j-1,l)+uR2(k+1,j,l))/4;
		    	drdzinf=(uT1(k,j,l+1)-uT1(k,j,l)).*difzdz;
		    	drdzsup=(uT1(k,j,l)-uT1(k,j,l-1)).*difzdz;
			dp=g.*(eta2(k+1,j,1)-eta2(k,j,1))./(dtheta*R);
			coriolis = fco*umedR*rho1(k,j,l);
			atrito_fundo = rfric*uR1(k,j,l);
			difusao=(drdzinf-drdzsup)/dz(k,j);
			%adveccao = uR2(k,j,l)*((uT1(k+1,j,l)-uT1(k,j,l))/dr)...
			%+(uT1(k,j,l)/R)*((uT1(k+1,j,l)-uT1(k,j,l))/dtheta)-(uT1(k,j,l)*uR1(k,j,l))/R...
			%+w2(k,j,l)*((uT1(k+1,j,l)-uT1(k,j,l))./dz(k,j));

		    	forc=coriolis-dp+difusao-atrito_fundo;%+adveccao;

		    	uT2(k,j,l)=uT1(k,j,l)*rho1(k,j,l)+forc*dt;

		    end
            end
end
end


% Condicoes de contorno com extrapolacao linear
for j=1:nr
    eta2(1,j,1)=(2*eta2(2,j,1)-eta2(3,j,1))*kmare(1,j);
    eta2(ntheta,j,1)=(2*eta2(ntheta-1,j,1)-eta2(ntheta-2,j,1))*kmare(ntheta,j);
    w2(1,j,2:nz-1)=(2*w2(2,j,2:nz-1)-w2(3,j,2:nz-1))*kmare(1,j);
    w2(ntheta,j,2:nz-1)=(2*w2(ntheta-1,j,2:nz-1)-w2(ntheta-2,j,2:nz-1))*kmare(ntheta,j);
    uR2(1,j,2:nz-1)=(2*uR2(2,j,2:nz-1)-uR2(3,j,2:nz-1))*kmarR(1,j);
    %uT2(1,j,2:nz-1)=(2*uT2(2,j,2:nz-1)-uT2(3,j,2:nz-1))*kmarT(1,j);
    uR2(ntheta,j,2:nz-1)=(2*uR2(ntheta-1,j,2:nz-1)-uR2(ntheta-2,j,2:nz-1))*kmarR(ntheta,j);
    %uT2(ntheta,j,2:nz-1)=(2*uT2(ntheta-1,j,2:nz-1)-uT2(ntheta-2,j,2:nz-1))*kmarT(ntheta,j);
end
for k=1:ntheta
    eta2(k,1,1)=(2*eta2(k,2,1)-eta2(k,3,1))*kmare(k,1);
    eta2(k,nr,1)=(2*eta2(k,nr-1,1)-eta2(k,nr-2,1))*kmare(k,nr);
    w2(k,1,2:nz-1)=(2*w2(k,2,2:nz-1)-w2(k,3,2:nz-1))*kmare(k,1);
    w2(k,nr,2:nz-1)=(2*w2(k,nr-1,2:nz-1)-w2(k,nr-2,2:nz-1))*kmare(k,nr);
    uR2(k,1,2:nz-1)=(2*uR2(k,2,2:nz-1)-uR2(k,3,2:nz-1))*kmarR(k,1);
    uR2(k,nr,2:nz-1)=(2*uR2(k,nr-1,2:nz-1)-uR2(k,nr-2,2:nz-1))*kmarR(k,nr);
    if ntheta>k
    uT2(k,1,2:nz-1)=(2*uT2(k,2,2:nz-1)-uT2(k,3,2:nz-1))*kmarT(k,1);
    uT2(k,nr,2:nz-1)=(2*uT2(k,nr-1,2:nz-1)-uT2(k,nr-2,2:nz-1))*kmarT(k,nr);
    end
end
    
% Renovando as variaveis no tempo
w1=w2;
uR1=uR2;
uT1=uT2;
eta1=eta2;
%rho1 = 1000*(1-((T2-288.9114)/(5089288*(T-68.129630)))*((T2-3.9863).^2))	
n
max(max(uR1))
end
   
