clear;
clf;
global g;

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
% mise en place
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
a=-1;
b=1; %bornes du domaine
g=9.81;

% maillage en espace
J=20; %nombre de mailles en espace
dx=(b-a)/J; %taille du pas d espace

x=a+[0:J]*dx; %definition du maillage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Topographies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fond plat
%ZB=zeros([0:J]);

% fond plat avec bosse parabolique
%ZB=zeros([0:J]);
%for i=1:J+1 % une boucle definit directement un vecteur colonne
%	if (x(i)>0.5 & x(i)<1.5)
%		ZB(i)=-(x(i)-0.5)*(x(i)-1.5);
%	end
%end

% fond plat avec bosse carrée
%ZB=zeros([0:J]);
%for i=1:J+1 % une boucle definit directement un vecteur colonne
%if (x(i)>0.98 & x(i)<1.02)
%	ZB(i)=ZB(i)*(1.+.05);
%end
%end

% pente constante (négative)
%S = -.05;
%ZB=S*(x-b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
%Condition Initiale (doit fournir un vecteur *colonne*)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
% Reversed dam break
ul=0;ur=0;hl=1;hr=2;
for i=1:J+1 % une boucle definit directement un vecteur colonne
	if (x(i)<0)
		U0(i,1)=hl;
		U0(i,2)=ul;
		fPrimeInit(i)=Fprime(U0(i,:));
		Froude(i)=Fr(U0(i,:));
	else
		U0(i,1)=hr;
		U0(i,2)=ur;
		fPrimeInit(i)=Fprime(U0(i,:));
		Froude(i)=Fr(U0(i,:));
	end
end

% lac au repos
% for i=1:J+1 % une boucle definit directement un vecteur colonne
% 	U0(i,1)=1.;	%hauteur h+z = cte
% 	U0(i,2)=0.;	% au repos
% 	fPrimeInit(i)=Fprime(U0(i,:));
% 	Froude(i)=Fr(U0(i,:));
% end

% "caillou"
% for i=1:J+1 % une boucle definit directement un vecteur colonne
% 	U0(i,1)=1.;	%hauteur
% 	U0(i,2)=1.;	%fluvial
% %	U0(i,2)=5.;	%torrentiel
% 	if (x(i)>0.98 & x(i)<1.02)
% 		U0(i,1)=U0(i,1)*(1.+.05);
% 	end
% 	fPrimeInit(i)=Fprime(U0(i,:));
% 	Froude(i)=Fr(U0(i,:));
% end

%u0=sin(2*%pi*x);  %definit un vecteur ligne
%u0=u0'; 	% pour passer en colonne


%plot2d(x,ZB,7);
%halt;

% vélocité max initiale pour le calcul de la CFL
c=max(abs(fPrimeInit));

sigma=.5; %CFL  sigma=0.1;
dt=sigma*dx/c; %taille du pas de temps
lambda=dt/dx;

% initialisation boucle en temps
dimT=.1; %temps de simulation
N=dimT/dt; %nombre d iterations en temps
nTime=0;
T=0;

% initialisation de la solution
u=U0;
u1=U0;
u_prime=zeros(J+1,2);
u_prime1=zeros(J+1,2);
ue=U0;
%uFrot=U0;
%uFrot_prime=zeros(J+1,2);

% boucle en temps
while T < dimT
    minmod=zeros(length(u),2);
    for i = 2:J
        minmod(i,:) = HancockSlopes(u(i-1,:),u(i,:),u(i+1,:),dx);
    end
% boucle en espace
  for j=2:J
    u_prime(j,:) = u(j,:) - lambda*(FluxHLL(u(j+1,:)-minmod(j+1,:)*dx/2-lambda/2*(F(u(j+1,:)+minmod(j+1,:)*dx/2)-F(u(j+1,:)-minmod(j+1,:)*dx/2)),u(j,:)+minmod(j,:)*dx/2-lambda/2*(F(u(j,:)+minmod(j,:)*dx/2)-F(u(j,:)-minmod(j,:)*dx/2)))-FluxHLL(u(j,:)-minmod(j,:)*dx/2-lambda/2*(F(u(j,:)+minmod(j,:)*dx/2)-F(u(j,:)-minmod(j,:)*dx/2)),u(j-1,:)+minmod(j-1,:)*dx/2-lambda/2*(F(u(j-1,:)+minmod(j-1,:)*dx/2)-F(u(j-1,:)-minmod(j-1,:)*dx/2))));
    u_prime1(j,:) = u1(j,:) - lambda*(FluxHLL(u1(j+1,:),u1(j,:))-FluxHLL(u1(j,:),u1(j-1,:)));
    ue(j,:)=Esol(x(j),T+dt,g,ul,ur,hl,hr);
% topographie (explicite)
%    u_prime(j,2)= u_prime(j,2) - dt*g*u(j,1)*(ZB(j+1)-ZB(j-1))/(2*dx);
    Froude(j)=Fr(u_prime(j,:));
    fPrime(j)=Fprime(u_prime(j,:));
    Froude1(j)=Fr(u_prime1(j,:));
    fPrime1(j)=Fprime(u_prime1(j,:));
  end
% fin espace
%conditions artificielles d'extrapolation constante
  u_prime(1,:)=u_prime(2,:);
  u_prime(J+1,:)=u_prime(J,:);
  u_prime1(1,:)=u_prime1(2,:);
  u_prime1(J+1,:)=u_prime1(J,:);
  ue(1,:)=ue(2,:);
  ue(J+1,:)=ue(J,:);

% évolution en temps
  c=max(max(abs(fPrime)),max(abs(fPrime1)));
  dt=sigma*dx/c; %taille du pas de temps
  lambda=dt/dx; 

  u=u_prime;
  u1=u_prime1;

  nTime=nTime+1;
  T=T+dt;

% fin temps
hold off
plot(x,u1(:,1),'red')
hold on
plot(x,u1(:,2),'green')
plot(x,u(:,1),'blue')
plot(x,u(:,2),'magenta')
plot(x,ue(:,1),'black')
plot(x,ue(:,2),'cyan')
title('Shallow water equations with HLL flux & MUSCL - Hancock reconstruction')
legend('height (HLL)','momentum (HLL)','height (MUSCL-Hancock)','momentum (MUSCL-Hancock)','height (exact)','momentum (exact)')
drawnow
end
% afficher la donnée initiale
figure
plot(x,U0(:,1),'r');
hold on
plot(x,U0(:,2))
title('Initial conditions')
legend('height','momentum')
toc