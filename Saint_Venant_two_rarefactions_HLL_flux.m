%*************************************************************
% 9 octobre 2012
% Schema Rusanov pour la resolution de système de Saint-Venant avec ou sans topo
%*************************************************************
clear;
clf;
global g;

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
% mise en place
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
a=0;
b=2; %bornes du domaine
g=9.81;

% maillage en espace
J=100; %nombre de mailles en espace
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
% Two outgoing rarefactions
for i=1:J+1 % une boucle definit directement un vecteur colonne
	if (x(i)<1)
		U0(i,1)=1.;
		U0(i,2)=-1.;
		fPrimeInit(i)=Fprime(U0(i,:));
		Froude(i)=Fr(U0(i,:));
	else
		U0(i,1)=1.0;
		U0(i,2)=1.;
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
u_prime=zeros(J+1,2);
%uFrot=U0;
%uFrot_prime=zeros(J+1,2);

% boucle en temps
while T < dimT
% boucle en espace
  for j=2:J
    u_prime(j,:) = u(j,:) - lambda*(FluxHLL(u(j+1,:),u(j,:))-FluxHLL(u(j,:),u(j-1,:)));
% topographie (explicite)
%    u_prime(j,2)= u_prime(j,2) - dt*g*u(j,1)*(ZB(j+1)-ZB(j-1))/(2*dx);
    Froude(j)=Fr(u_prime(j,:));
    fPrime(j)=Fprime(u_prime(j,:));
  end
% fin espace
%conditions artificielles d'extrapolation constante
  u_prime(1,:)=u_prime(2,:);
  u_prime(J+1,:)=u_prime(J,:);

% évolution en temps
  c=max(abs(fPrime));
  dt=sigma*dx/c; %taille du pas de temps
  lambda=dt/dx; 

  u=u_prime;

  nTime=nTime+1;
  T=T+dt;

% fin temps
hold off
plot(x,u(:,1),'r')
hold on
plot(x,u(:,2))
title('Shallow water equations with HLL flux')
legend('height','momentum')
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