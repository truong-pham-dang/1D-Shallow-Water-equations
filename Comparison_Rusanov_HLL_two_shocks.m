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
J=40; %nombre de mailles en espace
dx=(b-a)/J; %taille du pas d espace

x=a+[0:J]*dx; %definition du maillage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
%Condition Initiale (doit fournir un vecteur *colonne*)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
% Two outgoing shocks
hl=1;hr=1;
ul=3;ur=-3;
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

% boucle en temps
while T < dimT
% boucle en espace
  for j=2:J
    u_prime(j,:) = u(j,:) - lambda*(FluxRus(u(j+1,:),u(j,:))-FluxRus(u(j,:),u(j-1,:)));
    u_prime1(j,:) = u1(j,:) - lambda*(FluxHLL(u1(j+1,:),u1(j,:))-FluxHLL(u1(j,:),u1(j-1,:)));
    ue(j,:)=Esol(x(j),T+dt,g,ul,ur,hl,hr);
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
plot(x,u(:,1),'red')
hold on
plot(x,u(:,2),'green')
plot(x,u1(:,1),'blue')
plot(x,u1(:,2),'magenta')
plot(x,ue(:,1),'black')
plot(x,ue(:,2),'cyan')
axis([-1 1 -4 4])
title('Shallow water equations with Rusanov & HLL flux')
legend('height (Rusanov)','momentum (Rusanov)','height (HLL)','momentum (HLL)','height (exact)','momentum (exact)')
legend('Location','Southwest')
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