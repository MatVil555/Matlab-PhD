clear all
close all


%constant
epsilon = 8.85e-12;
mu = 12.56e-7;

%value that modify the analysis
Nz = 4; %number of division of the cable
Nt = 25; %number of time division
deltaZ = 1e-6; %spatial step
V0 = 5; %amplitude of the voltage generator
freq = 1e13; %frequency of the voltage generator
periods = 0.05; %number of periods of the voltage generator to be taken into account

%do not modify
Nv = Nz-2; %number of unknown voltage point
Ni = Nz-1; %number of unknown current point

if (Nz<3) then
	%insufficient number of points
    exit;
end;

lengthC = deltaZ*Nz; %length of the cable (never used)
deltaT = periods/(freq*Nt) %length of a slice of time
T = deltaT * Nt; %total time for the analysis

C = sqrt(epsilon*mu)/50*deltaZ; %capacity
L = 50*sqrt(epsilon*mu)*deltaZ; %inductance
L=1.67E-15;
C=6.67E-15
R = 50; %load resistence


%problem in the form A*X = B

%fill the matrix with zeros
A = zeros(Nz*2-3,Nz*2-3);

%construncting the matrix of the equations
A(1,Nz-1:Nz) = [1, -1];
for i = 1:(Nv-1)
    A(i+1,i) = 1;
    A(i+1,Nv+1+i:Nv+2+i) = [-deltaT/C,  deltaT/C];
end
A(Nv+1,1) = deltaT/L;
A(Nv+1,Nv+2) = 1;
for i = 1:(Ni-2)
    A(Nv+1+i,Nv+2+i) = 1;
    A(Nv+1+i,i:i+1) = [-deltaT/L,  deltaT/L];
end
A(Nv+Ni,Nv) = -1/R;
A(Nv+Ni,Ni+Nv) = 1;

%A = [1,         0,          0,  -deltaT/C,  deltaT/C;
%    0,          0,          1,  -1,         0;
%    deltaT/L,   0,          0,  1,          0;
%    -deltaT/L,  deltaT/L,   0,  0,          1;
%    0,          -1/R,       0,  0,          1];

A
pause

%array of the unknown
X = zeros(Nv+Ni,Nt); %(v2;v3;...vN;i0;i1;i2;...;i(N-1))

V = zeros(Nz,Nt);
I = zeros(Nz,Nt);
%V1 = V0*ones(1,Nt);

for t=1:Nt
V1(t) = V0*sin(2*pi*freq*t*deltaT); %shape of the voltage generator
end

for t = 1:(Nt-1)
    %array of the constant
	B(1,t) = (C/deltaT)*(V1(t+1)-V1(t));
    B(2:Nv,t) = X(1:Nv-1,t);
    B(Nv+1,t) = X(Nv+2,t)+deltaT/L*V1(t+1);
    B(Nv+2:Nv+Ni-1,t) = X(Nv+3:Nv+Ni,t);
    B(Nv+Ni,t) = 0;
    %B2(:,t) = [X(1,t);  (C/deltaT)*(V1(t+1)-V1(t));  X(4,t)+deltaT/L*V1(t+1);  X(5,t);  0];
    X(:,t+1) = A^(-1) * B(:,t);

    solucio=X(:,t+1)
    terme=B(:,t)    
    pause
end

%x & y grid to make plots
[x,y] = meshgrid ([1:1:Nt], [1:1:Nz]);

%transferring the information from the array X to the array I & V
%in the first two points the voltage is the one of the generator
V(1,:) = (V1');
V(2,:) = (V1');
V(3:Nz,:) = X(1:(Nz-2),:);

I(1:(Nz-1),:) = X(Nv+1:(Nv+Nz-1),:);
I(Nz,:) = I(Nz-1,:);

figure
    xlabel('time [s]');
    ylabel('voltage [V]');
    semilogy(x*deltaT,V(Nz,:),'r',x*deltaT,V(1,:),'b');

load PUNTERPSPICE.TXT
t=PUNTERPSPICE(:,1);
vin=PUNTERPSPICE(:,2);
iin=PUNTERPSPICE(:,3);
vout=PUNTERPSPICE(:,4);
iout=PUNTERPSPICE(:,5);

hold on
semilogy(t,vin,'b-o',t,vout,'r-x')
    
    