%**************************************************************************
%*              Numeric Solution to Schrodinger Equation                  *
%**************************************************************************     

clear all
close all

%***************** Constants definition  **********************************

Dx=5E-10;           % Spatial step: Width of cells in x axis
Dt=0.05E-16;         % Time step
% Dx=5E-10;           % Spatial step: Width of cells in x axis
% Dt=0.05E-16;         % Time step
L=10e-07;
%T=12000*Dt;         % Duration of experiment
%T=100000*Dt;         % Duration of experiment
T=4e-13;

m=0.041*(9.109e-31);  % Particle mass
h=1.054e-34;        % Planck's constant over 2Pi
qe=1.6e-19;         % Electron charge
m_1=m;

xpoints=floor(L/Dx);% Number of longitude elements in the box
tpoints=floor(T/Dt);% Number of lapses of time in the experiment
tstart1=tpoints/2;
tstart2=tstart1+tpoints/5;
tkill=floor(tpoints/1)+20000;
tpas=floor(tpoints/40);          % Interval de temps per dibuixar
tpasprob=floor(tpoints/200);          % Interval de temps per dibuixar 
ttt=floor(tpoints/20);   
t_transition=tpoints/20;

dispersionx=0.5e-7;  % Spatial dispersion of the wave packet
xcen=-3.5*dispersionx;
Einj=0.12*qe;      % 0.017 first level 0.07 second level
E_final=0.03*qe;%084*qe;   % 0.084 first velevel 0.306 second level

scat_step=1;
 
 k=sqrt(2*m_1*Einj/h^2); 
 k_final=sqrt(2*m_1*E_final/h^2)       %final k of the scattered wavefunction
 k_offset=-(k-k_final);          % + means subtracting

 Energy=(Einj-E_final);      % Minimum energy of wave F2 vanish in the end(g1)
pasene=Energy*4/qe;

%lambda=h*3e8;
%lambda=h*sqrt(2*m*(Energy)/h^2);                  %Auxiliar values
lambda=k_offset;

barlen=2e-9;
barwell=10e-9;
barini=L/2-barwell/2-barlen;
barvalue=0.2*qe;

% % % % C1=1i*h*Dt/(m*Dx^2); % Auxiliar calculations - constants
% % % % C2=-1i*2*Dt/h;
% % % % C4=-2*Dt/(2*m*Dx);
% % % % C5=-qe*1i*2*Dt/h/(qe*2*m);
% [ C1 C2 C4 C5]
% stop
    C1=1i*h*Dt/(m*Dx^2); % Auxiliar calculations - constants
 C2=-qe*1i*2*Dt/h;
    C4=-2*Dt/(2*m*Dx);
    C5=C2/(qe*2*m);
    HB=h
%******************** Potentials *******************%

for j=1:xpoints
    V(j) =0;
end

for j=floor(barini/Dx):floor((barini+barlen)/Dx)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    V(j) = V(j)+barvalue;
end

for j=floor((barini+8*barlen)/Dx):floor((barini+9*barlen)/Dx)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    V(j) = V(j)+barvalue;
end

for j=1:xpoints
    if V(j)>barvalue
    V(j) = barvalue;
    end
end

for j=1:xpoints    
         x(j)=j*Dx-(barini+barlen+barwell/2);
end



t=0;                                      %Initial time = 0    
a=dispersionx;
thetaa = 1/2*atan(2*h*t/(m*a^2));
phi = -thetaa-((h*k^2)*t/(2*m));

for j=1:xpoints                           %Set up gaussian wavefunction at t=0 in F1.
        F1(j,1)=((2*a^2/pi)^(1/4))*(exp(i*phi)/((a^4+(4*h^2*t^2)/m^2)^(1/4)))*exp(i*k*(x(j)-xcen))*exp(-(x(j)-xcen-h*k*t/m)^2/(a^2+2i*h*t/m)); % Gaussian Wave
end

t=Dt;                                      %Initial time = Dt
thetaa = 1/2*atan(2*h*t/(m*a^2));        
phi = -thetaa-((h*k^2)*t/(2*m));
for j=1:xpoints                            %Set up gaussian wavefunction at t=0 in F1.
        F1(j,2)=((2*a^2/pi)^(1/4))*(exp(i*phi)/((a^4+(4*h^2*t^2)/m^2)^(1/4)))*exp(i*k*(x(j)-xcen))*exp(-(x(j)-xcen-h*k*t/m)^2/(a^2+2i*h*t/m)); % Gaussian Wave
end


figure(1)
plot(x,abs(F1(:,1)))
hold on
plot(abs(F1(:,2)))

k_steps=120;% Number of k values
ine=0;
clear S
Emaxim=0.4; %don't keep it too small
kmaxim=sqrt(2*m_1*qe*(Emaxim)/(h^2));  %k vector of the first guess
dk=kmaxim/k_steps;

offset_scat=floor(k_offset/dk);
kom=linspace(-kmaxim,kmaxim,2*k_steps);
for i=1:2*k_steps
Eom(i)=h^2*(kom(i))^2/(2*qe*m(1));
end
n=3;
nt=1;
while (n<tpoints)  %evolution after the scattering 
    
    if (n<=tstart1)
       theta(n)=0; 
    elseif ((n>tstart1)&(n<tstart2))
            theta(n)=(n-tstart1)/(tstart2-tstart1);%*lambda;
    else
       theta(n)=1; 
    end    
    
    F1(1,3)=0;       
    F1(xpoints,3)=0;
    F2(1,3)=0;       
    F2(xpoints,3)=0;

    for j=2:(xpoints-1)   %Schrodinger numeric solution for all points but the ends
         F1(j,3)=F1(j,1)+C1*(F1(j-1,2)-2*F1(j,2)+F1(j+1,2))+C2*V(j)/qe*F1(j,2)+C4*theta(n)*lambda*h*(F1(j+1,2)-F1(j-1,2))+C5*lambda^2*h*h*theta(n)^2*F1(j,2); 
    end                                        %Main F3 shifts to auxiliar F2
        
        F1(:,1)=F1(:,2);                                %Auxiliar F2 shifts to auxiliar F1
        F1(:,2)=F1(:,3); 
        
         if mod(n,tpasprob)==0 
        figure(105)
       hold on
       aux22=abs(F1(:,3)).^2;
       aux44=zeros(xpoints,1);
       aux44(:)=n*Dt;
       plot3(aux44(1:xpoints),x(1:xpoints),aux22(1:xpoints),'b-') 
         end
        
        
     if mod(n,tpas)==0 
         tic
       temps(nt)=n*Dt;                  %Create time vector
       toc
       clear S
       tic
      [S, K, U, K_in, U_in, WF_norma]=WF_mom(F1(:,2), x, V, theta(n)*lambda, Dx, L, k_steps, Emaxim);
       toc
  %S_tot(:,:,nt)=S;
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%post-processing of WF function
c=figure(1111);
tempo=tpas*Dt*nt*1e15;
mesh(x(1:10:end),kom,real(S(1:end-1,1:10:end)))
title(['Time: ' num2str(tempo) ' fs'])
xlabel('Position x (m)')
zlabel('Probability')
ylabel('Wavevector k (m^-1)')
view(0,90);
axis([-5e-7 5e-7 -6e8 6e8])

figure(2222)
mesh(x(1:10:end),kom,imag(S(1:end-1,1:10:end)))

for ii=1:k_steps*2
k_comp(ii)=sum(S(ii,:))*Dx;
% figure(999)
% plot(x,abs(S(ii,:)))
end

for ii=1:xpoints
position(ii)=sum(S(:,ii))*dk;
end

cc=figure(1110);

subplot(3,1,1)
plot(x,abs(position))
xlabel('Position x (m)')
ylabel('Probability')
subplot(3,1,2)
plot(kom,abs(k_comp))
xlabel('Wavevector k (m^-1)')
ylabel('Probability')
len_Eom=length(Eom);
Eom_prim=[-Eom(1:floor(len_Eom/2)) Eom(floor(len_Eom/2)+1:end)];
    
subplot(3,1,3)
plot(Eom_prim,abs(k_comp))
xlabel('Energy (eV)')
ylabel('Probability')


for ii=1:k_steps*2
k_comp(ii)=sum(S(ii,:))*Dx;
% figure(999)
% plot(x,abs(S(ii,:)))
end

for ii=floor(barini/Dx):floor((barini+barwell+2*barlen)/Dx)
position(ii)=sum(S(:,ii))*dk;
end

bb=figure(90);

subplot(3,1,1)
plot(x,abs(position))
xlabel('Position x (m)')
ylabel('Probability')
subplot(3,1,2)
plot(kom,abs(k_comp))
xlabel('Wavevector k (m^-1)')
ylabel('Probability')
len_Eom=length(Eom);
Eom_prim=[-Eom(1:floor(len_Eom/2)) Eom(floor(len_Eom/2)+1:end)];
    
subplot(3,1,3)
plot(Eom_prim,abs(k_comp))
xlabel('Energy (eV)')
ylabel('Probability')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% save figures
       tic
filename=['WF_' num2str(nt) '.png'];
saveas(c,filename)

filename=['WF_integrated_' num2str(nt) '.png'];
saveas(cc,filename)

filename=['WF_integrated_in_' num2str(nt) '.png'];
saveas(bb,filename)
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% just look at the evolution
figure(105)
       hold on
       aux2=abs(F1(:,3)).^2;
       aux4=zeros(xpoints,1);
       aux4(:)=n*Dt;
       plot3(aux4(1:xpoints),x(1:xpoints),aux2(1:xpoints),'b-') 
       aux2=abs(F2(:,3)).^2;
       plot3(aux4(1:xpoints),x(1:xpoints),aux2(1:xpoints),'r-') 

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       nt=nt+1;
       
       
    end   
       n=n+1; 
end
