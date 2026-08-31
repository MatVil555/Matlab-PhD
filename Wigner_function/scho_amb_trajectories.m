%**************************************************************************
%*              Numeric Solution to Schrodinger Equation                  *
%**************************************************************************     

clear all
close all
 
%***************** Constants definition  **********************************

pas=30;             % Nombre de particules
numlim=1000;        % Limt pel calcul de les trajectories  

Dx=5E-10;           % Spatial step: Width of cells in x axis
Dt=0.05E-16;         % Time step

L=4E-07;            % Box length

T=10000*Dt;         % Duration of experiment
m=0.15*(9.109e-31);  % Particle mass

h=1.054e-34;        % Planck's constant over 2Pi
qe=1.6e-19;         % Electron charge

xpoints=floor(L/Dx);% Number of longitude elements in the box
tpoints=floor(T/Dt);% Number of lapses of time in the experiment
tstart1=tpoints/4;  % Temps per activar l'scatering
tstart2=tpoints/4+tpoints/20; % temps per desactivar l'scattering
%mac_tstrat1=tstart1; %before this moment for sure no scattering;

tpas=floor(tpoints/10);          % Interval de temps per dibuixar 

xcen=1.5E-07;         % Central position of the wave packet
dispersionx=75*Dx;  % Spatial dispersion of the wave packet
Energy=0.025*qe;      % Minimum energy of wave F2 vanish in the end(g1)
pasene=Energy*4/qe;

%lambda=h*3e8;
lambda=h*sqrt(2*m*(Energy)/h^2);                  %Auxiliar values

k=sqrt(2*m*Energy/h^2);                  %Auxiliar values

Energy/qe
h^2*(k+lambda/h)^2/(2*m*qe)


wosc=2*3.1415*1e13; % Oscilatory frequency  
aosc=-qe*0.016;     % Oscilatory amplitude  V(t)=aosc*cos(wosc*t)

barlen=0.8e-9;
barini=L/2-3*barlen;
barvalue=0.3*qe*0;

C1=i*h*Dt/(m*Dx^2); % Auxiliar calculations - constants
C2=-qe*i*2*Dt/h;
C4=-2*Dt/(2*m*Dx);
C5=C2/(qe*2*m);
[ C1 C2 C4 C5]

%******************** Potentials *******************%

barpoten=aosc*qe*cos(wosc*0);

for j=1:xpoints
    V(j) = -barpoten;
end

for j=floor(barini/Dx):floor((barini+barlen)/Dx)
    V(j) = V(j)+barvalue;
end

for j=floor((barini+5*barlen)/Dx):floor((barini+6*barlen)/Dx)
    V(j) = V(j)+barvalue;
end
 
for j=1:xpoints    
         x(j)=j*Dx;
end


% Set up F1

t=0;                                      %Initial time = 0    
a=dispersionx;
k=sqrt(2*m*Energy/h^2);                  %Auxiliar values
thetaa = 1/2*atan(2*h*t/(m*a^2));
phi = -thetaa-((h*k^2)*t/(2*m));

for j=1:xpoints                           %Set up gaussian wavefunction at t=0 in F1.
        F1(j,1)=((2*a^2/pi)^(1/4))*(exp(i*phi)/((a^4+(4*h^2*t^2)/m^2)^(1/4)))*exp(i*k*(x(j)-xcen))*exp(-(x(j)-xcen-h*k*t/m)^2/(a^2+2i*h*t/m)); % Gaussian Wave
 end

norma1(1)=(F1(:,1)'*F1(:,1))*Dx;           %Set up norm = 1

t=Dt;                                      %Initial time = Dt
thetaa = 1/2*atan(2*h*t/(m*a^2));        
phi = -thetaa-((h*k^2)*t/(2*m));
for j=1:xpoints                            %Set up gaussian wavefunction at t=0 in F1.
        F1(j,2)=((2*a^2/pi)^(1/4))*(exp(i*phi)/((a^4+(4*h^2*t^2)/m^2)^(1/4)))*exp(i*k*(x(j)-xcen))*exp(-(x(j)-xcen-h*k*t/m)^2/(a^2+2i*h*t/m)); % Gaussian Wave
end


norma1(2)=(F1(:,2)'*F1(:,2))*Dx;            %Set up norm = 1


% CALCULATION OF THE SYSTEM EVOLUTION

control=1;
nt=1;
temps(1)=0;
temps(2)=Dt;

norma=(F1(:,2)'*F1(:,2))*Dx;

n=3;

%%%%%%%%%%%%% Calcul de la distribucio de moments  ************
        for indE=1:400
        Ene(indE)=pasene*(indE-200)/200;    
        kvec(indE)=sign(Ene(indE))*sqrt(2*m*qe*abs(Ene(indE))/h^2);     
        ak(indE)=0;    
        for j=2:(xpoints-1)  
           ak(indE)=ak(indE)+F1(j,2)*exp(-1i*kvec(indE)*x(j))*Dx; 
        end
        end
        figure(500)
        plot(Ene,abs(ak).^2,'red')
        hold on


for j=1:pas
position(j,n)=random('normal',xcen,a/2);
inp(j,n)=floor(position(j,n)./Dx)+1;
end

figure(600)
[hy,hx]=hist(position(:,n),150);
plotyy(hx,hy,x,abs(F1(:,2)).^2)

tempstra(1:tpoints)=[1:tpoints]*Dt;
 C1=i*h*Dt/(m*Dx^2); %Calculation of the enewrgy
      %energia
   Cte=-h^2/(2*m*Dx*qe); % Auxiliar calculations - constants
   
   
 for ip=1:pas
     n=3;
while (n<tpoints)                         
    if (n<=tstart1)
       theta(n)=0; 
    elseif ((n>tstart1)&(n<tstart2))
       theta(n)=(n-tstart1)/(tstart2-tstart1)*lambda;
    else
       theta(n)=1*lambda; 
    end    
    %Calculation of new wavefunction F1 and F2   
    F1(1,3)=0;       
    F1(xpoints,3)=0;

    for j=2:(xpoints-1)   %Schrodinger numeric solution for all points but the ends
        F1(j,3)=F1(j,1)+C1*(F1(j-1,2)-2*F1(j,2)+F1(j+1,2))+C2*V(j)/qe*F1(j,2)+C4*theta(n)*(F1(j+1,2)-F1(j-1,2))+C5*theta(n)^2*F1(j,2); 
        velo1(j)=h/m/Dx*imag((F1(j+1,3)-F1(j-1,3))/F1(j,3))+lambda*theta(n)/m;
    end

%     if (n>7000) 
% %     C4
% %     C5
% %     theta(n)
%     
%         
%     end
    %norma=(F1(:,3)'*F1(:,3))*Dx;
    
    
    % Auxiliar calculations - constants

%     for j=2:xpoints-1
%      
%     end
    velo1(1)=velo1(2);
    velo1(xpoints)=velo1(xpoints-1);
 
    for j=1:xpoints
     if (abs(F1(j,3))<0.1) velo1(j)=0; end
    end
 
    %%%%%%  Trajectories quantiques
    
       

        control=0;
        tempo_restante=Dt;

        while control==0

        inp(ip,n)=floor(position(ip,n)/Dx)+1;
        if (inp(ip,n)<3) 
           inp(ip,n)=3; 
        end    
        if (inp(ip,n)>xpoints-3) 
           inp(ip,n)=xpoints-3; 
        end
        vbohm=h/m*imag((F1(inp(ip,n)+1,3)-F1(inp(ip,n)-1,3))/(2*Dx*F1(inp(ip,n),3)))+theta(n)/m;

        if vbohm>0
            dt1=abs(((inp(ip,n)+1)*Dx-position(ip,n))/vbohm);
            if dt1<Dt/numlim 
                dt1=abs(Dx/(vbohm*numlim));
            end
        else
            dt1=abs((position(ip,n)-inp(ip,n)*Dx)/vbohm); 
            if dt1<Dt/numlim 
                dt1=abs(Dx/(vbohm*numlim));
            end     
        end

        tempo=min(Dt,min(dt1,tempo_restante));

        position(ip,n)=position(ip,n)+vbohm*tempo;
        
        

        tempo_restante = tempo_restante -tempo;

        if tempo_restante < Dt/numlim
           control=1;
           position(ip,n+1)=position(ip,n);
           inp(ip,n+1)=inp(ip,n);
        end
        end %while Dt
         %particules
    
    
    %%%%%%  swicht punts i contador
    
    
     for j=1:xpoints
        F1(j,1)=F1(j,2);                                %Auxiliar F2 shifts to auxiliar F1
        F1(j,2)=F1(j,3);                                %Main F3 shifts to auxiliar F2
      end                           
  
    % Result plot 
    if mod(n,tpas)==0                                 
        
        temps(nt)=n*Dt;                  %Create time vector
        %Calculation of the norm -- Step advance of auxiliar functions F1 and F2

       % norma1(nt)=(F1(:,3)'*F1(:,3))*Dx;
 
        %moment
        moment1(nt)=0;
        moment2(nt)=0;
        for j=2:(xpoints-1)  
           moment1(nt)=moment1(nt)-1i*h*conj(F1(j,2))*(F1(j+1,2)-F1(j-1,2))/2; 
        end
        moment2(nt)=lambda*theta(n);

  
        
%         energia0(nt)=0;
%         energia1(nt)=0;
%         energia2(nt)=0;
%         energia3(nt)=0;
%         for j=2:(xpoints-1)  
%            energia0(nt)=energia0(nt)+conj(F1(j,2))*1i*h*(F1(j,3)-F1(j,1))*Dx/(Dt*qe); 
%            energia1(nt)=energia1(nt)+Cte*conj(F1(j,2))*(F1(j-1,2)-2*F1(j,2)+F1(j+1,2)); 
%            energia2(nt)=energia2(nt)+conj(F1(j,2))*F1(j,2)*V(j)*Dx/qe; 
%         end
%         energia3(nt)=lambda^2/(2*m)/qe*theta(n)^2+lambda*moment1(nt)/m/qe*theta(n);
                
%     figure(1)
%     subplot(2,2,1)
%     plot(x,V/qe)
%     xlabel('Position (m)')
%     ylabel('Potential Energy (eV)') 
%     
%     subplot(2,2,3)
%     plotyy(x,velo1,x,abs(F1(:,3)).^2)
%     xlabel('Position (m)')
%     ylabel('velocity (m/s)') 
%    
%     subplot(2,2,2)
%     plot(temps, norma1)
%     xlabel('Time(s)')
%     ylabel('Norm') 
%    
%     subplot(2,2,4)
%     plotyy(temps,(energia1+energia2+energia3),temps,real((moment1+moment2)/m))
%     xlabel('Time (s)')
%     ylabel('Kinetic energy (eV)') 
% 
      aux=abs(F1(:,3)).^2;
%      
      figure(100)
      hold on
      aux4=zeros(xpoints,1);
      aux4(:)=n*Dt;
      plot3(aux4(1:xpoints),x(1:xpoints),aux(1:xpoints),'b-')  
      hold on
%        

       
 
      nt=nt+1;
    end
    
    n=n+1;
   
end
toc
 end


         for indE=1:400
        Ene(indE)=pasene*(indE-200)/200;    
        kvec(indE)=sign(Ene(indE))*sqrt(2*m*qe*abs(Ene(indE))/h^2);     
        ak(indE)=0;    
        for j=2:(xpoints-1)  
           ak(indE)=ak(indE)+F1(j,2)*exp(-1i*kvec(indE)*x(j))*Dx; 
        end
        end
        figure(501)
        plot(Ene,abs(ak).^2,'red')
        hold on
 
  for index=1:pas
 aux2=zeros(tpoints,1);
 figure(100)
 hold on
 plot3(tempstra(3:10:tpoints),position(index,3:10:tpoints),aux2(3:10:tpoints),'r-')
 end
 stop
 
% % figure(300)
% % plotyy(temps,real((moment1+moment2)/m),temps,(energia1+energia2+energia3))
% % title('Moment and energy with phonon')
% % 
% % figure(400)
% % plot(temps,real((moment1+moment2)/m),'r-o',temps,imag((moment1+moment2)/m),'b-x')
% % title('Real(red)  and imaginary(blue)  moment with phonon')
% % 
% % 
% % figure(200)
% % plot(temps,(energia1+energia2+energia3),'b-')
% % hold on
% % plot(temps,(energia0),'g-')
% plot(temps,(energia1),'p-')
% plot(temps,(energia2),'m-')
% plot(temps,(energia3),'k-')




%%%%%%%%%%%%% Calcul de la distribucio de moments  ************
        for indE=1:400
        Ene(indE)=pasene*(indE-200)/200;    
        kvec(indE)=sign(Ene(indE))*sqrt(2*m*qe*abs(Ene(indE))/h^2);     
        ak(indE)=0;    
        for j=2:(xpoints-1)  
           ak(indE)=ak(indE)+F1(j,2)*exp(-1i*kvec(indE)*x(j))*Dx; 
        end
        end
        figure(500)
        plot(Ene,abs(ak).^2,'blue')
        hold on
        
 
 
 stop
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% Set up F1

t=0;                                      %Initial time = 0    
a=dispersionx;
k=sqrt(2*m*Energy/h^2);                  %Auxiliar values
thetaa = 1/2*atan(2*h*t/(m*a^2));
phi = -thetaa-((h*k^2)*t/(2*m));
for j=1:xpoints                           %Set up gaussian wavefunction at t=0 in F1.
        F1(j,1)=((2*a^2/pi)^(1/4))*(exp(i*phi)/((a^4+(4*h^2*t^2)/m^2)^(1/4)))*exp(i*k*(x(j)-xcen))*exp(-(x(j)-xcen-h*k*t/m)^2/(a^2+2i*h*t/m)); % Gaussian Wave
 end

norma1(1)=(F1(:,1)'*F1(:,1))*Dx;           %Set up norm = 1

t=Dt;                                      %Initial time = Dt
thetaa = 1/2*atan(2*h*t/(m*a^2));        
phi = -thetaa-((h*k^2)*t/(2*m));
for j=1:xpoints                            %Set up gaussian wavefunction at t=0 in F1.
        F1(j,2)=((2*a^2/pi)^(1/4))*(exp(i*phi)/((a^4+(4*h^2*t^2)/m^2)^(1/4)))*exp(i*k*(x(j)-xcen))*exp(-(x(j)-xcen-h*k*t/m)^2/(a^2+2i*h*t/m)); % Gaussian Wave
end

norma1(2)=(F1(:,2)'*F1(:,2))*Dx;            %Set up norm = 1


% CALCULATION OF THE SYSTEM EVOLUTION

control=1;
nt=1;
temps(1)=0;
temps(2)=Dt;

norma=(F1(:,2)'*F1(:,2))*Dx;

n=3;


tempstra(1:tpoints)=[1:(tpoints)]*Dt;

while (n<tpoints)                         
    

    %Calculation of new wavefunction F1 and F2   
    F1(1,3)=0;       
    F1(xpoints,3)=0;

    for j=2:(xpoints-1)   %Schrodinger numeric solution for all points but the ends
        F1(j,3)=F1(j,1)+C1*(F1(j-1,2)-2*F1(j,2)+F1(j+1,2))+C2*V(j)/qe*F1(j,2);  
    end
    
    norma=(F1(:,3)'*F1(:,3))*Dx;
    
    %Calculation of the enewrgy
    C1=i*h*Dt/(m*Dx^2); % Auxiliar calculations - constants

    for j=2:xpoints-1
     velo1(j)=h/m/Dx*imag((F1(j+1,3)-F1(j-1,3))/F1(j,3));
    end
    velo1(1)=velo1(2);
    velo1(xpoints)=velo1(xpoints-1);
    
    for j=1:xpoints
     if (abs(F1(j,3))<0.1) velo1(j)=0; end
    end
    
    
    %%%%%%  Trajectories quantiques
    
        for ip=1:pas

        control=0;
        tempo_restante=Dt;

        while control==0

        inp(ip,n)=floor(position(ip,n)/Dx)+1;
        if (inp(ip,n)<3) 
           inp(ip,n)=3; 
        end    
        if (inp(ip,n)>xpoints-3) 
           inp(ip,n)=xpoints-3; 
        end
        vbohm=h/m*imag((F1(inp(ip,n)+1,3)-F1(inp(ip,n)-1,3))/(2*Dx*F1(inp(ip,n),3)));

        if vbohm>0
            dt1=abs(((inp(ip,n)+1)*Dx-position(ip,n))/vbohm);
            if dt1<Dt/numlim 
                dt1=abs(Dx/(vbohm*numlim));
            end
        else
            dt1=abs((position(ip,n)-inp(ip,n)*Dx)/vbohm); 
            if dt1<Dt/numlim 
                dt1=abs(Dx/(vbohm*numlim));
            end     
        end

        tempo=min(Dt,min(dt1,tempo_restante));

        position(ip,n)=position(ip,n)+vbohm*tempo;

        tempo_restante = tempo_restante -tempo;

        if tempo_restante < Dt/numlim
           control=1;
           position(ip,n+1)=position(ip,n);
           inp(ip,n+1)=inp(ip,n);
        end
        end %while Dt
        end %particules
    
    
    %%%%%%  swicht punts i contador
    
    
     for j=1:xpoints
        F1(j,1)=F1(j,2);                                %Auxiliar F2 shifts to auxiliar F1
        F1(j,2)=F1(j,3);                                %Main F3 shifts to auxiliar F2
      end                           
  
    % Result plot 
    if mod(n,tpas)==0                                 
        
        temps(nt)=n*Dt;                  %Create time vector
        %Calculation of the norm -- Step advance of auxiliar functions F1 and F2

        norma1(nt)=(F1(:,3)'*F1(:,3))*Dx;
 
        %moment
        moment1(nt)=0;
        moment2(nt)=0;
        for j=2:(xpoints-1)  
           moment1(nt)=moment1(nt)-1i*h*conj(F1(j,2))*(F1(j,2)-F1(j-1,2)); 
        end
        moment2(nt)=0;
        
        %energia
        Cte=-h^2/(2*m*Dx*qe); % Auxiliar calculations - constants
        
        energia1(nt)=0;
        energia2(nt)=0;
        energia3(nt)=0;
        for j=2:(xpoints-1)  
           energia1(nt)=energia1(nt)+Cte*conj(F1(j,2))*(F1(j-1,2)-2*F1(j,2)+F1(j+1,2)); 
           energia2(nt)=energia2(nt)+conj(F1(j,2))*F1(j,2)*V(j)*Dx/qe; 
           energia3(nt)=0.0; 
        end
        
    figure(1)
    subplot(2,2,1)
    plot(x,V/qe)
    xlabel('Position (m)')
    ylabel('Potential Energy (eV)') 
    
    subplot(2,2,3)
    plotyy(x,velo1,x,abs(F1(:,3)).^2)
    xlabel('Position (m)')
    ylabel('velocity (m/s)') 
   
    subplot(2,2,2)
    plot(temps, norma1)
    xlabel('Time(s)')
    ylabel('Norm') 
   
    subplot(2,2,4)
    plotyy(temps,(energia1+energia2+energia3),temps,real((moment1+moment2)/m))
    xlabel('Time (s)')
    ylabel('Kinetic energy (eV)') 
   
    
     aux=abs(F1(:,3)).^2;
    
     figure(100)
     hold on
     aux4=zeros(xpoints,1);
     aux4(:)=n*Dt;
     plot3(aux4(1:10:xpoints),x(1:10:xpoints),aux(1:10:xpoints),'r-')
         
        nt=nt+1;
    end
    
    n=n+1;
end

figure(200)
plot(temps,(energia1+energia2+energia3),'r-')
hold on
xlabel('Time(s)')
ylabel('Energy(eV)')


 for index=1:pas
 figure(100)
 plot3(tempstra(3:10:tpoints),position(index,3:10:tpoints),aux2(3:10:tpoints),'r-')
 hold on
 end
 
 aux4(:)=0;
 plot3(aux4(1:xpoints),x(1:xpoints),5*V(1:xpoints)*max(aux)/qe,'g-')

 grid on
 ylabel('Position (m)')
 xlabel('Time(s)')
 zlabel('Conditional wave-fucntion modulus')
 


 %%%%%%%%%%%%% Calcul de la distribucio de moments  ************
        for indE=1:400
        Ene(indE)=pasene*(indE-200)/200;    
        kvec(indE)=sign(Ene(indE))*sqrt(2*m*qe*abs(Ene(indE))/h^2);     
        ak(indE)=0;    
        for j=2:(xpoints-1)  
           ak(indE)=ak(indE)+F1(j,2)*exp(-1i*kvec(indE)*x(j))*Dx; 
        end
        end
        figure(500)
        plot(Ene,abs(ak).^2,'black')
        hold on
