 %**************************************************************************
%*              Numeric Solution to Schrodinger Equation                  *
%**************************************************************************     
% stop
 clear all
 close all

 
%***************** Constants definition  **********************************
pas=0;
Dx=5E-10;           % Spatial step: Width of cells in x axis
Dt=0.05E-16;         % Time step
L_active=40e-07;
L_abs=0.0E-07;
L=L_active+2*L_abs;            % Box length
n_s=1;
%T=12000*Dt;         % Duration of experiment
T=50000*Dt;         % Duration of experiment
T=5e-12;
tot_time=T;
m=0.041*(9.109e-31);  % Particle mass

h=1.054e-34;        % Planck's constant over 2Pi
qe=1.6e-19;         % Electron charge
%alpha
m_1=m;

xpoints=floor(L/Dx);% Number of longitude elements in the box
tpoints=floor(T/Dt);% Number of lapses of time in the experiment
tpoints1=tpoints/1;
tkill=floor(tpoints/1)+20000;
tpas=floor(tpoints/2000);          % Interval de temps per dibuixar
ddtt=tpas*Dt;
tpasprob=floor(tpoints/200);          % Interval de temps per dibuixar 
ttt=floor(tpoints/20);   
om=(0.128-0.07)*qe/h;
om=(0.25-0.064)*qe/h;
Dq=0.5*sqrt(m_1)*Dx;
maxq=1000;
L_q=maxq*Dq;

num_sampling=1;

ampl=0.02;
om_AC=3.33e12
om_AC=logspace(11,13,1);
om_AC=1e12*pi*2; %om=2*pi*f=2*pi/T
ampl=[0.001 0.01 0.02 0.05 0.2]*0;
for ii=1:1
for jj=3:3

    clear D_tot_real
    clear D_real_large
    clear D_tot
    clear D
    close all
% [ii jj]
ampl_n=ampl(jj);
om_AC_n=om_AC(ii);

for p=0:num_sampling-1
for n=1:tpoints
V_AC(n,p+1)=ampl_n*sin(om_AC_n*n*Dt+p*2*pi/num_sampling);
end
end
figure(22)
plot(V_AC)

numlim=1000;
tempstra(1:tpoints)=[1:tpoints]*Dt;
%D(:,:)=zeros(200,1);

dispersionx=250*Dx;  % Spatial dispersion of the wave packet
% dispersionx=35*Dx;  % Spatial dispersion of the wave packet
% dispersionx=45*Dx;
xcen=L/2-3.5*dispersionx

Einj=0.058*qe;      % 0.017 first level 0.07 second level
E_final=0.084*qe;%084*qe;   % 0.084 first velevel 0.306 second level


 %alpha(j)=5e-1*qe/L;

k=sqrt(2*m_1*Einj/h^2); 
k_final=sqrt(2*m_1*E_final/h^2)       %final k of the scattered wavefunction
k_offset=+(k-k_final);          % + means subtracting

%k=sqrt(2*m*Energy/h^2);                  %Auxiliar values

Einj/qe

barlen=2e-9;
barwell=10e-9;
barini=L/2-barwell/2-barlen;
barvalue=0.5*qe;
barvalue2=0.5*qe;
C1=1i*h*Dt/(m*Dx^2); % Auxiliar calculations - constants
C2=-1i*2*Dt/h;

%******************** Potentials *******************%

barpoten=0;

for j=1:xpoints
    V_or(j) = -barpoten;
end

% for j=floor((barini-barlen)/Dx):floor((barini)/Dx)
%     V(j) = V(j)+barvalue*(floor((barini-barlen)/Dx)-j)/(floor((barini-barlen)/Dx)-floor((barini)/Dx));%%%
% end

for j=floor(barini/Dx):floor((barini+1*barlen)/Dx)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    V_or(j) = V_or(j)+barvalue;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%engineering of the well
%%%%%%%%%%square
% for j=floor((barini+5*barlen)/Dx):floor((barini+9*barlen)/Dx)%%%3.1THz% con 24nm%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%     V_or(j) = V_or(j)+barvalue/6;
% end

% for j=floor((barini+4*barlen)/Dx):floor((barini+10*barlen)/Dx)%% 2.5THz con 24nm 
%     V_or(j) = V_or(j)+barvalue/5;
% end

for j=floor((barini+5*barlen)/Dx):floor((barini+11*barlen)/Dx)%% 1.7THz con 28nm 
    V_or(j) = V_or(j)+barvalue/5;
end

% for j=floor((barini+5.5*barlen)/Dx):floor((barini+10.5*barlen)/Dx)%% 1.7THz con 28nm 
%     V_or(j) = V_or(j)+barvalue/6;
% end

%%%%%%%%%%triangle
% for j=floor((barini+4*barlen)/Dx):floor((barini+7*barlen)/Dx)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%     V_or(j) = V_or(j)+barvalue/10*(j*Dx-(barini+4*barlen))/((barini+7*barlen)-(barini+4*barlen));
% end
% 
% for
% j=floor((barini+7*barlen)/Dx)+1:floor((barini+10*barlen)/Dx)%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%3.6THz con con 24nm
%     V_or(j) = V_or(j)+barvalue/10*((barini+10*barlen)-j*Dx)/(barini+10*barlen-(barini+7*barlen));
% end


%%%%%%%%%%%%%%%parabola
% for j=floor((barini+4*barlen)/Dx):floor((barini+10*barlen)/Dx)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%     V_or(j) = V_or(j)+barvalue/5*(j*Dx-(barini+4*barlen))/((barini+10*barlen)-(barini+4*barlen));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=floor((barini+15*barlen)/Dx):floor((barini+16*barlen)/Dx)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    V_or(j) = V_or(j)+barvalue;
end

% for j=floor((barini+9.5*barlen)/Dx):floor((barini+10.5*barlen)/Dx)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%     V(j) = V(j)+barvalue3*0;
% end

for j=1:xpoints    
         x(j)=j*Dx;%-(barini+barlen+barwell/2);
end

% Potential for the quantum-well alone
ccc=0;  %ccc=30 or 50, gives me transition from 1st-2nd (combination) to 2nd-3rd!!!!
% 60, 45they become the same, the red is higher than the blue
% 30  they oscillate with frequency 1.5 THz important: the frequency
% depends on the injectione enrgy. If E is lower, f is higher.



%stop

% decrease deltaE
for j=1:xpoints    
         if x(j)>L/2+barlen+barwell/2
            xp(j)=x(floor((barini+barwell+2*barlen)/Dx)); 
            xp(j)=0;
%             alpha(j)=0;
%             alpha(j)=1e2*qe/L;
         elseif x(j)<L/2-(barlen+barwell/2)
            xp(j)=x(floor((barini)/Dx));
            xp(j)=0;
%             alpha(j)=0;
%             alpha(j)=1e2*qe/L;
         else if x(j)>=L/2-(barlen+barwell/2) && x(j)<L/2-barwell/2
            xp(j)=x(j)-L/2;
%             alpha(j)=ccc*qe/L*((x(j)+barlen+barwell/2)/(-barwell/2+barlen+barwell/2))
         else if x(j)<=L/2+(barlen+barwell/2) && x(j)>L/2+barwell/2
            xp(j)=x(j)-L/2;
%             alpha(j)=ccc*qe/L*((x(j)-(barlen+barwell/2))/(barwell/2-(barlen+barwell/2)))
         else
            xp(j)=x(j)-L/2; 
         end
        % xp(j)=x(j);
%          xp(j)=x(j);
end
end
alpha(j)=ccc*qe/L;
end

figure(1)
plot(x,V_or/qe,'r-o')
% figure(2)
% plot(xp)
% figure(3)
% plot(alpha)
% figure(4)
% plot(alpha.*xp,'-o')


% Set up F1

t=0;                                      %Initial time = 0    
a=dispersionx;
%k=sqrt(2*m*Energy/h^2);                  %Auxiliar values
thetaa = 1/2*atan(2*h*t/(m*a^2));
phi = -thetaa-((h*k^2)*t/(2*m));
p=1;
for j=1:xpoints                           %Set up gaussian wavefunction at t=0 in F1.
        F1(j,1)=((2*a^2/pi)^(1/4))*(exp(i*phi)/((a^4+(4*h^2*t^2)/m^2)^(1/4)))*exp(i*k*(x(j)-xcen))*exp(-(x(j)-xcen-h*k*t/m)^2/(a^2+2i*h*t/m)); % Gaussian Wave
        %F2(j,1)=((2*a^2/pi)^(1/4))*(exp(i*phi)/((a^4+(4*h^2*t^2)/m^2)^(1/4)))*exp(i*k*(x(j)-xcen))*exp(-(x(j)-xcen-h*k*t/m)^2/(a^2+2i*h*t/m)); % Gaussian Wave
        %F1(j,1)=0;
        F2(j,1)=0;
end

norma1(1)=(F1(:,1)'*F1(:,1))*Dx;           %Set up norm = 1
norma2(1)=(F2(:,1)'*F2(:,1))*Dx;           %Set up norm = 1



t=Dt;                                      %Initial time = Dt
thetaa = 1/2*atan(2*h*t/(m*a^2));        
phi = -thetaa-((h*k^2)*t/(2*m));
for j=1:xpoints                            %Set up gaussian wavefunction at t=0 in F1.
        F1(j,2)=((2*a^2/pi)^(1/4))*(exp(i*phi)/((a^4+(4*h^2*t^2)/m^2)^(1/4)))*exp(i*k*(x(j)-xcen))*exp(-(x(j)-xcen-h*k*t/m)^2/(a^2+2i*h*t/m)); % Gaussian Wave
        %F2(j,2)=((2*a^2/pi)^(1/4))*(exp(i*phi)/((a^4+(4*h^2*t^2)/m^2)^(1/4)))*exp(i*k*(x(j)-xcen))*exp(-(x(j)-xcen-h*k*t/m)^2/(a^2+2i*h*t/m)); % Gaussian Wave
        F2(j,2)=0;
        %F1(j,2)=0;
end

for p=1:num_sampling
    F1(:,1,p)=F1(:,1);
    F1(:,2,p)=F1(:,2);
    F2(:,1,p)=F2(:,1);
    F2(:,2,p)=F2(:,2);
end

% F_ini(:,1)=F1(:,1);
% F_ini(:,2)=F1(:,2);

xpoints_abs=floor(L_abs/Dx);
norma1(2)=(F1(:,2)'*F1(:,2))*Dx;            %Set up norm = 1
norma2(2)=(F2(:,2)'*F2(:,2))*Dx;            %Set up norm = 1
par_abs=50
par_abs2=200


% CALCULATION OF THE SYSTEM EVOLUTION

nt=1;
n=3;
for yy=1:1%calculation of the T and scattering 
k_steps=2000;
Emaxim=0.15;
kmaxim=sqrt(2*m_1*qe*(Emaxim)/(h^2));  %k vector of the first guess
dk=kmaxim/k_steps;
% dE=(dk^2)*(h^2)/2/m_1/qe
% dE=Emaxim/(k_steps+2);
% ddkk=sqrt(2*m_1*qe*(dE)/(h^2));
offset_scat=floor(k_offset/dk)*0;

A=(10*Dx^2)/12;
AA=Dx^2/12;

ine=0;

%%%%%%%%%%%%%%%%%%%%%
parpar=1;
dx=Dx/parpar;
par=floor(Dx/dx);
xpoints1=floor(L/dx);
pos_l=3;
pos_r=xpoints1-3;
% if par==1
 Va(:)=V_or(:)
% xx(:)=x(:)
% else
% % for iii=1:xpoints1    
% %     xx(iii)=iii*dx;
% % end
for iii=1:xpoints1    
         xx(iii)=iii*dx-(barini+barlen+barwell/2);
end


% % % figure(11)
% % % plot(xx,Va)



tic
 for kright=-kmaxim:dk:0
  %for kright=0:0   
   ine=ine+1;
   E=h^2/(2*m_1)*kright^2+abs(Va(xpoints1));  % The negative sign of Eind is an artifical convention to refer to right injection. E has the correct physical sign, but it has E=0 at the bottom of the condution band at the source-left. 
   Eaxis=-E/qe;
   

   %CALCULO TEORICO DE LA FUxpointsCIOxpoints DE OxpointsDA 
   if E>Va(1) 
     kleft=-sqrt((2*m_1*(E-Va(1)))/(h^2));% 1/m %this is the k vector in the drain, so m=m_1 

     for j=1:xpoints1
           FO(j)=cos(kleft*(xx(j)-xx(pos_l)))+1i*sin(kleft*(xx(j)-xx(pos_l))); %first guess sin or con xpointsumerov 
           
     end

    FOxpoints(1,ine)=FO(1);
    FOxpoints(2,ine)=FO(2);
       for j=3:xpoints1
        %FOxpoints(j,ine)=FOxpoints(j-1,ine)*((2*m(j-1)*Dx^2*(V(j-1)-E)/(h^2))+1+(m(j-1)/m(j-2)))-FOxpoints(j-2,ine)*(m(j-1)/m(j-2)); 
        FOxpoints(j,ine)=FOxpoints(j-1,ine)*((2*m*dx^2*(Va(j-1)-E)/(h^2))+2)-FOxpoints(j-2,ine);
       end
      
      DER=(-FOxpoints(pos_r+2,ine)+8*FOxpoints(pos_r+1,ine)-8*FOxpoints(pos_r-1,ine)+FOxpoints(pos_r-2,ine))/(12*dx); 
      FOxpoints1RE=real(FOxpoints(pos_r,ine));
      DERRE=real(DER);
      FOxpoints1IMAG=imag(FOxpoints(pos_r,ine));
      DERIMAG=imag(DER);

      T1=(FOxpoints1RE+DERIMAG/kright)^2+(FOxpoints1IMAG-DERRE/kright)^2;
      T=4*kleft/(T1*kright); %in function of space
      if T<1.1 || barvalue==0
     Tran(ine)=T;
      
      norma=sqrt(T/(2*pi)*kright/kleft);  %%%%% (ORIOLS) RIGHT IxpointsCIDExpointsT
      for j=1:xpoints1                           %%%%% (ORIOLS) 
            FOxpoints(j,ine)=FOxpoints(j,ine)*norma;                %%%%% (ORIOLS)  
      end                                 %%%%% (ORIOLS) 
       
      else if T>1.1 && barvalue>0 
      Tran(ine)=0;

      norma=0;                            %%%%% (ORIOLS) RIGHT IxpointsCIDExpointsT
      for j=1:xpoints1                           %%%%% (ORIOLS) 
            FOxpoints(j,ine)=FOxpoints(j,ine)*norma;                %%%%% (ORIOLS)  
      end                                 %%%%% (ORIOLS)    
      
         end
     end
      Ene(ine)=Eaxis;
   
   end
 end
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
 for kleft=0:dk:kmaxim    %%%%% (ORIOLS) LEFT IxpointsCIDExpointsT 
     
   ine=ine+1;
   E=h^2/(2*m_1)*kleft^2+Va(1); %in J % The negative sign of Eind is an artifical convention to refer to right injection. E has the correct physical sign, but it has E=0 at the bottom of the condution band at the source-left.                
   E/qe;
   Eaxis=E/qe;
  
   if E>Va(xpoints1)
       
   %CALCULO TEORICO DE LA FUxpointsCIOxpoints DE OxpointsDA
   kright=sqrt(2*m_1*(E-Va(xpoints1))/(h^2));   
 for j=1:xpoints1
       FO(j)=cos(kright*(xx(j)-xx(pos_r)))+1i*sin(kright*(xx(j)-xx(pos_r))); %first guess sin or con xpointsumerov
 end
 FOxpoints(xpoints1,ine)=FO(xpoints1);
 FOxpoints(xpoints1-1,ine)=FO(xpoints1-1);
 for j=xpoints1-2:-1:1
    FOxpoints(j,ine)=FOxpoints(j+1,ine)*((2*m*dx^2*(Va(j+1)-E)/(h^2))+2)-FOxpoints(j+2,ine);
 end
 
 DER=(-FOxpoints(pos_l+2,ine)+8*FOxpoints(pos_l+1,ine)-8*FOxpoints(pos_l-1,ine)+FOxpoints(pos_l-2,ine))/(12*dx); 

 FOxpoints1RE=real(FOxpoints(pos_l,ine));
 DERRE=real(DER);
 FOxpoints1IMAG=imag(FOxpoints(pos_l,ine));
 DERIMAG=imag(DER);
 
 T1=(FOxpoints1RE+DERIMAG/kleft)^2+(FOxpoints1IMAG-DERRE/kleft)^2;  
  ine;
 T=4*kright/(T1*kleft); %in function of space
 if T<1.1 || barvalue==0
 Tran(ine)=T;
  
 norma=sqrt(T/(2*pi)*kleft/kright);  
 for j=1:xpoints1                           
      FOxpoints(j,ine)=FOxpoints(j,ine)*norma;                  
 end                                 
  
 else if T>1.1 && barvalue>0 
 Tran(ine)=0;
  
 norma=0;                           
 for j=1:xpoints1                            
      FOxpoints(j,ine)=FOxpoints(j,ine)*norma;                  
 end                                
  
     end
  end
 Ene(ine)=Eaxis;
  
end
  
 end
 toc

 if offset_scat>0
        lim1=1;
        lim2=ine-offset_scat;
    else
        lim1=1-offset_scat;
        lim2=ine;       
    end
 
for ine2=1:ine %decomposition of the components of the wavefunction
ak(ine2,1)=0;
for j=1:xpoints                      
      ak(ine2,1)=ak(ine2,1)+conj(FOxpoints(j,ine2))*F1(j,1)*Dx;
      ak2(ine2,1)=ak(ine2,1)+conj(FOxpoints(j,ine2))*F2(j,1)*Dx;
end      
 ine2;
end

figure(66)
plot(Ene,(Tran))
%plot(Ene,log(1))
ylabel('T')
xlabel('Energias (eV)')

% figure(60)
% plot((Tran))
% %plot(Ene,log(1))
% ylabel('T')
% xlabel('Energias (eV)')

figure(6)
subplot(3,1,1)
plot(Ene,log(Tran))
%plot(Ene,log(1))
ylabel('T')
xlabel('Energias (eV)')



old=abs(ak2(1:ine,1)).^2/max(abs(ak2(1:ine,1)).^2);
  %new=abs(ak(lim1+offset_scat:lim2+offset_scat,1)).^2/max(abs(ak(1:ine,1)).^2)
  plo2=abs(ak2(lim1:lim2,1)).^2/max(abs(ak2(1:ine,1)).^2)
  
subplot(3,1,2)
%plot(Ene,abs(ak2(1:ine,1)).^2/max(abs(ak2(1:ine,1)).^2),'-b')
%hold on
plot(Ene(lim1:lim2),log(plo2),'-r')
ylabel('ak')
xlabel('Energias (eV)')
for kk=1:length(plo2)
plo3(kk)=log(Tran(kk)*plo2(kk));
end

subplot(3,1,3)
plot(Ene(lim1:lim2),plo3,'-r')
ylabel('ak*T')
xlabel('Energias (eV)')

end
% % % %  par_gauss_q=2*Dq*0.5*1e-24;
% % % % q=linspace(-L_q/2,L_q/2,maxq);
q=linspace(0,L_q,maxq);
par_gauss_q=(L_q/10)*0.2*1e-24;


a_q=sqrt(par_gauss_q);
%a_q=2*sqrt(par_gauss_q);

for j=1:pas
%     tic
position_x(j,n)=random('normal',xcen,a/2);
inp(j)=floor(position_x(j,n)./Dx)+1;
position_q(j,n)=random('normal',L_q/2,a_q);
inp_q(j)=floor(position_q(j,n)./Dq)+1;
% toc
end

%%%%%%%%%%%%%%%%%%%%%% photon states harmonic oscillator


for j=1:maxq  
    harm_1(j)=(pi)^(-1/4)*exp(0.5*(-(q(j)-(L_q/2)).^2)/par_gauss_q);
    harm_2(j)=(1/sqrt(2))*(pi)^(-1/4)*(2*(q(j)-(L_q/2)))*exp((-(q(j)-(L_q/2)).^2)/(par_gauss_q));
end
integr_1=sum(harm_1.^2)*Dq
integr_2=sum(harm_2.^2)*Dq

harm_1=harm_1.*sqrt((1/integr_1));
harm_2=harm_2.*sqrt((1/integr_2));

count_ini=zeros(1,maxq)
for ip=1:pas
    for i=1:maxq
    if position_q(ip,3)>=i*Dq && position_q(ip,3)<(i+1)*Dq
        count_ini(i)=count_ini(i)+1;
    end
    end
end

% figure(7)
% plot(q,harm_1)
% 
% figure(71)
% plot(harm_1)
% hold on
% bar(count_ini*(max(harm_1)/max(count_ini)))
% 
% figure(77)
% plot(q,harm_2)

integr_norm_1=sum(harm_1.^2)*Dq
integr_norm_2=sum(harm_2.^2)*Dq

%stop

count_ini_x=zeros(1,xpoints)
for ip=1:pas
    for i=1:xpoints
    if position_x(ip,3)>=i*Dx && position_x(ip,3)<(i+1)*Dx
        count_ini_x(i)=count_ini_x(i)+1;
    end
    end
end

% figure(81)
% 
% plot(abs(F1(:,1)).^2)
% hold on
% bar(count_ini_x*(max(abs(F1(:,1)).^2)/max(count_ini_x)))
% 
% figure(77)
% plot(q,harm_2)
% % % V_or(:)=V(:);
flag_ini=0;
for j=1:xpoints
    if V_or(j)~=0 && flag_ini==0;
      barini_n=j+2; 
      barini_n_real=j;
      flag_ini=1;
    end
end

flag_end=0;
for j=xpoints:-1:1
if V_or(j)~=0 && flag_end==0
      barend_n=j-2;  
      barend_n_real=j;
      flag_end=1;
end
end
%stop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (n<tpoints1)  %evolution after the scattering 
    flag_nt=0;
    for p=1:num_sampling %boocle on the phases
    %update the potential profile
    for j=1:xpoints
        if j<barini_n_real-15
        V(j,p)=V_or(j);
        else if j<=barend_n_real+15
        V(j,p)=V_or(j)+V_AC(n,p)*qe*(j-(barini_n_real-15))/(barend_n_real-barini_n_real+30);        
        else
        V(j,p)=V_or(j)+V_AC(n,p)*qe;  
            end
        end
    end
    
%     figure(11)
%     plot(V(:,:))

    
    %Calculation of new wavefunction F1 and F2   
    F1(1,3,p)=0;       
    F1(xpoints,3,p)=0;
% 
%     F2(1,3,p)=0;       
%     F2(xpoints,3,p)=0;
%  tic
%     for j=2:(xpoints-1)   %Schrodinger numeric solution for all points but the ends
%      if V(j)==0
%          F1(j,3,p)=F1(j,1,p)+C1*(F1(j-1,2,p)-2*F1(j,2,p)+F1(j+1,2,p));
%      else
%         F1(j,3,p)=F1(j,1,p)+C1*(F1(j-1,2,p)-2*F1(j,2,p)+F1(j+1,2,p))+C2*V(j,p)*F1(j,2,p);%-1i*2*Dt*0.5*om*F1(j,2,p)+C2*alpha(j)*xp(j)*F2(j,2,p);
%        % F2(j,3,p)=F2(j,1,p);%+C1*(F2(j-1,2,p)-2*F2(j,2,p)+F2(j+1,2,p))+C2*V(j,p)*F2(j,2,p)-1i*2*Dt*1.5*om*F2(j,2,p)+C2*alpha(j)*xp(j)*F1(j,2,p); %%%+C1*(F2(j-1,2)-2*F2(j,2)+F2(j+1,2))
%      end
%     end
%  toc
 %tic

        F1(2:xpoints-1,3,p)=F1(2:xpoints-1,1,p)+C1*(F1(1:xpoints-2,2,p)-2*F1(2:xpoints-1,2,p)+F1(3:xpoints,2,p))+C2*V(2:xpoints-1,p).*F1(2:xpoints-1,2,p);%-1i*2*Dt*0.5*om*F1(j,2,p)+C2*alpha(j)*xp(j)*F2(j,2,p);
        F2(2:xpoints-1,3,p)=F2(2:xpoints-1,1,p)+C1*(F2(1:xpoints-2,2,p)-2*F2(2:xpoints-1,2,p)+F2(3:xpoints,2,p))+C2*V(2:xpoints-1,p).*F2(2:xpoints-1,2,p);%-1i*2*Dt*0.5*om*F1(j,2,p)+C2*alpha(j)*xp(j)*F2(j,2,p);
        

 %toc


          F2(:,1,p)=F2(:,2,p); 
          F1(:,1,p)=F1(:,2,p); 
                                           %Auxiliar F2 shifts to auxiliar F1
        F2(:,2,p)=F2(:,3,p);                                %Main F3 shifts to auxiliar F2                            %Auxiliar F2 shifts to auxiliar F1
        F1(:,2,p)=F1(:,3,p); 

        %F(:,:)=F1(:,3)*harm_1(:)'+F2(:,3)*harm_2(:)';

    if pas>0
    for ip=1:pas

        control=0;
        tempo_restante=Dt;

        while control==0
        % index in space
        inp(ip)=floor(position_x(ip,n)/Dx)+1;
        if (inp(ip)<3) 
           inp(ip)=3; 
        end    
        if (inp(ip)>xpoints-3) 
           inp(ip)=xpoints-3; 
        end
        % index in q
        inp_q(ip)=floor(position_q(ip,n)/Dq)+1;
        if (inp_q(ip)<3) 
           inp_q(ip)=3; 
        end    
        if (inp_q(ip)>maxq-3) 
           inp_q(ip)=maxq-3; 
        end
        
%%%%%%%%%%%%%%5%speed in space
        vbohm=imag((F(inp(ip)+1,inp_q(ip))-F(inp(ip)-1,inp_q(ip)))/(2*Dx*F(inp(ip),inp_q(ip))));%+lambda*theta(n);
        vbohm=vbohm*h/m_1;

        if vbohm>0
            dt1=abs(((inp(ip)+1)*Dx-position_x(ip,n))/vbohm);
            if dt1<Dt/numlim 
                dt1=abs(Dx/(vbohm*numlim));
            end
        else
            dt1=abs((position_x(ip,n)-(inp(ip)-1)*Dx)/vbohm); 
            if dt1<Dt/numlim 
                dt1=abs(Dx/(vbohm*numlim));
            end     
        end

        %%%%%%%%%%speed in q
         vbohm_q=imag((F(inp(ip),inp_q(ip)+1)-F(inp(ip),inp_q(ip)-1))/(2*Dq*F(inp(ip),inp_q(ip))));%+lambda*theta(n);
        vbohm_q=vbohm_q*h;
       % V_bohm_q_mem(ip,n)=vbohm_q;
        if vbohm_q>0
            dt2=abs(((inp_q(ip)+1)*Dq-position_q(ip,n))/vbohm_q);
            if dt2<Dt/numlim 
                dt2=abs(Dq/(vbohm_q*numlim));
            end
        else
            dt2=abs((position_q(ip,n)-(inp_q(ip)-1)*Dq)/vbohm_q); 
            if dt2<Dt/numlim 
                dt2=abs(Dq/(vbohm_q*numlim));
            end     
        end
 
 %%%%%%%%%%%%%%%%%calculate and update new positions       
        tempo=min(Dt,min(dt1,min(dt2,tempo_restante)));

        position_x(ip,n)=position_x(ip,n)+vbohm*tempo;
        position_q(ip,n)=position_q(ip,n)+vbohm_q*tempo;
        
        tempo_restante = tempo_restante -tempo;

        if tempo_restante < Dt/numlim
           control=1;
           position_x(ip,n+1)=position_x(ip,n);
           position_q(ip,n+1)=position_q(ip,n);
        end
        end %while Dt
        end %particules
    end
        
    
    
    % Result plot 
    if mod(n,tpas)==0    
        tic
    F_t_1(:,p)= F1(barini_n:barend_n,3,p);
    %F_t_2(:,p)= F2(barini_n:barend_n,3,p);
          lenlen1=length(F_t_1(:,p));
%          lenlen2=length(F_t_2(:,p));
%          if lenlen1~=lenlen2
%              stop
%          else
             lenlen=lenlen1;
%          end
        for j=2:lenlen-1
        DER_F_t_1(j,p)=(F_t_1(j+1,p)-F_t_1(j-1,p));
%        DER_F_t_2(j,p)=(F_t_2(j+1,p)-F_t_2(j-1,p));
        %DER_F_t(j)=(-F_t(j+2)+8*F_t(j+1)-8*F_t(j-1)+F_t(j-2))/(12*Dx); 
        end
        DER_F_t_1(:,p)=DER_F_t_1(:,p)/(2*Dx);
 %       DER_F_t_2(:,p)=DER_F_t_2(:,p)/(2*Dx);        
        for j=2:lenlen-1
        J_1(j,p)=real(conj(F_t_1(j,p)).*(1/1i).*DER_F_t_1(j,p));%/abs(F_t(j))^2
%        J_2(j,p)=real(conj(F_t_2(j,p)).*(1/1i).*DER_F_t_2(j,p));%/abs(F_t(j))^2
        end
        J_1(:,p)=(h/m).*J_1(:,p);
%        J_2(:,p)=(h/m).*J_2(:,p);
        J(:,p)=J_1(:,p);%+J_2(:,p)
%         show_1=max(abs(F_t).^2)/max(V(barini_n:barend_n));
%         show_2=max(J)/max(V(barini_n:barend_n));

        
%         if nt>=2
        D_prim(nt,p)=sum(J(:,p));%*Dx
        D(nt,p)=D_prim(nt,p)/(lenlen);%/Dx
%         if mod(n,tpas*10)==0
        figure(30)
        if p==1
        hold off
        end
        plot(tempstra(1:tpas:n),D(1:nt,p))
        hold on
        toc
%         end
%         end

        
        
       temps(nt)=n*Dt;                  %Create time vector
       %Calculation of the norm -- Step advance of auxiliar functions F1 and F2

      
      
if p==floor(num_sampling/2) || num_sampling==1

for ine2=1:ine %decomposition of the components of the wavefunction
ak(ine2,1)=0;
ak2(ine2,1)=0;
for j=barini_n:barend_n                   
      ak(ine2,1)=ak(ine2,1)+conj(FOxpoints(j*par,ine2))*F1(j,1)*Dx;
     % ak2(ine2,1)=ak2(ine2,1)+conj(FOxpoints(j*par,ine2))*F2(j,1)*Dx;
end      
  
    
% for j=floor((barini+barlen)/Dx):floor((barini+4*barlen)/Dx)                   
%       ak_1(ine2,1)=ak_1(ine2,1)+conj(FOxpoints(j*par,ine2))*F1(j,1)*Dx;
%       ak2_1(ine2,1)=ak2_1(ine2,1)+conj(FOxpoints(j*par,ine2))*F2(j,1)*Dx;
% end      
  

%     for j=floor((barini+5*barlen)/Dx):floor((barini+10*barlen)/Dx)                   
%       ak_2(ine2,1)=ak_2(ine2,1)+conj(FOxpoints(j*par,ine2))*F1(j,1)*Dx;
%       ak2_2(ine2,1)=ak2_2(ine2,1)+conj(FOxpoints(j*par,ine2))*F2(j,1)*Dx;
% end      
end   


prob_half(nt)=sum(abs(F1(floor((barini+barlen)/Dx):floor((barini+5*barlen)/Dx),1).^2));%+abs(F2(floor((barini+barlen)/Dx):floor((barini+5*barlen)/Dx),1).^2))*Dx;

Energy_1(:,nt)=abs(ak(1:ine,1)).^2;
%Energy_2(:,nt)=abs(ak2(1:ine,1)).^2;
toc
% Energy_1_1(:,nt)=abs(ak_1(1:ine,1)).^2;
% Energy_2_1(:,nt)=abs(ak2_1(1:ine,1)).^2;
% Energy_1_2(:,nt)=abs(ak_2(1:ine,1)).^2;
% Energy_2_2(:,nt)=abs(ak2_2(1:ine,1)).^2;
       
% toc
       end
       if flag_nt==0
       nt=nt+1
       flag_nt=1;
       end
    end   


        % Result plot 
    if mod(n,tpasprob)==0                                 
     n;

       figure(1001)
       hold on
       aux=abs(F1(:,3)).^2;
 %      phi_1(:,nt)=aux;
       aux4=zeros(xpoints,1);
       aux4(:)=n*Dt;
       plot3(aux4(1:xpoints),x(1:xpoints),aux(1:xpoints),'b-') 
       aux=abs(F2(:,3)).^2;
%       phi_2(:,nt)=aux;
       plot3(aux4(1:xpoints),x(1:xpoints),aux(1:xpoints),'r-') 
        axis([0 n*Dt barini_n*Dx-L/30 barend_n*Dx+L/30 0 5e6]);
        view(45,45)
%        phi_e(:,nt)=phi_1(:,nt)+phi_2(:,nt);    
    end

end
    n=n+1;
    
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(8)
hold on
plot(x,abs(F1(:,1)).^2)
ccc
steps=1;


stop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%post processing D%%%%%%%%%%%%%%%%%%
if num_sampling>1
num_phases=100;
else
  num_phases=num_sampling;
end
%D_real(:,n)=zeros(tpoints,tpoints);
for p=0:num_phases-1
   phase(p+1)=p/num_sampling; %periodo é 2*pi
   t(p+1)=phase(p+1)/om_AC*pi*2;
   t_l=t;
end
t=floor(t/Dt);


%t_p=
t_pas_2=floor(5e-16/Dt)

period_AC=1/om_AC_n*pi*2;
tpoints_period_AC=floor(period_AC/Dt);
time_t_pas_2=t_pas_2*Dt;

t_x=linspace(0,t(end)*Dt,floor(tpoints_period_AC*floor(num_phases/num_sampling)/t_pas_2)+1);

clear plus_f
clear minus_f
clear par_plus
clear par_minus

nntt=1

for n=1:t_pas_2:tpoints_period_AC*floor(num_phases/num_sampling)
    if num_phases>1
    for p=1:num_phases

        if p<=num_sampling
        if n>t(p) && n<=t(p+1)
        plus_f(nntt)=p+1;
        minus_f(nntt)=p;
        par_plus(nntt)=((t(p)-n)/(t(p)-t(p+1)));
        par_minus(nntt)=((n-t(p+1))/(t(p)-t(p+1)));
        %D_tot(:,n)=par_minus*D(:,minus_f)+par_plus*D(:,plus_f);
        end
        else if ceil(p/num_sampling)>1
        p2=p-num_sampling*floor(p/num_sampling);
        if n>t(p) && n<=t(p+1)
        plus_f(nntt)=p2+1;
        minus_f(nntt)=p2;
        par_plus(nntt)=((t(p)-n)/(t(p)-t(p+1)));
        par_minus(nntt)=((n-t(p+1))/(t(p)-t(p+1)));
        end
            end
        end
    end
    else
       plus_f(nntt)=1 ;
       minus_f(nntt)=1;
       par_plus(nntt)=1;
       par_minus(nntt)=0;
    end
 nntt=nntt+1;   
end
nntt=1;
for n=1:t_pas_2:tpoints_period_AC*floor(num_phases/num_sampling)
%    if  plus_f(n)>8
%        plus_f(n)=plus_f(n)-8;
%    end
      if  minus_f(nntt)==0
%           stop
       minus_f(nntt)=num_sampling;
      end
         if  plus_f(nntt)==num_sampling+1
%           stop
       plus_f(nntt)=1;
         end
   nntt=nntt+1;   
end

t_x=linspace(0,t(end)*Dt,length(minus_f));

figure(22)
plot(t_x,minus_f,'b')
hold on
plot(t_x,plus_f,'r')

% % t_pas_2=ceil(tpoints/length(D(:,1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_real=1;
for n=1:floor((tpoints_period_AC*floor(num_phases/num_sampling))/(t_pas_2))
D_real(:,n_real)=par_minus(n).*D(:,minus_f(n))+par_plus(n).*D(:,plus_f(n));
    %n_real is related to the injection time
n_real=n_real+1;
end


figure(33)
mesh(D_real(:,1:50:end))

clear D_real_large

nt=0;
tic
for n=1:10:length(D_real(1,:))
nt=nt+1;
 D_real_large(nt,:)=[zeros(n,1)', D_real(:,n)', zeros(length(D_real(1,:))-n-1,1)'];
end
toc

figure(4)
mesh(D_real_large(1:10:end,1:10:end))


D_tot=zeros(length(D_real(1,:)),1);

for n=1:length(D_real_large(1,:))
    D_tot(n)=sum(D_real_large(:,n))*ddtt;
end
figure(5)
plot(D_tot)
hold on

% stop
% clear D_real_large
% clear par_minus
% clear par_plus
% 
% D_tot_real=D_tot(1:end-304)-D_tot(305:end);
% %D_tot_real=D_tot-D_tot;
% figure(7)
% plot(D_tot_real)
% i = 10;
% % fname = 
% % nom=sprintf('alpha=0_f_%dampl_%d.mat', om_AC(ii), ampl(jj));
% % save(nom,'D_tot_real','D_tot')

%%%%%%%%%%%%%%%%%%%%%%%%%%save files

D_ampl(jj,:)=D_tot;

end
end

for jj=1:length(D_ampl(:,1))
D_SD(jj,:)=D_ampl(jj,1:end-floor(tpoints_period_AC/2/t_pas_2))-D_ampl(jj,floor(tpoints_period_AC/2/t_pas_2)+1:end);
end
figure(9)
plot(t_x,D_SD(3,1:length(t_x)))
figure(88)
plot(D_ampl(jj,1:end-floor(tpoints_period_AC/2/t_pas_2)))
hold on
plot(D_ampl(jj,floor(tpoints_period_AC/2/t_pas_2)+1:end))
lenlen2=length(Energy_1(1,:));

lenlen2=length(Energy_1(1,:));
%%%%%%%%%%%%%calcolo probability

for jj=1:lenlen2
prob_Ene1_1(jj)=sum(Energy_1(2000:3000,jj));%*Emaxim/k_steps;
prob_Ene1_2(jj)=sum(Energy_2(2000:3000,jj));%*Emaxim/k_steps;

prob_Energy_2_1(jj)=sum(Energy_1(3000:end,jj));%*Emaxim/k_steps;
prob_Energy_2_2(jj)=sum(Energy_2(3000:end,jj));%*Emaxim/k_steps;

% prob_Energy_3_1(jj)=sum(Energy_1(8300:end,jj));%*Emaxim/k_steps;
% prob_Energy_3_2(jj)=sum(Energy_2(8300:end,jj));%*Emaxim/k_steps;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Prob_norm(jj)=prob_Ene1_1(jj)+prob_Ene1_2(jj)+prob_Energy_2_1(jj)+prob_Energy_2_2(jj);
end

figure(13)
plot(temps,prob_Ene1_1)
hold on
plot(temps,prob_Ene1_2)

figure(14)
plot(temps,prob_Energy_2_1)
hold on
plot(temps,prob_Energy_2_2)

figure(15)
plot(temps,prob_Ene1_1.*(1./Prob_norm))
hold on
plot(temps,prob_Ene1_2.*(1./Prob_norm))

figure(16)
plot(temps,prob_Energy_2_1.*(1./Prob_norm))
hold on
plot(temps,prob_Energy_2_2.*(1./Prob_norm))

% % % figure(15)
% % % plot(temps,prob_Energy_3_1)
% % % hold on
% % % plot(temps,prob_Energy_3_2)

% figure(131)
% plot(temps,prob_Ene1_1_QW1)
% hold on
% plot(temps,prob_Ene1_2_QW1)
% 
% figure(141)
% plot(temps,prob_Ene2_1_QW1)
% hold on
% plot(temps,prob_Ene2_2_QW1)

% % % figure(151)
% % % plot(temps,prob_Ene3_1_QW1)
% % % hold on
% % % plot(temps,prob_Ene3_2_QW1)


% figure(132)
% plot(temps,prob_Ene1_1_QW2)
% hold on
% plot(temps,prob_Ene1_2_QW2)
% 
% figure(142)
% plot(temps,prob_Ene2_1_QW2)
% hold on
% plot(temps,prob_Ene2_2_QW2)


% % % % figure(2211)
% % % % plot(temps,prob_Enet_QW1)
% % % % hold on
% % % % plot(temps,prob_Enet_QW2)


% 
% figure(15)
% plot(temps,prob_Ene1_1_norm)
% hold on
% plot(temps,prob_Ene1_2_norm)
% 
% figure(16)
% plot(temps,prob_Energy_2_1_norm)
% hold on
% plot(temps,prob_Energy_2_2_norm)

figure(17)
plot(temps,prob_half)


%exit
stop







