%**************************************************************************
%*              Numeric Solution to Schrodinger Equation                  *
%**************************************************************************     
% the frequency of the oscillating electron si not necessarely the same as
% the photon frequency!
clear all
close all
 
%***************** Constants definition  **********************************

pas=0;             % Nombre de particules
numlim=1000;        % Limt pel calcul de les trajectories  

Dx=5E-10;           % Spatial step: Width of cells in x axis
Dt=0.05E-16;         % Time step

L=30E-07;            % Box length

%T=1000*Dt;         % Duration of experiment
T=2.2e-12;
m=0.041*(9.109e-31);  % Particle mass
%m=0.5*(9.109e-31);
h=1.054e-34;        % Planck's constant over 2Pi
qe=1.6e-19;         % Electron charge

xpoints=floor(L/Dx);% Number of longitude elements in the box
tpoints=floor(T/Dt);% Number of lapses of time in the experiment
tstart1=floor(2e-13/Dt);  % Temps per activar l'scatering
tstart2=tstart1+floor(1e-13/Dt); % temps per desactivar l'scattering
%mac_tstrat1=tstart1; %before this moment for sure no scattering;

tpas=floor(tpoints/2000);          % Interval de temps per dibuixar 
tpas2=floor(tpoints/2000);

xcen=1.5E-07;         % Central position of the wave packet
dispersionx=1e-7;  % Spatial dispersion of the wave packet
xcen=L/2-dispersionx*3.5;         % Central position of the wave packet
Energy=0.105*qe;      % Minimum energy of wave F2 vanish in the end(g1)
pasene=0.45;
k_steps=2000;

barlen=2e-9;
barini=L/2-3*barlen;
barvalue=0.5*qe;

num_sampling=1;


ddtt=tpas*Dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

om_AC=1e12*pi*2; %om=2*pi*f=2*pi/T
ampl=0.05;

ampl_n=ampl*0;
om_AC_n=om_AC*0;

for p=0:num_sampling-1
for n=1:tpoints
V_AC(n,p+1)=ampl_n*sin(om_AC_n*n*Dt+p*2*pi/num_sampling);
end
end
figure(22)
plot(V_AC)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C1=i*h*Dt/(m*Dx^2); % Auxiliar calculations - constants
C2=-qe*i*2*Dt/h;
C4=-2*Dt/(2*m*Dx);
C5=C2/(qe*2*m);
%[ C1 C2 C4 C5]



%******************** Potentials *******************%

% % % barpoten=aosc*qe*cos(wosc*0);

for j=1:xpoints
    V_or(j) = 0;
end

for j=floor(barini/Dx):floor((barini+barlen)/Dx)
    V_or(j) = V_or(j)+barvalue;
end

for j=floor((barini+9*barlen)/Dx):floor((barini+10*barlen)/Dx)
    V_or(j) = V_or(j)+barvalue;
end
 
for j=1:xpoints    
         x(j)=j*Dx;
end
Emaxim=pasene; %don't keep it too small
kmaxim=sqrt(2*m*qe*(Emaxim)/(h^2));  %k vector of the first guess
dk=kmaxim/k_steps;
tic
[Tran, Ene_old, FOxpoints]=matrix(V_or, Emaxim, k_steps, x, L, Dx);
toc
Ene=Ene_old(2:end-2)
figure(66)
plot(Ene,Tran(4:end))




flag_ini=0;
for j=floor(xpoints/2):-1:1
    if V_or(j)~=0 && flag_ini==0;
    barini_in_n=j;
    flag_ini=1;
    end
end
barini_out_n=barini_in_n-floor(barlen/Dx);

flag_end=0;
for j=floor(xpoints/2):xpoints
    if V_or(j)~=0 && flag_end==0;
    barend_in_n=j-1;
    flag_end=1;
    end
end
barend_out_n=barend_in_n+floor(barlen/Dx);

L_well=8*barlen;
%calc the dipole
psi_1(:)=zeros(1,xpoints);
psi_2(:)=zeros(1,xpoints);
for  j=barini_in_n:barend_in_n
psi_1(j)=sqrt(2/L_well).*sin(1*pi.*(x(j)-barini_in_n*Dx)./(barend_in_n*Dx-barini_in_n*Dx));    %E_1 = 238.5 meV
psi_2(j)=sqrt(2/L_well).*sin(2*pi.*(x(j)-barini_in_n*Dx)./(barend_in_n*Dx-barini_in_n*Dx));     
end

figure(1)
plot(x,psi_1)
hold on
plot(x,psi_2)
%for  j=barini_in_n:barend_in_n
dip12  = sum(psi_1.*x.*psi_2)*Dx;            %dip =  -... m in meters 
dip12=dip12;                             %dip =  -... nm
%end

%%%%%%%%%%%%%%%%%%%estiamte alpha
LL=1e-5;
ccc=7.4;
for j=1:xpoints
cost=ccc*qe/LL; %L_well & ccc=20: 2e-10 J/m ---> 1.2500 eV/nm ---> 1250 meV/nm
                               %L & ccc=20: 8e-13 J/m ---> 0.005 eV/nm ---> 5 meV/nm
                               %L & ccc=7.4: 2.96e-13 J/m ---> 0.0019 eV/nm ---> 1.9 meV/nm
alpha(j)=abs(dip12)*cost;                     
end
                %cost=J/m and dip=m
alpha12=alpha;  %alpha=J

%freq of oscillation

f_R=alpha(1)/(h);
om=f_R*2*pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% amplitude of oscillation
stop
E_f=7;    %electric field (meV/nm)  for a well of around 10 nm, I have 0.1V
L_nm=L_well/1e-9;
E_lll=E_f*1e-3*L_nm;   
%lambda1=h*sqrt(2*m*(E_lll*qe)/h^2);
lambda=h*sqrt(2*m*(alpha(1))/h^2);
% stop
k=sqrt(2*m*Energy/h^2);                  %Auxiliar values

h^2*(k+lambda/h)^2/(2*m*qe);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V_or(:)=V(:);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flag_ini=0;
for j=1:xpoints
    if V_or(j)>1e-25 && flag_ini==0;
      barini_n=j+2; 
      barini_n_real=j;
      flag_ini=1;
    end
end

flag_end=0;
for j=xpoints:-1:1
if V_or(j)>1e-25 && flag_end==0
      barend_n=j-2;  
      barend_n_real=j;
      flag_end=1;
end
end

% Set up F1

t=0;                                     %Initial time = 0    
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

for p=1:num_sampling
    F1(:,1,p)=F1(:,1);
    F1(:,2,p)=F1(:,2);
end

norma1(2)=(F1(:,2)'*F1(:,2))*Dx;            %Set up norm = 1


% CALCULATION OF THE SYSTEM EVOLUTION

control=1;
nt=1;
temps(1)=0;
temps(2)=Dt;

norma=(F1(:,2)'*F1(:,2))*Dx;

n=3;

%%%%%%%%%% Calcul de la distribucio de moments  ************
        for indE=1:2*k_steps
 %        Ene(indE)=pasene*(indE-k_steps)/k_steps;    
        kvec(indE)=sign(Ene(indE))*sqrt(2*m*qe*abs(Ene(indE))/h^2);     
        ak(indE)=0;    
        for j=2:(xpoints-1)  
           ak(indE)=ak(indE)+F1(j,2)*exp(-1i*kvec(indE)*x(j))*Dx; 
        end
        end
        figure(500)
        plot(Ene,abs(ak).^2,'red')
        hold on
figure(7)
subplot(2,1,1)
plot(Ene,log(Tran(1:length(Ene))))
hold on
subplot(2,1,2)
plot(Ene,log(abs(ak).^2))

for j=1:pas
position(j,n)=random('normal',xcen,a/2);
inp(j,n)=floor(position(j,n)./Dx)+1;
end
if pas>1
figure(600)
[hy,hx]=hist(position(:,n),150);
plotyy(hx,hy,x,abs(F1(:,2)).^2)
end
tempstra(1:tpoints)=[1:tpoints]*Dt;
 C1=i*h*Dt/(m*Dx^2); %Calculation of the enewrgy
      %energia
   for par=1:1
       n=0;
       nt=1;
   %om=0.25*par*(Energy/qe-0.028)*qe/h; %the frequency of oscillation of theta then will be om/2/pi
   clear D
   clear D_prim
   
  %om=(Energy/qe-0.028)*qe/h;

% % % %    figure(111)
% % % %    plot(x_cm)
% % % %    stop
        V(:,:)=V_or(:)/qe;
        tic
      for n=1:tpoints
              if n<tstart1
      theta(n)=0;
      else if n>=tstart1 && n<tstart2
      theta(n)=lambda*cos(om*n*Dt)*(tstart1-n)/(tstart1-tstart2); 
      else
      theta(n)=lambda*cos(om*n*Dt); %I thried more but it is not wokring, that's it, I can only see the nu oscillation, non the Rabi.
      end
      end
        
      E_theta(n)=h^2*theta(n)^2/2/m/h^2/qe;     %the energy related to the electric field
      tttheta(n)=theta(n)^2;
        end
        toc
       
     n=3;
while (n<tpoints)   
     par_j(:)=zeros(xpoints,1);
    flag_nt=0;
      Ctheta4=C4*theta(n);
      %Ctheta5=C5*tttheta(n);

     %for p=1:num_sampling 

       for j=barini_n_real-15:barini_n_real+15
        V(j,:)=V_or(j)/qe+V_AC(n,:)*qe*(j-(barini_n_real-15))/(barend_n_real-barini_n_real+30)/qe;        
       end
       for j=barini_n_real+16:xpoints
               V(j,:)=V_or(j)/qe+V_AC(n,:)*qe;
       end
%        tic
%        V_C2(:,:)=;
%        toc
     %end

      
      
    %Calculation of new wavefunction F1 and F2   
    F1(1,3,:)=0;       
    F1(xpoints,3,:)=0;


%   tic
 
  F1(2:(xpoints-1),3,:)=F1(2:xpoints-1,1,:)+C1*(F1(1:xpoints-2,2,:)-2*F1(2:xpoints-1,2,:)+F1(3:xpoints,2,:))+C2*V(2:xpoints-1,:).*F1(2:xpoints-1,2,:);%+Ctheta4*(F1(3:(xpoints),2,:)-F1(1:(xpoints-2),2,:));%+Ctheta5*F1(2:(xpoints-1),2,p); 
   %AUX(:)=(F1(3:(xpoints),2,p)-F1(1:(xpoints-2),2,p))
  for j=barini_in_n:barend_in_n   %Schrodinger numeric solution for all points but the ends
      F1(j,3,:)=F1(j,3,:)+Ctheta4*(F1(j+1,2,:)-F1(j-1,2,:));%+Ctheta5*F1(j,2,:);
      par_j(j)=1;
  end
  
  for j=barini_out_n:barini_in_n   %Schrodinger numeric solution for all points but the ends
      par_j(j)=(barini_out_n-j)/(barini_out_n-barini_in_n);
      F1(j,3,:)=F1(j,3,:)+par_j(j)*Ctheta4*(F1(j+1,2,:)-F1(j-1,2,:));%+Ctheta5*F1(j,2,:);
  end

  for j=barend_in_n:barend_out_n   %Schrodinger numeric solution for all points but the ends
      par_j(j)=(barend_out_n-j)/(barend_out_n-barend_in_n);
      F1(j,3,:)=F1(j,3,:)+par_j(j)*Ctheta4*(F1(j+1,2,:)-F1(j-1,2,:));%+Ctheta5*F1(j,2,:);
  end

%   toc  

     for p=1:num_sampling 


        if mod(n,tpas)==0                                 
        
        temps(nt)=n*Dt;            
      %phi_e(:,nt)=abs(F1(:,3)).^2;

        %Create time vector
        %Calculation of the norm -- Step advance of auxiliar functions F1 and F2
% 
%        % norma1(nt)=(F1(:,3)'*F1(:,3))*Dx;
%  
%         %moment
%         momentum1(nt)=0;
%         momentum2(nt)=0;
%         for j=2:(xpoints-1)  
%            momentum1(nt)=momentum1(nt)-1i*h*conj(F1(j,2))*(F1(j+1,2)-F1(j-1,2))/2; 
%         end
%         momentum2(nt)=lambda*theta(n);

if mod(n,tpas2)==0 
      aux=abs(F1(:,3)).^2;
      phi_e(:,nt)=aux;
      if p==floor(num_sampling/2) || num_sampling==1
      figure(100)
      
      hold on

      
      aux4=zeros(xpoints,1);
      aux4(:)=n*Dt;
      plot3(aux4(1:xpoints),x(1:xpoints),aux(1:xpoints),'b-')  
      axis([0 n*Dt barini_n*Dx-L/20 barend_n*Dx+L/20 0 5e6]);
      view(45,45)
      hold on
%       
      end
      
% %         for indE=1:2*k_steps
% % %        Ene(indE)=pasene*(indE-k_steps)/k_steps;    
% %        % kvec(indE,nt)=sign(Ene(indE))*sqrt(2*m*qe*abs(Ene(indE))/h^2);     
% %         ak(indE,nt)=0;    
% % %         for j=2:(xpoints-1)  
% %         for j=barini_n:barend_n
% %            ak(indE,nt)=ak(indE,nt)+F1(j,2)*conj(FOxpoints(j,indE))*Dx; 
% %           % ak(indE,nt)=ak(indE,nt)+F1(j,2)*exp(-1i*kvec(indE,nt)*x(j))*Dx; 
% %         end
% %         end
% %  
% % 
% %        figure(6)
% %       hold on
% %       aux6=zeros(indE,1);
% %       aux6(:)=n*Dt;
% %       plot3(aux6(1:indE),Ene,abs(ak(:,nt)).^2,'b-')  
% %       hold on
% % 
% %       Energy_1(:,nt)=abs(ak(:,nt)).^2;
      end
      
      
      %%%%%%%%%%%%%%%%%estimation D
    F_t_1(:,p)= F1(barini_in_n:barend_in_n,3,p);
         lenlen=length(F_t_1(:,p));
        for j=2:lenlen-1
        DER_F_t_1(j,p)=(F_t_1(j+1,p)-F_t_1(j-1,p));
        end
        DER_F_t_1(:,p)=DER_F_t_1(:,p)/(2*Dx);       
        for j=2:lenlen-1
        J_1(j,p)=real(conj(F_t_1(j,p)).*(1/1i).*DER_F_t_1(j,p));%/abs(F_t(j))^2
        end
        J(:,p)=(h/m).*J_1(:,p);

        
        D_prim(nt,p)=sum(J(:,p));%*Dx
        D(nt,p)=D_prim(nt,p)/(lenlen);%/Dx
        figure(30)
        if p==1
        hold off
        end
        plot(temps(1:nt),D(1:nt,p))
        hold on
      
      
       if flag_nt==0
       nt=nt+1
       flag_nt=1;
       end

%       nt=nt+1
        end
    
     end
     %toc
%       tic
        F1(:,1,:)=F1(:,2,:);                                %Auxiliar F2 shifts to auxiliar F1
        F1(:,2,:)=F1(:,3,:);                                %Main F3 shifts to auxiliar F2
%        toc
   n=n+1;

end
D_par(par,:)=D(:,1);

   end
t_theta=linspace(0,T,tpoints)
figure(1)
plot(t_theta(1:end-1),theta)


%%%%%%%%%%%%%%%%%%%%%%%post processing for D

num_phases=100;
%D_real(:,n)=zeros(tpoints,tpoints);
for p=0:num_phases-1
   phase(p+1)=p/num_sampling; %periodo é 2*pi
   t(p+1)=phase(p+1)/om_AC*pi*2;
   t_l=t;
end
t=floor(t/Dt);


%t_p=
t_pas_2=floor(1e-16/Dt)

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

clear D_real_large_square

nt=0;
tic
for n=1:10:length(D_real(1,:))
nt=nt+1;
 D_real_large_square(nt,:)=[zeros(nt,1)', D_real(1:10:end,n)', zeros(floor(length(D_real(1,:))/10)-nt-1,1)'];
end
toc
% % 
% % figure(4)
% % mesh(D_real_large(1:10:end,1:10:end))

figure(44)
mesh(D_real_large_square(1:end,1:end))

D_tot_sq=zeros(length(D_real(1,:)),1);

for n=1:length(D_real_large_square(1,:))
    D_tot_sq(n)=sum(D_real_large_square(:,n))*ddtt;
end
% % figure(5)
% % plot(D_tot)
% % hold on

figure(55)
plot(D_tot_sq)
%%%%%%%%%%%%%%%%%%%%%%% see the eenrgy components

lenlen=length(Energy_1(1,:));

for jj=1:lenlen
prob_Ene_lev_1(jj)=sum(Energy_1(2000:2500,jj));%*Emaxim/k_steps;
prob_Ene_lev_2(jj)=sum(Energy_1(2500:2800,jj));%*Emaxim/k_steps;
Prob_norm(jj)=prob_Ene_lev_1(jj)+prob_Ene_lev_2(jj);
end

figure(14)
plot(temps(1:length(prob_Ene_lev_1)),prob_Ene_lev_1)
hold on
plot(temps(1:length(prob_Ene_lev_1)),prob_Ene_lev_2)

figure(16)
plot(temps(1:length(prob_Ene_lev_1)),prob_Ene_lev_1./Prob_norm,'r-')

figure(17)
plot(temps(1:length(prob_Ene_lev_1)),prob_Ene_lev_2./Prob_norm)

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
 plot3(tempstra(3:50:tpoints),position(index,3:50:tpoints),aux2(3:50:tpoints),'r-')
 end
 stop
 