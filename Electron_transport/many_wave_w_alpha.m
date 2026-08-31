%**************************************************************************
%*              Numeric Solution to Schrodinger Equation                  *
%**************************************************************************     
% stop
 clear all
 close all

 


%ampl=0.02;
% om_AC=3.33e12
% om_AC=logspace(11,13,1);


%f=[2 5 7.5 10 12.5 15 17.5 20 22.5 25 27.5 30 40]*1e11;
f=[3 5 7.5 10 12.5 15 17.5 20 22.5 25 27.5 30 40]*1e11;
for ii=1:length(f);


    
clearvars -except ii
% clearvars -except ampl
% clearvars -except f
% clearvars -except nt
% clearvars -except jj
    close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=[3 5 7.5 10 12.5 15 17.5 20 22.5 25 27.5 30 40]*1e11;
%f=[5 7.5 10 12.5 15 17.5 20 25 30 40 50]*1e11;
    nt=0;
ampl=0.001;
for jj=1:1
%***************** Constants definition  **********************************
rr=1;
pas=0;
Dx=6E-10;           % Spatial step: Width of cells in x axis
Dt=0.05E-16;         % Time step
L_active=14e-07;
% L_active=7e-07;
L_abs=2E-07;
L=L_active+2*L_abs;            % Box length
n_s=1;
%T=12000*Dt;         % Duration of experiment
T=50000*Dt;         % Duration of experiment
T=1.2e-12;
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
tpas=floor(tpoints/3000);          % Interval de temps per dibuixar
ddtt=tpas*Dt;
tpasprob=floor(tpoints/1);          % Interval de temps per dibuixar 
ttt=floor(tpoints/20);   
om=(0.128-0.07)*qe/h;
om=(0.113-0.028)*qe/h;
Dq=0.5*sqrt(m_1)*Dx;
maxq=1000;
L_q=maxq*Dq;

num_sampling=16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ii jj]
ampl_n=ampl(jj);
om_AC=f(ii)*pi*2; %om=2*pi*f=2*pi/T
 om_AC_n=om_AC(rr);

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

dispersionx=1e-7;  % Spatial dispersion of the wave packet
% dispersionx=0.5e-7;  % Spatial dispersion of the wave packet
%dispersionx=0.5e-7;  % Spatial dispersion of the wave packet
% dispersionx=35*Dx;  % Spatial dispersion of the wave packet
% dispersionx=45*Dx;
xcen=L/2-3.5*dispersionx

Einj=0.132*qe;      % 0.017 first level 0.07 second level
%Einj=(ii*0.005+0.01)*qe
%Einj=0.07*qe;      % 0.017 first level 0.07 second level
E_final=0.084*qe;%084*qe;   % 0.084 first velevel 0.306 second level


 %alpha(j)=5e-1*qe/L;

k=sqrt(2*m_1*Einj/h^2); 
k_final=sqrt(2*m_1*E_final/h^2) ;      %final k of the scattered wavefunction
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

for j=floor(barini/Dx):floor((barini+barlen)/Dx)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    V_or(j) = V_or(j)+barvalue;
end


for j=floor((barini+4*barlen)/Dx):floor((barini+5*barlen)/Dx)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    V_or(j) = V_or(j)+barvalue2;
end


for j=floor((barini+9*barlen)/Dx):floor((barini+10*barlen)/Dx)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    V_or(j) = V_or(j)+barvalue2;
end

% for j=floor((barini+9.5*barlen)/Dx):floor((barini+10.5*barlen)/Dx)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%     V(j) = V(j)+barvalue3*0;
% end

for j=1:xpoints    
         x(j)=j*Dx;%-(barini+barlen+barwell/2);
end

% Potential for the quantum-well alone
ccc=20; %eV!!! %ccc=30 or 50, gives me transition from 1st-2nd (combination) to 2nd-3rd!!!!
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
L_well=16e-9;

end
flag_ini=0;
for j=floor(xpoints/2):-1:1
    if V_or(j)~=0 && flag_ini==0;
    barini_in_n=j;
    flag_ini=1;
    end
end

flag_end=0;
for j=floor(xpoints/2):xpoints
    if V_or(j)~=0 && flag_end==0;
    barend_in_n=j-1;
    flag_end=1;
    end
end


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
dip12=dip12*1e9;                             %dip =  -... nm
%end
%alpha12=alpha*abs(dip12); %this is energy: meV/m*m
ccc=7.4;
ccc=0;
LL=16e-7;
for j=1:xpoints
cost=ccc*qe/L; %L_well & ccc=20: 2e-10 J/m ---> 1.2500 eV/nm ---> 1250 mEV/nm
cost=ccc*qe/LL;
                               %L & ccc=20: 8e-13 J/m ---> 0.005 eV/nm ---> 5 meV/nm
                               %L & ccc=7.4: 2.96e-13 J/m ---> 0.0019 eV/nm ---> 1.9 meV/nm
alpha(j)=abs(dip12)*cost;                     
end
                        %cost=J/m and dip=nm
alpha12=alpha;          % alpha12=J/m*nm

%per capire
om_R=alpha12(1)*1e-9/h

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
k_steps=1000;

tic
Emaxim=0.45;
kmaxim=sqrt(2*m_1*qe*(Emaxim)/(h^2));  %k vector of the first guess
dk=kmaxim/k_steps;
[Tran, Ene_old, FOxpoints]=matrix(V_or, Emaxim, k_steps, x, L, Dx);
Ene=Ene_old(2:end-2)
figure(66)
plot(Ene,Tran(4:end))
offset_scat=0;
toc

ine=length(Ene)-1;
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
plot(Ene,(Tran(1:length(Ene))))
%plot(Ene,log(1))
ylabel('T')
xlabel('Energias (eV)')

figure(6)
subplot(3,1,1)
plot(Ene,log(Tran(1:length(Ene))))
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

clear FOxpoints

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
xpoints_abs=floor(L_abs/Dx);
% norma1(2)=(F1(:,2)'*F1(:,2))*Dx;            %Set up norm = 1
% norma2(2)=(F2(:,2)'*F2(:,2))*Dx;            %Set up norm = 1
par_abs=50
par_abs2=40
for j=2:xpoints
if j>2 && j<xpoints_abs
    k_abs(j)=(1-exp(-(j-2)/par_abs2));
else if j>xpoints-xpoints_abs && j<=xpoints
    k_abs(j)=(1-exp((j-(xpoints))/par_abs2));
    else
        k_abs(j)=1;
    end
    end


% if j>2 && j<floor((barini)/Dx)
%     k_abs2(j)=(1-exp(-(j-2)/par_abs2));
% else if j>floor((barini+barwell+2*barlen)/Dx) && j<xpoints
%     k_abs2(j)=(1-exp((j-(xpoints))/par_abs2));
%     else
%         k_abs2(j)=1;
%     end
% end
%     k_abs(j)=1;
%     k_abs2(j)=1;
end
        
        figure(999)
        plot(k_abs);
%         hold on
%         plot(k_abs2);
A1=1i*2*Dt*0.5*om;
A2=1i*2*Dt*1.5*om;
C2_a=C2*alpha(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (n<tpoints1)  %evolution after the scattering 
%     for j=1:xpoints
%     if j<L_abs


%     else if j>L-L_abs
%           F1(j,1,:)=F1(j,1,:)
%           F1(j,2,:)=F1(j,2,:)
%           F1(j,3,:)=F1(j,3,:)
%           F2(j,1,:)=F2(j,1,:)
%           F2(j,2,:)=F2(j,2,:)
%           F2(j,3,:)=F2(j,3,:)
%          
%         end
%     end
        
    
    flag_nt=0;
    for p=1:num_sampling %boocle on the phases
    %update the potential profile
    for j=1:xpoints
        if j<barini_n_real-15
        V(j,1,p)=V_or(j);
        else if j<=barend_n_real+15
        V(j,1,p)=V_or(j)+V_AC(n,p)*qe*(j-(barini_n_real-15))/(barend_n_real-barini_n_real+30);        
        else
        V(j,1,p)=V_or(j)+V_AC(n,p)*qe;  
            end
        end
    end
    end
%     figure(11)
%     plot(V(:,:))

    
    %Calculation of new wavefunction F1 and F2   
    F1(1,3,:)=0;       
    F1(xpoints,3,:)=0;

    F2(1,3,:)=0;       
    F2(xpoints,3,:)=0;

%     for j=2:(xpoints-1)   %Schrodinger numeric solution for all points but the ends
%    
%         
%         F1(j,3,p)=F1(j,1,p)+C1*(F1(j-1,2,p)-2*F1(j,2,p)+F1(j+1,2,p))+C2*V(j,p)*F1(j,2,p)-1i*2*Dt*0.5*om*F1(j,2,p)+C2*alpha(j)*xp(j)*F2(j,2,p);
%         F2(j,3,p)=F2(j,1,p)+C1*(F2(j-1,2,p)-2*F2(j,2,p)+F2(j+1,2,p))+C2*V(j,p)*F2(j,2,p)-1i*2*Dt*1.5*om*F2(j,2,p)+C2*alpha(j)*xp(j)*F1(j,2,p); %%%+C1*(F2(j-1,2)-2*F2(j,2)+F2(j+1,2))
% 
%         
%     end
        F1(2:(xpoints-1),3,:)=F1(2:(xpoints-1),1,:)+C1*(F1(1:(xpoints-2),2,:)-2*F1(2:(xpoints-1),2,:)+F1(3:(xpoints),2,:))+C2*V(2:(xpoints-1),1,:).*F1(2:(xpoints-1),2,:)-A1*F1(2:(xpoints-1),2,:)+C2_a.*xp(2:(xpoints-1))'.*F2(2:(xpoints-1),2,:);
        F2(2:(xpoints-1),3,:)=F2(2:(xpoints-1),1,:)+C1*(F2(1:(xpoints-2),2,:)-2*F2(2:(xpoints-1),2,:)+F2(3:(xpoints),2,:))+C2*V(2:(xpoints-1),1,:).*F2(2:(xpoints-1),2,:)-A2*F2(2:(xpoints-1),2,p)+C2_a.*xp(2:(xpoints-1))'.*F1(2:(xpoints-1),2,:); %%%+C1*(F2(j-1,2)-2*F2(j,2)+F2(j+1,2))

          F2(:,1,:)=F2(:,2,:).*k_abs(:); 
          F1(:,1,:)=F1(:,2,:).*k_abs(:); 
                                           %Auxiliar F2 shifts to auxiliar F1
          F2(:,2,:)=F2(:,3,:);                                %Main F3 shifts to auxiliar F2                            %Auxiliar F2 shifts to auxiliar F1
          F1(:,2,:)=F1(:,3,:); 
% 
%           F1(:,1,:)=F1(:,1,:);
%           F2(:,1,:)=F2(:,1,:);

        
        %F(:,:)=F1(:,3)*harm_1(:)'+F2(:,3)*harm_2(:)';
 for p=1:num_sampling %boocle on the phases
        
    % Result plot 
    if mod(n,tpas)==0                                 
    F_t_1(:,p)= F1(barini_n:barend_n,3,p);
    F_t_2(:,p)= F2(barini_n:barend_n,3,p);
         lenlen1=length(F_t_1(:,p));
         lenlen2=length(F_t_2(:,p));
         if lenlen1~=lenlen2
             stop
         else
             lenlen=lenlen1;
         end
        for j=2:lenlen-1
        DER_F_t_1(j,p)=(F_t_1(j+1,p)-F_t_1(j-1,p));
        DER_F_t_2(j,p)=(F_t_2(j+1,p)-F_t_2(j-1,p));
        %DER_F_t(j)=(-F_t(j+2)+8*F_t(j+1)-8*F_t(j-1)+F_t(j-2))/(12*Dx); 
        end
        DER_F_t_1(:,p)=DER_F_t_1(:,p)/(2*Dx);
        DER_F_t_2(:,p)=DER_F_t_2(:,p)/(2*Dx);        
        for j=2:lenlen-1
        J_1(j,p)=real(conj(F_t_1(j,p)).*(1/1i).*DER_F_t_1(j,p));%/abs(F_t(j))^2
        J_2(j,p)=real(conj(F_t_2(j,p)).*(1/1i).*DER_F_t_2(j,p));%/abs(F_t(j))^2
        end
        J_1(:,p)=(h/m).*J_1(:,p);
        J_2(:,p)=(h/m).*J_2(:,p);
        J(:,p)=J_1(:,p)+J_2(:,p);

        
%         if nt>=2
        D_prim(nt,p)=sum(J(:,p));%*Dx
        D(nt,p)=D_prim(nt,p)/(lenlen);%/Dx
%         if mod(n,tpas*10)==0
        figure(30)
        if p==1
        hold off
        end
        plot(D(1:nt,p)')
        hold on
        
%         end
%         end

        
        
       temps(nt)=n*Dt;                  %Create time vector
       %Calculation of the norm -- Step advance of auxiliar functions F1 and F2

      
       if flag_nt==0 && p==num_sampling
       nt=nt+1;
       [nt ii]
       flag_nt=1;
       end
    end   

    
        % Result plot 
    if mod(n,tpasprob)==0                                 
     n;

       figure(1001)
       hold on
       aux=abs(F1(:,3)).^2;
       aux4=zeros(xpoints,1);
       aux4(:)=n*Dt;
%        plot3(aux4(1:xpoints),x(1:xpoints),aux(1:xpoints),'b-') 
       plot3(aux4(1:xpoints),x(1:xpoints),aux(1:xpoints),'r-') 
        axis([0 tot_time 0 10e-7 0 1e7]);
        axis([0 tot_time barini_n*Dx-L/20 barend_n*Dx+L/20 0 0.1e6]);
        view(45,45) 
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%post processing D%%%%%%%%%%%%%%%%%%

n_periods=5;
num_phases=n_periods*num_sampling;
%D_real(:,n)=zeros(tpoints,tpoints);
for p=0:num_phases
   phase(p+1)=2*pi*p/num_sampling; %periodo é 2*pi
   time_phase(p+1)=phase(p+1)/om_AC_n;
   t_l=time_phase;
end

 time_phase=floor(time_phase/ddtt);
% time_phase=ddtt;

time_t_pas_2=ddtt; %reduce sampling to post-process
t_pas_2=floor(time_t_pas_2/ddtt)

period_AC=1/om_AC_n*pi*2;
tpoints_period_AC=floor(period_AC/ddtt);
% time_t_pas_2=t_pas_2*Dt;

t_x=linspace(0,time_phase(end)*ddtt,floor(tpoints_period_AC*n_periods/t_pas_2)+1);

clear plus_f
clear minus_f
clear par_plus
clear par_minus

nntt=1
for n=1:t_pas_2:tpoints_period_AC*n_periods
    for p=1:num_phases
        if p<=num_sampling
        if n>time_phase(p) && n<=time_phase(p+1)
        plus_f(nntt)=p+1;
        minus_f(nntt)=p;
        par_plus(nntt)=((time_phase(p)-n)/(time_phase(p)-time_phase(p+1)));
        par_minus(nntt)=((n-time_phase(p+1))/(time_phase(p)-time_phase(p+1)));
        %D_tot(:,n)=par_minus*D(:,minus_f)+par_plus*D(:,plus_f);
        end
        else if ceil(p/num_sampling)>1
        p2=p-num_sampling*floor(p/num_sampling);
        if n>time_phase(p) && n<=time_phase(p+1)
        plus_f(nntt)=p2+1;
        minus_f(nntt)=p2;
        par_plus(nntt)=((time_phase(p)-n)/(time_phase(p)-time_phase(p+1)));
        par_minus(nntt)=((n-time_phase(p+1))/(time_phase(p)-time_phase(p+1)));
        end
            end
        end
    end
 nntt=nntt+1;   
end
nntt=1;
for n=1:t_pas_2:tpoints_period_AC*n_periods-1
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

t_x_2=linspace(0,time_phase(end)*ddtt,length(minus_f));

figure(22)
plot(t_x_2,minus_f,'b')
hold on
plot(t_x_2,plus_f,'r')

% % t_pas_2=ceil(tpoints/length(D(:,1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_real=1;
tic
for n=1:1:floor(tpoints_period_AC*n_periods/(t_pas_2))-1
D_real(:,n_real)=par_minus(n).*D(:,minus_f(n))+par_plus(n).*D(:,plus_f(n));
    %n_real is related to the injection time
n_real=n_real+1;
end
toc

figure(33)
mesh(D_real(:,1:100:end))

clear D_real_large

nt=0;
tic
for n=1:10:length(D_real(1,:))
nt=nt+1;
 D_real_large(nt,:)=[zeros(n,1)', D_real(:,n)', zeros(length(D_real(1,:))-n-1,1)'];
end
 toc
% 
% figure(4)
% mesh(D_real_large(1:1:end,1:10:end))
% 
% close 4
D_tot=zeros(length(D_real(1,:)),1);

for n=1:length(D_real_large(1,:))
    D_tot(n)=sum(D_real_large(:,n))*10*t_pas_2*ddtt; %J is calculated not eevry Dt, but every ddtt=tpas*Dt
end
figure(55)
plot(t_x_2,D_tot(1:length(t_x_2)))
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
end
% D_ampl(ii)=D_tot(floor(end/2));


    clear D_tot_real
    clear D_real_large
    clear D_real
    
    nom=['Res_alpha=0_ampla_0001_freq_E=E_2' num2str(f(ii)) '.mat']
    save(nom)
    
D_Ene(ii,:)=D_tot(:);

end


stop

% for jj=1:length(D_tot(:,1))
D_SD(:)=D_tot(1:end-floor(tpoints_period_AC/2/t_pas_2),1)-D_tot(floor(tpoints_period_AC/2/t_pas_2)+1:end,1);
% end
figure(10)
plot(t_x(1,1:length(D_SD)),D_SD)
figure(88)
plot(D_ampl(jj,1:end-floor(tpoints_period_AC/2/t_pas_2)))
hold on
plot(D_ampl(jj,floor(tpoints_period_AC/2/t_pas_2)+1:end))
lenlen2=length(Energy_1(1,:));

lenlen2=length(Energy_1(1,:));
%%%%%%%%%%%%%calcolo probability

for jj=1:lenlen2
prob_Ene1_1(jj)=sum(Energy_1(2000:2800,jj));%*Emaxim/k_steps;
prob_Ene1_2(jj)=sum(Energy_2(2000:2800,jj));%*Emaxim/k_steps;

prob_Energy_2_1(jj)=sum(Energy_1(2800:3200,jj));%*Emaxim/k_steps;
prob_Energy_2_2(jj)=sum(Energy_2(2800:3200,jj));%*Emaxim/k_steps;

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







