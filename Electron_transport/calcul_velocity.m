close all
clear all

pas=10;

tpoints=1500000;
T_total=6e-12;
Dt=T_total/tpoints;

xpoints=4000;
Lx=18e-09;
Dx=Lx/xpoints;

m=0.1*(9.109e-31);  % Particle mass
Dq=3*sqrt(m)*Dx;
qpoints=4000;
Lq=qpoints*Dq;

k_a=0;
k_b=1;
k_c=1;
k_d=0;

hbar=1.054e-34;
qe=1.6e-19;
numlim=1000;
tempstra(1:tpoints)=[1:tpoints]*Dt;

E1=1^2*(pi^2)*hbar^2/(2*m*Lx^2);
E2=2^2*(pi^2)*hbar^2/(2*m*Lx^2);
omega=(E2-E1)/hbar;
alpha=0.0005*qe;

Ea=E1+0.5*hbar*omega;
Ed=E2+1.5*hbar*omega;
Ec=E1+1.5*hbar*omega;
Eb=E2+0.5*hbar*omega;

omegaa=alpha/hbar;
omegat=sqrt(omegaa^2+omega^2);

cteaa=0.5;
cteab=-(omega-omegat)/(omega+omegat)*cteaa;

cteba=0.5;
ctebb=0.5;

for j=1:xpoints    
     x(j)=j*Dx-Lx/2;
     inf_1(j)=sqrt(2/Lx)*cos(pi*x(j)/Lx);
     inf_2(j)=sqrt(2/Lx)*sin(2*pi*x(j)/Lx);
end

for j=1:qpoints
    q(j)=j*Dq-Lq/2;
    harm_1(j)=(omega/(pi*hbar))^(1/4)*exp(-q(j)^2*omega/(2*hbar));
    harm_2(j)=(omega/(4*pi*hbar))^(1/4)*(omega/hbar)^(1/2)*exp(-q(j)^2*omega/(2*hbar))*2*q(j);
end
xcen=0;


t(1)=Dt;
j=1;
ap(j)=(cteaa*exp(1i*omegat*t(j))+cteab*exp(-1i*omegat*t(j)))*exp(-1i*(Ed+Ea)/(2*hbar)*t(j));
dp(j)=((omega-omegat)/omegaa*cteaa*exp(1i*omegat*t(j))+(omega+omegat)/omegaa*cteab*exp(-1i*omegat*t(j)))*exp(-1i*(3*Ed-Ea)/(2*hbar)*t(j));
bp(j)=(cteba*exp(1i*omegaa*t(j))+ctebb*exp(-1i*omegaa*t(j)))*exp(-1i*(Ec)/(hbar)*t(j));
cp(j)=(-cteba*exp(1i*omegaa*t(j))+ctebb*exp(-1i*omegaa*t(j)))*exp(-1i*(Ec)/(hbar)*t(j));
i1h1(1)=ap(j)*k_a;
i2h1(1)=bp(j)*k_b;
i1h2(1)=cp(j)*k_c;
i2h2(1)=dp(j)*k_d;

tic
% for jj=1:xpoints
%    for jjj=1:qpoints
% F1(jj,jjj)=i1h1(1)*inf_1(jj)*harm_1(jjj)+i1h2(1)*inf_1(jj)*harm_2(jjj)+i2h1(1)*inf_2(jj)*harm_1(jjj)+i2h2(1)*inf_2(jj)*harm_2(jjj);
%    end
% end

F1(:,:)=i1h1(1)*inf_1(:)*harm_1(:)'+i1h2(1)*inf_1(:)*harm_2(:)'+i2h1(1)*inf_2(:)*harm_1(:)'+i2h2(1)*inf_2(:)*harm_2(:)';

toc
figure(3)
mesh((-floor(qpoints/2-1):20:floor(qpoints/2))*Dq,(-floor(xpoints/2-1):20:floor(xpoints/2))*Dx,abs(F1(1:20:xpoints,1:20:qpoints)).^2)
% % % position_x(j,1)=-(Lx-10*Dx)/2+(Lx-10*Dx)*(jj-1)/(pas1-1);
% % % position_q(j,1)=-Lq/3+2*Lq/3*(jjj-1)/(pas2-1);


% NORM(j)=sum(sum(abs(F1(1:xpoints,1:qpoints)).^2))*Dx*Dq



% pd=fitdist(abs(F1(:,1).^2),'Normal');

position_x_ini(1:pas,1) = Dx*datasample(-floor(xpoints/2-1):floor(xpoints/2),pas,'Weights',abs(F1(:,1).^2))
position_q_ini(1:pas,1) = Dq*datasample(-floor(qpoints/2-1):floor(qpoints/2),pas,'Weights',abs(F1(1,:).^2))

for j=1:pas
% % % position_q(j,1)=random('normal',L_q/2,a_q);
% % % toc
inp_x_ini(j,1)=floor(position_x_ini(j,1)./Dx)+1;
inp_q_ini(j,1)=floor(position_q_ini(j,1)./Dq)+1;
end

aux=zeros(1,pas)+1e32;

hold on
figure(3)
for i=1:pas
    plot3(position_q_ini,position_x_ini,aux,'r*')
end





% stop
%    end
% end

%%%%%%%%%%%%%%%%%%%%%% photon states harmonic oscillator

figure(1)
plot(q,harm_1)
hold on
plot(q,harm_2)

figure(2)
plot(x,inf_1)
hold on
plot(x,inf_2)
pause(0.1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j=1
   position_x(:,j)=position_x_ini(:,1);
   position_q(:,j)=position_q_ini(:,1);
   inp(:,j)=inp_x_ini(:,1);
   inp_q(:,j)=inp_q_ini(:,1);

for j=1:tpoints

% mean(inp)
% pause(0.5)

t(j)=(j-1)*Dt;

ap(j)=(cteaa*exp(1i*omegat*t(j))+cteab*exp(-1i*omegat*t(j)))*exp(-1i*(Ed+Ea)/(2*hbar)*t(j));
dp(j)=((omega-omegat)/omegaa*cteaa*exp(1i*omegat*t(j))+(omega+omegat)/omegaa*cteab*exp(-1i*omegat*t(j)))*exp(-1i*(3*Ed-Ea)/(2*hbar)*t(j));

bp(j)=(cteba*exp(1i*omegaa*t(j))+ctebb*exp(-1i*omegaa*t(j)))*exp(-1i*(Ec)/(hbar)*t(j));
cp(j)=(-cteba*exp(1i*omegaa*t(j))+ctebb*exp(-1i*omegaa*t(j)))*exp(-1i*(Ec)/(hbar)*t(j));

%%%%% Variables internes
i1h1(j)=ap(j)*exp(-1i*Ea*t(j)/hbar)*k_a;
i2h1(j)=bp(j)*exp(-1i*Eb*t(j)/hbar)*k_b;
i1h2(j)=cp(j)*exp(-1i*Ec*t(j)/hbar)*k_c;
i2h2(j)=dp(j)*exp(-1i*Ed*t(j)/hbar)*k_d;
%%%%%
for ip=1:pas

        control=0;
        tempo_restante=Dt;

        while control==0
        inp(ip)=floor((position_x(ip,j)+Lx/2)/Dx)+1;
        if (inp(ip)<1) 
           inp(ip)=1;
           position_x(ip,j)=x(1);
        end    
        if (inp(ip)>xpoints) 
           inp(ip)=xpoints; 
           position_x(ip,j)=x(xpoints);
        end
        inp_q(ip)=floor((position_q(ip,j)+Lq/2)/Dq)+1;
        if (inp_q(ip)<1) 
           inp_q(ip)=1; 
           position_q(ip,j)=q(1);
        end    
        if (inp_q(ip)>qpoints) 
           inp_q(ip)=qpoints; 
           position_q(ip,j)=q(qpoints);
        end
        
%%%%%%%%%%%%%%%speed in space
jj=inp(ip)+1;
jjj=inp_q(ip);
if (jj>xpoints)
Fxp1=0;   
else    
Fxp1=i1h1(j)*inf_1(jj)*harm_1(jjj)+i1h2(j)*inf_1(jj)*harm_2(jjj)+i2h1(j)*inf_2(jj)*harm_1(jjj)+i2h2(j)*inf_2(jj)*harm_2(jjj);
end
jj=inp(ip)-1;
jjj=inp_q(ip);
if (jj<1)
Fxm1=0;
else
Fxm1=i1h1(j)*inf_1(jj)*harm_1(jjj)+i1h2(j)*inf_1(jj)*harm_2(jjj)+i2h1(j)*inf_2(jj)*harm_1(jjj)+i2h2(j)*inf_2(jj)*harm_2(jjj);
end
jj=inp(ip);
jjj=inp_q(ip);
Fx=i1h1(j)*inf_1(jj)*harm_1(jjj)+i1h2(j)*inf_1(jj)*harm_2(jjj)+i2h1(j)*inf_2(jj)*harm_1(jjj)+i2h2(j)*inf_2(jj)*harm_2(jjj);

vbohm=imag((Fxp1-Fxm1)/(2*Dx*Fx));
vbohm=vbohm*hbar/m;

if vbohm>0
      dt1=abs(((inp(ip)+1)*Dx-position_x(ip,j))/vbohm);
      if dt1<Dt/numlim 
         dt1=abs(Dx/(vbohm*numlim));
      end
else
      dt1=abs((position_x(ip,j)-(inp(ip)-1)*Dx)/vbohm); 
      if dt1<Dt/numlim 
         dt1=abs(Dx/(vbohm*numlim));
      end     
end

%%%%%%%%%% speed in q
jj=inp(ip);
jjj=inp_q(ip)+1;
if (jjj>xpoints)
Fqp1=0;
else
Fqp1=i1h1(j)*inf_1(jj)*harm_1(jjj)+i1h2(j)*inf_1(jj)*harm_2(jjj)+i2h1(j)*inf_2(jj)*harm_1(jjj)+i2h2(j)*inf_2(jj)*harm_2(jjj);;
end
jj=inp(ip);
jjj=inp_q(ip)-1;
if (jjj<1)
Fqm1=0;
else
Fqm1=i1h1(j)*inf_1(jj)*harm_1(jjj)+i1h2(j)*inf_1(jj)*harm_2(jjj)+i2h1(j)*inf_2(jj)*harm_1(jjj)+i2h2(j)*inf_2(jj)*harm_2(jjj);
end
jj=inp(ip);
jjj=inp_q(ip);
Fq=i1h1(j)*inf_1(jj)*harm_1(jjj)+i1h2(j)*inf_1(jj)*harm_2(jjj)+i2h1(j)*inf_2(jj)*harm_1(jjj)+i2h2(j)*inf_2(jj)*harm_2(jjj);

vbohm_q=imag((Fqp1-Fqm1)/(2*Dq*Fq));
vbohm_q=vbohm_q*hbar;

if vbohm_q>0
            dt2=abs(((inp_q(ip)+1)*Dq-position_q(ip,j))/vbohm_q);
            if dt2<Dt/numlim 
                dt2=abs(Dq/(vbohm_q*numlim));
            end
        else
            dt2=abs((position_q(ip,j)-(inp_q(ip)-1)*Dq)/vbohm_q); 
            if dt2<Dt/numlim 
                dt2=abs(Dq/(vbohm_q*numlim));
            end     
        end
 
 %%%%%%%%%%%%%%%%%calculate and update new positions       
        tempo=min(Dt,min(dt1,min(dt2,tempo_restante)));

        position_x(ip,j)=position_x(ip,j)+vbohm*tempo;
        position_q(ip,j)=position_q(ip,j)+vbohm_q*tempo;
       
        tempo_restante = tempo_restante -tempo;
        if tempo_restante < Dt/numlim
           control=1;
           position_x(ip,j+1)=position_x(ip,j);
           position_q(ip,j+1)=position_q(ip,j);
%          inp(ip,n+1)=inp(ip,n);
%            inp_q(ip,n+1)=inp_q(ip,n);
        end
        end %while Dt
        end %particules
        
   halfcharge(j)=0;    
   for jjjj=1:pas
       if position_x(jjjj,j)>0
           halfcharge(j)=halfcharge(j)+1;
       end 
   end  
   j
end


figure(10)
for ip=1:pas
    plot(position_x(ip,1:tpoints),t)
    hold on
end
axis([-Lx/2, Lx/2, 0 t(end)])
plot(Lx*cos(omega*t)/2,t)

figure(11)
for ip=1:pas
    plot(position_q(ip,1:tpoints),t)
    hold on
end
axis([-Lq/2, Lq/2, 0 t(end)])
plot(Lq*cos(omega*t)/2,t)

figure(6)
for ip=1:pas
plot3(position_x(ip,1:tpoints),position_q(ip,1:tpoints),t)
hold on
end
figure
plot(t,abs(i1h1).^2,t,abs(i1h2).^2,t,abs(i2h1).^2,t,abs(i2h2).^2)

figure(7)
plot(t,halfcharge*1.6*1e-19,t,cos(omega*t)*pas/2)

figure(8)
plot(t,halfcharge*1.6*1e-19)