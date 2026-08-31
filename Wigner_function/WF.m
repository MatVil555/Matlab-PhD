 function [SS, K_ini, U_ini, K_ini_in, U_ini_in, WF_norma_ini]=WF(FFF, x, V)
 
 
Dx=5E-10;    
L=10e-07;
xpoints=floor(L/Dx);% Number of longitude elements in the box
m=0.041*(9.109e-31);  % Particle mass
h=1.054e-34;        % Planck's constant over 2Pi
qe=1.6e-19;         % Electron charge
m_1=m;
barlen=2e-9;
barwell=10e-9;
barini=L/2-barwell/2-barlen;
barvalue=0.2*qe;
 k_steps=120;% Number of k values
ine=0;
clear S
Emaxim=0.4; %don't keep it too small
kmaxim=sqrt(2*m_1*qe*(Emaxim)/(h^2));  %k vector of the first guess
dk=kmaxim/k_steps;

kom=linspace(-kmaxim,kmaxim,2*k_steps);
for i=1:2*k_steps
Eom(i)=h^2*(kom(i))^2/(2*qe*m(1));
end



S=zeros(2*k_steps+1,xpoints);
N=floor(xpoints/2);
F1_l=[zeros(floor(xpoints/2),1) ; FFF ; zeros(floor(xpoints/2),1)];

for j=-xpoints+1:xpoints-1
    for i=1:xpoints
    rho_l(i,j+xpoints)=conj(F1_l(i+round(j/2)+N))*(F1_l(i-round(j/2)+N));
    end
    i;
end


N=floor(length(x)/2);
for j=1:xpoints+2*N    
         x_l(j)=j*Dx-(barini+barlen+barwell/2)-N*Dx;
end

hh=0;
x_l_lim=2*max(x);
x_l_n=size(rho_l(1,:));
x_l(:)=linspace(-x_l_lim,x_l_lim,x_l_n(2));

for hh=1:k_steps*2-1
    aux(:)=exp(-1i*kom(hh+1)*x_l(1:end));
for g=1:xpoints
     S(hh+1,g)=rho_l(g,:)*aux(:);
%      S(hh+1,g)=rho_l(g,:)*exp(1i*kom(hh+1)*x_l(1:end))';
end
end

SS=Dx*S*(1/2/pi);
WF_norma_ini=sum(sum(abs(SS)))*Dx*dk;
K_ini=0;
U_ini=0;
K_ini_in=0;
U_ini_in=0;
for hh=1:k_steps*2
    for i=1:xpoints
K_ini=K_ini+real(SS(hh,i))*h^2*kom(hh)^2*Dx*dk/2/m;
U_ini=U_ini+real(SS(hh,i))*V(i)*Dx*dk;
    end
end
for hh=1:k_steps*2
    for i=floor(barini/Dx):floor((barini+10*barlen)/Dx)
K_ini_in=K_ini_in+real(SS(hh,i))*h^2*kom(hh)^2*Dx*dk/2/m;
U_ini_in=U_ini_in+real(SS(hh,i))*V(i)*Dx*dk;
    end
end
E_ini_in=K_ini_in+U_ini_in;
E_ini=K_ini+U_ini;

S_ini=S;


 end
%clear rho_l