function [Tran, Ene]=WF(V, Emaxim, k_steps, x, L, Dx)
ine=1;
m=0.041*(9.109e-31);  % Particle mass
h=1.054e-34;        % Planck's constant over 2Pi
qe=1.6e-19;         % Electron charge
m_1=m;
Dt=0.05E-16;         % Time step
kmaxim=sqrt(2*m_1*qe*(Emaxim)/(h^2));  %k vector of the first guess
dk=kmaxim/k_steps;
xpoints=floor(L/Dx);% Number of longitude elements in the box
xpoints1=xpoints;
Va=V;
parpar=1;
dx=Dx/parpar;
par=floor(Dx/dx);
% xpoints1=floor(L/dx);
pos_l=3;
pos_r=xpoints1-3;
% if par==1
 Va(:)=V(:)
barlen=2e-9;
barwell=10e-9;
barini=L/2-barwell/2-barlen;
barvalue=0.3*qe;

i_barini=floor(barini/Dx);
i_barend=floor((L/2+barwell/2+barlen)/Dx);
xx=x;


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

