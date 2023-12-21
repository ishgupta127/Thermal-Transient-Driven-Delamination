clear all
%close all
clc
prompt = 'Enter 1 for dirichlet, 2 for euler';
bctype = input(prompt);
prompt = 'Enter 1 for plane stress, 2 for plane strain,3 for biaxial state';
x = input(prompt);

nlayer = 2; %Number of material layers
tboun=[2700 300]; %Initial temperature of the boundaries, a steady-state profile satisfying these is used as initial temperature profile
tspan = [0 0.01933e-3];% Time span for simulation
   
tau=0.5; %cooling rate, refer paper "Critical cooling rates to avoid transient-driven cracking in thermal barrier coating (TBC) systems"
tamb=300; %Ambient temperature

nel=[300 200]; %Number of elements in each layer

ho=[1e-4 (1e-3)-(1e-4)];  %Thickness of each layer
pcp=[4e6 4e6]; %rho*cp value of each layer

k=[385 385]; %conductivity of each layer

e=[110e9 110e9];   %Young's modulus of each layer

nu=[0.3 0.3];   %Poisson's ratio of each layer

alpha=[24.8e-6 24.8e-6];  %TEC of each layer

esize=ho./nel; %Element size in each layer

if x==1
    ebar=e;
    c=ones(1,nlayer);
elseif x==2
    ebar=e./(1-nu.^2);
    c=1.+nu;
elseif x==3
    ebar=e./(1-nu);
    c=ones(1,nlayer);      
end

eltot=sum(nel);
M=zeros(eltot+1,eltot+1);
K=zeros(eltot+1,eltot+1);
F=zeros(eltot+1,1);

q=1;
for layer=1:1:nlayer
    
    pcpl=pcp(layer);
    esizel=esize(layer);
    kl=k(layer);

for index=q:1:sum(nel(1:layer))
    
    M(index:index+1,index:index+1)=(pcpl*esizel*(1/6)*[2 1;1 2])+M(index:index+1,index:index+1);
    K(index:index+1,index:index+1)=((kl/esizel)*[1 -1;-1 1])+K(index:index+1,index:index+1);
    F(index:index+1,1)=[0;0]+F(index:index+1,1);
end

q=sum(nel(1:layer))+1;
end

if bctype==1

%T=sym('T', [1 nlayer+1]);
T(1)=tboun(1);
T(2)=(((tboun(2)-tboun(1))/1e-3)*ho(1))+tboun(1);
T(3)=tboun(2);

% for index=1:1:nlayer
%     term(index)=(k(index)*(T(index)-T(index+1)))/ho(index); 
% end
% for num=1:1:nlayer-1
%     eqn(num)=term(num+1)-term(num);
% end
% S = solve(eqn==0);
% 
% if nlayer~=2
% C = struct2cell(S);
% end
% 
% index=2;
% for i=1:1:nlayer-1
%     if nlayer~=2
%     T(index)=C{i,1};
%     else
%       T(index)=S;
%     end
%     
%     index=index+1;
% end

q=1;
for layer=1:1:nlayer
    
for index=q:1:sum(nel(1:layer))+1
if layer~=1
position(index)=sum(ho(1:layer-1))+esize(layer)*(index-q);
else
    position(index)=esize(layer)*(index-q);
end
end

for flag=1:1:length(position)
    y0(flag)=tboun(1)+((tboun(2)-tboun(1))/(1e-3))*position(flag);
end
% figure
% plot(position,y0)

q=sum(nel(1:layer))+1;
end
ti=y0;
y0(1)=[];
y0(sum(nel))=[];

[t,y] = ode15s(@(t,y) odefcnew(t,y,bctype,M,K,F,nel,tau,tamb,tboun), tspan, y0);

if tspan(end)<=0.0129e-3
   t1=tboun(1)+((2351-tboun(1))*tspan(end))/(0.0129e-3);
    t2=(tboun(2));
elseif tspan(end)>0.0129e-3 && tspan(end)<=0.01933e-3
t1=2351+(((3280-2351))/(6.43e-6))*(tspan(end)-0.0129e-3);
    t2=(tboun(2));
else
    t1=3280+(((2700-3280))/(3.225e-6))*(tspan(end)-0.01933e-3);
    t2=(tboun(2));
end

tfield=[t1 y(end,:) t2];

end

figure
plot(position,tfield)

hintpos(1)=0;
for int=2:1:nlayer+1
    hintpos(int)=hintpos(int-1)+ho(int-1);
end

%Below part of the code is for calculations of the stresses and SERR 
u1(1)=0;
last=eltot+1;
main=0;
for node=1:1:last
    main=main+1;
   
    for layer=1:1:nlayer
        if position(node)>hintpos(layer) && position(node)<hintpos(layer+1)
        clayer=layer;
        clayer2=layer;
        break;
        elseif position(node)==hintpos(layer+1)
            clayer=layer;
            clayer2=layer+1;
          break;
        else
            clayer2=1;
            
        end
    end
     if node~=1
        a11=0;
        a12=0;
        a21=0;
        a22=0;
        b1=0;
        b2=0;
for tl=1:1:clayer
    
                alphal=alpha(tl);
                ebarl=ebar(tl);
                cl=c(tl);
                esizel=esize(tl);
       
       if tl==1
           ini=1;
       else
           ini=sum(nel(1:tl-1))+1;
       end
       
       if tl==clayer
           fin=node-1;
       else
           fin=sum(nel(1:tl));
       end
       
       
       for num=ini:1:fin
                       yt=position(num+1);
            yb=position(num);
            tht=alphal*(tfield(num+1)-ti(num+1));
            thb=alphal*(tfield(num)-ti(num));
            a11=a11+(ebarl*esizel);
            a12=a12-((ebarl*((yt)^2-(yb)^2)/2));
            a21=-a12;
            a22=a22-((ebarl*((yt)^3-(yb)^3)/3));
            b1=b1+(1/2)*((cl*ebarl*esizel*(tht+thb)));
            b2=b2+((cl*ebarl*esizel*((esizel*(1/3)*(tht-thb))+((1/2)*(tht*yb+thb*yt)))));
           
       end   
end

          A=[a11 a12;a21 a22]; 
          b=[b1; b2];
          para = A\b;
          
       sig1=[];
       sigz1=[];
       se1=[];
       x1=[];
          flag=1;
for tl=1:1:clayer
    
                alphal=alpha(tl);
                ebarl=ebar(tl);
                cl=c(tl);
                nul=nu(tl);
                el=e(tl);
       
       if tl==1
           ini=1;
       else
           ini=sum(nel(1:tl-1))+1;
       end
       
       if tl==clayer
           fin=node;
       else
           fin=sum(nel(1:tl))+1;
       end

       for index=ini:1:fin
            
                 y=position(index);
        th=alphal*(tfield(index)-ti(index));
           
    sig1(flag)=ebarl*(para(1)-(para(2)*y)-(cl*th));
    if x==2
    sigz1(flag)=nul*sig1(flag)-el*(th);
    end
    se1(flag)=(((sig1(flag))^2+(sigz1(flag))^2-(2*nul*(sig1(flag))*(sigz1(flag))))/(2*el));
    x1(flag)=y;
    flag=flag+1;
           
       end 
 
       
       u1(main)=trapz(x1,se1);
end

    end


   %%other side of the crack
   
   if node~=eltot+1
   
        a11=0;
        a12=0;
        a21=0;
        a22=0;
        b1=0;
        b2=0;
for tl=clayer2:1:nlayer

    
                alphal=alpha(tl);
                ebarl=ebar(tl);
                cl=c(tl);
                esizel=esize(tl);
       
       if tl==clayer2
           ini=node;
       else
           ini=sum(nel(1:tl-1))+1;
       end
       
           fin=sum(nel(1:tl));
       
       
       for num=ini:1:fin
                       yt=position(num+1);
            yb=position(num);
            tht=alphal*(tfield(num+1)-ti(num+1));
            thb=alphal*(tfield(num)-ti(num));
            a11=a11+(ebarl*esizel);
            a12=a12-((ebarl*((yt)^2-(yb)^2)/2));
            a21=-a12;
            a22=a22-((ebarl*((yt)^3-(yb)^3)/3));
            b1=b1+(1/2)*((cl*ebarl*esizel*(tht+thb)));
            b2=b2+((cl*ebarl*esizel*((esizel*(1/3)*(tht-thb))+((1/2)*(tht*yb+thb*yt)))));
           
       end   
end

          A=[a11 a12;a21 a22]; 
          b=[b1; b2];
          para = A\b;
          
       sig2=[];
       sigz2=[];
       se2=[];
       x2=[];
       flag=1;
for tl=clayer2:1:nlayer
    
                alphal=alpha(tl);
                ebarl=ebar(tl);
                cl=c(tl);
                nul=nu(tl);
                el=e(tl);
       
       if tl==clayer2
           ini=node;
       else
           ini=sum(nel(1:tl-1))+1;
       end
       
           fin=sum(nel(1:tl))+1;
       
       for index=ini:1:fin
            
                 y=position(index);
        th=alphal*(tfield(index)-ti(index));
           
    sig2(flag)=ebarl*(para(1)-(para(2)*y)-(cl*th));
    if x==2
    sigz2(flag)=nul*sig2(flag)-el*(th);
    end
    se2(flag)=(((sig2(flag))^2+(sigz2(flag))^2-(2*nul*(sig2(flag))*(sigz2(flag))))/(2*el));
    x2(flag)=y;
    flag=flag+1;
          
       end 
%        figure
%        plot(x2,sig2/10^6)
       
end
u2(main)=trapz(x2,se2);
   end 
    pos(main)=position(node);
end

if last==eltot+1
    u2(main)=0;
end
   
    for main=1:1:length(u1)
   setot(main)=u1(main)+u2(main); 
   
    end
  
for i=1:1:length(u1)
    serr(i)=setot(1)-setot(i);
end
%     
% figure
% plot(pos,u1)
% figure
% plot(pos,u2)
% figure
% plot(pos,setot)
figure
plot(pos,serr)
      figure
       plot(x1,sig1/10^6)