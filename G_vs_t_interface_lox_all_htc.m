clear all
%close all
%Code written by Isha Gupta: last updated on 10/22/2021
clc
prompt = 'Enter 1 for plane stress, 2 for plane strain,3 for biaxial state';
x = input(prompt);

%loxt1=1; % time span for first flusing with liquid oxygen at 90.15 K
lox0=0.15;
hox1=150;% take-off burn
st1=3;
lox1=0.15;
vac1=90;%MECO
lox11=0.15;
hox2=30;
st2=6;
lox2=0.15;
vac2=210;
lox22=0.15;
hox3=20;
st3=3;
lox3=0.15;
vac3=40;
lox33=0.15;
hox4=32;
st4=3;
lox4=0.15;

templox=90;
temphox=770;

%sc=1 is for bulk heating
%sc=2 is for flash heating

sigres1=0; % Initializing stress distribution in the part above the crack
sigres2=0; % Initializing stress distribution in the part below the crack
temp=300; %Initializing temperature field

%% 
%If there is a residual stress from the formation process, calculate it separately, define it as a new quantity here and add it to the stress
%expressions
%%
tic
%%Choose the engine operation stage here
inc=0;
for laststep=1:1:4
%%
for sc=1:1:laststep
    
    tdur=[lox0 hox1 st1 lox1 vac1 lox11 hox2 st2 lox2 vac2 lox22 hox3 st3 lox3 vac3 lox33 hox4 st4 lox4];

    if sc<laststep
tspan = [0 tdur(sc)];% Time span for simulation
    else
        tspan=[0 tdur(sc)];
    end


nlayer=5; %Number of material layers
nel=[100 100 300 100 100]; %Number of elements in each layer
tboun=[300 300]; %Initial temperature of the boundaries, a steady-state profile satisfying these is used as initial temperature profile


tamb=300; %Ambient temperature


%We are only solving for the symmetric half, third entry is half-blade
%thickness.

ho=[100e-6 50e-6 10e-3 50e-6 100e-6];  %Thickness of each layer%Thickness of the blade can vary from 3-30 mm
pcp=[4000000 3933800 3600000 3933800 4000000]; %rho*cp value of each layer

cond=[7.4 70.4 11.4 70.4 7.4]; %conductivity of each layer

e=[83e9 205e9 205e9 205e9 83e9];   %Young's modulus of each layer

nu=[0.3 0.3 0.3 0.3 0.3];   %Poisson's ratio of each layer

alpha=[9.9e-6 12.6e-6 13e-6 12.6e-6 9.9e-6];  %TEC of each layer

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
    kl=cond(layer);

for index=q:1:sum(nel(1:layer))
    
    M(index:index+1,index:index+1)=(pcpl*esizel*(1/6)*[2 1;1 2])+M(index:index+1,index:index+1);
    K(index:index+1,index:index+1)=((kl/esizel)*[1 -1;-1 1])+K(index:index+1,index:index+1);
    F(index:index+1,1)=[0;0]+F(index:index+1,1);
end

q=sum(nel(1:layer))+1;
end

%%
%This part is for obtaining temperature profile

T=sym('T', [1 nlayer+1]);
T(1)=tboun(1);
T(nlayer+1)=tboun(2);
for index=1:1:nlayer
    term(index)=(cond(index)*(T(index)-T(index+1)))/ho(index); 
end
for num=1:1:nlayer-1
    eqn(num)=term(num+1)-term(num);
end
S = solve(eqn==0);

if nlayer~=2
C = struct2cell(S);
end

index=2;
for i=1:1:nlayer-1
    if nlayer~=2
    T(index)=C{i,1};
    else
      T(index)=S;
    end
    
    index=index+1;
end
if sc==1
q=1;
for layer=1:1:nlayer
    
for index=q:1:sum(nel(1:layer))+1
y0(index)=double(T(layer)+((T(layer+1)-T(layer))/nel(layer))*(index-q));
if layer~=1
position(index)=sum(ho(1:layer-1))+esize(layer)*(index-q);
else
    position(index)=esize(layer)*(index-q);
end
end

q=sum(nel(1:layer))+1;
end
ti=y0;
tfield1=[];
else
 y0=tfield1;
end

tic
options = odeset('RelTol',1e-4);
[t,y] = ode15s(@(t,y) odefcnew2_lox_allhtc(t,y,sc,M,K,F,tdur),tspan,y0,options);
toc
    
    tfield1=y(end,:);    
what1=(t+sum(tdur(1:sc-1)));
what2=(y(:,1));
%     figure
% plot(t+sum(tdur(1:sc-1)),y(:,1))

end

% figure
% plot(position,tfield1)
% 
for bec=1:1:length(t)

    tfield=y(bec,:);    
        
% what1=(position(1:sum(nel(1:2))+1)*10^6)';
% what2=(tfield(1:sum(nel(1:2))+1))';
% % % figure
% % % plot(position(1:sum(nel(1:2))+nel(3)/2+1)*1000,ti(1:sum(nel(1:2))+nel(3)/2+1))
% % hold on
% figure
% plot(what1,what2)


hintpos(1)=0;
for int=2:1:nlayer+1
    hintpos(int)=hintpos(int-1)+ho(int-1);
end

%Below part of the code is for calculations of the stresses and SERR 
u1(1)=0;
%last=eltot+1;
last=sum(nel(1:2))+1;
main=0;

%%
% sum(nel(2)+1) gives SERR at electroless Ni/IN718 interface, change this to sum(nel(1)+1)
nodeset=[eltot+1 sum(nel(1)+1)];
%%
for specs=1:1:2
    
    node=nodeset(specs);

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
            
                 ypos=position(index);
        th=alphal*(tfield(index)-ti(index));
           
    sig1(flag)=(ebarl*(para(1)-(para(2)*ypos)-(cl*th)));
    if x==2
    sigz1(flag)=nul*sig1(flag)-el*(th);
    end
    se1(flag)=(((sig1(flag))^2+(sigz1(flag))^2-(2*nul*(sig1(flag))*(sigz1(flag))))/(2*el));
    x1(flag)=ypos;
    flag=flag+1;
           
       end 
se1u=se1(1:length(x1));     
       u1(main)=trapz(x1,se1u);
end

% if sc==laststep
% what3=(x1(1:sum(nel(1:2))+1)*10^6)';
% what4=(sig1((1:sum(nel(1:2))+1)))'/10^6;
% 
% 
%        
% end

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
            
                 ypos=position(index);
        th=alphal*(tfield(index)-ti(index));
           
    sig2(flag)=ebarl*((para(1)-(para(2)*ypos)-(cl*th)));
    if x==2
    sigz2(flag)=nul*sig2(flag)-el*(th);
    end
    se2(flag)=(((sig2(flag))^2+(sigz2(flag))^2-(2*nul*(sig2(flag))*(sigz2(flag))))/(2*el));
    x2(flag)=ypos;
    flag=flag+1;
          
end 

       
end
se2u=se2(1:length(x2));
%        figure
%        plot(x2,sig2/10^6)
u2(main)=trapz(x2,se2u);
   end 
    pos(main)=position(node);

end

toc

if last==eltot+1
    u2(main)=0;
end

temp=temp+(tfield-ti);
% figure
% plot(pos,u1)
% figure
% plot(pos,u2)

 for main=1:1:length(u1)
   setot(main)=u1(main)+u2(main);   
 end
  
for i=1:1:length(u1)
    serr(i)=setot(1)-setot(i);
end
if laststep~=1
serrint(bec+inc)=serr(2);
time(bec+inc)=t(bec)+sum(tdur(1:laststep-1));
else
serrint(bec+inc)=serr(2);
time(bec+inc)=t(bec);
end
tempsurf(bec+inc)=y(bec,1);
end
% 
%   figure
%   plot(what3,what4)
%    
% figure
% plot(pos,u1)
% figure
% plot(pos,u2)
% figure
% plot(pos,setot)
% figure
% plot(pos,serr)
% figure
% what5=(pos(1:nel(1)+nel(2)+1))'*10^6;
% what6=serr(1:nel(1)+nel(2)+1)';
% plot(what5,what6)
% figure
% plot(time,serrint)
% time=time';
% serrint=serrint';

inc=inc+bec;
   
% hold on
% plot(time,serrint)
time=time';
serrint=serrint';
end

% figure
% plot(pos,u1)
% figure
% plot(pos,u2)
% figure
% plot(pos,serr)

figure
plot(time,serrint)