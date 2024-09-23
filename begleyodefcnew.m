function dydt = begleyodefcnew(t,y,bctype,M,K,F,nel,tau,tamb,tboun)

if bctype==1
    
    temp1=tamb+(tboun(1)-tamb)*exp(-t/tau);
    temp2=tamb+(tboun(2)-tamb)*exp(-t/tau);
    
    dtemp1=(tboun(1)-tamb)*(-1/tau)*exp(-t/tau);
    dtemp2=(tboun(2)-tamb)*(-1/tau)*exp(-t/tau);

    
    M(1,:)=[];
    K(1,:)=[];
    F(1,:)=[];

    M(sum(nel),:)=[];
    K(sum(nel),:)=[];
    F(sum(nel),:)=[];

F=-(M(:,1)*dtemp1)-(K(:,1)*temp1)-(M(:,sum(nel)+1)*dtemp2)-K(:,sum(nel)+1)*temp2; 

M(:,1)=[];
K(:,1)=[];
M(:,sum(nel))=[];
K(:,sum(nel))=[];

end
dydt=(M)\(F-K*y);

end
