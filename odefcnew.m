function dydt = odefcnew(t,y,bctype,M,K,F,nel,tau,tamb,tboun)

if bctype==1
    
    if t<=0.0129e-3
    temp1=tboun(1)+((2351-tboun(1))*t)/(0.0129e-3);
    temp2=(tboun(2));
    
    dtemp1=((2351-tboun(1)))/(0.0129e-3);
    dtemp2=0;

    elseif t>0.0129e-3 && t<=0.01933e-3
    temp1=2351+(((3280-2351))/(6.43e-6))*(t-0.0129e-3);
    temp2=(tboun(2));
    
    dtemp1=(((3280-2351))/(6.43e-6));
    dtemp2=0; 
    else
          temp1=3280+(((2700-3280))/(3.225e-6))*(t-0.01933e-3);
    temp2=(tboun(2));
    
    dtemp1=(((2700-3280))/(3.225e-6));
    dtemp2=0;   
    end
    
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
