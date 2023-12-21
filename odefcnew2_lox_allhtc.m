function dydt = odefcnew2_lox_allhtc(t,y,sc,M,K,F,tdur)

templox=90;
temphox=770;
temphox2=420;
 hox=631765.730766605; %Heat transfer coefficient
% hfilm=130;
 hlox=41000;

 if sc==3 || sc==8 || sc==13 || sc==18

     temphoxe=temphox+((temphox2-temphox)/tdur(sc))*t;

    K(1,1)=K(1,1)+(hox);
    
    F(1,:)=F(1,:)+(hox*temphoxe);
    
        K(end,end)=K(end,end)+(hox);
    
    F(end,:)=F(end,:)+(hox*temphoxe); 


 elseif sc==1 || sc==4 || sc==6|| sc==9|| sc==11 || sc==14|| sc==16 || sc==19

                 K(1,1)=K(1,1)+(hlox);
    
    F(1,:)=F(1,:)+(hlox*templox);
    
        K(end,end)=K(end,end)+(hlox);
    
    F(end,:)=F(end,:)+(hlox*templox);

elseif sc==2 || sc==7 || sc==12 || sc==17
    
    K(1,1)=K(1,1)+(hox);
    
    F(1,:)=F(1,:)+(hox*temphox);
    
        K(end,end)=K(end,end)+(hox);
    
    F(end,:)=F(end,:)+(hox*temphox); 
    
end
dydt=(M)\(F-K*y);

end
