clear

dt=0.1;
rep=1;
TF=300;
NC=[1,1];
MU=[1,1,-1,1/60,0.5];
X0=[50,0,1,0,0];
X=zeros(5,round(TF/dt)+1,rep);


% for j=1:rep
%     X(:,1,j)=X0;
%     for i=1:round(TF/dt*1)-1
%         simulo l'aggiunta di un core
%         if(i==50)
%             X(4,i,j)=3;
%         end
%         Xi=delayQN(X(:,i,j),MU,NC,dt,1,dt);
%         X(:,i+1,j) = Xi(:,end,1);
%     end
% end

j=1;
X(:,1,j)=X0;
for i=1:round(TF/dt*1)-1

    if(i==50)
        X(4,i,j)=3;
    end

    [t,y,ssTR,ssRT]=delayQN_ODE(X(:,i,j),MU,NC,TF,dt);
    X(:,i+1,j) = y(end,:);
end

xtime=linspace(0,TF,round(TF/dt)+1);
plot(xtime,mean(X,3)')