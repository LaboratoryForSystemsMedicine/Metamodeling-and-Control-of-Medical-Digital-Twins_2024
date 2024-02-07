
Par=[2041.47268252353;5.00396479173130e-05;1.58271706006361e-06;2.89933858118167e-05;2.01974664064434e-05;0.294006511024592;-0.243605436640806;1.03399329930110;0.656475505289393;-0.153503134191659;0.0397709188045785;-0.0293737138451394;0.126726352542406;0.912919840877928;1.30595144065291;1.47935348102370;0.930658675595974;-0.0395399531877818;0.0748025259974428;-0.178759987291515;0.118259892028700;0.626936103552563;1.00043598523427;0.890719469533799]



lst_ODE_model(Par,1)
drawnow

% options = optimset('Display','iter','MaxFunEvals',3000,'MaxIter',3000);
% new_par=fminsearch(@lst_ODE_model,Par,options,0)
% lst_ODE_model(new_par,1)


return

function z=lst_ODE_model(Par,g)

Par(1:6)=abs(Par(1:6));


ok=0;lst1=0;lst2=0;lstl1=0;
% Par([11 12 ])=0;

load('Dados01')
    
Tspan=50:1000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6);     
[t,y] = ode45(@ODE_SWG2,Tspan,[Ya(Tspan(1),:)],options,Par);
if t(end)<Tspan(end)
    ok=1e12
end
Yss=y(end,:);
if g==1
    newcolors = [0.25 0.80 0.54;0.7 0.7 0.7;0.83 0.14 0.14];
         
    
    figure(1)
    colororder(newcolors)
    plot(t,y, 'LineWidth', 1.5)
    
    hold on
    plot(Ta,Ya,'--', 'LineWidth', 1.5)
    hold off
    legend('grass','sheep','wolves')

end

dadosS=interp1(t,y,Ta);
lst=(dadosS-Ya).^2;
lst(isnan(lst))=0;
lst1=sum(sum(lst));
%%
load('Dados02')

Tspan=1051:2000;
     
[t,y] = ode45(@ODE_SWG2,Tspan,[Ya(Tspan(1),:)],options,Par);

if t(end)<Tspan(end)
    ok=1e12
end

if g==1
    figure(2)
    colororder(newcolors)
    plot(Ta(1000:end),Ya(1000:end,:),'--', 'LineWidth', 1.5)
    hold on
    plot(t,y, 'LineWidth', 1.5)
    hold off
    legend('grass','sheep','wolves')
    xlim([1000 2000])

end
dadosS=interp1(t,y,Ta);
lst=(dadosS-Ya).^2;
lst(isnan(lst))=0;
lst2=sum(sum(lst));
%%

if g==1 disp([lst1 lst2]);end

z=lst1+lst2+ok;
end
function dxdt=ODE_SWG2(t,x,p) %This is the ODE function
K=255*255;

Pa=reshape(p,[6 4]);
Pr=nan(6,4);
Pr(:,1)=Pa(:,1);
Pr(:,2)=x(1).^Pa(:,2);
Pr(:,3)=x(2).^Pa(:,3);
Pr(:,4)=x(3).^Pa(:,4);

Pr=prod(Pr,2);

dxdt=[Pr(1)-Pr(2)
      Pr(3)-Pr(4)
      Pr(5)-Pr(6)];

end