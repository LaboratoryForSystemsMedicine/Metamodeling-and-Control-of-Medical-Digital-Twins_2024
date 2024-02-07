
Par=[0.271534422167105;5.06439269612949e-06;3.54767246927118e-05;1.06349863925472e-06;0.0100168833738356;8.57919266464095e-06;1.79504478308065e-05;0.0751134696364733]


lst_ODE_model(Par,1)
drawnow; pause(.1)

% options = optimset('Display','iter','MaxFunEvals',3000,'MaxIter',3000);
% new_par=fminsearch(@lst_ODE_model,Par,options,0)
% lst_ODE_model(new_par,1)



return

function z=lst_ODE_model(Par,g)
Par=abs(Par);
ok=0;

load('Dados01')
Tspan=50:1000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6);     
[t,y] = ode45(@ODE_SWG,Tspan,[Ya(Tspan(1),:)],options,Par);
if t(end)<Tspan(end)
    ok=1e12;
end

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
     
[t,y] = ode45(@ODE_SWG,Tspan,[Ya(Tspan(1),:)],options,Par);
if t(end)<Tspan(end)
    ok=1e12;
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
    
%     y(end,:)
end
dadosS=interp1(t,y,Ta);
lst=(dadosS-Ya).^2;
lst(isnan(lst))=0;
lst2=sum(sum(lst));




z=lst1+lst2+ok;
end
function dxdt=ODE_SWG(t,x,p) %This is the ODE function

  
dxdt=[p(1)*x(1)-p(2)*x(1)^2-p(3)*x(1)*x(2)
      p(4)*x(1)*x(2)-p(5)*x(2)-p(6)*x(2)*x(3)
      p(7)*x(2)*x(3)-p(8)*x(3)];

end
