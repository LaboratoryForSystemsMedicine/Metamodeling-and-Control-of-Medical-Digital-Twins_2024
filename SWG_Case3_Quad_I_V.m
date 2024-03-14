
Par=[-0.129035738251663;0.00632601829325315;0.00373433796409941;-0.597124203090022;0.00220652636250568;0.0488143549423052;0.0874831708463893;-0.0345319439969162;0.000330096930337369;8.09500415602465e-07;6.97706847500802e-07;1.02716995346397e-07;5.84723925196616e-05;-3.02585228817813e-06;2.14409470324196e-07;5.75787779880825e-05;1.39698681877772e-05;3.46802383792615e-07;6.50792419160402e-06;5.72736636683956e-06;7.81969326867586e-07;-0.000144635260567241;-3.56782396476998e-05;3.61976464804188e-05;-1.25552327204832e-05;-6.78848173698218e-06;3.08180235028640e-06]


old_lst=lst_ODE_model(Par,1)
drawnow 



% options = optimset('Display','iter','MaxFunEvals',5000,'MaxIter',5000);
% new_par=fminsearch(@lst_ODE_model,Par,options,0)
% new_lst=lst_ODE_model(new_par,1)


return

function z=lst_ODE_model(Par,g)
ok=0;lst1=0;lst2=0;


load('Dados01')
SS=mean(Ya(900:1000,:));

tic
Tspan=50:1000;

options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events',@odeEvents);     
[t,y] = ode45(@ODE_SWG,Tspan,[Ya(Tspan(1),:) ],options,Par,SS);
if t(end)<Tspan(end)
    ok=ok+1e12
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

load('Dados02')
Tspan=1051:2000;
tic
[t,y] = ode45(@ODE_SWG,Tspan,[Ya(Tspan(1),:)],options,Par,SS);
if t(end)<Tspan(end)
    ok=ok+1e12
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
if length(t)>10
    dadosS=interp1(t,y,Ta);
    lst=(dadosS-Ya).^2;
    lst(isnan(lst))=0;
    lst2=sum(sum(lst));
else
    lst2=1e12
end

%%
load('DadosConGrass2')
load('DadosConGrass2s')
Tspan=500:1500;
tic
options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events',@odeEvents2);         
[t,y] = ode45(@ODE_SWG_wC2,Tspan,Ya(Tspan(1),:),options,Par,SS,[0.02 0 0]);
if t(end)<Tspan(end)
    ok=ok+1e12
end

if g==1
    figure(4)
    colororder(newcolors)
    plot(t,y, 'LineWidth', 1.5)
    hold on
    plot(Ta(500:1500),Ya(500:1500,:),'--', 'LineWidth', 1.5)
    plot(Ta(500:1500),Ya(500:1500,:)+3*Yas(500:1500,:),'--', 'LineWidth', 1)
    plot(Ta(500:1500),Ya(500:1500,:)-3*Yas(500:1500,:),'--', 'LineWidth', 1)
    
    hold off
    legend('grass','sheep','wolves')
    xlim([500 1500])
end


dadosS=interp1(t,y,Ta);
lst=(dadosS-Ya).^2;
lst(isnan(lst))=0;
lst4=sum(sum(lst));

%%
load('DadosConSheep2')
load('DadosConSheep2s')
Tspan=500:1500;
     
[t,y] = ode45(@ODE_SWG_wC2,Tspan,Ya(Tspan(1),:),options,Par,SS,[0 0.02 0]);
if t(end)<Tspan(end)
    ok=ok+1e12
end
if sum(sum(y<0))>0
    ok=ok+sum(sum(y<0))*1e8;
end

if g==1
    figure(5)
    colororder(newcolors)
    plot(t,y, 'LineWidth', 1.5)
    hold on
    plot(Ta(500:1500),Ya(500:1500,:),'--', 'LineWidth', 1.5)
    plot(Ta(500:1500),Ya(500:1500,:)+3*Yas(500:1500,:),'--', 'LineWidth', 1)
    plot(Ta(500:1500),Ya(500:1500,:)-3*Yas(500:1500,:),'--', 'LineWidth', 1)
    
    hold off
    legend('grass','sheep','wolves')
    xlim([500 1500])
end


dadosS=interp1(t,y,Ta);
lst=(dadosS-Ya).^2;
lst(isnan(lst))=0;
lst5=sum(sum(lst));

%%
load('DadosConWolves1.5.mat')
load('DadosConWolves1.5s.mat')
Tspan=500:1500;
     
[t,y] = ode45(@ODE_SWG_wC2,Tspan,Ya(Tspan(1),:),options,Par,SS,[0 0 0.015]);
if t(end)<Tspan(end)
    ok=ok+1e12
end

if g==1
    figure(6)
    colororder(newcolors)
    plot(t,y, 'LineWidth', 1.5)
    hold on
    plot(Ta(500:1500),Ya(500:1500,:),'--', 'LineWidth', 1.5)
    plot(Ta(500:1500),Ya(500:1500,:)+3*Yas(500:1500,:),'--', 'LineWidth', 1)
    plot(Ta(500:1500),Ya(500:1500,:)-3*Yas(500:1500,:),'--', 'LineWidth', 1)
    
    hold off
    legend('grass','sheep','wolves')
    xlim([500 1500])
end



dadosS=interp1(t,y,Ta);
lst=(dadosS-Ya).^2;
lst(isnan(lst))=0;
lst6=sum(sum(lst));
%%
% Code to exclude other fixed points. The DE should only have the main SS as a
% fixed point, since other scondary fixed points may likely be unsteable.
% this way we are looking for the best fitting system that only has one
% fixed point given by the SS of the ABM.

J=reshape(Par(1:9),3,3);
H=reshape(Par(10:27),3,6);

syms x y z
eqns = [J(1,:)*[x;y;z]+H(1,:)*[x^2;y^2;z^2;x*y;y*z;x*z] == 0, J(2,:)*[x;y;z]+H(2,:)*[x^2;y^2;z^2;x*y;y*z;x*z] == 0, J(3,:)*[x;y;z]+H(3,:)*[x^2;y^2;z^2;x*y;y*z;x*z] == 0 ];
vars = [x y z];
eqnr = [x > 0-SS(1), x < 60000-SS(1), y > 0-SS(2), y < 10000-SS(2), z > 0-SS(3), z < 7500-SS(3) ];
[solx, soly, solz] = solve([eqns eqnr],vars);
sol=vpa([solx, soly, solz]);
sol2=eval(sol+SS);
limites=sol2(:,1)>0&sol2(:,2)>0&sol2(:,3)>0&sol2(:,1)<60000&sol2(:,2)<10000&sol2(:,3)<7500;
nzeros=length(limites);

lst7=(nzeros-1)*1e10;


if g==1 disp([lst1 lst2 lst4 lst5 lst6 lst7]);end

z=lst1+lst2+lst4+lst5+lst6+lst7+ok;
end
function dxdt=ODE_SWG(t,x,p,ss) %Hes
K=255*255;

Pr=reshape(p(1:9),3,3);
H=reshape(p(10:27),3,6);

X=(x-ss');
X2=[(x(1)-ss(1))*(x(1)-ss(1))
   (x(2)-ss(2))*(x(2)-ss(2))
   (x(3)-ss(3))*(x(3)-ss(3))
   (x(1)-ss(1))*(x(2)-ss(2))
   (x(2)-ss(2))*(x(3)-ss(3))
   (x(1)-ss(1))*(x(3)-ss(3))];

dxdt=Pr*X+H*X2;



end
function dxdt=ODE_SWG_wC2(t,x,p,ss,uu) %Hes

K=255*255;

u=[0 0 0];
if t>1000 u=uu; end

Pr=reshape(p(1:9),3,3);
H=reshape(p(10:27),3,6);

X=(x-ss');
X2=[(x(1)-ss(1))*(x(1)-ss(1))
    (x(2)-ss(2))*(x(2)-ss(2))
    (x(3)-ss(3))*(x(3)-ss(3))
    (x(1)-ss(1))*(x(2)-ss(2))
    (x(2)-ss(2))*(x(3)-ss(3))
    (x(1)-ss(1))*(x(3)-ss(3))];

dxdt=Pr*X+H*X2;

dxdt=dxdt-u'.*x;

end
function [value,isterminal,direction] = odeEvents(t,y,Par,SSa)
x=toc-10;
if x<0
    x=0;
end
value = x;     % Detect height = 0
isterminal = 1;   % Stop the integration
direction = 0;   
end
function [value,isterminal,direction] = odeEvents2(t,y,Par,SSa,uu)
x=toc-60;
if x<0
    x=0;
end
value = x;     % Detect height = 0
isterminal = 1;   % Stop the integration
direction = 0;   
end
