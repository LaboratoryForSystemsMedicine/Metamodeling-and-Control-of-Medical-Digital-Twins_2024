
Par=[-0.157480962626959;0.00704751987606036;-0.000611887284221500;-0.598248145268082;8.24147957006338e-06;0.0148722971881940;-9.02894813734507e-06;-0.0704180185709875;7.76515493728275e-07;3.16465826229703e-06;1.47130639137662e-07;-4.40097827285390e-08;-7.09888481748566e-06;1.16900188311027e-06;2.63799816823151e-07;-3.59864109201160e-06;-1.82459502749679e-06;-1.22787874461126e-06;-8.36738732291953e-06;3.62968247715953e-06;3.72131159944533e-07;-1.30994768705092e-05;-2.42854373282788e-06;1.31232704683145e-06;1.02251531487329e-06;4.26066768464608e-06;-6.74103923509882e-07]


old_lst=lst_ODE_model(Par,1)
drawnow; pause(.1)


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

options = odeset('RelTol',1e-6,'AbsTol',1e-6);     
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
%%
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


if g==1 disp([lst1 lst2 lst7 ok]);end

z=lst1+lst2+lst7+ok;
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