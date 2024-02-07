
Par=[707.898797713976;2.74290709712722e-05;1.18449032179785e-06;7.63614372391576;2.47724294281499e-06;1.08954552559593e-05;0.129450032776724;-0.131322619148285;1.01002513278854;0.937497531445507;-0.940942480785047;-0.00313911396964099;0.0379856779311007;0.00216129790049945;0.0652059903326178;0.915235668006856;0.988680344035039;0.680196173029810;1.43529013066231;0.998693159011258;0.00973677896764006;0.138430231017208;-0.0624222334509402;0.0585716046376420;0.702043076249022;0.660600296309528;1.03520222002640;0.942070121122811]



lst_ODE_model(Par,1)
drawnow; pause(.1)


return

function z=lst_ODE_model(Par,g)
Par(1:7)=abs(Par(1:7));
Par(8:end)=Par(8:end).*(abs(Par(8:end))>0.03);
ok=0;lst1=0;lst2=0;


load('Dados01')
Tspan=50:1000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6);     
[t,y] = ode45(@ODE_SWG,Tspan,[Ya(Tspan(1),:) ],options,Par);
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
     
[t,y] = ode45(@ODE_SWG,Tspan,[Ya(Tspan(1),:)],options,Par);
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
lstl1=(reshape(Par,[7 4]));
lstl1=((lstl1(:,2:4))*1000).^2; 
lstl1(isnan(lstl1))=0;
lstl1=sum(sum(lstl1))*10;

if g==1 disp([lst1 lst2 lstl1 ok]);end

z=lst1+lst2+lstl1+ok;
end
function dxdt=ODE_SWG(t,x,p) %This is the ODE function

K=255*255;


Pn=reshape(1:28,[7 4]);
Pr(:,1)=p(Pn(:,1));
Pr(:,2:4)=repmat([x(1) x(2) x(3)],7,1).^p(Pn(:,2:4));
Pr=prod(Pr,2);

dxdt=[Pr(1)-Pr(2)
      Pr(3)-Pr(4)-Pr(5)
      Pr(6)-Pr(7)];

end
