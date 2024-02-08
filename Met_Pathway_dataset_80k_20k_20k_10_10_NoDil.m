%function z=Met_Pathway_matlab()
%v2 implementing the changes from the ODE:
% ParA=[1327/500; 7547; 22350; 1588];
% ParE=[1340/500*1.5; 7794/2; 20830];
% ParI=[1317/500; 17630/2; 21350];
% ParO=[1295/500; 7053; 4459/500; 615; 8388/4; 24700; 98920];
clear

%% External parameters
par.SubS0=80000;
par.SubP0=20000;
par.SubQ0=20000;
par.SubR0=10;
par.SubT0=10;

par.EnzA0=200;
par.EnzE0=200;
par.EnzI0=200;
par.EnzO0=200;

par.ComplAS0=0;
par.ComplAP0=0;
par.ComplAR0=0;

par.ComplEP0=0;
par.ComplEQ0=0;

par.ComplIQ0=0;
par.ComplIR0=0;

par.ComplOP0=0;
par.ComplOT0=0;
par.ComplOR0=00;
par.ComplORP0=0;
par.ComplORT0=0;

par.worldx=100;
par.worldy=100;
par.maxmol=10000+sum([par.SubS0 par.SubP0 par.SubQ0 par.SubR0 par.SubT0])+100;
par.maxenz=300;

par.maxtime=100000;
par.graf=0;
par.addnum=1;  %5
par.molmovement=1;
par.enzmovement=0.1;

par.bindingAS=4;
par.bindingAP=2;
par.releaseAS=2;
par.releaseAP=4;
par.catAS=2;
par.bindingAR=8;
par.releaseAR=1;

% ParE=[1340/500*1.5; 7794/2; 20830];
par.bindingEP=4*2;
par.bindingEQ=2;
par.releaseEP=2;
par.releaseEQ=4;
par.catEP=2*1.5;

% ParI=[1317/500; 17630/2; 21350];
par.bindingIQ=1*2;
par.bindingIR=2;
par.releaseIQ=0.5;
par.releaseIR=4;
par.catIQ=2;

% ParO=[1295/500; 7053; 4459/500; 615; 8388/4; 24700; 98920];
par.bindingOP=4;
par.bindingOT=2;
par.releaseOP=2;
par.releaseOT=4;
par.catOP=2;

par.bindingORP=8*4;
par.bindingORT=1;
par.releaseORP=1;
par.releaseORT=8;
par.catORP=10;                          %Aqui

par.bindingOR=16;
par.releaseOR=1;

%internal parameters
mapacores=[1 1 1 ;hsv(26)];


tic1=tic;
% Noutput=[];
for k=1:par.addnum
    
    
    %% Setup
    
    
    %% create-SubS
    %                         1    2 3
    SubS=nan(par.maxmol,4); %alive, x,y
    if par.SubS0>0
        SubS(1:par.SubS0,1)=1;
        X = lhsdesign(par.SubS0,2);
        SubS(1:par.SubS0,2)=X(1:par.SubS0,1)*par.worldx+1;
        SubS(1:par.SubS0,3)=X(1:par.SubS0,2)*par.worldy+1;  %
        SubS(:,4)=floor(SubS(:,2)-1)+1+(floor(SubS(:,3)-1))*par.worldx;
    end
    
    %% create-SubP
    %                         1    2 3
    SubP=nan(par.maxmol,4); %alive, x,y
    if par.SubP0>0
        SubP(1:par.SubP0,1)=1;
        X = lhsdesign(par.SubP0,2);
        SubP(1:par.SubP0,2)=X(1:par.SubP0,1)*par.worldx+1;
        SubP(1:par.SubP0,3)=X(1:par.SubP0,2)*par.worldy+1;  %
        SubP(:,4)=floor(SubP(:,2)-1)+1+(floor(SubP(:,3)-1))*par.worldx;
    end
    
    %% create-SubQ
    %                         1    2 3
    SubQ=nan(par.maxmol,4); %alive, x,y
    if par.SubQ0>0
        SubQ(1:par.SubQ0,1)=1;
        X = lhsdesign(par.SubQ0,2);
        SubQ(1:par.SubQ0,2)=X(1:par.SubQ0,1)*par.worldx+1;
        SubQ(1:par.SubQ0,3)=X(1:par.SubQ0,2)*par.worldy+1;  %
        SubQ(:,4)=floor(SubQ(:,2)-1)+1+(floor(SubQ(:,3)-1))*par.worldx;
    end
    
    %% create-SubT
    %                         1    2 3
    SubT=nan(par.maxmol,4); %alive, x,y
    if par.SubT0>0
        SubT(1:par.SubT0,1)=1;
        X = lhsdesign(par.SubT0,2);
        SubT(1:par.SubT0,2)=X(1:par.SubT0,1)*par.worldx+1;
        SubT(1:par.SubT0,3)=X(1:par.SubT0,2)*par.worldy+1;  %
        SubT(:,4)=floor(SubT(:,2)-1)+1+(floor(SubT(:,3)-1))*par.worldx;
    end
    
    %% create-SubR
    %                         1    2 3
    SubR=nan(par.maxmol,4); %alive, x,y
    if par.SubR0>0
        SubR(1:par.SubR0,1)=1;
        X = lhsdesign(par.SubR0,2);
        SubR(1:par.SubR0,2)=X(1:par.SubR0,1)*par.worldx+1;
        SubR(1:par.SubR0,3)=X(1:par.SubR0,2)*par.worldy+1;  %
        SubR(:,4)=floor(SubR(:,2)-1)+1+(floor(SubR(:,3)-1))*par.worldx;
    end
    
    %% create-EnzA
    EnzA=nan(par.maxenz,4); %alive, x,y
    if par.EnzA0>0
        EnzA(1:par.EnzA0,1)=1;
        X = lhsdesign(par.EnzA0,2);
        EnzA(1:par.EnzA0,2)=X(1:par.EnzA0,1)*par.worldx+1;
        EnzA(1:par.EnzA0,3)=X(1:par.EnzA0,2)*par.worldy+1;
        EnzA(:,4)=floor(EnzA(:,2)-1)+1+(floor(EnzA(:,3)-1))*par.worldx;
    end
    
    %% create-EnzE
    EnzE=nan(par.maxenz,4); %alive, x,y
    if par.EnzE0>0
        EnzE(1:par.EnzE0,1)=1;
        X = lhsdesign(par.EnzE0,2);
        EnzE(1:par.EnzE0,2)=X(1:par.EnzE0,1)*par.worldx+1;
        EnzE(1:par.EnzE0,3)=X(1:par.EnzE0,2)*par.worldy+1;
        EnzE(:,4)=floor(EnzE(:,2)-1)+1+(floor(EnzE(:,3)-1))*par.worldx;
    end
    
    %% create-EnzI
    EnzI=nan(par.maxenz,4); %alive, x,y
    if par.EnzI0>0
        EnzI(1:par.EnzI0,1)=1;
        X = lhsdesign(par.EnzI0,2);
        EnzI(1:par.EnzI0,2)=X(1:par.EnzI0,1)*par.worldx+1;
        EnzI(1:par.EnzI0,3)=X(1:par.EnzI0,2)*par.worldy+1;
        EnzI(:,4)=floor(EnzI(:,2)-1)+1+(floor(EnzI(:,3)-1))*par.worldx;
    end
    
    %% create-EnzO
    EnzO=nan(par.maxenz,4); %alive, x,y
    if par.EnzO0>0
        EnzO(1:par.EnzO0,1)=1;
        X = lhsdesign(par.EnzO0,2);
        EnzO(1:par.EnzO0,2)=X(1:par.EnzO0,1)*par.worldx+1;
        EnzO(1:par.EnzO0,3)=X(1:par.EnzO0,2)*par.worldy+1;
        EnzO(:,4)=floor(EnzO(:,2)-1)+1+(floor(EnzO(:,3)-1))*par.worldx;
    end
    
    %% create-ComplOR
    ComplOR=nan(par.maxenz,4); %alive, x,y
    if par.ComplOR0>0
        ComplOR(1:par.ComplOR0,1)=1;
        X = lhsdesign(par.ComplOR0,2);
        ComplOR(1:par.ComplOR0,2)=X(1:par.ComplOR0,1)*par.worldx+1;
        ComplOR(1:par.ComplOR0,3)=X(1:par.ComplOR0,2)*par.worldy+1;
        ComplOR(:,4)=floor(ComplOR(:,2)-1)+1+(floor(ComplOR(:,3)-1))*par.worldx;
    end
    
    %% create empties
    %                         1    2 3
    
    ComplAS=nan(par.maxenz,4); %alive, x,y
    ComplAP=nan(par.maxenz,4); %alive, x,y
    ComplAR=nan(par.maxenz,4); %alive, x,y
    
    ComplEP=nan(par.maxenz,4); %alive, x,y
    ComplEQ=nan(par.maxenz,4); %alive, x,y
    
    ComplIQ=nan(par.maxenz,4); %alive, x,y
    ComplIR=nan(par.maxenz,4); %alive, x,y
    
    ComplOP=nan(par.maxenz,4); %alive, x,y
    ComplOT=nan(par.maxenz,4); %alive, x,y
    
    ComplORP=nan(par.maxenz,4); %alive, x,y
    ComplORT=nan(par.maxenz,4); %alive, x,y
    
    NonInt=0;
    
    time=0;
    if isfinite(par.maxtime)
        totaloutput=nan(par.maxtime,22);
        output=nan(par.maxtime,6);
    else
        totaloutput=[];
        output=[];
    end
    tic2=tic;
    while time<par.maxtime
        %         time
        if isfinite(par.maxtime)
            totaloutput(time+1,:)=[time, sum(SubS(:,1)>0),sum(ComplAS(:,1)>0), sum(EnzA(:,1)>0),sum(ComplAR(:,1)>0),sum(ComplAP(:,1)>0),sum(SubP(:,1)>0),sum(ComplEP(:,1)>0),sum(EnzE(:,1)>0),sum(ComplEQ(:,1)>0),sum(SubQ(:,1)>0),sum(ComplIQ(:,1)>0),sum(EnzI(:,1)>0),sum(ComplIR(:,1)>0),sum(SubR(:,1)>0),sum(ComplOP(:,1)>0),sum(EnzO(:,1)>0),sum(ComplOT(:,1)>0),sum(SubT(:,1)>0),sum(ComplORP(:,1)>0),sum(ComplOR(:,1)>0),sum(ComplORT(:,1)>0)];
            output(time+1,:)=[time,sum(SubS(:,1)>0),sum(SubP(:,1)>0),sum(SubQ(:,1)>0),sum(SubT(:,1)>0),sum(SubR(:,1)>0) ];
            %             output(time+1,:)=[time, sum(SubS(:,1)>0),sum(SubP(:,1)>0), sum(SubQ(:,1)>0),0,0];
            %             output(time+1,:)=[time, sum(SubS(:,1)>0),sum(ComplAS(:,1)>0), sum(EnzA(:,1)>0),sum(ComplAR(:,1)>0),sum(ComplAP(:,1)>0),sum(SubP(:,1)>0),sum(ComplEP(:,1)>0),sum(EnzE(:,1)>0),sum(ComplEQ(:,1)>0),sum(SubQ(:,1)>0),sum(ComplIQ(:,1)>0),sum(EnzI(:,1)>0),sum(ComplIR(:,1)>0),sum(SubR(:,1)>0),sum(ComplOP(:,1)>0),sum(EnzO(:,1)>0),sum(ComplOT(:,1)>0),sum(SubT(:,1)>0),sum(ComplORP(:,1)>0),sum(ComplOR(:,1)>0),sum(ComplORT(:,1)>0)];
            %             output(time+1,:)=[time, sum(SubS(:,1)>0),sum(SubP(:,1)>0),sum(SubQ(:,1)>0),sum(SubT(:,1)>0),sum(SubR(:,1)>0)];
        else
            totaloutput=[totaloutput; time, sum(SubS(:,1)>0),sum(ComplAS(:,1)>0), sum(EnzA(:,1)>0),sum(ComplAR(:,1)>0),sum(ComplAP(:,1)>0),sum(SubP(:,1)>0),sum(ComplEP(:,1)>0),sum(EnzE(:,1)>0),sum(ComplEQ(:,1)>0),sum(SubQ(:,1)>0),sum(ComplIQ(:,1)>0),sum(EnzI(:,1)>0),sum(ComplIR(:,1)>0),sum(SubR(:,1)>0),sum(ComplOP(:,1)>0),sum(EnzO(:,1)>0),sum(ComplOT(:,1)>0),sum(SubT(:,1)>0),sum(ComplORP(:,1)>0),sum(ComplOR(:,1)>0),sum(ComplORT(:,1)>0)];
            output=[output; time, sum(SubS(:,1)>0),sum(ComplAS(:,1)>0), sum(EnzA(:,1)>0),sum(ComplAP(:,1)>0),sum(SubP(:,1)>0),sum(ComplEP(:,1)>0), sum(EnzE(:,1)>0),sum(ComplEQ(:,1)>0),sum(SubQ(:,1)>0)];
        end
        
        
        %% Plotting
        if par.graf==1
            figure(1)
            
            
            mapa=zeros(par.worldx*par.worldy,1);
            mapa(EnzA(EnzA(:,1)>0,4))=2;
            mapa(EnzE(EnzE(:,1)>0,4))=3;
            mapa(EnzI(EnzI(:,1)>0,4))=4;
            mapa(EnzO(EnzO(:,1)>0,4))=5;
            
            mapa(SubS(SubS(:,1)>0,4))=6;
            mapa(SubP(SubP(:,1)>0,4))=7;
            mapa(SubQ(SubQ(:,1)>0,4))=8;
            mapa(SubR(SubR(:,1)>0,4))=9;
            mapa(SubT(SubT(:,1)>0,4))=10;
            
            mapa(ComplAS(ComplAS(:,1)>0,4))=11;
            mapa(ComplAP(ComplAP(:,1)>0,4))=12;
            mapa(ComplAR(ComplAR(:,1)>0,4))=13;
            
            mapa(ComplEP(ComplEP(:,1)>0,4))=14;
            mapa(ComplEQ(ComplEQ(:,1)>0,4))=15;
            mapa(ComplIQ(ComplIQ(:,1)>0,4))=16;
            mapa(ComplIR(ComplIR(:,1)>0,4))=17;
            
            mapa(ComplOP(ComplOP(:,1)>0,4))=18;
            mapa(ComplOT(ComplOT(:,1)>0,4))=19;
            mapa(ComplOR(ComplOR(:,1)>0,4))=20;
            mapa(ComplORP(ComplORP(:,1)>0,4))=21;
            mapa(ComplORT(ComplORT(:,1)>0,4))=22;
            if exist('KillSub','var')
                mapa(KillSub)=25;
            end
            
            mapa=reshape(mapa,par.worldx,par.worldy);
            mapa(par.worldx,par.worldy)=25;
            
            mapa=transpose(mapa);
            
            %             heatmap(mapa)
            image(mapa)
            colormap(mapacores)
            title(num2str(time))
            %             drawnow
            
            
        end
        
        %% Move
        %SubS
        if sum(SubS(:,1)>0)>0
            % move
            X=rand(length(SubS),1)*2*pi;
            SubS(:,2)=SubS(:,2)+par.molmovement*cos(X);
            SubS(:,3)=SubS(:,3)-par.molmovement*sin(X);
            
            SubS((SubS(:,2)<1),2)=SubS((SubS(:,2)<1),2)+par.worldx;
            SubS((SubS(:,2)>par.worldx+1),2)=SubS((SubS(:,2)>par.worldx+1),2)-par.worldx;
            
            SubS((SubS(:,3)<1),3)=SubS((SubS(:,3)<1),3)+par.worldy;
            SubS((SubS(:,3)>par.worldy+1),3)=SubS((SubS(:,3)>par.worldy+1),3)-par.worldy;
            
            SubS(:,4)=floor(SubS(:,2)-1)+1+(floor(SubS(:,3)-1))*par.worldx;
        end
        %SubP
        if sum(SubP(:,1)>0)>0
            % move
            X=rand(length(SubP),1)*2*pi;
            SubP(:,2)=SubP(:,2)+par.molmovement*cos(X);
            SubP(:,3)=SubP(:,3)-par.molmovement*sin(X);
            
            SubP((SubP(:,2)<1),2)=SubP((SubP(:,2)<1),2)+par.worldx;
            SubP((SubP(:,2)>par.worldx+1),2)=SubP((SubP(:,2)>par.worldx+1),2)-par.worldx;
            
            SubP((SubP(:,3)<1),3)=SubP((SubP(:,3)<1),3)+par.worldy;
            SubP((SubP(:,3)>par.worldy+1),3)=SubP((SubP(:,3)>par.worldy+1),3)-par.worldy;
            
            SubP(:,4)=floor(SubP(:,2)-1)+1+(floor(SubP(:,3)-1))*par.worldx;
        end
        %SubQ
        if sum(SubQ(:,1)>0)>0
            % move
            X=rand(length(SubQ),1)*2*pi;
            SubQ(:,2)=SubQ(:,2)+par.molmovement*cos(X);
            SubQ(:,3)=SubQ(:,3)-par.molmovement*sin(X);
            
            SubQ((SubQ(:,2)<1),2)=SubQ((SubQ(:,2)<1),2)+par.worldx;
            SubQ((SubQ(:,2)>par.worldx+1),2)=SubQ((SubQ(:,2)>par.worldx+1),2)-par.worldx;
            
            SubQ((SubQ(:,3)<1),3)=SubQ((SubQ(:,3)<1),3)+par.worldy;
            SubQ((SubQ(:,3)>par.worldy+1),3)=SubQ((SubQ(:,3)>par.worldy+1),3)-par.worldy;
            
            SubQ(:,4)=floor(SubQ(:,2)-1)+1+(floor(SubQ(:,3)-1))*par.worldx;
        end
        %SubR
        if sum(SubR(:,1)>0)>0
            % move
            X=rand(length(SubR),1)*2*pi;
            SubR(:,2)=SubR(:,2)+par.molmovement*cos(X);
            SubR(:,3)=SubR(:,3)-par.molmovement*sin(X);
            
            SubR((SubR(:,2)<1),2)=SubR((SubR(:,2)<1),2)+par.worldx;
            SubR((SubR(:,2)>par.worldx+1),2)=SubR((SubR(:,2)>par.worldx+1),2)-par.worldx;
            
            SubR((SubR(:,3)<1),3)=SubR((SubR(:,3)<1),3)+par.worldy;
            SubR((SubR(:,3)>par.worldy+1),3)=SubR((SubR(:,3)>par.worldy+1),3)-par.worldy;
            
            SubR(:,4)=floor(SubR(:,2)-1)+1+(floor(SubR(:,3)-1))*par.worldx;
        end
        %SubT
        if sum(SubT(:,1)>0)>0
            % move
            X=rand(length(SubT),1)*2*pi;
            SubT(:,2)=SubT(:,2)+par.molmovement*cos(X);
            SubT(:,3)=SubT(:,3)-par.molmovement*sin(X);
            
            SubT((SubT(:,2)<1),2)=SubT((SubT(:,2)<1),2)+par.worldx;
            SubT((SubT(:,2)>par.worldx+1),2)=SubT((SubT(:,2)>par.worldx+1),2)-par.worldx;
            
            SubT((SubT(:,3)<1),3)=SubT((SubT(:,3)<1),3)+par.worldy;
            SubT((SubT(:,3)>par.worldy+1),3)=SubT((SubT(:,3)>par.worldy+1),3)-par.worldy;
            
            SubT(:,4)=floor(SubT(:,2)-1)+1+(floor(SubT(:,3)-1))*par.worldx;
        end
        
        %EnzA
        if sum(EnzA(:,1)>0)>0
            % move
            X=rand(length(EnzA),1)*2*pi;
            EnzA(:,2)=EnzA(:,2)+par.enzmovement*cos(X);
            EnzA(:,3)=EnzA(:,3)-par.enzmovement*sin(X);
            
            EnzA((EnzA(:,2)<1),2)=EnzA((EnzA(:,2)<1),2)+par.worldx;
            EnzA((EnzA(:,2)>par.worldx+1),2)=EnzA((EnzA(:,2)>par.worldx+1),2)-par.worldx;
            
            EnzA((EnzA(:,3)<1),3)=EnzA((EnzA(:,3)<1),3)+par.worldy;
            EnzA((EnzA(:,3)>par.worldy+1),3)=EnzA((EnzA(:,3)>par.worldy+1),3)-par.worldy;
            
            EnzA(:,4)=floor(EnzA(:,2)-1)+1+(floor(EnzA(:,3)-1))*par.worldx;
        end
        %EnzE
        if sum(EnzE(:,1)>0)>0
            % move
            X=rand(length(EnzE),1)*2*pi;
            EnzE(:,2)=EnzE(:,2)+par.enzmovement*cos(X);
            EnzE(:,3)=EnzE(:,3)-par.enzmovement*sin(X);
            
            EnzE((EnzE(:,2)<1),2)=EnzE((EnzE(:,2)<1),2)+par.worldx;
            EnzE((EnzE(:,2)>par.worldx+1),2)=EnzE((EnzE(:,2)>par.worldx+1),2)-par.worldx;
            
            EnzE((EnzE(:,3)<1),3)=EnzE((EnzE(:,3)<1),3)+par.worldy;
            EnzE((EnzE(:,3)>par.worldy+1),3)=EnzE((EnzE(:,3)>par.worldy+1),3)-par.worldy;
            
            EnzE(:,4)=floor(EnzE(:,2)-1)+1+(floor(EnzE(:,3)-1))*par.worldx;
        end
        %EnzI
        if sum(EnzI(:,1)>0)>0
            % move
            X=rand(length(EnzI),1)*2*pi;
            EnzI(:,2)=EnzI(:,2)+par.enzmovement*cos(X);
            EnzI(:,3)=EnzI(:,3)-par.enzmovement*sin(X);
            
            EnzI((EnzI(:,2)<1),2)=EnzI((EnzI(:,2)<1),2)+par.worldx;
            EnzI((EnzI(:,2)>par.worldx+1),2)=EnzI((EnzI(:,2)>par.worldx+1),2)-par.worldx;
            
            EnzI((EnzI(:,3)<1),3)=EnzI((EnzI(:,3)<1),3)+par.worldy;
            EnzI((EnzI(:,3)>par.worldy+1),3)=EnzI((EnzI(:,3)>par.worldy+1),3)-par.worldy;
            
            EnzI(:,4)=floor(EnzI(:,2)-1)+1+(floor(EnzI(:,3)-1))*par.worldx;
        end
        %EnzO
        if sum(EnzO(:,1)>0)>0
            % move
            X=rand(length(EnzO),1)*2*pi;
            EnzO(:,2)=EnzO(:,2)+par.enzmovement*cos(X);
            EnzO(:,3)=EnzO(:,3)-par.enzmovement*sin(X);
            
            EnzO((EnzO(:,2)<1),2)=EnzO((EnzO(:,2)<1),2)+par.worldx;
            EnzO((EnzO(:,2)>par.worldx+1),2)=EnzO((EnzO(:,2)>par.worldx+1),2)-par.worldx;
            
            EnzO((EnzO(:,3)<1),3)=EnzO((EnzO(:,3)<1),3)+par.worldy;
            EnzO((EnzO(:,3)>par.worldy+1),3)=EnzO((EnzO(:,3)>par.worldy+1),3)-par.worldy;
            
            EnzO(:,4)=floor(EnzO(:,2)-1)+1+(floor(EnzO(:,3)-1))*par.worldx;
        end
        
        %ComplAS
        if sum(ComplAS(:,1)>0)>0
            % move
            X=rand(length(ComplAS),1)*2*pi;
            ComplAS(:,2)=ComplAS(:,2)+par.enzmovement*cos(X);
            ComplAS(:,3)=ComplAS(:,3)-par.enzmovement*sin(X);
            
            ComplAS((ComplAS(:,2)<1),2)=ComplAS((ComplAS(:,2)<1),2)+par.worldx;
            ComplAS((ComplAS(:,2)>par.worldx+1),2)=ComplAS((ComplAS(:,2)>par.worldx+1),2)-par.worldx;
            
            ComplAS((ComplAS(:,3)<1),3)=ComplAS((ComplAS(:,3)<1),3)+par.worldy;
            ComplAS((ComplAS(:,3)>par.worldy+1),3)=ComplAS((ComplAS(:,3)>par.worldy+1),3)-par.worldy;
            
            ComplAS(:,4)=floor(ComplAS(:,2)-1)+1+(floor(ComplAS(:,3)-1))*par.worldx;
        end
        %ComplAP
        if sum(ComplAP(:,1)>0)>0
            % move
            X=rand(length(ComplAP),1)*2*pi;
            ComplAP(:,2)=ComplAP(:,2)+par.enzmovement*cos(X);
            ComplAP(:,3)=ComplAP(:,3)-par.enzmovement*sin(X);
            
            ComplAP((ComplAP(:,2)<1),2)=ComplAP((ComplAP(:,2)<1),2)+par.worldx;
            ComplAP((ComplAP(:,2)>par.worldx+1),2)=ComplAP((ComplAP(:,2)>par.worldx+1),2)-par.worldx;
            
            ComplAP((ComplAP(:,3)<1),3)=ComplAP((ComplAP(:,3)<1),3)+par.worldy;
            ComplAP((ComplAP(:,3)>par.worldy+1),3)=ComplAP((ComplAP(:,3)>par.worldy+1),3)-par.worldy;
            
            ComplAP(:,4)=floor(ComplAP(:,2)-1)+1+(floor(ComplAP(:,3)-1))*par.worldx;
        end
        %ComplAR
        if sum(ComplAR(:,1)>0)>0
            % move
            X=rand(length(ComplAR),1)*2*pi;
            ComplAR(:,2)=ComplAR(:,2)+par.enzmovement*cos(X);
            ComplAR(:,3)=ComplAR(:,3)-par.enzmovement*sin(X);
            
            ComplAR((ComplAR(:,2)<1),2)=ComplAR((ComplAR(:,2)<1),2)+par.worldx;
            ComplAR((ComplAR(:,2)>par.worldx+1),2)=ComplAR((ComplAR(:,2)>par.worldx+1),2)-par.worldx;
            
            ComplAR((ComplAR(:,3)<1),3)=ComplAR((ComplAR(:,3)<1),3)+par.worldy;
            ComplAR((ComplAR(:,3)>par.worldy+1),3)=ComplAR((ComplAR(:,3)>par.worldy+1),3)-par.worldy;
            
            ComplAR(:,4)=floor(ComplAR(:,2)-1)+1+(floor(ComplAR(:,3)-1))*par.worldx;
        end
        
        %ComplEQ
        if sum(ComplEQ(:,1)>0)>0
            % move
            X=rand(length(ComplEQ),1)*2*pi;
            ComplEQ(:,2)=ComplEQ(:,2)+par.enzmovement*cos(X);
            ComplEQ(:,3)=ComplEQ(:,3)-par.enzmovement*sin(X);
            
            ComplEQ((ComplEQ(:,2)<1),2)=ComplEQ((ComplEQ(:,2)<1),2)+par.worldx;
            ComplEQ((ComplEQ(:,2)>par.worldx+1),2)=ComplEQ((ComplEQ(:,2)>par.worldx+1),2)-par.worldx;
            
            ComplEQ((ComplEQ(:,3)<1),3)=ComplEQ((ComplEQ(:,3)<1),3)+par.worldy;
            ComplEQ((ComplEQ(:,3)>par.worldy+1),3)=ComplEQ((ComplEQ(:,3)>par.worldy+1),3)-par.worldy;
            
            ComplEQ(:,4)=floor(ComplEQ(:,2)-1)+1+(floor(ComplEQ(:,3)-1))*par.worldx;
        end
        %ComplEP
        if sum(ComplEP(:,1)>0)>0
            % move
            X=rand(length(ComplEP),1)*2*pi;
            ComplEP(:,2)=ComplEP(:,2)+par.enzmovement*cos(X);
            ComplEP(:,3)=ComplEP(:,3)-par.enzmovement*sin(X);
            
            ComplEP((ComplEP(:,2)<1),2)=ComplEP((ComplEP(:,2)<1),2)+par.worldx;
            ComplEP((ComplEP(:,2)>par.worldx+1),2)=ComplEP((ComplEP(:,2)>par.worldx+1),2)-par.worldx;
            
            ComplEP((ComplEP(:,3)<1),3)=ComplEP((ComplEP(:,3)<1),3)+par.worldy;
            ComplEP((ComplEP(:,3)>par.worldy+1),3)=ComplEP((ComplEP(:,3)>par.worldy+1),3)-par.worldy;
            
            ComplEP(:,4)=floor(ComplEP(:,2)-1)+1+(floor(ComplEP(:,3)-1))*par.worldx;
        end
        %ComplIQ
        if sum(ComplIQ(:,1)>0)>0
            % move
            X=rand(length(ComplIQ),1)*2*pi;
            ComplIQ(:,2)=ComplIQ(:,2)+par.enzmovement*cos(X);
            ComplIQ(:,3)=ComplIQ(:,3)-par.enzmovement*sin(X);
            
            ComplIQ((ComplIQ(:,2)<1),2)=ComplIQ((ComplIQ(:,2)<1),2)+par.worldx;
            ComplIQ((ComplIQ(:,2)>par.worldx+1),2)=ComplIQ((ComplIQ(:,2)>par.worldx+1),2)-par.worldx;
            
            ComplIQ((ComplIQ(:,3)<1),3)=ComplIQ((ComplIQ(:,3)<1),3)+par.worldy;
            ComplIQ((ComplIQ(:,3)>par.worldy+1),3)=ComplIQ((ComplIQ(:,3)>par.worldy+1),3)-par.worldy;
            
            ComplIQ(:,4)=floor(ComplIQ(:,2)-1)+1+(floor(ComplIQ(:,3)-1))*par.worldx;
        end
        %ComplIR
        if sum(ComplIR(:,1)>0)>0
            % move
            X=rand(length(ComplIR),1)*2*pi;
            ComplIR(:,2)=ComplIR(:,2)+par.enzmovement*cos(X);
            ComplIR(:,3)=ComplIR(:,3)-par.enzmovement*sin(X);
            
            ComplIR((ComplIR(:,2)<1),2)=ComplIR((ComplIR(:,2)<1),2)+par.worldx;
            ComplIR((ComplIR(:,2)>par.worldx+1),2)=ComplIR((ComplIR(:,2)>par.worldx+1),2)-par.worldx;
            
            ComplIR((ComplIR(:,3)<1),3)=ComplIR((ComplIR(:,3)<1),3)+par.worldy;
            ComplIR((ComplIR(:,3)>par.worldy+1),3)=ComplIR((ComplIR(:,3)>par.worldy+1),3)-par.worldy;
            
            ComplIR(:,4)=floor(ComplIR(:,2)-1)+1+(floor(ComplIR(:,3)-1))*par.worldx;
        end
        
        %ComplOP
        if sum(ComplOP(:,1)>0)>0
            % move
            X=rand(length(ComplOP),1)*2*pi;
            ComplOP(:,2)=ComplOP(:,2)+par.enzmovement*cos(X);
            ComplOP(:,3)=ComplOP(:,3)-par.enzmovement*sin(X);
            
            ComplOP((ComplOP(:,2)<1),2)=ComplOP((ComplOP(:,2)<1),2)+par.worldx;
            ComplOP((ComplOP(:,2)>par.worldx+1),2)=ComplOP((ComplOP(:,2)>par.worldx+1),2)-par.worldx;
            
            ComplOP((ComplOP(:,3)<1),3)=ComplOP((ComplOP(:,3)<1),3)+par.worldy;
            ComplOP((ComplOP(:,3)>par.worldy+1),3)=ComplOP((ComplOP(:,3)>par.worldy+1),3)-par.worldy;
            
            ComplOP(:,4)=floor(ComplOP(:,2)-1)+1+(floor(ComplOP(:,3)-1))*par.worldx;
        end
        %ComplOT
        if sum(ComplOT(:,1)>0)>0
            % move
            X=rand(length(ComplOT),1)*2*pi;
            ComplOT(:,2)=ComplOT(:,2)+par.enzmovement*cos(X);
            ComplOT(:,3)=ComplOT(:,3)-par.enzmovement*sin(X);
            
            ComplOT((ComplOT(:,2)<1),2)=ComplOT((ComplOT(:,2)<1),2)+par.worldx;
            ComplOT((ComplOT(:,2)>par.worldx+1),2)=ComplOT((ComplOT(:,2)>par.worldx+1),2)-par.worldx;
            
            ComplOT((ComplOT(:,3)<1),3)=ComplOT((ComplOT(:,3)<1),3)+par.worldy;
            ComplOT((ComplOT(:,3)>par.worldy+1),3)=ComplOT((ComplOT(:,3)>par.worldy+1),3)-par.worldy;
            
            ComplOT(:,4)=floor(ComplOT(:,2)-1)+1+(floor(ComplOT(:,3)-1))*par.worldx;
        end
        %ComplOR
        if sum(ComplOR(:,1)>0)>0
            % move
            X=rand(length(ComplOR),1)*2*pi;
            ComplOR(:,2)=ComplOR(:,2)+par.enzmovement*cos(X);
            ComplOR(:,3)=ComplOR(:,3)-par.enzmovement*sin(X);
            
            ComplOR((ComplOR(:,2)<1),2)=ComplOR((ComplOR(:,2)<1),2)+par.worldx;
            ComplOR((ComplOR(:,2)>par.worldx+1),2)=ComplOR((ComplOR(:,2)>par.worldx+1),2)-par.worldx;
            
            ComplOR((ComplOR(:,3)<1),3)=ComplOR((ComplOR(:,3)<1),3)+par.worldy;
            ComplOR((ComplOR(:,3)>par.worldy+1),3)=ComplOR((ComplOR(:,3)>par.worldy+1),3)-par.worldy;
            
            ComplOR(:,4)=floor(ComplOR(:,2)-1)+1+(floor(ComplOR(:,3)-1))*par.worldx;
        end
        %ComplORP
        if sum(ComplORP(:,1)>0)>0
            % move
            X=rand(length(ComplORP),1)*2*pi;
            ComplORP(:,2)=ComplORP(:,2)+par.enzmovement*cos(X);
            ComplORP(:,3)=ComplORP(:,3)-par.enzmovement*sin(X);
            
            ComplORP((ComplORP(:,2)<1),2)=ComplORP((ComplORP(:,2)<1),2)+par.worldx;
            ComplORP((ComplORP(:,2)>par.worldx+1),2)=ComplORP((ComplORP(:,2)>par.worldx+1),2)-par.worldx;
            
            ComplORP((ComplORP(:,3)<1),3)=ComplORP((ComplORP(:,3)<1),3)+par.worldy;
            ComplORP((ComplORP(:,3)>par.worldy+1),3)=ComplORP((ComplORP(:,3)>par.worldy+1),3)-par.worldy;
            
            ComplORP(:,4)=floor(ComplORP(:,2)-1)+1+(floor(ComplORP(:,3)-1))*par.worldx;
        end
        %ComplORT
        if sum(ComplORT(:,1)>0)>0
            % move
            X=rand(length(ComplORT),1)*2*pi;
            ComplORT(:,2)=ComplORT(:,2)+par.enzmovement*cos(X);
            ComplORT(:,3)=ComplORT(:,3)-par.enzmovement*sin(X);
            
            ComplORT((ComplORT(:,2)<1),2)=ComplORT((ComplORT(:,2)<1),2)+par.worldx;
            ComplORT((ComplORT(:,2)>par.worldx+1),2)=ComplORT((ComplORT(:,2)>par.worldx+1),2)-par.worldx;
            
            ComplORT((ComplORT(:,3)<1),3)=ComplORT((ComplORT(:,3)<1),3)+par.worldy;
            ComplORT((ComplORT(:,3)>par.worldy+1),3)=ComplORT((ComplORT(:,3)>par.worldy+1),3)-par.worldy;
            
            ComplORT(:,4)=floor(ComplORT(:,2)-1)+1+(floor(ComplORT(:,3)-1))*par.worldx;
        end
        
        %ComplEP ComplEQ ComplIQ ComplIR
        %ComplOP ComplOT ComplOR ComplORP ComplORT
        
        %% Complex formation
        % EnzA
        
        activeEnzA=find(EnzA(:,1)>0);
        sequence=randperm(length(activeEnzA));
        for m=1:length(activeEnzA)
            mm=activeEnzA(sequence(m));
            SubShere=find(SubS(:,4)==EnzA(mm,4));
            SubShere(:,2)=1;
            SubPhere=find(SubP(:,4)==EnzA(mm,4));
            SubPhere(:,2)=2;
            SubRhere=find(SubR(:,4)==EnzA(mm,4));
            SubRhere(:,2)=3;
            Molhere=[SubShere;SubPhere;SubRhere];
            if ~isempty(Molhere)
                seque2=randperm(size(Molhere,1));
                kk=[par.bindingAS par.bindingAP par.bindingAR];
                for m2=1:size(Molhere,1)
                    mm2=seque2(m2);
                    if rand*100<kk(Molhere(mm2,2))
                        if EnzA(mm,1)==1
                            
                            switch Molhere(mm2,2)
                                case 1
                                    SubS(Molhere(mm2,1),1:4)=nan;
                                    emptyComplAS=isnan(ComplAS(:,1));
                                    emptyComplAS=find(emptyComplAS>0);
                                    ComplAS(emptyComplAS(1),1)=-1;
                                    ComplAS(emptyComplAS(1),2:4)=EnzA(mm,2:4);
                                case 2
                                    SubP(Molhere(mm2,1),1:4)=nan;
                                    emptyComplAP=isnan(ComplAP(:,1));
                                    emptyComplAP=find(emptyComplAP>0);
                                    ComplAP(emptyComplAP(1),1)=-1;
                                    ComplAP(emptyComplAP(1),2:4)=EnzA(mm,2:4);
                                case 3
                                    SubR(Molhere(mm2,1),1:4)=nan;
                                    emptyComplAR=isnan(ComplAR(:,1));
                                    emptyComplAR=find(emptyComplAR>0);
                                    ComplAR(emptyComplAR(1),1)=-1;
                                    ComplAR(emptyComplAR(1),2:4)=EnzA(mm,2:4);
                            end
                            EnzA(mm,1:4)=nan;
                        end
                    end
                    
                end
            end
            
        end
        
        % EnzE
        activeEnzE=find(EnzE(:,1)>0);
        sequence=randperm(length(activeEnzE));
        for m=1:length(activeEnzE)
            mm=activeEnzE(sequence(m));
            SubPhere=find(SubP(:,4)==EnzE(mm,4));
            SubPhere(:,2)=1;
            SubQhere=find(SubQ(:,4)==EnzE(mm,4));
            SubQhere(:,2)=2;
            Molhere=[SubPhere;SubQhere];
            if ~isempty(Molhere)
                seque2=randperm(size(Molhere,1));
                kk=[par.bindingEP par.bindingEQ];
                for m2=1:size(Molhere,1)
                    mm2=seque2(m2);
                    if rand*100<kk(Molhere(mm2,2))
                        if EnzE(mm,1)==1
                            switch Molhere(mm2,2)
                                case 1
                                    SubP(Molhere(mm2,1),1:4)=nan;
                                    emptyComplEP=isnan(ComplEP(:,1));
                                    emptyComplEP=find(emptyComplEP>0);
                                    ComplEP(emptyComplEP(1),1)=-1;
                                    ComplEP(emptyComplEP(1),2:4)=EnzE(mm,2:4);
                                case 2
                                    SubQ(Molhere(mm2,1),1:4)=nan;
                                    emptyComplEQ=isnan(ComplEQ(:,1));
                                    emptyComplEQ=find(emptyComplEQ>0);
                                    ComplEQ(emptyComplEQ(1),1)=-1;
                                    ComplEQ(emptyComplEQ(1),2:4)=EnzE(mm,2:4);
                            end
                            EnzE(mm,1:4)=nan;
                        end
                    end
                    
                end
            end
            
        end
        
        % EnzI
        activeEnzI=find(EnzI(:,1)>0);
        sequence=randperm(length(activeEnzI));
        for m=1:length(activeEnzI)
            mm=activeEnzI(sequence(m));
            SubQhere=find(SubQ(:,4)==EnzI(mm,4));
            SubQhere(:,2)=1;
            SubRhere=find(SubR(:,4)==EnzI(mm,4));
            SubRhere(:,2)=2;
            Molhere=[SubQhere;SubRhere];
            if ~isempty(Molhere)
                seque2=randperm(size(Molhere,1));
                kk=[par.bindingIQ par.bindingIR];
                for m2=1:size(Molhere,1)
                    mm2=seque2(m2);
                    if rand*100<kk(Molhere(mm2,2))
                        if EnzI(mm,1)==1
                            
                            switch Molhere(mm2,2)
                                case 1
                                    SubQ(Molhere(mm2,1),1:4)=nan;
                                    emptyComplIQ=isnan(ComplIQ(:,1));
                                    emptyComplIQ=find(emptyComplIQ>0);
                                    ComplIQ(emptyComplIQ(1),1)=-1;
                                    ComplIQ(emptyComplIQ(1),2:4)=EnzI(mm,2:4);
                                case 2
                                    SubR(Molhere(mm2,1),1:4)=nan;
                                    emptyComplIR=isnan(ComplIR(:,1));
                                    emptyComplIR=find(emptyComplIR>0);
                                    ComplIR(emptyComplIR(1),1)=-1;
                                    ComplIR(emptyComplIR(1),2:4)=EnzI(mm,2:4);
                            end
                            EnzI(mm,1:4)=nan;
                        end
                    end
                    
                end
            end
            
        end
        
        % EnzO updated
        
        activeEnzO=find(EnzO(:,1)>0);
        sequence=randperm(length(activeEnzO));
        for m=1:length(activeEnzO)
            mm=activeEnzO(sequence(m));
            SubPhere=find(SubP(:,4)==EnzO(mm,4));
            SubPhere(:,2)=1;
            SubThere=find(SubT(:,4)==EnzO(mm,4));
            SubThere(:,2)=2;
            SubRhere=find(SubR(:,4)==EnzO(mm,4));
            SubRhere(:,2)=3;
            Molhere=[SubPhere;SubThere;SubRhere];
            if ~isempty(Molhere)
                seque2=randperm(size(Molhere,1));
                kk=[par.bindingOP par.bindingOT par.bindingOR];
                for m2=1:size(Molhere,1)
                    mm2=seque2(m2);
                    if rand*100<kk(Molhere(mm2,2))
                        if EnzO(mm,1)==1
                            
                            switch Molhere(mm2,2)
                                case 1
                                    SubP(Molhere(mm2,1),1:4)=nan;
                                    emptyComplOP=isnan(ComplOP(:,1));
                                    emptyComplOP=find(emptyComplOP>0);
                                    ComplOP(emptyComplOP(1),1)=-1;
                                    ComplOP(emptyComplOP(1),2:4)=EnzO(mm,2:4);
                                case 2
                                    SubT(Molhere(mm2,1),1:4)=nan;
                                    emptyComplOT=isnan(ComplOT(:,1));
                                    emptyComplOT=find(emptyComplOT>0);
                                    ComplOT(emptyComplOT(1),1)=-1;
                                    ComplOT(emptyComplOT(1),2:4)=EnzO(mm,2:4);
                                case 3
                                    SubR(Molhere(mm2,1),1:4)=nan;
                                    emptyComplOR=isnan(ComplOR(:,1));
                                    emptyComplOR=find(emptyComplOR>0);
                                    ComplOR(emptyComplOR(1),1)=-1;
                                    ComplOR(emptyComplOR(1),2:4)=EnzO(mm,2:4);
                            end
                            EnzO(mm,1:4)=nan;
                        end
                    end
                    
                end
            end
            
        end
        
        % ComplOR
        
        activeComplOR=find(ComplOR(:,1)>0);
        sequence=randperm(length(activeComplOR));
        for m=1:length(activeComplOR)
            mm=activeComplOR(sequence(m));
            SubPhere=find(SubP(:,4)==ComplOR(mm,4));
            SubPhere(:,2)=1;
            SubThere=find(SubT(:,4)==ComplOR(mm,4));
            SubThere(:,2)=2;
            Molhere=[SubPhere;SubThere];
            if ~isempty(Molhere)
                seque2=randperm(size(Molhere,1));
                kk=[par.bindingORP par.bindingORT];
                for m2=1:size(Molhere,1)
                    mm2=seque2(m2);
                    if rand*100<kk(Molhere(mm2,2))
                        if ComplOR(mm,1)==1
                            
                            switch Molhere(mm2,2)
                                case 1
                                    SubP(Molhere(mm2,1),1:4)=nan;
                                    emptyComplORP=isnan(ComplORP(:,1));
                                    emptyComplORP=find(emptyComplORP>0);
                                    ComplORP(emptyComplORP(1),1)=-1;
                                    ComplORP(emptyComplORP(1),2:4)=ComplOR(mm,2:4);
                                case 2
                                    SubT(Molhere(mm2,1),1:4)=nan;
                                    emptyComplORT=isnan(ComplORT(:,1));
                                    emptyComplORT=find(emptyComplORT>0);
                                    ComplORT(emptyComplORT(1),1)=-1;
                                    ComplORT(emptyComplORT(1),2:4)=ComplOR(mm,2:4);
                            end
                            ComplOR(mm,1:4)=nan;
                        end
                    end
                    
                end
            end
        end
        
        % ComplOP +SubR --> ORP
        activeComplOP=find(ComplOP(:,1)>0);
        sequence=randperm(length(activeComplOP));
        for m=1:length(activeComplOP)
            mm=activeComplOP(sequence(m));
            SubRhere=find(SubR(:,4)==ComplOP(mm,4));
            SubRhere(:,2)=1;
            Molhere=[SubRhere];
            if ~isempty(Molhere)
                seque2=randperm(size(Molhere,1));
                kk=[par.bindingOR];
                for m2=1:size(Molhere,1)
                    mm2=seque2(m2);
                    if rand*100<kk(Molhere(mm2,2))
                        if ComplOP(mm,1)==1
                            
                            switch Molhere(mm2,2)
                                case 1
                                    SubR(Molhere(mm2,1),1:4)=nan;
                                    emptyComplORP=isnan(ComplORP(:,1));
                                    emptyComplORP=find(emptyComplORP>0);
                                    ComplORP(emptyComplORP(1),1)=-1;
                                    ComplORP(emptyComplORP(1),2:4)=ComplOP(mm,2:4);
                                    
                            end
                            ComplOP(mm,1:4)=nan;
                        end
                    end
                    
                end
            end
        end
        
        % ComplOT +SubR --> ORT
        activeComplOT=find(ComplOT(:,1)>0);
        sequence=randperm(length(activeComplOT));
        for m=1:length(activeComplOT)
            mm=activeComplOT(sequence(m));
            SubRhere=find(SubR(:,4)==ComplOT(mm,4));
            SubRhere(:,2)=1;
            Molhere=[SubRhere];
            if ~isempty(Molhere)
                seque2=randperm(size(Molhere,1));
                kk=[par.bindingOR];
                for m2=1:size(Molhere,1)
                    mm2=seque2(m2);
                    if rand*100<kk(Molhere(mm2,2))
                        if ComplOT(mm,1)==1
                            
                            switch Molhere(mm2,2)
                                case 1
                                    SubR(Molhere(mm2,1),1:4)=nan;
                                    emptyComplORT=isnan(ComplORT(:,1));
                                    emptyComplORT=find(emptyComplORT>0);
                                    ComplORT(emptyComplORT(1),1)=-1;
                                    ComplORT(emptyComplORT(1),2:4)=ComplOT(mm,2:4);
                                    
                            end
                            ComplOT(mm,1:4)=nan;
                        end
                    end
                    
                end
            end
        end
        
        %% Complex separation
        % ComplAS
        
        activeComplAS=find(ComplAS(:,1)>0);
        sequence=randperm(length(activeComplAS));
        for m=1:length(activeComplAS)
            mm=activeComplAS(sequence(m));
            dice=rand*100;
            if dice<par.releaseAS+par.catAS
                %cat or release
                if dice<par.releaseAS
                    %release ComplAS --> SubS + EnzA
                    emptyEnzA=isnan(EnzA(:,1));
                    emptyEnzA=find(emptyEnzA>0);
                    EnzA(emptyEnzA(1),1)=-1;
                    EnzA(emptyEnzA(1),2:4)=ComplAS(mm,2:4);
                    
                    emptySubS=isnan(SubS(:,1));
                    emptySubS=find(emptySubS>0);
                    SubS(emptySubS(1),1)=-1;
                    SubS(emptySubS(1),2:4)=ComplAS(mm,2:4);
                    
                    ComplAS(mm,1:4)=NaN;
                else
                    %cat ComplAS --> ComplAP
                    emptyComplAP=isnan(ComplAP(:,1));
                    emptyComplAP=find(emptyComplAP>0);
                    ComplAP(emptyComplAP(1),1)=-1;
                    ComplAP(emptyComplAP(1),2:4)=ComplAS(mm,2:4);
                    
                    ComplAS(mm,1:4)=NaN;
                end
            end
        end
        
        % ComplAP
        activeComplAP=find(ComplAP(:,1)>0);
        sequence=randperm(length(activeComplAP));
        for m=1:length(activeComplAP)
            mm=activeComplAP(sequence(m));
            dice=rand*100;
            if dice<par.releaseAP
                %release ComplAP --> SubP + EnzA
                emptyEnzA=isnan(EnzA(:,1));
                emptyEnzA=find(emptyEnzA>0);
                EnzA(emptyEnzA(1),1)=-1;
                EnzA(emptyEnzA(1),2:4)=ComplAP(mm,2:4);
                
                emptySubP=isnan(SubP(:,1));
                emptySubP=find(emptySubP>0);
                SubP(emptySubP(1),1)=-1;
                SubP(emptySubP(1),2:4)=ComplAP(mm,2:4);
                
                ComplAP(mm,1:4)=NaN;
            end
        end
        % ComplAR
        activeComplAR=find(ComplAR(:,1)>0);
        sequence=randperm(length(activeComplAR));
        for m=1:length(activeComplAR)
            mm=activeComplAR(sequence(m));
            dice=rand*100;
            if dice<par.releaseAR
                %release ComplAR --> SubR + EnzA
                emptyEnzA=isnan(EnzA(:,1));
                emptyEnzA=find(emptyEnzA>0);
                EnzA(emptyEnzA(1),1)=-1;
                EnzA(emptyEnzA(1),2:4)=ComplAR(mm,2:4);
                
                emptySubR=isnan(SubR(:,1));
                emptySubR=find(emptySubR>0);
                SubR(emptySubR(1),1)=-1;
                SubR(emptySubR(1),2:4)=ComplAR(mm,2:4);
                
                ComplAR(mm,1:4)=NaN;
            end
        end
        
        % ComplEP
        
        activeComplEP=find(ComplEP(:,1)>0);
        sequence=randperm(length(activeComplEP));
        for m=1:length(activeComplEP)
            mm=activeComplEP(sequence(m));
            dice=rand*100;
            if dice<par.releaseEP+par.catEP
                %cat or release
                if dice<par.releaseEP
                    %release ComplEP --> SubP + EnzE
                    emptyEnzE=isnan(EnzE(:,1));
                    emptyEnzE=find(emptyEnzE>0);
                    EnzE(emptyEnzE(1),1)=-1;
                    EnzE(emptyEnzE(1),2:4)=ComplEP(mm,2:4);
                    
                    emptySubP=isnan(SubP(:,1));
                    emptySubP=find(emptySubP>0);
                    SubP(emptySubP(1),1)=-1;
                    SubP(emptySubP(1),2:4)=ComplEP(mm,2:4);
                    
                    ComplEP(mm,1:4)=NaN;
                else
                    %cat ComplEP --> ComplEQ
                    emptyComplEQ=isnan(ComplEQ(:,1));
                    emptyComplEQ=find(emptyComplEQ>0);
                    ComplEQ(emptyComplEQ(1),1)=-1;
                    ComplEQ(emptyComplEQ(1),2:4)=ComplEP(mm,2:4);
                    
                    ComplEP(mm,1:4)=NaN;
                end
            end
        end
        
        % ComplEQ
        activeComplEQ=find(ComplEQ(:,1)>0);
        sequence=randperm(length(activeComplEQ));
        for m=1:length(activeComplEQ)
            mm=activeComplEQ(sequence(m));
            dice=rand*100;
            if dice<par.releaseEQ
                %release ComplEQ --> SubQ + EnzE
                emptyEnzE=isnan(EnzE(:,1));
                emptyEnzE=find(emptyEnzE>0);
                EnzE(emptyEnzE(1),1)=-1;
                EnzE(emptyEnzE(1),2:4)=ComplEQ(mm,2:4);
                
                emptySubQ=isnan(SubQ(:,1));
                emptySubQ=find(emptySubQ>0);
                SubQ(emptySubQ(1),1)=-1;
                SubQ(emptySubQ(1),2:4)=ComplEQ(mm,2:4);
                
                ComplEQ(mm,1:4)=NaN;
            end
        end
        
        % ComplIQ
        
        activeComplIQ=find(ComplIQ(:,1)>0);
        sequence=randperm(length(activeComplIQ));
        for m=1:length(activeComplIQ)
            mm=activeComplIQ(sequence(m));
            dice=rand*100;
            if dice<par.releaseIQ+par.catIQ
                %cat or release
                if dice<par.releaseIQ
                    %release ComplIQ --> SubQ + EnzI
                    emptyEnzI=isnan(EnzI(:,1));
                    emptyEnzI=find(emptyEnzI>0);
                    EnzI(emptyEnzI(1),1)=-1;
                    EnzI(emptyEnzI(1),2:4)=ComplIQ(mm,2:4);
                    
                    emptySubQ=isnan(SubQ(:,1));
                    emptySubQ=find(emptySubQ>0);
                    SubQ(emptySubQ(1),1)=-1;
                    SubQ(emptySubQ(1),2:4)=ComplIQ(mm,2:4);
                    
                    ComplIQ(mm,1:4)=NaN;
                else
                    %cat ComplIQ --> ComplIR
                    emptyComplIR=isnan(ComplIR(:,1));
                    emptyComplIR=find(emptyComplIR>0);
                    ComplIR(emptyComplIR(1),1)=-1;
                    ComplIR(emptyComplIR(1),2:4)=ComplIQ(mm,2:4);
                    
                    ComplIQ(mm,1:4)=NaN;
                end
            end
        end
        
        % ComplIR
        activeComplIR=find(ComplIR(:,1)>0);
        sequence=randperm(length(activeComplIR));
        for m=1:length(activeComplIR)
            mm=activeComplIR(sequence(m));
            dice=rand*100;
            if dice<par.releaseIR
                %release ComplIR --> SubR + EnzI
                emptyEnzI=isnan(EnzI(:,1));
                emptyEnzI=find(emptyEnzI>0);
                EnzI(emptyEnzI(1),1)=-1;
                EnzI(emptyEnzI(1),2:4)=ComplIR(mm,2:4);
                
                emptySubR=isnan(SubR(:,1));
                emptySubR=find(emptySubR>0);
                SubR(emptySubR(1),1)=-1;
                SubR(emptySubR(1),2:4)=ComplIR(mm,2:4);
                
                ComplIR(mm,1:4)=NaN;
            end
        end
        
        % ComplOP
        
        activeComplOP=find(ComplOP(:,1)>0);
        sequence=randperm(length(activeComplOP));
        for m=1:length(activeComplOP)
            mm=activeComplOP(sequence(m));
            dice=rand*100;
            if dice<par.releaseOP+par.catOP
                %cat or release
                if dice<par.releaseOP
                    %release ComplOP --> SubP + EnzO
                    emptyEnzO=isnan(EnzO(:,1));
                    emptyEnzO=find(emptyEnzO>0);
                    EnzO(emptyEnzO(1),1)=-1;
                    EnzO(emptyEnzO(1),2:4)=ComplOP(mm,2:4);
                    
                    emptySubP=isnan(SubP(:,1));
                    emptySubP=find(emptySubP>0);
                    SubP(emptySubP(1),1)=-1;
                    SubP(emptySubP(1),2:4)=ComplOP(mm,2:4);
                    
                    ComplOP(mm,1:4)=NaN;
                else
                    %cat ComplOP --> ComplOT
                    emptyComplOT=isnan(ComplOT(:,1));
                    emptyComplOT=find(emptyComplOT>0);
                    ComplOT(emptyComplOT(1),1)=-1;
                    ComplOT(emptyComplOT(1),2:4)=ComplOP(mm,2:4);
                    
                    ComplOP(mm,1:4)=NaN;
                end
            end
        end
        
        % ComplOT
        activeComplOT=find(ComplOT(:,1)>0);
        sequence=randperm(length(activeComplOT));
        for m=1:length(activeComplOT)
            mm=activeComplOT(sequence(m));
            dice=rand*100;
            if dice<par.releaseOT
                %release ComplOT --> SubT + EnzO
                emptyEnzO=isnan(EnzO(:,1));
                emptyEnzO=find(emptyEnzO>0);
                EnzO(emptyEnzO(1),1)=-1;
                EnzO(emptyEnzO(1),2:4)=ComplOT(mm,2:4);
                
                emptySubT=isnan(SubT(:,1));
                emptySubT=find(emptySubT>0);
                SubT(emptySubT(1),1)=-1;
                SubT(emptySubT(1),2:4)=ComplOT(mm,2:4);
                
                ComplOT(mm,1:4)=NaN;
            end
        end
        
        % ComplORP    updated
        
        activeComplORP=find(ComplORP(:,1)>0);
        sequence=randperm(length(activeComplORP));
        for m=1:length(activeComplORP)
            mm=activeComplORP(sequence(m));
            dice=rand*100;
            if (dice<par.releaseORP+par.catORP+par.releaseOR)&&(dice>=par.releaseORP+par.catORP)
                %release ComplORP --> SubR + ComplOP
                emptyComplOP=isnan(ComplOP(:,1));
                emptyComplOP=find(emptyComplOP>0);
                ComplOP(emptyComplOP(1),1)=-1;
                ComplOP(emptyComplOP(1),2:4)=ComplORP(mm,2:4);
                
                emptySubR=isnan(SubR(:,1));
                emptySubR=find(emptySubR>0);
                SubR(emptySubR(1),1)=-1;
                SubR(emptySubR(1),2:4)=ComplORP(mm,2:4);
                
                ComplORP(mm,1:4)=NaN;
            end
            if dice<par.releaseORP+par.catORP
                %cat or release
                if dice<par.releaseORP
                    %release ComplORP --> SubP + ComplOR
                    emptyComplOR=isnan(ComplOR(:,1));
                    emptyComplOR=find(emptyComplOR>0);
                    ComplOR(emptyComplOR(1),1)=-1;
                    ComplOR(emptyComplOR(1),2:4)=ComplORP(mm,2:4);
                    
                    emptySubP=isnan(SubP(:,1));
                    emptySubP=find(emptySubP>0);
                    SubP(emptySubP(1),1)=-1;
                    SubP(emptySubP(1),2:4)=ComplORP(mm,2:4);
                    
                    ComplORP(mm,1:4)=NaN;
                else
                    %cat ComplORP --> ComplORT
                    emptyComplORT=isnan(ComplORT(:,1));
                    emptyComplORT=find(emptyComplORT>0);
                    ComplORT(emptyComplORT(1),1)=-1;
                    ComplORT(emptyComplORT(1),2:4)=ComplORP(mm,2:4);
                    
                    ComplORP(mm,1:4)=NaN;
                end
            end
        end
        
        % ComplORT
        activeComplORT=find(ComplORT(:,1)>0);
        sequence=randperm(length(activeComplORT));
        for m=1:length(activeComplORT)
            mm=activeComplORT(sequence(m));
            dice=rand*100;
            if dice<par.releaseORT+par.releaseOR
                %release T or release R
                if dice<par.releaseORT
                    %release ComplORT --> SubT + ComplOR
                    emptyComplOR=isnan(ComplOR(:,1));
                    emptyComplOR=find(emptyComplOR>0);
                    ComplOR(emptyComplOR(1),1)=-1;
                    ComplOR(emptyComplOR(1),2:4)=ComplORT(mm,2:4);
                    
                    emptySubT=isnan(SubT(:,1));
                    emptySubT=find(emptySubT>0);
                    SubT(emptySubT(1),1)=-1;
                    SubT(emptySubT(1),2:4)=ComplORT(mm,2:4);
                    
                    ComplORT(mm,1:4)=NaN;
                else
                    %release ComplORT --> SubR + ComplOT
                    emptyComplOT=isnan(ComplOT(:,1));
                    emptyComplOT=find(emptyComplOT>0);
                    ComplOT(emptyComplOT(1),1)=-1;
                    ComplOT(emptyComplOT(1),2:4)=ComplORT(mm,2:4);
                    
                    emptySubR=isnan(SubR(:,1));
                    emptySubR=find(emptySubR>0);
                    SubR(emptySubR(1),1)=-1;
                    SubR(emptySubR(1),2:4)=ComplORT(mm,2:4);
                    
                    ComplORT(mm,1:4)=NaN;
                end
            end
        end
        
        % ComplOR
        activeComplOR=find(ComplOR(:,1)>0);
        sequence=randperm(length(activeComplOR));
        for m=1:length(activeComplOR)
            mm=activeComplOR(sequence(m));
            dice=rand*100;
            if dice<par.releaseOR
                %release ComplOR --> SubR + EnzO
                emptyEnzO=isnan(EnzO(:,1));
                emptyEnzO=find(emptyEnzO>0);
                EnzO(emptyEnzO(1),1)=-1;
                EnzO(emptyEnzO(1),2:4)=ComplOR(mm,2:4);
                
                emptySubR=isnan(SubR(:,1));
                emptySubR=find(emptySubR>0);
                SubR(emptySubR(1),1)=-1;
                SubR(emptySubR(1),2:4)=ComplOR(mm,2:4);
                
                ComplOR(mm,1:4)=NaN;
            end
        end
        
        %%
        % Adding SubS
        
        Feed=10*0;  %0.4
        
        Feed=Feed+NonInt;
        
        NonInt=Feed-floor(Feed);
        Feed=floor(Feed);
        
        if Feed>0
            NewS=zeros(Feed,4);
            X = lhsdesign(Feed,2);
            NewS(:,1)=-1;
            NewS(:,2)=X(:,1)*par.worldx+1;
            NewS(:,3)=X(:,2)*par.worldy+1;
            NewS(:,4)=floor(NewS(:,2)-1)+1+(floor(NewS(:,3)-1))*par.worldx;
            emptySubS=isnan(SubS(:,1));
            emptySubS=find(emptySubS>0);
            SubS(emptySubS(1:Feed),:)=NewS;
        end
        
        %%
        %Dilusion
        
        %        D=1;
        %        if D>0
        %            D=floor(D*100);
        %            ini=par.worldx*par.worldy*0.75-par.worldx/2+1-floor(D/2);
        %            fim=ini+D-1;
        %            KillSub=ini:fim;
        %            for m=1:D
        %                SubS(find(SubS(:,4)==KillSub(m)),:)=nan;
        %                SubP(find(SubP(:,4)==KillSub(m)),:)=nan;
        %                SubQ(find(SubQ(:,4)==KillSub(m)),:)=nan;
        %                SubR(find(SubR(:,4)==KillSub(m)),:)=nan;
        %                SubT(find(SubT(:,4)==KillSub(m)),:)=nan;
        %            end
        %        end
        
        D=0.01*0;
        if D>0
            
            NumS=sum(SubS(:,1)>0);
            NumP=sum(SubP(:,1)>0);
            NumQ=sum(SubQ(:,1)>0);
            NumR=sum(SubR(:,1)>0);
            NumT=sum(SubT(:,1)>0);
             if NumS>0
                 ActiveSubS=find(SubS(:,1)>0);
                 for m=1:NumS
                     if rand<D
                         SubS(ActiveSubS(m),1:4)=NaN;
                     end
                 end
             end
            if NumP>0
                 ActiveSubP=find(SubP(:,1)>0);
                 for m=1:NumP
                     if rand<D
                         SubP(ActiveSubP(m),1:4)=NaN;
                     end
                 end
            end
             if NumQ>0
                 ActiveSubQ=find(SubQ(:,1)>0);
                 for m=1:NumQ
                     if rand<D
                         SubQ(ActiveSubQ(m),1:4)=NaN;
                     end
                 end
             end
            if NumR>0
                 ActiveSubR=find(SubR(:,1)>0);
                 for m=1:NumR
                     if rand<D
                         SubR(ActiveSubR(m),1:4)=NaN;
                     end
                 end
            end
             if NumT>0
                 ActiveSubT=find(SubT(:,1)>0);
                 for m=1:NumT
                     if rand<D
                         SubT(ActiveSubT(m),1:4)=NaN;
                     end
                 end
             end
        end
        
        
        %%
        
        SubS(SubS(:,1)==-1,1)=1;
        ComplAS(ComplAS(:,1)==-1,1)=1;
        EnzA(EnzA(:,1)==-1,1)=1;
        ComplAP(ComplAP(:,1)==-1,1)=1;
        SubP(SubP(:,1)==-1,1)=1;
        
        SubQ(SubQ(:,1)==-1,1)=1;
        SubR(SubR(:,1)==-1,1)=1;
        SubT(SubT(:,1)==-1,1)=1;
        
        EnzE(EnzE(:,1)==-1,1)=1;
        EnzI(EnzI(:,1)==-1,1)=1;
        EnzO(EnzO(:,1)==-1,1)=1;
        
        ComplAR(ComplAR(:,1)==-1,1)=1;
        ComplEP(ComplEP(:,1)==-1,1)=1;
        ComplEQ(ComplEQ(:,1)==-1,1)=1;
        
        ComplIQ(ComplIQ(:,1)==-1,1)=1;
        ComplIR(ComplIR(:,1)==-1,1)=1;
        
        ComplOP(ComplOP(:,1)==-1,1)=1;
        ComplOT(ComplOT(:,1)==-1,1)=1;
        ComplOR(ComplOR(:,1)==-1,1)=1;
        ComplORP(ComplORP(:,1)==-1,1)=1;
        ComplORT(ComplORT(:,1)==-1,1)=1;
        
        time=time+1;
    end
    
    disp( strcat('Time per iteration:',num2str(toc(tic2)/par.maxtime),'s Simulation time:',num2str(toc(tic2)),'s or ',num2str(toc(tic2)/60),'min'))
    if par.addnum>1
        clear deadsheep
        clear deadwolves
        clear emptysheep
        clear emptywolves
        clear livesheep
        clear livewolves
        clear m
        clear mapa
        clear mm
        clear sequence
        clear seque2
        clear sheephere
        clear tobefeed
        clear grass
        clear green
        clear sheep
        clear wolves
        clear grasstemp
    end
    
    if exist('Noutput','var')
        Noutput=Noutput+output;
        NTotaloutput=NTotaloutput+totaloutput;
    else
        NTotaloutput=totaloutput;
        Noutput=output;
    end
    
    
end
disp([toc(tic1)/par.addnum/60 toc(tic1)/60])
z=Noutput/k;
figure(5)
plot(z(:,1),z(:,2:end))
% legend('S','AS','A','AR','AP','P','EP','E','EQ','Q','IQ','I','IR','R','OP','O','OT','T','ORP','OR','ORT')%
legend('S','P','Q','T','R')%
% legend('S','AS','A','AP','P')%

format shortG; disp(z(600,2:end)-z(100,2:end)); format short


