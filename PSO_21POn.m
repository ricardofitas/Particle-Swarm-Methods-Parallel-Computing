% application of Particle Swarm Optimisation (PSO) - inertia weight and
% restriction factor - Eberhart et al., 2000 - Parallel On. Computing
% Ricardo Fitas, 2020

%% initial code
try
    myCluster = parcluster('local');
    parpool(myCluster.NumWorkers);
   
end
main_PSOst = evalin('base','main_PSOst');
if main_PSOst.selectfunc == 1
    1
    main_PSOst = Part1(main_PSOst);
end
if main_PSOst.selectfunc == 2
    2
    main_PSOst = Part2(main_PSOst);
end
if main_PSOst.selectfunc == 3
    3
    main_PSOst = Part3(main_PSOst);
end
assignin('base','main_PSOst',main_PSOst)

%% functions
function main_PSOst = Part1(main_PSOst)
%% handle of variables
X_min = main_PSOst.X_min;
X_max = main_PSOst.X_max;
P = main_PSOst.P;
D = main_PSOst.D;

%% code
main_PSOst.X= round(rand(D,P).*(X_max-X_min)+X_min);
main_PSOst.xpbest = main_PSOst.X;
end

function main_PSOst = Part2(main_PSOst)
%% handle of variables
K = main_PSOst.K;
X = main_PSOst.X;
P = main_PSOst.P;
XG = main_PSOst.XG;
YG = main_PSOst.YG;
ZG = main_PSOst.ZG;
f1 = main_PSOst.f1;
user_push1 = evalin('base','user_push1');
detaux = str2func(main_PSOst.KDSstr{1});
%[Scherwinkel,XT,YT,ZT]=detaux(K, X(:,1), XG, YG, ZG);
%% code
parfor ii = 1:P
        assignin('base','user_push1','user_push1');
        [Scherwinkel,XT,YT,ZT]=detaux(K, X(:,ii), XG, YG, ZG);
        ii
        f1(ii) = max(max(abs(Scherwinkel)));
        Scherwinkelmat(:,:,ii) = Scherwinkel;
        XTmat(:,:,ii) =  XT;
        YTmat(:,:,ii) =  YT;
        ZTmat(:,:,ii) =  ZT;
end

%% return of handle
[main_PSOst.gbest,k] = min(f1);
main_PSOst.xgbest = X(:,k);
main_PSOst.XT = XTmat(:,:,k);
main_PSOst.YT = YTmat(:,:,k);
main_PSOst.ZT = ZTmat(:,:,k);
main_PSOst.Scherwinkel = Scherwinkelmat(:,:,k);
main_PSOst.f1 = f1;
main_PSOst.XTmat = XTmat;
main_PSOst.YTmat = YTmat;
main_PSOst.ZTmat = ZTmat;
main_PSOst.Scherwinkelmat = Scherwinkelmat;
end

function main_PSOst = Part3(main_PSOst)
%% handle of variables
K = main_PSOst.K;
X = main_PSOst.X;
P = main_PSOst.P;
XG = main_PSOst.XG;
YG = main_PSOst.YG;
ZG = main_PSOst.ZG;
f1 = main_PSOst.f1;
X_min = main_PSOst.X_min;
X_max = main_PSOst.X_max;
P = main_PSOst.P;
D = main_PSOst.D;
iter = main_PSOst.iter;
v = main_PSOst.v;
main_PSOst.omega = main_PSOst.o_max - iter*...
    (main_PSOst.o_max-main_PSOst.o_min)/(main_PSOst.maxiter);
c1 = main_PSOst.c1;
c2 = main_PSOst.c2;

%% code
K_K= 2/abs(2-c1-c2-sqrt((c1+c2)^2-4*(c1+c2)));
pbest = main_PSOst.pbest;
xpbest = main_PSOst.xpbest;
assignin('base','main_PSOst','main_PSOst')
detaux = str2func(main_PSOst.KDSstr{1});
f2 = str2func('run');
user_push1 = evalin('base','user_push1');
gbest = main_PSOst.gbest;
xgbest = main_PSOst.xgbest;
commPSO.gbest = gbest;
commPSO.xgbest = xgbest;
commPSO.flag = 0;
commPSO.XT = main_PSOst.XT;
commPSO.YT = main_PSOst.YT;
commPSO.ZT = main_PSOst.ZT;
commPSO.Scherwinkel = main_PSOst.Scherwinkel;

nPSO = parcluster('local').NumWorkers;

spmd
    labSend(commPSO,[1:(labindex-1), (labindex+1):numlabs]);
    if labindex == 1
        gbest = min(main_PSOst.gbest);
    end
    for ii = (labindex-1):(nPSO-1):P
        if labindex ~= 1
            commPSO = labReceive(1);
            xgbest = commPSO.xgbest;
            gbest = commPSO.gbest;
            assignin('base','user_push1','user_push1');
            [Scherwinkel,XT,YT,ZT]=detaux(K, X(:,ii), XG, YG, ZG);
            Scherwinkelmat(:,:,ii) = Scherwinkel;
            XTmat(:,:,ii) =  XT;
            YTmat(:,:,ii) =  YT;
            ZTmat(:,:,ii) =  ZT;
            f1(ii) = max(max(abs(Scherwinkel)));
            pbest(ii)=min(f1(ii), pbest(ii));
    
        gbest2 = min(gbest,pbest(ii));
        if gbest>gbest2
            gbest = gbest2
            string1 = "Found a new gbest at iteration "+ string(iter)+ ": "+ string(gbest);
            message2 ="----- Notification from MATLAB: main_PSOsts PC ------" + newline + string1 + newline + "Kind regards,"+ newline + "Ricardo Fitas";
            message = "MATLAB: Notification of new gbest";
            f2('sendemail.m');
            xgbest = X(:,ii);
            commPSO.xgbest = xgbest;
            commPSO.gbest = gbest;
            commPSO.flag =1;
            commPSO.Scherwinkel = Scherwinkelmat(:,:,ii);
            commPSO.XT = XTmat(:,:,ii);
            commPSO.YT = YTmat(:,:,ii);
            commPSO.ZT = ZTmat(:,:,ii);
            
        end 
        labSend(commPSO,1);
        r1 = rand(1,1);
        r2 = rand(1,1);
        v(:,ii) = K*(v(:,ii)+c1.*r1.*(xpbest(:,ii)-X(:,ii)) + c2.*r2.*(xgbest-X(:,ii)));
        X(:,ii) = X(:,ii) + v(:,ii);
        X(:,ii) = round(min(X_max, max(X_min,X(:,ii))));
        else
            for jjj = 2:numlabs
                data32 = labReceive(jjj);
                if data32.flag == 1
                    gbest2 = min(gbest,data32.gbest);
                    if gbest2 < gbest
                        gbest = gbest2;
                        commPSO = data32;
                        commPSO.flag = 0;
                    end
                end
            end
            labSend(commPSO,2:numlabs);
        end
    end
    
end

%% return of handle
comp1 = commPSO{1};
main_PSOst.gbest = comp1.gbest;
main_PSOst.xgbest = comp1.xgbest;
main_PSOst.XT= comp1.XT;
main_PSOst.YT= comp1.YT;
main_PSOst.ZT= comp1.ZT;
main_PSOst.Scherwinkel =  comp1.Scherwinkel;
nPSO = (nPSO-1);
for jjj = 2:numel(commPSO)
    comp = X{jjj};
    main_PSOst.X(:,(jjj-1):nPSO:P)= comp(:,(jjj-1):nPSO:P);
    comp = v{jjj};
    main_PSOst.v(:,(jjj-1):nPSO:P)= comp(:,(jjj-1):nPSO:P);
    comp = pbest{jjj};
    main_PSOst.pbest(:,(jjj-1):nPSO:P)= comp(:,(jjj-1):nPSO:P);
%     xpbest
%     comp = xpbest{jjj};
%     main_PSOst.xpbest(:,(jjj-1):nPSO:P)= comp(:,(jjj-1):nPSO:P);
	try
        comp = XTmat{jjj};
        main_PSOst.XTmat(:,:,(jjj-1):nPSO:P)= comp(:,:,(jjj-1):nPSO:P);
        comp = YTmat{jjj};
        main_PSOst.YTmat(:,:,(jjj-1):nPSO:P)= comp(:,:,(jjj-1):nPSO:P);
        comp = ZTmat{jjj};
        main_PSOst.ZTmat(:,:,(jjj-1):nPSO:P)= comp(:,:,(jjj-1):nPSO:P);
        comp =Scherwinkelmat{jjj};
        main_PSOst.Scherwinkelmat(:,:,(jjj-1):nPSO:P)= comp(:,:,(jjj-1):nPSO:P);
    end
    
end
main_PSOst.xpbest = xpbest;


end

