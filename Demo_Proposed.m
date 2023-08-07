
%% Image clustering Test
clear, clc
fullpath = mfilename('fullpath');
[path,name,ext] = fileparts(fullpath);
cd(path)
%% Experiment Settings
addpath(genpath(pwd));
% load data
database = {'ORL','Umist','MNIST','COIL20','YALEB','COIL40'};
numdatas =  length(database);

ProjectionType = 0;          % data projection type    
NormalizationType = 2;         % data normalization type               

for dataindex = 6:6%:numdatas   % different databases
    DataName = database{dataindex}
    load(database{dataindex});

    X = fea;
    L = gnd;

    if min(unique(L)) == 0
        L = L + 1;
    end
    nbcluster = max(unique(L));
   

    % projection
    switch ProjectionType
        case 0
            X = X;
        case 1
            X =  DataProjection(X, dim);
    end
 
    % normalization
    switch NormalizationType
        case 0
            X = X;
        case 1
            X = mexNormalize(X);
        case 2 
            if max(max(X)) > 1
                X = X./repmat((255)*ones(1,size(X,2)),size(X,1),1);
            end
    end

    % parameters 
    alpha = [0.01, 0.1, 0.5, 1, 2, 5, 10, 50];
    beta =  [0.0001, 0.001, 0.01, 0.1, 0.5, 1, 2, 5, 10, 50];
%     alpha = [1e-5,1e-4,1e-3,5e-3,0.01,0.05,0.1,0.5];
%     beta  = [1e-5,1e-4,1e-3,5e-3,0.01,0.05,0.1,0.5];
%     eta_list = 0:0.1:1;
    eta_list = [1,0.95:-0.1:0.05,0]
    for e =12:-1:10%9:-1:1
        Timecell = []
        acc_array = []
        nmi_array = []
         for a=1:length(alpha)
            for b=1:length(beta)
                tic
                [H, C] =  sc_agf0(X, alpha(a), beta(b), eta_list(e));        
                toc

                Timecell(a,b) = toc
                C1 = C;
                [idx,~] = clu_ncut(C1,nbcluster);
                acc_array(a, b) = compacc(idx',L)
                nmi_array(a, b) = nmi(L, idx')
                
              
                % save features and coefficient matrices
                folder ="./Results/" + database{dataindex};
                if exist(folder,'dir')~=7
                    mkdir(folder);
                end
                cd(folder)
                subfolder = num2str(eta_list(e),3);
                if ~exist(folder+"/"+subfolder,'dir')
                    mkdir(subfolder);
                end
                cd(subfolder)
                tfilename = num2str(a) + "_" + num2str(b) ;
                save(tfilename, "H", "C")  
                cd ..
                cd ..
                cd ..
            end
        end
    
        % save accuracy and nmi
        cd(folder)
        if eta_list(e)==0
            strname = '0_0';
        else
            if eta_list(e)==1
                 strname = '1_0';
            else
                str = num2str(eta_list(e),3);
                sstr = split(str,'.');
                strname = sstr{1} + "_" + sstr{2};
            end
        end
        filename = DataName + "_" + num2str(ProjectionType)+ "_"+ num2str(NormalizationType) + "_" + strname;
        save(filename, "acc_array", "nmi_array", "Timecell")
        cd ..
        cd ..
    end
    
end