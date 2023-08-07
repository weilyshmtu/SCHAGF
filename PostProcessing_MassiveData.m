% To obtain the results of TAGCSC

clear
clc
fullpath = mfilename('fullpath');
[path,name,ext] = fileparts(fullpath);
cd(path)
addpath '.\Databases\Face Image Databases'
addpath '.\Databases\MNIST'
addpath '.\Databases\COIL'
addpath '.\Measurements'
addpath '.\Algorithms'
addpath(genpath(pwd));
database = {'ORL','Umist','COIL20','MNIST','COIL40','YALEB'};
alpha = [0.01, 0.1, 0.5, 1, 2, 5, 10, 50];
beta =  [0.0001, 0.001, 0.01, 0.1, 0.5, 1, 2, 5, 10, 50];
para_rows = length(alpha);
para_cols = length(beta);

for dataindex = 5:5%length(database)

    % load data
    load(database{dataindex});
    X = fea;
    L = gnd;
    nbcluster = max(unique(L));
    
    NeighborSize = 3:15;
    post_acc_results = cell(1, length(NeighborSize));
    post_nmi_results = cell(1, length(NeighborSize));

    acc_array = 0;
    nmi_array = 0;

    tdir = ".\Results\" + database{dataindex};
    cd(tdir)
    eta_list = [0,0.05:0.1:0.95, 1];
    for e = 6:length(eta_list)
        subtdir = num2str(eta_list(e),3);
        cd(subtdir)
        for ni = 1:length(NeighborSize)
            acc_array = min(acc_array, 0);
            nmi_array = min(nmi_array, 0);
            for r=1:para_rows
                for c=1:para_cols
                    tfilename = num2str(r) + "_" + num2str(c);
                    load(tfilename)
                    tC = C;
                    [Z] = refinecoefficient(tC, NeighborSize(ni));
        
                    [idx,~] = clu_ncut(Z,nbcluster);
                    acc_array(r, c) = compacc(idx',L)
                    nmi_array(r, c) = nmi(L, idx')
                end
            end
            post_acc_results{ni} = acc_array
            post_nmi_results{ni} = nmi_array
        end
        for j=1:length(post_acc_results)
            max_acc_array(j) = max(max(post_acc_results{j}));
            max_nmi_array(j) = max(max(post_nmi_results{j}));
        end
        max_acc_array
        max_nmi_array
        %     
        cd ..
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
        filename = "Post_Processing_" +database{dataindex} + "_" + strname ;
        save(filename, "post_acc_results", "post_nmi_results", "max_acc_array", "max_nmi_array")
        
        
    end
end


