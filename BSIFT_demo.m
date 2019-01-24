% This is implementation of "BSIFT: Boosting SIFT Using Principal Component
% Analysis. 22nd Iranian Conference on Electrical Engineering, 20-22 May 
% 2014, Tehran, Iran. M. Fotouhi, S. E. Mirsadeghi, S. Kasaei, K. Faez.
% Please cite our paper if you use this code for your research.

clc ; clear all ; close all
% BSIFT_demo
TT1 = zeros(5,5);
TT2 = zeros(1,5);

% Loading image files
imgPath = 'leuven\'; dCell = dir([imgPath '*.ppm']);
dCell2 = dir([imgPath '*.txt']);
disp('Loading image and Homography transform files');

L = length(dCell);
H = cell(1,length(dCell2));
imgd = cell(1,length(dCell));
img_stg = imgd;

for d = 1:L
    imgd{d} = imread([imgPath dCell(d).name]);
    img_stg{d} = im2single(rgb2gray(imgd{d}));
    if d <= length(dCell2)
        tmp_PCA = load([imgPath dCell2(d).name]);
        H{d} = tmp_PCA;
    end
end

%% Pre-defined Values
% Threshold for matching keypoints and different PCA-dimnesion
tr_bbf = 1.25;
tr_pca = 1.00;
lambda = 5.25;
PCA_energy = 0.5:0.1:0.9;
% PCA_energy = 0.8;
% Distance for considering mathed keypoint as correct matches with respect
% to Homography transform
dist_pixel = 2;

% Pre-Alocation
f = cell(1,L);
d = cell(1,L);

% Variable for storing elapsed time for each method
t_bbf = zeros(1,L-1);
PerKeyTime_BBF = zeros(1,L-1);
t_pca = zeros(length(PCA_energy),L-1);
PerKeyTime_PCA = zeros(length(PCA_energy),L-1);
IdBBF = cell (1,L-1);
NumT  = zeros(1,L-1);
Num_BBF = zeros(2,L-1);
IdPCA = cell(length(PCA_energy),L-1);
di = zeros(L-1,length(PCA_energy));
Num_PCA = cell(1,length(PCA_energy));

%% Compute SIFT Keypoints and Descriptors

% Refrence Image
tic; [f{1},d{1}] = vl_sift(img_stg{1},'EdgeThresh',40);t_sift1 = toc;
t_sift1 = t_sift1/numel(img_stg{1});

for i = 2:L
    clc
    disp('------------  BBF Algorithm  -----------')
    tic;[f{i},d{i}] = vl_sift(img_stg{i},'EdgeThresh',40);t_sift2 = toc;
    t_sift2 = t_sift2/numel(img_stg{i});
    
    fprintf('Matching Image 1 to Image %d , Threshold = %.2f \n',i,tr_bbf)
    % Match keypoints based on BBF algorithm proposed by UBC
    tic;[IdBBF{1,i-1},~] = vl_ubcmatch(d{1},d{i},tr_bbf);t_bbf(i-1) = toc;
    PerKeyTime_BBF(1,i-1) = t_bbf(i-1)/(size(d{1},2)*size(d{i},2));
    % Keypoint Location in the reference image
    Ptmp_bbf = f{1}(1:2,IdBBF{1,i-1}(1,:));
    
    P1_bbf = f{1}(1:2,:);
    P2_bbf = f{i}(1:2,:);
    
    % Transfer keypoint of reference imafe to the next image by
    % Homography transform between two view.
    % Normalized transformed keypoint by third element
    Correct_bbf = TransformPoint(Ptmp_bbf,H{1,i-1});
    P3_bbf = TransformPoint(P1_bbf,H{1,i-1});
    % Delete out of image transformed keypoints which have negative location
    P3_bbf(:,sum(P3_bbf <= 0) > 0) = [];
    minNum = min(size(P2_bbf,2),size(P3_bbf,2));
    ref = minNum == size(P2_bbf,2);
    indx = zeros(1,minNum);
    for z = 1:minNum
        if ref == true
            Key = kron(P2_bbf(:,z),ones(1,size(P3_bbf,2)));
            isExist = double(sum(sum(abs(P3_bbf - Key) < dist_pixel) == 2) == 1);
            NumT(i-1) = NumT(i-1) + isExist;
        else
            Key = kron(P3_bbf(:,z),ones(1,size(P2_bbf,2)));
            isExist = double(sum(sum( abs(P2_bbf - Key) < dist_pixel) == 2) == 1);
            NumT(i-1) = NumT(i-1) + isExist;
        end
    end
    % Matched keypoint by BBF-UBC algorithm
    Matched_BBF = f{i}(1:2,IdBBF{1,i-1}(2,:));
    % Compute number of mistmatched keypoint if they have
    % distance more than "dist_pixel" variable
    Num_wrong = sum(sum(abs(Correct_bbf - Matched_BBF) < dist_pixel) ~= 2);
    Num_total = size(Matched_BBF,2);
    Num_BBF(:,i-1) = [Num_total Num_wrong];
end

%% Our Method - PCA in all feature vectors

for Index = 1:length(PCA_energy)
    for i = 2:L
        clc;disp('------------  PCA Algorithm  -----------')
        Dcat = [d{1} d{i}];
        % PCA of Concatenated descriptors
        tic;[Y,di(i-1,Index)] = comp_pca(double(Dcat),PCA_energy(Index));tmp_pca = toc;
        a1 = Y(:,1:size(d{1},2));
        a2 = Y(:,size(d{1},2)+1:end);
        fprintf('Matching Image 1 to Image %d ,PCA_Dim = %d , Threshold = %.2f \n',i,PCA_energy(Index),tr_pca)
        timer1 = tic;
        [IdPCA{Index,i-1},~] = vl_ubcmatch(a1,a2,tr_pca);
        
        D1 = d{1}(:,IdPCA{Index,i-1}(1,:));
        D2 = d{i}(:,IdPCA{Index,i-1}(2,:));
        
        M = SKLD(D1,D2,lambda);
        KL_tr = 0.075;
        Idx = (M > KL_tr)';
        
        IdPCA{Index,i-1} = IdPCA{Index,i-1}(:,Idx);
        
        t_pca(Index,i-1) = tmp_pca + toc(timer1);
        PerKeyTime_PCA(Index,i-1) = t_pca(Index,i-1)/(size(d{1},2)*size(d{i},2));
        P_tmp_pca = (f{1}(1:2,IdPCA{Index,i-1}(1,:)));
        P_pca = [P_tmp_pca;ones(1,size(P_tmp_pca,2))];
        P_pca = H{1,i-1} * P_pca;
        
        Correct_Pca = [P_pca(1,:)./P_pca(3,:) ; P_pca(2,:)./P_pca(3,:)];
        Matched_Pca = f{i}(1:2,IdPCA{Index,i-1}(2,:));
        
        Num_wrong = sum(sum(abs(Correct_Pca-Matched_Pca)< dist_pixel ) ~=2 );
        Num_total = size(Matched_Pca,2);
        Num_PCA{Index}(:,i-1) = [Num_total Num_wrong];
    end
end
%% Plot Results
% PerKeyTime_PCA = TT1/KK;
% PerKeyTime_BBF = TT2/KK;

dbbf = (100*Num_BBF(2,:)./Num_BBF(1,:));
for i = 1:length(Num_PCA)
    dpca(i,:) = (100*Num_PCA{i}(2,:)./Num_PCA{i}(1,:));
end

%%

incr = [dbbf;dpca];
TimePerKey = [PerKeyTime_BBF;PerKeyTime_PCA]*1e+6;
Co = kron(TimePerKey(1,:)',ones(1,size(TimePerKey,1)-1))';
Speed_Up = TimePerKey;
Speed_Up(1,:) = NaN;
Speed_Up(2:end,:) = Co./TimePerKey(2:end,:);

plot_flag = 1;
if plot_flag == 1
    figure
    plot(incr','-s','LineWidth',2.5);grid on
    title('Quality Measure','FontName','Times New Roman','FontWeight','bold','FontSize',12.5);grid on
    ylabel('Correct Rate (%)');xlabel('Image Pairs')
    set(gca,'XTick',1:1:5)
    set(gca,'XTickLabel',{'I1-->I2','I1-->I3','I1-->I4','I1-->I5','I1-->I6'})
    legend('SIFT','BSIFT-50%','BSIFT-60%','BSIFT-70%','BSIFT-80%','BSIFT-90%')
    axis([1 5 0 100])
    
    figure
    plot(TimePerKey','-s','LineWidth',2.5);grid on
    title('Elapsed Time','FontName','Times New Roman','FontWeight','bold','FontSize',12.5);grid on
    ylabel('Elapsed Time (\musec)');xlabel('Image Pairs')
    set(gca,'XTick',1:1:5)
    set(gca,'XTickLabel',{'I1-->I2','I1-->I3','I1-->I4','I1-->I5','I1-->I6'})
    legend('SIFT','BSIFT-50%','BSIFT-60%','BSIFT-70%','BSIFT-80%','BSIFT-90%')
    
    figure
    plot(Speed_Up','-s','LineWidth',2.5);grid on
    title('Speed Comparison','FontName','Times New Roman','FontWeight','bold','FontSize',12.5);grid on
    ylabel('Speed Up Ratio');xlabel('Image Pairs')
    set(gca,'XTick',1:1:5)
    set(gca,'XTickLabel',{'I1-->I2','I1-->I3','I1-->I4','I1-->I5','I1-->I6'})
    legend('','BSIFT-50%','BSIFT-60%','BSIFT-70%','BSIFT-80%','BSIFT-90%')
    
end
tmp_PCA = [];
for m = 1:size(Num_PCA,2)
    tm = [Num_PCA{1,m}(1,:);Num_PCA{1,m}(1,:)-Num_PCA{1,m}(2,:)];
    tmp_PCA = [tmp_PCA;tm];
end
tmp_BBF=[Num_BBF(1,:);Num_BBF(1,:)-Num_BBF(2,:)];

Pr_bbf = (tmp_BBF(2,:)./tmp_BBF(1,:));
for i = 1:length(Num_PCA)
    ttp = [Num_PCA{1,i}(1,:);Num_PCA{1,i}(1,:)-Num_PCA{1,i}(2,:)];
    Pr_pca(i,:) = (ttp(2,:)./ttp(1,:));
end

Pr = [Pr_bbf;Pr_pca];
BSIFT.Speed_Up = Speed_Up;
BSIFT.TimePerKey = TimePerKey;
BSIFT.Precision = Pr;

BSIFT.PCAMatches = tmp_PCA;
BSIFT.BBFMatches = Num_BBF;

%% Visualize Keypoints

MarkerStyle{1} = '-s';
MarkerStyle{2} = '-o';
MarkerStyle{3} = '-d';
MarkerStyle{4} = '-p';
MarkerStyle{5} = '-h';
MarkerStyle{6} = '-x';

MarkerColor{1} = [1 0 0];
MarkerColor{2} = [0 1 0];
MarkerColor{3} = [0 0 1];
MarkerColor{4} = [0.8 0.8 0];
MarkerColor{5} = [0 0.8 0.8];
MarkerColor{6} = [0.8 0 0.8];
figure
for j = 2:size(BSIFT.Speed_Up,1)
    plot(BSIFT.Speed_Up(j,:),MarkerStyle{j},'LineWidth',2.5,'color',MarkerColor{j},'MarkerSize',10)
    hold on
end
ylabel('Speed Up Ratio','FontName','Times New Roman','FontWeight','bold','FontSize',12.5);
xlabel('Image Pairs'   ,'FontName','Times New Roman','FontWeight','bold','FontSize',12.5);
set(gca,'XTick',1:1:5)
set(gca,'XTickLabel',{'I1-->I2','I1-->I3','I1-->I4','I1-->I5','I1-->I6'})
title('Speed Comparison','FontName','Times New Roman','FontWeight','bold','FontSize',12.5);grid on
hteg = legend('BSIFT-50%','BSIFT-60%','BSIFT-70%','BSIFT-80%','BSIFT-90%','location','SouthWest');
set(hteg,'FontName','Times New Roman','FontWeight','bold','FontSize',10)