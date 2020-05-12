    clear all;
    clc;
    close all;
    ROIs=csvimport('ROI_names.csv','noHeader',true);
   
    %%%% YOUNG POPULATION %%%
    load('avg_FC_young.mat','avg_FC_young'); %average of all young FC thresholded at 10 percent
    thresh_FC_norm= weight_conversion(avg_FC_young,'normalize'); %binarized adjacency matrix

    thresh_FC=threshold_proportional(avg_FC_young,1);
    adjFC= weight_conversion(thresh_FC,'binarize'); %binarized adjacency matrix
    L=weight_conversion(adjFC,'lengths');
    
%%% Degree Distribution %%%
    figure;
    deg_w=degrees_und(adjFC);
    [deg_sort,I]=sort(deg_w,'descend'); 
    ROIs_sorted=ROIs(I);
    b=bar(deg_sort);
    deg_young=deg_w;
    y=mean(deg_sort)+std(deg_sort);
    title('Degree distribution-Younger Group');
    set(0,'defaulttextinterpreter','none');
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );
    setbarcolor(b,deg_sort,y);
    deg_Ind=deg_w>y;
%     
% %%% Strength Distribution %%%
%    
    figure;
    str=strengths_und(thresh_FC_norm);
    [str_sort,S]=sort(str,'descend');
    ROIs_sorted=ROIs(S);
    s=bar(str_sort);
    str_young=str;
    y=mean(str_sort)+std(str_sort);
    title('Strength distribution-Younger Group');
    set(0,'defaulttextinterpreter','none')
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );
    setbarcolor(s,str_sort,y);
    Ind=str>y;
%   
%  %%% Betweenness Centrality %%%
%    
    hold on; 
    figure;
    [EBC BC]=edge_betweenness_wei(L);
    [BC_sort,B]=sort(BC,'descend');
    BC_young=BC;
    ROIs_sorted=ROIs(B);
    bc=bar(BC_sort);
    y=mean(BC_sort)+std(BC_sort);

    title('Betweenness Centrality-Younger Group');
    set(0,'defaulttextinterpreter','none')
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );
%     setbarcolor(bc,BC_sort,deg_Ind(B));
    setbarcolor(bc,BC_sort,y);

    
     %%%% Eigenvector Centrality %%%
    hold on; 
    figure;
    v = eigenvector_centrality_und(thresh_FC_norm);
    [eig_sort,e]=sort(v,'descend');
    ROIs_sorted=ROIs(e);
    eig=bar(eig_sort);
    y=mean(eig_sort)+std(eig_sort);
    v_young=v;
    title('Eigenvector Centrality-Younger Group');
    set(0,'defaulttextinterpreter','none')
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );
%     setbarcolor(eig,eig_sort,deg_Ind(e));
    setbarcolor(eig,eig_sort,y);

    %%% closeness centrality %%%% 
%     hold on; 
%     figure;
%     clos = closeness(thresh_FC);
%     [clos_sort,c]=sort(clos,'descend');
%     ROIs_sorted=ROIs(c);
%     cl=bar(clos_sort);
%     y=mean(clos_sort)+std(clos_sort);
% 
%     title('Closeness Centrality-Weighted');
%     set(0,'defaulttextinterpreter','none')
%     set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
%     rotateXLabels( gca(), 90 );
%     setbarcolor(cl,clos_sort,Ind(c));   
    
  %% Community Structure%%
    Q_young=0;
   
         for i=1:500
             [Ci0 Q0]= modularity_louvain_und_sign(thresh_FC_norm);
             Q0=round(Q0*1e4)/1e4;
             Q_list(i)=Q0;
             Ci_iter(i,:)=Ci0;
                         
         end
         [Q0,i] =max(Q_list);
          Ci0=Ci_iter(i,:);
          while Q0>Q_young
             [Ci Q_young] = modularity_finetune_und_sign(thresh_FC_norm,'pos',Ci0);  
             Q_young=round(Q_young*1e4)/1e4;
             end

         
    
     [Ci_sorted,Ci_I]=sort(Ci);
     figure;
     imagesc(thresh_FC_norm(Ci_I,Ci_I));
     hold on;   % hold on to overlay community visualization
    [X,Y,INDSORT] = grid_communities(Ci); % call function
    imagesc(thresh_FC_norm(INDSORT,INDSORT));           % plot ordered adjacency matrix
    hold on;                                 % hold on to overlay community visualization
    plot(X,Y,'y','linewidth',3);
    title('Community Structure-Younger');

    Ci_young=Ci_I;
    
      Z_score=module_degree_zscore(thresh_FC_norm,Ci,0);
        [Z_sort,Z]=sort(Z_score);
        figure;
        ROIs_sorted=ROIs(Z);
        z=bar(Z_sort);
         title('Module-degree-Z-score-Younger');
         set(0,'defaulttextinterpreter','none')
         set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
        rotateXLabels( gca(), 90 );%     [M,Q] = community_louvain(thresh_SC,0.3,Ci0,'potts');
        setbarcolor(z,Z_sort,1.5);
    
    %%% Participation Coefficient %%%
       hold on; 
    figure;
    PC=participation_coef(thresh_FC_norm,Ci,0);
    [PC_sort,P]=sort(PC,'descend');
    PC_young=PC;
    ROIs_sorted=ROIs(P);
    p= bar(PC_sort);
    title('Participation Coefficient-Younger');
    set(0,'defaulttextinterpreter','none')
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );
    y=mean(PC_sort)+std(PC_sort);
%     setbarcolor(p,PC_sort,deg_Ind(P));
    setbarcolor(p,PC_sort,y);
    
    adjFC_young=adjFC;
    thresh_FC_young=thresh_FC;
    
%%% OLD POPULATION %%% 
 
   load('avg_FC_old.mat','avg_FC_old'); %average of all young FC thresholded at 10 percent

    thresh_FC_norm= weight_conversion(avg_FC_old,'normalize'); %binarized adjacency matrix
    thresh_FC=threshold_proportional(avg_FC_old,1);
    adjFC= weight_conversion(thresh_FC,'binarize'); %binarized adjacency matrix
    L=weight_conversion(adjFC,'lengths');
    
%%% Degree Distribution %%%
    figure;
    deg_w=degrees_und(adjFC);
    [deg_sort,I]=sort(deg_w,'descend'); 
    ROIs_sorted=ROIs(I);
    b=bar(deg_sort);
    y=mean(deg_sort)+std(deg_sort);
    title('Degree distribution-Older Group');
    set(0,'defaulttextinterpreter','none');
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );
    setbarcolor(b,deg_sort,y);
    deg_Ind=deg_w>y;
%     
% %%% Strength Distribution %%%
%    
    figure;
    str=strengths_und(thresh_FC);
    [str_sort,S]=sort(str,'descend');
    ROIs_sorted=ROIs(S);
    s=bar(str_sort);
    y=mean(str_sort)+std(str_sort);
    title('Strength distribution-Older Group');
    set(0,'defaulttextinterpreter','none')
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );
    setbarcolor(s,str_sort,y);
    Ind=str>y;
%   
%  %%% Betweenness Centrality %%%
%    
    hold on; 
    figure;
    [EBC BC]=edge_betweenness_wei(L);
    [BC_sort,B]=sort(BC,'descend');
    ROIs_sorted=ROIs(B);
    bc=bar(BC_sort);
    y=mean(BC_sort)+std(BC_sort);

    title('Betweenness Centrality-Older Group');
    set(0,'defaulttextinterpreter','none')
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );
%     setbarcolor(bc,BC_sort,deg_Ind(B));
    setbarcolor(bc,BC_sort,y);

    
     %%%% Eigenvector Centrality %%%
    hold on; 
    figure;
    v = eigenvector_centrality_und(thresh_FC);
    [eig_sort,e]=sort(v,'descend');
    ROIs_sorted=ROIs(e);
    eig=bar(eig_sort);
    y=mean(eig_sort)+std(eig_sort);

    title('Eigenvector Centrality-Older Group');
    set(0,'defaulttextinterpreter','none')
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );
%     setbarcolor(eig,eig_sort,deg_Ind(e));
    setbarcolor(eig,eig_sort,y);

    %%%% closeness centrality %%%% 
%     hold on; 
%     figure;
%     clos = closeness(thresh_FC);
%     [clos_sort,c]=sort(clos,'descend');
%     ROIs_sorted=ROIs(c);
%     cl=bar(clos_sort);
%     y=mean(clos_sort)+std(clos_sort);
% 
%     title('Closeness Centrality-Weighted');
%     set(0,'defaulttextinterpreter','none')
%     set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
%     rotateXLabels( gca(), 90 );
%     setbarcolor(cl,clos_sort,Ind(c));   
    
  %%% Community Structure%%
    Q_old=0;
   
         for i=1:500
             [Ci0 Q0]= modularity_louvain_und_sign(thresh_FC);
             Q0=round(Q0*1e4)/1e4;
             Q_list(i)=Q0;
             Ci_iter(i,:)=Ci0;
                         
         end
         [Q0,i] =max(Q_list);
          Ci0=Ci_iter(i,:);
          while Q0>Q_old
             [Ci Q_old] = modularity_finetune_und_sign(thresh_FC,'smp',Ci0);  
             Q_old=round(Q_old*1e4)/1e4;
          end
    
     [Ci_sorted,Ci_I]=sort(Ci);
     figure;
     imagesc(thresh_FC(Ci_I,Ci_I));
     hold on;   % hold on to overlay community visualization
    [X,Y,INDSORT] = grid_communities(Ci); % call function
    imagesc(thresh_FC(INDSORT,INDSORT));           % plot ordered adjacency matrix
    hold on;                                 % hold on to overlay community visualization
    plot(X,Y,'y','linewidth',3); 
    title('Community Structure-Older');
    
    
      Z_score=module_degree_zscore(thresh_FC,Ci,0);
        [Z_sort,Z]=sort(Z_score);
        figure;
        ROIs_sorted=ROIs(Z);
        z=bar(Z_sort);
         title('Module-degree-Z-score-Older');
         set(0,'defaulttextinterpreter','none')
         set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
        rotateXLabels( gca(), 90 );%     [M,Q] = community_louvain(thresh_SC,0.3,Ci0,'potts');
        setbarcolor(z,Z_sort,1.5);
    
    %%% Participation Coefficient %%%
       hold on; 
    figure;
    PC=participation_coef(thresh_FC,Ci,0);
    [PC_sort,P]=sort(PC,'descend');
    ROIs_sorted=ROIs(P);
    p= bar(PC_sort);
    title('Participation Coefficient-Older Group');
    set(0,'defaulttextinterpreter','none')
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );
    y=mean(PC_sort)+std(PC_sort);
%     setbarcolor(p,PC_sort,deg_Ind(P));
    setbarcolor(p,PC_sort,y);