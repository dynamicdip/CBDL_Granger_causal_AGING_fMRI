    clear all;
    clc;
    close all;
    
    %%%% YOUNG %%%%%%%%%%%%
    
    load('avg_SC_young.mat','avg_SC_young');
    ROIs=csvimport('ROI_names.csv','noHeader',true);
    thresh_SC=threshold_proportional(avg_SC_young,1); % 30 percent of strongest weights preserved
    adjSC= weight_conversion(thresh_SC,'binarize'); %binarized adjacency matrix
    thresh_SC= weight_conversion(thresh_SC,'normalize');
   
    distSC= distance_bin(adjSC);% shortest path length between each pair
    L=weight_conversion(thresh_SC,'lengths') ;%weight to length matrix
    
    [distwSC,short_ed]= distance_wei(L); %weighted shortest path length
    
   % t=visualize(adjSC);
    mycolor=[0 0 0;0 0 1;1 0 0];
%%% Degree Distribution (Degree Centrality)%%%
    figure;
    deg=degrees_und(adjSC);
    [deg_sort,I]=sort(deg,'descend'); 
    ROIs_sorted=ROIs(I);
    b=bar(deg_sort);
    title('Degree distribution');
    y=mean(deg_sort)+std(deg_sort);
    set(0,'defaulttextinterpreter','none')
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );
    setbarcolor(b,deg_sort,y);
    Ind_deg=deg>y; 
    deg_young=deg;
    
     %%% Strength Distribution %%%
   hold on;
    figure;
    str=strengths_und(thresh_SC);
    [str_sort,S]=sort(str,'descend');
    ROIs_sorted=ROIs(S);
   s= bar(str_sort);
    title('Strength distribution');
     y=mean(str_sort)+std(str_sort);
    set(0,'defaulttextinterpreter','none')
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );
     setbarcolor(s,str_sort,y);
   Ind=str>y;
    str_young=str;
    
    %%% Betweenness %%%
    hold on; 
    figure;
    BC = betweenness_wei(L);
    [BC_sort,B]=sort(BC,'descend');
    ROIs_sorted=ROIs(B);
    bc=bar(BC_sort);
    y=mean(BC_sort)+std(BC_sort);

    title('Betweenness Centrality-Weighted');
    set(0,'defaulttextinterpreter','none')
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );
    %setbarcolor(bc,BC_sort,Ind(B));
    setbarcolor(bc,BC_sort,y);
    BC_young=BC;

    
    %%%% Eigenvector Centrality %%%
    hold on; 
    figure;
    v = eigenvector_centrality_und(thresh_SC);
    [eig_sort,e]=sort(v,'descend');
    ROIs_sorted=ROIs(e);
    eig=bar(eig_sort);
    y=mean(eig_sort)+std(eig_sort);

    title('Eigenvector Centrality-Weighted');
    set(0,'defaulttextinterpreter','none')
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );
   % setbarcolor(eig,eig_sort,Ind(e));
    setbarcolor(eig,eig_sort,y);

    v_young=v;

    %%%% closeness centrality %%%% 
%      hold on; 
%     figure;
%     clos = closeness(thresh_SC);
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
    
    %%community structure%%
    Q_young=0;
    for i=1:500
             [Ci0 Q0]= modularity_louvain_und_sign(thresh_SC);
             Q0=round(Q0*1e4)/1e4;
             Q_list(i)=Q0;
             Ci_iter(i,:)=Ci0;
                         
         end
         [Q0,i] =max(Q_list);
          Ci=Ci_iter(i,:);
          while Q0>Q_young
             [Ci Q_young] = modularity_finetune_und_sign(thresh_SC,'gja',Ci);  
             Q_young=round(Q_young*1e4)/1e4;
            end
  [Ci Q_young]= community_louvain(thresh_SC,1,Ci);
   
   [Ci_sorted,Ci_I]=sort(Ci);
   figure;
   imagesc(thresh_SC(Ci_I,Ci_I));
   hold on;                                 % hold on to overlay community visualization

   [X,Y,INDSORT] = grid_communities(Ci); % call function
   imagesc(thresh_SC(INDSORT,INDSORT));           % plot ordered adjacency matrix
   hold on;                                 % hold on to overlay community visualization
    plot(X,Y,'y','linewidth',3);
   
    Z_score=module_degree_zscore(thresh_SC,Ci,0);
    [Z_sort,Z]=sort(Z_score);
    figure;
    ROIs_sorted=ROIs(Z);
    z=bar(Z_sort);
    title('Module-degree-Z-score');
    set(0,'defaulttextinterpreter','none')
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );%     [M,Q] = community_louvain(thresh_SC,0.3,Ci0,'potts');
    setbarcolor(z,Z_sort,2);
        Ci_young=Ci;
%     %% k-core analysis %%
%     [coreness,kn]=kcoreness_centrality_bu(adjSC);
%     [core_sort,K]=sort(coreness);
%     figure;
%     ROIs_sorted=ROIs(K);
%     core=bar(core_sort);
%     title('k-core decomposition');
%     set(0,'defaulttextinterpreter','none')
%     set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
%     rotateXLabels( gca(), 90 );
%     setbarcolor(core,core_sort,Ind_deg(K));
     
    %%Participation Coefficient %%
    
     hold on; 
    figure;
    PC=participation_coef(thresh_SC,Ci,0);
    [PC_sort,P]=sort(PC,'descend');
    ROIs_sorted=ROIs(P);
    p= bar(PC_sort);
    title('Participation Coefficient');
    set(0,'defaulttextinterpreter','none')
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );
    y=mean(PC_sort)+std(PC_sort);
%     setbarcolor(p,PC_sort,Ind(P));
      setbarcolor(p,PC_sort,y);

    PC_young=PC;
    
  %%%% Overlapping PC %%%%
  M=link_communities(adjSC);
  O_PC=sum(M,1);
  [O_PC_sort,OP]=sort(O_PC,'descend');
  ROIs_sorted=ROIs(OP);
  p= bar(O_PC_sort);
    title('Overlapping Participation Coefficient');
    set(0,'defaulttextinterpreter','none')
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );
    y=mean(O_PC_sort)+std(O_PC_sort);
%     setbarcolor(p,PC_sort,Ind(P));
   setbarcolor(p,O_PC_sort,y);
   O_PC_young=O_PC; 
    
  %%%% OLD %%%%%
  
    
    load('avg_SC_old.mat','avg_SC_old');
    ROIs=csvimport('ROI_names.csv','noHeader',true);
    thresh_SC=threshold_proportional(avg_SC_old,1); % 30 percent of strongest weights preserved
    adjSC= weight_conversion(thresh_SC,'binarize'); %binarized adjacency matrix
    thresh_SC= weight_conversion(thresh_SC,'normalize');
   
    distSC= distance_bin(adjSC);% shortest path length between each pair
    L=weight_conversion(thresh_SC,'lengths') ;%weight to length matrix
    
    [distwSC,short_ed]= distance_wei(L); %weighted shortest path length
    
   % t=visualize(adjSC);
    mycolor=[0 0 0;0 0 1;1 0 0];
%%% Degree Distribution (Degree Centrality)%%%
    figure;
    deg=degrees_und(adjSC);
    [deg_sort,I]=sort(deg,'descend'); 
    ROIs_sorted=ROIs(I);
    b=bar(deg_sort);
    title('Degree distribution');
    y=mean(deg_sort)+std(deg_sort);
    set(0,'defaulttextinterpreter','none')
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );
    setbarcolor(b,deg_sort,y);
    Ind_deg=deg>y;
    deg_w=deg;
    
    %%% Strength Distribution %%%
   hold on;
    figure;
    str=strengths_und(thresh_SC);
    [str_sort,S]=sort(str,'descend');
    ROIs_sorted=ROIs(S);
   s= bar(str_sort);
    title('Strength distribution');
     y=mean(str_sort)+std(str_sort);
    set(0,'defaulttextinterpreter','none')
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );
     setbarcolor(s,str_sort,y);
     Ind=str>y;

    %%% Betweenness %%%
    hold on; 
    figure;
    BC = betweenness_wei(L);
    [BC_sort,B]=sort(BC,'descend');
    ROIs_sorted=ROIs(B);
    bc=bar(BC_sort);
    y=mean(BC_sort)+std(BC_sort);

    title('Betweenness Centrality-Weighted');
    set(0,'defaulttextinterpreter','none')
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );
%     setbarcolor(bc,BC_sort,Ind(B));
    setbarcolor(bc,BC_sort,y);

    
    %%%% Eigenvector Centrality %%%
    hold on; 
    figure;
    v = eigenvector_centrality_und(thresh_SC);
    [eig_sort,e]=sort(v,'descend');
    ROIs_sorted=ROIs(e);
    eig=bar(eig_sort);
    y=mean(eig_sort)+std(eig_sort);

    title('Eigenvector Centrality-Weighted');
    set(0,'defaulttextinterpreter','none')
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );
%     setbarcolor(eig,eig_sort,Ind(e));
    setbarcolor(eig,eig_sort,y);


    %%%% closeness centrality %%%% 
%      hold on; 
%     figure;
%     clos = closeness(thresh_SC);
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
    
    %%community structure%%
    Q_old=0;
    for i=1:500
             [Ci0 Q0]= modularity_louvain_und_sign(thresh_SC);
             Q0=round(Q0*1e4)/1e4;
             Q_list(i)=Q0;
             Ci_iter(i,:)=Ci0;
                         
         end
         [Q0,i] =max(Q_list);
          Ci=Ci_iter(i,:);
          while Q0>Q_old
             [Ci Q_old] = modularity_finetune_und_sign(thresh_SC,'gja',Ci);  
             Q_old=round(Q_old*1e4)/1e4;
             end
   
   
      [Ci_sorted,Ci_I]=sort(Ci);
   figure;
   imagesc(thresh_SC(Ci_I,Ci_I));
   hold on;                                 % hold on to overlay community visualization

   [X,Y,INDSORT] = grid_communities(Ci); % call function
   imagesc(thresh_SC(INDSORT,INDSORT));           % plot ordered adjacency matrix
   hold on;                                 % hold on to overlay community visualization
    plot(X,Y,'y','linewidth',3);
   

   Z_score=module_degree_zscore(thresh_SC,Ci,0);
    [Z_sort,Z]=sort(Z_score);
    figure;
    ROIs_sorted=ROIs(Z);
    z=bar(Z_sort);
    title('Module-degree-Z-score');
    set(0,'defaulttextinterpreter','none')
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );%     [M,Q] = community_louvain(thresh_SC,0.3,Ci0,'potts');
    setbarcolor(z,Z_sort,2);
    
    %% k-core analysis %%
%     [coreness,kn]=kcoreness_centrality_bu(adjSC);
%     [core_sort,K]=sort(coreness);
%     figure;
%     ROIs_sorted=ROIs(K);
%     core=bar(core_sort);
%     title('S-core decomposition');
%     set(0,'defaulttextinterpreter','none')
%     set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
%     rotateXLabels( gca(), 90 );
%     setbarcolor(core,core_sort,Ind_deg(K));
%     setbarcolor(core,core_sort,y);

    %% Participation Coefficient %%
    
     
    hold on; 
    figure;
    PC=participation_coef(thresh_SC,Ci,0);
    [PC_sort,P]=sort(PC,'descend');
    ROIs_sorted=ROIs(P);
    p= bar(PC_sort);
    title('Participation Coefficient');
    set(0,'defaulttextinterpreter','none')
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );
    y=mean(PC_sort)+std(PC_sort);
%     setbarcolor(p,P,Ind(P));  
     setbarcolor(p,PC_sort,y);  

       %%%% Overlapping PC %%%%
  M=link_communities(adjSC);
  O_PC=sum(M,1);
  [O_PC_sort,OP]=sort(O_PC,'descend');
  ROIs_sorted=ROIs(OP);
  p= bar(O_PC_sort);
    title('Overlapping Participation Coefficient');
    set(0,'defaulttextinterpreter','none')
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );
    y=mean(O_PC_sort)+std(O_PC_sort);
%     setbarcolor(p,PC_sort,Ind(P));
      setbarcolor(p,O_PC_sort,y);
     