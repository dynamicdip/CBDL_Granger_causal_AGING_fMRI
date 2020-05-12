 %%%%%%%%%%% YOUNG %%%%%%%%%%%%
  yBCnorm=(BC_young-min(BC_young))/std(BC_young);
 ydegnorm=(deg_young-min(deg_young))/std(deg_young);
 yvnorm=(v_young-min(v_young))/std(v_young);
 ystrnorm=(str_young-min(str_young))/std(str_young);
 ypcnorm = (PC_young -min(PC_young))/std(PC_young);
% yopcnorm = (O_PC_young -min(O_PC_young))/std(O_PC_young);
 
 for i=1:68 avg_rank(i)=yBCnorm(i)+ydegnorm(i)+yvnorm(i)+ystrnorm(i)+ypcnorm(i); end
    hold on; 
    figure;
    
    [avg_sort,a]=sort(avg_rank,'descend');
    ROIs_sorted=ROIs(a);
    avg_bar=bar(avg_sort);
    y=mean(avg_sort)+std(avg_sort);
    title('Average ranking-Younger Group');
    set(0,'defaulttextinterpreter','none')
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );
    setbarcolor(avg_bar,avg_sort,y);




%%%%%%%%% OLD %%%%%%%%%%%%%
 BCnorm=(BC-min(BC))/std(BC);
 degnorm=(deg_w-min(deg_w))/std(deg_w);
 vnorm=(v-min(v))/std(v);
 strnorm=(str-min(str))/std(str);
 pcnorm = (PC -min(PC))/std(PC);
 %opcnorm=(O_PC-min(O_PC))/std(O_PC);
 
 
 
for i=1:68 avg_rank(i)=BCnorm(i)+degnorm(i)+vnorm(i)+strnorm(i)+pcnorm(i); end
    hold on; 
    figure;
    
    [avg_sort,a]=sort(avg_rank,'descend');
    ROIs_sorted=ROIs(a);
    avg_bar=bar(avg_sort);
    y=mean(avg_sort)+std(avg_sort);
    title('Average ranking-Older Group');
    set(0,'defaulttextinterpreter','none')
    set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
    rotateXLabels( gca(), 90 );
    setbarcolor(avg_bar,avg_sort,y);