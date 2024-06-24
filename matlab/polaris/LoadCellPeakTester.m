folderPath='C:\Users\Sander\Documents\Imperial\matlab\Tests';
fileName=[folderPath,'\DSC_test_cut.csv'];
rawdata = xlsread(fileName);
LCvalues=(rawdata(:,2));

angles = 11;

rowpos = zeros(angles,4); % rows are different flex angles from 0-100 deg, columns are 1st,2nd,3rd repeats
Tension = zeros(angles,4);

val=0; %we want to target 0deg as full extension, but on occasion it doesn't quite reach this

% this 'for' loop finds the row position when the knee is extended back to zero deg, thus defining our three repeats 
oneortwo = 1;
row= find(round(LCvalues) <=1);

% row= find(round(FlexExt) ==0);
% row= find(round(FlexExt) <=-50);

for i = 2:length(row)
    if row(i)-row(i-1) > 20
        nextpeakstart(oneortwo)=row(i);
        oneortwo=oneortwo +1;
    end
end
disp(row)
disp(nextpeakstart)



    rowpos(1,1)= 1;
    rowpos(1,2)= nextpeakstart(1); 
    rowpos(1,3)= nextpeakstart(2);

for i=1:angles
    rowpos(i,1)= find((LCvalues(1:nextpeakstart(1))) <= (i-1)*-10,1,'first');
    rowpos(i,2)= find((LCvalues(nextpeakstart(1):nextpeakstart(2))) <= (i-1)*-10,1,'first')+nextpeakstart(1); 
    rowpos(i,3)= find((LCvalues(nextpeakstart(2):end)) <= (i-1)*-10,1,'first')+nextpeakstart(2);
end

% { for i=1:angles
%     a=abs(FlexExt - ((i-1)*-10));
%     rowpos(i,1)= find((a(1:nextpeakstart(1))) == min(a(1:nextpeakstart(1))),1,'first');
%      rowpos(i,2)= find((a(nextpeakstart(1):nextpeakstart(2))) == min(a(nextpeakstart(1):nextpeakstart(2))),1,'first')+nextpeakstart(1); 
%      rowpos(i,3)= find((a(nextpeakstart(2):end)) == min(a(nextpeakstart(2):end)),1,'first')+nextpeakstart(2);
%  end
%  
% for i=1
%    b=abs(FlexExt - val);
%     rowpos(i,1)= find((b(1:nextpeakstart(1))) == min(b(1:nextpeakstart(1))),1,'first');
%     rowpos(i,2)= find((b(nextpeakstart(1):nextpeakstart(2))) == min(b(nextpeakstart(1):nextpeakstart(2))),1,'first')+nextpeakstart(1); 
%     rowpos(i,3)= find((b(nextpeakstart(2):end)) == min(b(nextpeakstart(2):end)),1,'first')+nextpeakstart(2);
% end
% 
for i=angles+1
    rowpos(i,1)= find((LCvalues(1:nextpeakstart(1))) == min(LCvalues(1:nextpeakstart(1))),1,'first');
    rowpos(i,2)= find((LCvalues(nextpeakstart(1):nextpeakstart(2))) == min(LCvalues(nextpeakstart(1):nextpeakstart(2))),1,'first')+nextpeakstart(1); 
    rowpos(i,3)= find((LCvalues(nextpeakstart(2):end)) == min(LCvalues(nextpeakstart(2):end)),1,'first')+nextpeakstart(2);
end
    
for i=1:angles+1
    for j=1:3
       Tension(i,j)= LCvalues(rowpos(i,j));
    end
end

for i=1:angles+1
    Tension(i,4) = (Tension(i,1)+Tension(i,2)+Tension(i,3))/3;
end
  


% third=length(Flex)/3; %separate data into thirds to pature 3 repeats
% third=round(third); %round to nearest integer
% nextpeakstart = 1;
% 
% for i=1:11
%     rowpos(i,1)= find(Flex(1:third) <= (i-1)*-10,1,'first');
%     FlexValue(i,1)= Flex(rowpos(i,1));
%     VarValue(i,1)= Var(rowpos(i,1));
%     IntValue(i,1)= Int(rowpos(i,1));
%     AntValue(i,1)= Ant(rowpos(i,1));
%     
%     rowpos(i,2)= (find(Flex(third+1:2*third) <= (i-1)*-10,1,'first'))+third;
%     FlexValue(i,2)= Flex(rowpos(i,2));
%     VarValue(i,2)= Var(rowpos(i,2));
%     IntValue(i,2)= Int(rowpos(i,2));
%     AntValue(i,2)= Ant(rowpos(i,2));
%     
%     rowpos(i,3)= (find(Flex(2*third+1:end) <= (i-1)*-10,1,'first'))+2*third;
%     FlexValue(i,3)= Flex(rowpos(i,3));
%     VarValue(i,3)= Var(rowpos(i,3));
%     IntValue(i,3)= Int(rowpos(i,3));
%     AntValue(i,3)= Ant(rowpos(i,3));
% end



%end


