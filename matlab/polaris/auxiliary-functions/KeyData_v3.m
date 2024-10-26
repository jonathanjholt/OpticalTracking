function [nextpeakstart,rowpos,FlexExtValues,VarValValues,IntExtValues,AntPosValues,TensionValues] = KeyData_v3(VarVal,IntExt,FlexExt,AntPost,LCvalues)

angles = 11;

rowpos = zeros(angles,4); % rows are different flex angles from 0-100 deg, columns are 1st,2nd,3rd repeats
FlexExtValues = zeros(angles,4);
VarValValues = zeros(angles,4);
IntExtValues = zeros(angles,4);
AntPosValues = zeros(angles,4);
TensionValues= zeros(angles,4);
val=0; %we want to target 0deg as full extension, but on occasion it doesn't quite reach this

% this 'for' loop finds the row position when the knee is extended back to zero deg, thus defining our three repeats 
oneortwo = 1;
row= find(round(FlexExt) >=-10);
% row= find(round(FlexExt) ==0);
% row= find(round(FlexExt) <=-50);

for i = 2:length(row)
    if row(i)-row(i-1) > 20
        nextpeakstart(oneortwo)=row(i);
        oneortwo=oneortwo +1;
    end
end
%     rowpos(1,1)= 1;
%     rowpos(1,2)= nextpeakstart(1); 
%     rowpos(1,3)= nextpeakstart(2);

for i=1:angles
    rowpos(i,1)= find((FlexExt(1:nextpeakstart(1))) <= (i-1)*-10,1,'first');
    rowpos(i,2)= find((FlexExt(nextpeakstart(1):nextpeakstart(2))) <= (i-1)*-10,1,'first')+nextpeakstart(1); 
    rowpos(i,3)= find((FlexExt(nextpeakstart(2):end)) <= (i-1)*-10,1,'first')+nextpeakstart(2);
end
% for i=1:angles
%     a=abs(FlexExt - ((i-1)*-10));
%     rowpos(i,1)= find((a(1:nextpeakstart(1))) == min(a(1:nextpeakstart(1))),1,'first');
%     rowpos(i,2)= find((a(nextpeakstart(1):nextpeakstart(2))) == min(a(nextpeakstart(1):nextpeakstart(2))),1,'first')+nextpeakstart(1); 
%     rowpos(i,3)= find((a(nextpeakstart(2):end)) == min(a(nextpeakstart(2):end)),1,'first')+nextpeakstart(2);
% end
% 
for i=1
    b=abs(FlexExt - val);
    rowpos(i,1)= find((b(1:nextpeakstart(1))) == min(b(1:nextpeakstart(1))),1,'first');
    rowpos(i,2)= find((b(nextpeakstart(1):nextpeakstart(2))) == min(b(nextpeakstart(1):nextpeakstart(2))),1,'first')+nextpeakstart(1); 
    rowpos(i,3)= find((b(nextpeakstart(2):end)) == min(b(nextpeakstart(2):end)),1,'first')+nextpeakstart(2);
end

for i=angles+1
    rowpos(i,1)= find((FlexExt(1:nextpeakstart(1))) == min(FlexExt(1:nextpeakstart(1))),1,'first');
    rowpos(i,2)= find((FlexExt(nextpeakstart(1):nextpeakstart(2))) == min(FlexExt(nextpeakstart(1):nextpeakstart(2))),1,'first')+nextpeakstart(1); 
    rowpos(i,3)= find((FlexExt(nextpeakstart(2):end)) == min(FlexExt(nextpeakstart(2):end)),1,'first')+nextpeakstart(2);
end
    
for i=1:angles+1
    for j=1:3
        FlexExtValues(i,j)= FlexExt(rowpos(i,j));
        VarValValues(i,j)= VarVal(rowpos(i,j));
        IntExtValues(i,j)= IntExt(rowpos(i,j));
        AntPosValues(i,j)= AntPost(rowpos(i,j));
        %TensionValues(i,j)=LCvalues(rowpos(i,j));
    end
end

for i=1:angles+1
    FlexExtValues(i,4) = (FlexExtValues(i,1)+FlexExtValues(i,2)+FlexExtValues(i,3))/3;
    VarValValues(i,4) = (VarValValues(i,1)+VarValValues(i,2)+VarValValues(i,3))/3;
    IntExtValues(i,4) = (IntExtValues(i,1)+IntExtValues(i,2)+IntExtValues(i,3))/3;
    AntPosValues(i,4) = (AntPosValues(i,1)+AntPosValues(i,2)+AntPosValues(i,3))/3;
   % TensionValues(i,4) = (TensionValues(i,1)+TensionValues(i,2)+TensionValues(i,3))/3;
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



end


