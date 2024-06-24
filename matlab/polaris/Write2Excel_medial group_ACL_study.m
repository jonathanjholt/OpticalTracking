function [DataTable] = Write2Excel(FlexExt,VarVal,IntExt,AntPos,file2writeto,kneestate,loadtype,total_kneestates,LCvalues)

if strcmp(loadtype,'Neutral.csv')
    Columncounter=1;
elseif strcmp(loadtype,'Anterior.csv')
    Columncounter=2;
elseif strcmp(loadtype,'Posterior.csv')
    Columncounter=3;
elseif strcmp(loadtype,'External.csv')
    Columncounter=4;
elseif strcmp(loadtype,'Internal.csv')
    Columncounter=5;
elseif strcmp(loadtype,'Varus.csv')
    Columncounter=6;
elseif strcmp(loadtype,'Valgus.csv')
    Columncounter=7;
elseif strcmp(loadtype,'Anterior+ER.csv')
    Columncounter=8;
elseif strcmp(loadtype,'Pivot shift.csv')
    Columncounter=9;        
else
    Columncounter=10;
end
    
disp(Columncounter);
disp(loadtype)
Flexext_mean=FlexExt(1:11,4);
Varval_mean=VarVal(1:11,4);
Intext_mean=IntExt(1:11,4);
Antpos_mean=AntPos(1:11,4);
LCvalues_mean=LCvalues(1:11,4);

DataTable=[Flexext_mean,Varval_mean,Intext_mean,Antpos_mean,LCvalues_mean];

column_flexext=char(65+Columncounter);
column_antpos=char(65+Columncounter+10);
if Columncounter < 6
column_intext=char(65+Columncounter+20);
else 
column_intext=strcat(char(65),char(59+Columncounter));
end
column_varval=strcat(char(65),char(69+Columncounter));
column_LC=strcat(char(65),char(79+Columncounter));

filename=file2writeto;

if Columncounter==1 && kneestate==1
    for i=1:total_kneestates
        disp('Intact neutral')
        disp(i)
    writematrix(Flexext_mean,filename,'Sheet',i,'Range',strcat(column_flexext,'3'));
    writematrix(Antpos_mean,filename,'Sheet',i,'Range',strcat(column_antpos,'3'));
    writematrix(Intext_mean,filename,'Sheet',i,'Range',strcat(column_intext,'3'));
	writematrix(Varval_mean,filename,'Sheet',i,'Range',strcat(column_varval,'3'));
    end
elseif Columncounter==1
        disp('non intact neutral')
    writematrix(Flexext_mean,filename,'Sheet',kneestate,'Range',strcat(column_flexext,'31'));
    writematrix(Antpos_mean,filename,'Sheet',kneestate,'Range',strcat(column_antpos,'31'));
    writematrix(Intext_mean,filename,'Sheet',kneestate,'Range',strcat(column_intext,'31'));
    writematrix(Varval_mean,filename,'Sheet',kneestate,'Range',strcat(column_varval,'31'));
    writematrix(LCvalues_mean,filename,'Sheet',kneestate,'Range',strcat(column_LC,'3'));
else
    disp('non intact not neutral')
    writematrix(Flexext_mean,filename,'Sheet',kneestate,'Range',strcat(column_flexext,'3'));
    writematrix(Antpos_mean,filename,'Sheet',kneestate,'Range',strcat(column_antpos,'3'));
    writematrix(Intext_mean,filename,'Sheet',kneestate,'Range',strcat(column_intext,'3'));
    writematrix(Varval_mean,filename,'Sheet',kneestate,'Range',strcat(column_varval,'3'));
    writematrix(LCvalues_mean,filename,'Sheet',kneestate,'Range',strcat(column_LC,'3'));
end



 

