clear, close all ,clc
 load('NMCM35ADATA.mat')
name=string(celldata.cellname);
Namestring=string(fieldnames(celldata.data));
nomcap=celldata.datatable.Capacity';
cyclestring=string(fieldnames(celldata.cycledata));
cyclesets=extractBetween(cyclestring, "D", "C");
cyclesets= double(strrep(cyclesets,"_","."));
cyclesets(:,2)=double(extractBetween(extractAfter(cyclestring,"o"), "d", "D"))/100;
cyclesets(:,3)=double(extractBetween(cyclestring,"C","c"))+273.15;
chargecycle=extractBetween(extractAfter(cyclestring,"C"),"c","C");
cyclesets(:,4)=double(strrep(chargecycle,"_","."));
k=1;
indexcycle=0;
for i=1:length(cyclesets(:,1))
    for j=1:length(cyclesets(:,1))
        if cyclesets(i,3)==cyclesets(j,3) && cyclesets(i,1)~=cyclesets(j,1) && 0==sum(ismember(i,indexcycle))
            indexcycle(k)=i;
            Cycles(i).diff=diff(celldata.cycledata.(cyclestring{i}).Cycles);
            Cycles(i).throughput=celldata.cycledata.(cyclestring{i}).Cycles*nomcap*2;
            k=k+1;
        end
    end
end


totsoh=[];
totah=[];
totdis=[];
totchar=[];
for i=1:length(cyclestring)
   totah=[totah Cycles(i).throughput'];
   l=length(Cycles(i).throughput);
   totdis=[totdis ones(1,l)*cyclesets(i,1)];
   totchar=[totchar ones(1,l)*cyclesets(i,4)];
   totsoh=[totsoh celldata.cycledata.(cyclestring{i}).SoH'];
end
 g=fittype(@(A,B,Char,Ah,Cdis) (1-A.*(exp(B.*Cdis)+exp(B.*Char)).*Ah./2), dependent="SoH",independent=["Ah" "Cdis"],coefficients=["A" "B"],problem=["Char"]);
 X=[totah;totdis]';
 lower=[ 0 0 ];
 upper=[ 10 10];
 start=[ 0 4.3264];
 ft=fit(X,totsoh',g,'problem',totchar(1),'Lower',lower,'Upper',upper,'Start',start);
 plot(ft,X,totsoh');

