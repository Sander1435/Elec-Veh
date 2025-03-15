clc, clear, close all
% Tref=25;%reference temperature [C]
varNames = {'DoD','V'};
celldata.data.d0_2C = renamevars(readtable('NMC1discharge0_2C.csv','Delimiter','semi','DecimalSeparator','.'),["Var1","Var2"],["DoD","V"]);
celldata.data.d0_5C = renamevars(readtable('NMC1discharge0_5C.csv','Delimiter','semi','DecimalSeparator','.'),["Var1","Var2"],["DoD","V"]);
celldata.data.d1_0C = renamevars(readtable('NMC1discharge1_0C.csv','Delimiter','semi','DecimalSeparator','.'),["Var1","Var2"],["DoD","V"]);
celldata.data.d2_85C = renamevars(readtable('NMC1discharge2_85C.csv','Delimiter','semi','DecimalSeparator','.'),["Var1","Var2"],["DoD","V"]);
celldata.data.d0_5C45 = renamevars(readtable('NMC1discharge0_5C45Celcius.csv','Delimiter',',','DecimalSeparator','.'),["Var1","Var2"],["DoD","V"]);
celldata.data.d0_5C0 = renamevars(readtable('NMC1discharge0_5C0Celcius.csv','Delimiter','semi','DecimalSeparator','.'),["Var1","Var2"],["DoD","V"]);
celldata.curvepoints=renamevars(readtable('curvepoints2.csv','Delimiter','semi','DecimalSeparator','.'),["Var1","Var2"],["DoD","V"]);
celldata.cycledata.dod100D2_3C23c0_5C=renamevars(readtable('NMC1Cycle100DOD2_3C23Celcius.csv','Delimiter','semi','DecimalSeparator','.'),["Var1","Var2"],["Cycles","SoH"]);
celldata.cycledata.dod100D1_15C23c0_5C=renamevars(readtable('NMC1Cycle100DOD1_15C23Celcius.csv','Delimiter','semi','DecimalSeparator','.'),["Var1","Var2"],["Cycles","SoH"]);
celldata.cycledata.dod100D1_7C23c0_5C=renamevars(readtable('NMC1Cycle100DOD1_7C23Celcius.csv','Delimiter','semi','DecimalSeparator','.'),["Var1","Var2"],["Cycles","SoH"]);
% Capacity to DoD 
names=fieldnames(celldata.data);

%% Data from datasheet
name={};
unit={};
data=[];
% Capacity
name=[name,'Capacity'];
unit=[unit,"Ah"];
data=[data,3.5];
name=[name,'CapacityCurrent'];
unit=[unit,"A"];
data=[data,0.2*3.5];
name=[name,'CurvepointCurrent'];
unit=[unit,"A"];
data=[data,0.2*3.5];
% Voltage data
name=[name,'Nom_Voltage'];
unit=[unit,"V"];
data=[data,3.6];
name=[name,'Max_Voltage'];
unit=[unit,"V"];
data=[data,4.21];
name=[name,'Min_Voltage'];
unit=[unit,"V"];
data=[data,2.5];
% Current data
name=[name,'Max_Charge_Current'];
unit=[unit,"A"];
data=[data,1.7];
name=[name,'Max_disCharge_Current'];
unit=[unit,"A"];
data=[data,10];
% Temperature limits
name=[name,'Max_Charge_Temperature'];
unit=[unit,"K"];
data=[data,60+273.15];
name=[name,'Min_Charge_Temperature'];
unit=[unit,"K"];
data=[data,0+273.15];
name=[name,'Max_disCharge_Temperature'];
unit=[unit,"K"];
data=[data,60+273.15];
name=[name,'Min_disCharge_Temperature'];
unit=[unit,"K"];
data=[data,-40+273.15];
name=[name,'Reference_Temperature'];
unit=[unit,"K"];
data=[data,23+273.15];
% weight
name=[name,'Weight'];
unit=[unit,"kg"];
data=[data,0.048];
% dimensions
name=[name,'length'];
unit=[unit,"m"];
data=[data,0.0652];
name=[name,'diameter'];
unit=[unit,"m"];
data=[data,0.0186];
%% dod adjustment
DODorCAPACITY=2; %0 for DoD, 1 for Capacity from N-0, 2 for Capacity from 0-N

if DODorCAPACITY==1
     celldata.curvepoints.DoD=1-celldata.curvepoints.DoD./data(1);
for i=1:length(names)
   
celldata.data.(names{i}).DoD=1-celldata.data.(names{i}).DoD./data(1);
end
elseif DODorCAPACITY==2
 celldata.curvepoints.DoD=celldata.curvepoints.DoD./data(1);
for i=1:length(names)
celldata.data.(names{i}).DoD=celldata.data.(names{i}).DoD./data(1);
end
end

%% saving data
celldata.datatable=array2table(data);
celldata.datatable.Properties.VariableNames=name;
celldata.datatable.Properties.VariableUnits=unit;
celldata.cellname='NMCM35A';
celldata.chemistry="NMC";
celldata.realcellname='INR-18650-M35A NMC 3.5Ah';
clearvars unit data name Tref DODorCAPACITY capacityarray
save NMCM35ADATA.mat
