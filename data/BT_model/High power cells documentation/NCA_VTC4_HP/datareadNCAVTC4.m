clc, clear, close all
% Tref=25;%reference temperature [C]
varNames = {'DoD','V'};
celldata.data.d0_1C = renamevars(readtable('NCA1discharge0_1C.csv','Delimiter',',','DecimalSeparator','.'),["Var1","Var2"],["DoD","V"]);
celldata.data.d0_24C = renamevars(readtable('NCA1discharge0_24C.csv','Delimiter',',','DecimalSeparator','.'),["Var1","Var2"],["DoD","V"]);
celldata.data.d0_475C = renamevars(readtable('NCA1discharge0_475C.csv','Delimiter',',','DecimalSeparator','.'),["Var1","Var2"],["DoD","V"]);
celldata.data.d0_95C = renamevars(readtable('NCA1discharge0_95C.csv','Delimiter',',','DecimalSeparator','.'),["Var1","Var2"],["DoD","V"]);
celldata.data.d4_76C45 = renamevars(readtable('NCA1discharge4_76C45Celcius.csv','Delimiter',',','DecimalSeparator','.'),["Var1","Var2"],["DoD","V"]);
celldata.data.d4_76C60 = renamevars(readtable('NCA1discharge4_76C60Celcius.csv','Delimiter',',','DecimalSeparator','.'),["Var1","Var2"],["DoD","V"]);
celldata.curvepoints=renamevars(readtable('curvepoints1.csv','Delimiter','semi','DecimalSeparator','.'),["Var1","Var2"],["DoD","V"]);
celldata.cycledata.dod100D14_0C23c1_9C=renamevars(readtable('nca1Cycle100DOD14_0C23Celcius','Delimiter',',','DecimalSeparator','.'),["Var1","Var2"],["Cycles","SoH"]);
celldata.cycledata.dod100D4_76C23c1_9C=renamevars(readtable('nca1Cycle100DOD4_76C23Celcius','Delimiter',',','DecimalSeparator','.'),["Var1","Var2"],["Cycles","SoH"]);
% Capacity to DoD 
names=fieldnames(celldata.data);

%% Data from datasheet
name={};
unit={};
data=[];
% Capacity
name=[name,'Capacity'];
unit=[unit,"Ah"];
data=[data,2.1];
name=[name,'CapacityCurrent'];
unit=[unit,"A"];
data=[data,0.2*2.1];
name=[name,'CurvepointCurrent'];
unit=[unit,"A"];
data=[data,2];
% Voltage data
name=[name,'Nom_Voltage'];
unit=[unit,"V"];
data=[data,3.7];
name=[name,'Max_Voltage'];
unit=[unit,"V"];
data=[data,4.25];
name=[name,'Min_Voltage'];
unit=[unit,"V"];
data=[data,2.5];
% Current data
name=[name,'Max_Charge_Current'];
unit=[unit,"A"];
data=[data,12];
name=[name,'Max_disCharge_Current'];
unit=[unit,"A"];
data=[data,30];
% Temperature limits
name=[name,'Max_Charge_Temperature'];
unit=[unit,"K"];
data=[data,60+273.15];
name=[name,'Min_Charge_Temperature'];
unit=[unit,"K"];
data=[data,0+273.15];
name=[name,'Max_disCharge_Temperature'];
unit=[unit,"K"];
data=[data,80+273.15];
name=[name,'Min_disCharge_Temperature'];
unit=[unit,"K"];
data=[data,-20+273.15];
name=[name,'Reference_Temperature'];
unit=[unit,"K"];
data=[data,23+273.15];
% weight
name=[name,'Weight'];
unit=[unit,"kg"];
data=[data,0.045];
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
celldata.cellname='NCAVTC4';
celldata.chemistry="NCA";
celldata.realcellname='US18650VTC4 NCA 2.1Ah';
clearvars unit data name Tref DODorCAPACITY capacityarray
save NCAVTC4DATA.mat
