clc, clear, close all
% Tref=25;%reference temperature [C]
varNames = {'DoD','V'};
celldata.data.d1_0C = renamevars(readtable('lfp1discharge1_0C.csv','Delimiter','semi','DecimalSeparator',','),["Var1","Var2"],["DoD","V"]);
celldata.data.d2_0C = renamevars(readtable('lfp1discharge2_0C.csv','Delimiter','semi','DecimalSeparator',','),["Var1","Var2"],["DoD","V"]);
celldata.data.d4_0C = renamevars(readtable('lfp1discharge4_0C.csv','Delimiter','semi','DecimalSeparator',','),["Var1","Var2"],["DoD","V"]);
celldata.data.d6_0C = renamevars(readtable('lfp1discharge6_0C.csv','Delimiter','semi','DecimalSeparator',','),["Var1","Var2"],["DoD","V"]);
celldata.data.d8_0C = renamevars(readtable('lfp1discharge8_0C.csv','Delimiter','semi','DecimalSeparator',','),["Var1","Var2"],["DoD","V"]);
celldata.data.d10_0C = renamevars(readtable('lfp1discharge10_0C.csv','Delimiter','semi','DecimalSeparator',','),["Var1","Var2"],["DoD","V"]);
celldata.data.d2_5C45 = renamevars(readtable('lfp1discharge2_5C45Celcius.csv','Delimiter','semi','DecimalSeparator',','),["Var1","Var2"],["DoD","V"]);
celldata.data.d2_5C23 = renamevars(readtable('lfp1discharge2_5C23Celcius.csv','Delimiter','semi','DecimalSeparator',','),["Var1","Var2"],["DoD","V"]);
celldata.curvepoints=renamevars(readtable('curvepoints.csv','Delimiter','semi','DecimalSeparator','.'),["Var1","Var2"],["DoD","V"]);
celldata.cycledata.dod100D1_0C23c1_0C=renamevars(readtable('lfp1Cycle100DOD1_0C23Celcius.csv','Delimiter','semi','DecimalSeparator',','),["Var1","Var2"],["Cycles","SoH"]);
celldata.cycledata.dod100D1_0C45c1_2C=renamevars(readtable('lfp1Cycle100DOD1_0C45Celcius.csv','Delimiter','semi','DecimalSeparator',','),["Var1","Var2"],["Cycles","SoH"]);
celldata.cycledata.dod100D8_0C23c1_2C=renamevars(readtable('lfp1Cycle100DOD8_0C23Celcius.csv','Delimiter','semi','DecimalSeparator',','),["Var1","Var2"],["Cycles","SoH"]);
% Capacity to DoD 
names=fieldnames(celldata.data);
capacityarray=[2.56 2.56 2.56 2.56 2.56 2.56 2.56 2.5 2.5 2.5];

%% Data from datasheet
name={};
unit={};
data=[];
% Capacity
name=[name,'Capacity'];
unit=[unit,"Ah"];
data=[data,2.56];
name=[name,'CapacityCurrent'];
unit=[unit,"A"];
data=[data,0.5*2.56];
name=[name,'CurvepointCurrent'];
unit=[unit,"A"];
data=[data,2.56];
% Voltage data
name=[name,'Nom_Voltage'];
unit=[unit,"V"];
data=[data,3.3];
name=[name,'Max_Voltage'];
unit=[unit,"V"];
data=[data,3.6];
name=[name,'Min_Voltage'];
unit=[unit,"V"];
data=[data,2.0];
% Current data
name=[name,'Max_Charge_Current'];
unit=[unit,"A"];
data=[data,20.0];
name=[name,'Max_disCharge_Current'];
unit=[unit,"A"];
data=[data,120.0];
% Temperature limits
name=[name,'Max_Charge_Temperature'];
unit=[unit,"K"];
data=[data,60+273.15];
name=[name,'Min_Charge_Temperature'];
unit=[unit,"K"];
data=[data,-20+273.15];
name=[name,'Max_disCharge_Temperature'];
unit=[unit,"K"];
data=[data,60+273.15];
name=[name,'Min_disCharge_Temperature'];
unit=[unit,"K"];
data=[data,-30+273.15];
name=[name,'Reference_Temperature'];
unit=[unit,"K"];
data=[data,25+273.15];
% weight
name=[name,'Weight'];
unit=[unit,"kg"];
data=[data,0.076];

%% dod adjustment
DODorCAPACITY=2; %0 for DoD, 1 for Capacity from 0-N, 2 for Capacity from N-0

if DODorCAPACITY==1
     celldata.curvepoints.DoD=1-celldata.curvepoints.DoD./capacityarray(1);
for i=1:length(names)
   
celldata.data.(names{i}).DoD=1-celldata.data.(names{i}).DoD./capacityarray(i+1);
end
elseif DODorCAPACITY==2
 celldata.curvepoints.DoD=celldata.curvepoints.DoD./capacityarray(1);
for i=1:length(names)
celldata.data.(names{i}).DoD=celldata.data.(names{i}).DoD./capacityarray(i+1);
end
end

%% saving data
celldata.datatable=array2table(data);
celldata.datatable.Properties.VariableNames=name;
celldata.datatable.Properties.VariableUnits=unit;
celldata.cellname='LFPM1B';
celldata.chemistry="LFP";
celldata.realcellname='ANR26650M1B LFP 2.5Ah';
clearvars unit data name Tref DODorCAPACITY capacityarray
save LFPM1BDATA.mat
