clc, clear, close all
celldata.data.d3C = renamevars(readtable('nmcdischarge3_0C.csv','Delimiter','semi'),["x","y"],["DoD","V"]);
celldata.data.d3C.Temperature=(25+273.15)*(ones(height(celldata.data.d3C),1));
celldata.data.d2C = renamevars(readtable('nmcdischarge2_0C.csv','Delimiter','semi'),["x","y"],["DoD","V"]);
celldata.data.d2C.Temperature=(25+273.15)*(ones(height(celldata.data.d2C),1));
celldata.data.d1C = renamevars(readtable('nmcdischarge1_0C.csv','Delimiter','semi'),["x","y"],["DoD","V"]);
celldata.data.d1C.Temperature=(25+273.15)*(ones(height(celldata.data.d1C),1));
celldata.data.d0_5C = renamevars(readtable('nmcdischarge0_5C.csv','Delimiter','semi'),["x","y"],["DoD","V"]);
celldata.data.d0_5C.Temperature=(25+273.15)*(ones(height(celldata.data.d0_5C),1));
celldata.data.d0_2C = renamevars(readtable('nmcdischarge0_2C.csv','Delimiter','semi'),["x","y"],["DoD","V"]);
celldata.data.d0_2C.Temperature=(25+273.15)*(ones(height(celldata.data.d0_2C),1));
celldata.data.d0_5C55 = renamevars(readtable('nmcdischarge0_5C55Celcius.csv','Delimiter','semi'),["x","y"],["DoD","V"]);
celldata.data.d0_5C55.Temperature=(55+273.15)*(ones(height(celldata.data.d0_5C55),1));
celldata.data.d0_2Cmin20 = renamevars(readtable('nmcdischarge0_2Cmin20Celcius.csv','Delimiter','semi'),["x","y"],["DoD","V"]);
celldata.data.d0_2Cmin20.Temperature=(-20+273.15)*(ones(height(celldata.data.d0_2Cmin20),1));
celldata.curvepoints=renamevars(readtable('curvepoints.csv','Delimiter','semi'),["x","y"],["DoD","V"]);
%% Data from datasheet
name={};
unit={};
data=[];
% Capacity
name=[name,'Capacity'];
unit=[unit,"Ah"];
data=[data,2.6];
name=[name,'CapacityCurrent'];
unit=[unit,"A"];
data=[data,0.52];
% Voltage data
name=[name,'Nom_Voltage'];
unit=[unit,"V"];
data=[data,3.7];
name=[name,'Max_Voltage'];
unit=[unit,"V"];
data=[data,4.2];
name=[name,'Min_Voltage'];
unit=[unit,"V"];
data=[data,3.0];
% Current data
name=[name,'Max_Charge_Current'];
unit=[unit,"A"];
data=[data,1.3];
name=[name,'Max_disCharge_Current'];
unit=[unit,"A"];
data=[data,2.6];
% Temperature limits
name=[name,'Max_Charge_Temperature'];
unit=[unit,"K"];
data=[data,45+273.15];
name=[name,'Min_Charge_Temperature'];
unit=[unit,"K"];
data=[data,0+273.15];
name=[name,'Max_disCharge_Temperature'];
unit=[unit,"K"];
data=[data,60+273.15];
name=[name,'Min_disCharge_Temperature'];
unit=[unit,"K"];
data=[data,-20+273.15];
% weight
name=[name,'Weight'];
unit=[unit,"kg"];
data=[data,46.5e-3];
celldata.datatable=array2table(data);
celldata.datatable.Properties.VariableNames=name;
celldata.datatable.Properties.VariableUnits=unit;
celldata.cellname='NMC1';
celldata.realcellname='LIR18650 2600mAh';
clearvars unit data name 
save NMC1DATA.mat
