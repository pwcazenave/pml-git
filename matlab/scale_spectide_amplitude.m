% Open the spectide file and amend the values in the tidal amplitude to
% lower tides for the second type of model scenario.

clear
close all
clc

infile='/data/medusa/pica/models/FVCOM/runCO2_leak/input/configs/inputV7/co2_spectide.nc';

[A,B,C]=fileparts(infile);

% Open the file
ncid=netcdf.open(infile,'NC_NOWRITE');

% Get the variables
[numdims numvars numatts unlimdimID]=netcdf.inq(ncid);
for i=1:numvars
    [varnames{i},vartype{i},vardimid{i},varnatts{i}]=netcdf.inqVar(ncid,i-1);
end

% Get the ID of the variable we're interested in
for i=1:length(varnames)
    if strcmpi(varnames{i},'tide_Eamp')
        varID=i;
    end
end

% Get all the other data
for j=1:length(varnames)
    % Skip the tidal amplitude values
    if strcmpi(varnames{j},'tide_Eamp')~=1
        allData{j}=netcdf.getVar(ncid,netcdf.inqVarID(ncid,varnames{j}));
        varNamesOut{j}=varnames{j};
    end
end

% Get the current amplitudes
ncamp=netcdf.getVar(ncid,netcdf.inqVarID(ncid,varnames{varID}));

% Scale by some factor and put back in the allData cell array
tideFact=0.75;
allData{varID}=ncamp*tideFact;

netcdf.close(ncid)

% Create a new file
ncid=netcdf.create([A,'/',B,'_scaled',C],'NC_WRITE');
dimid=netcdf.defDim(ncid,'nobc',87);
my_varID=netcdf.defVar(ncid,varnames{1},'double',dimid);
netcdf.endDef(ncid)

netcdf.putVar(ncid,my_varID,allData{3});
my_varID=netcdf.defVar(ncid,varnames{3},'double',dimid);
netcdf.endDef(ncid)
netcdf.putVar(ncid,my_varID,allData{3});

% dimid=netcdf.defDim(ncid,'tidal_components',2);
% my_varID=netcdf.defVar(ncid,varnames{2},'double',dimid);
% netcdf.endDef(ncid)
% netcdf.putVar(ncid,my_varID,allData{2});


netcdf.close(ncid)


%%
% Copy the input file and use as the basis for the new file
copyfile(infile,[A,'/',B,'_scaled',C])
ncout=netcdf.open([A,'/',B,'_scaled',C],'NC_WRITE');

% Output the new variable to a new file 
dimid0=netcdf.defDim(ncout,'nobc',max(size(allData{1})));
% Do the 87s
my_varID=netcdf.getVar(ncout,0,'double');
netcdf.putVar(ncout,my_varID,allData{1});
netcdf.endDef(ncout)
my_varID2=netcdf.defVar(ncout,varnames{3},'double',dimid0);
netcdf.putVar(ncout,my_varID,allData{3});
netcdf.endDef(ncout)


netcdf.close(ncout)

%%
for k=1:length(varnames)
    % Check which variable we're using and adjust if we're at the amplitude
    dimid=netcdf.defDim(ncout,'my_dim',max(size(allData{k})));
    my_varID=netcdf.defVar(ncout,varnames{k},'double',dimid);
    netcdf.endDef(ncout)
    if strcmpi(varnames{k},'tide_Eamp')~=1
        netcdf.putVar(ncout,my_varID,allData{k});
        [varname2,xtype2,dimid2,natts2]=netcdf.inqVar(ncout,k-1)
    else
        netcdf.putVar(ncout,my_varID,ncampScale);
        [varname2,xtype2,dimid2,natts2]=netcdf.inqVar(ncout,k-1)
    end
    varnames{k}
    vartype{k}
    vardimid{k}
    varnatts{k}
end

netcdf.close(ncid)
netcdf.close(ncout)