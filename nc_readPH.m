%function [V G]=nc_readPH(name);
%Reads all the data in netcdf file using built-in matlab netcdf routines
%
%Peter Huybers
%Harvard
%2010

function [V G]=nc_readPH(name);

ncid = netcdf.open(name,'nc_nowrite');
[ndims nvars natts] = netcdf.inq(ncid);
for ct=1:nvars,
  [V(ct).varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,ct-1);
  V(ct).varid = netcdf.inqVarID(ncid,V(ct).varname);
  for ctn=1:numatts,
    V(ct).attname{ctn} = netcdf.inqAttName(ncid,V(ct).varid,ctn-1); % Get value of attribute.
    V(ct).attval{ctn} = netcdf.getAtt(ncid,V(ct).varid,V(ct).attname{ctn}); % Get name of global attribute
  end;
  temp=netcdf.getConstant('NC_GLOBAL');
  if temp~=-1, 
    V(ct).gattname = netcdf.inqAttName(ncid,temp,0);
    V(ct).gattval = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),V(ct).gattname); % Get value of global attribute.  
  end;
  V(ct).data = netcdf.getVar(ncid,V(ct).varid);
end;

if nargout==2,
  for ct=1:natts,
    G(ct).name = netcdf.inqAttName(ncid,netcdf.getConstant('NC_GLOBAL'),ct-1);
    G(ct).type = netcdf.inqAtt(ncid,netcdf.getConstant('NC_GLOBAL'),G(ct).name);
    G(ct).gattval = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),G(ct).name); % Get value of global attribute.  
  end;
end;

netcdf.close(ncid);
