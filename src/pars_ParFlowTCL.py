import numpy as np

def filter_blacklist(line):
    blacklist = ['#']
    for item in blacklist:
        if item in line:
            return False
    return True

def filter_whitelist(line):
    blacklist = ['pfset']
    for item in blacklist:
        if item in line:
            return True
    return False

def parse_TCLLine(line, prefix='pfset'):
    """ parses single lines from TCL-namelists and retuns key and value for dict

    This function parses single lines from TCL-namelists. Therefore it is
    assumed that each line starts with the a fix prefix (default: prefix=pfset)
    and that the variable name (eg. Solver.Nonlinear.VariableDz) does contain
    of a single string.

    returns tuple: (key, value) both with dtype=str
    key = the variable name (eg. Solver.Nonlinear.VariableDz)
    value = single value, empty string or list
    """
    # remove the prefix (usualy pfset)
    tmp = line.replace(prefix, '')
    # print(f'step1: {tmp}')
    # remove " and '
    # Does not care if there are no
    tmp = tmp.replace('"', '')
    tmp = tmp.replace("'", '')
    # print(f'step2: {tmp}')
    # remove leading and tailign blanks
    tmp = tmp.strip()
    # print(f'step3: {tmp}')
    # remove '\t' DO NEVER USE TAB IN PLAIN TEXT...
    tmp = tmp.replace('\t', '')
    # print(f'step3.5: {tmp}')
    # at this point the value should be left only
    # if this is a list parse it as list
    tmp = tmp.split(' ')
    # print(f'step4: {tmp}')
    # removing empty list entries which are from split()
    tmp = list(filter(None, tmp))
    # print(f'step5: {tmp}')

    # if value is not a list, do not return list
    if(len(tmp[1:])>1):
        #tmp[0] is allways the 'key' and the rest the 'value'
        return (tmp[0], tmp[1:])
    # if len(value) is 0 (happens if "" is in namelist) return ""
    elif(len(tmp[1:])==0):
        #tmp[0] is allways the 'key' and the rest the 'value'
        return (tmp[0], "")
    else:
        #tmp[0] is allways the 'key' and the rest the 'value'
        return (tmp[0], tmp[1])



file = './coup_oas.tcl'
with open(file,'r') as f:
    #return tcl content as list of lines
    tcl_content = f.read().splitlines()
    # do not keep empty lines
    # Note: filtered objects are iterables and not lists
    tcl_content = filter(None, tcl_content)
    # remove commented out lines
    tcl_content = filter(filter_blacklist, tcl_content)
    # keep lines started with 'pfset' only (constrain of parser function)
    tcl_content = filter(filter_whitelist, tcl_content)
    # convert to list
    tcl_content = list(tcl_content)

# parse eah line of TCL script and save as dict
PFNamelistDict = {}
PFNamelistDict.update({parse_TCLLine(line)[0]:parse_TCLLine(line)[1] for line in tcl_content})

# ------------------------------------------------------------------------------


# for key in PFNamelistDict:
#     print(key)
nz = int(PFNamelistDict['dzScale.nzListNumber'])
tcl_dz_keys = [f'Cell.{i}.dzScale.Value' for i in range(nz)]
dz_mult = [float(PFNamelistDict[tcl_dz_key]) for tcl_dz_key in tcl_dz_keys]
print('dz_mult')
print(dz_mult)

perm_patch_names = PFNamelistDict['Geom.Perm.Names']
tcl_perm_keys = [ ]
Perm = {perm_patch_name:float(PFNamelistDict[f'Geom.{perm_patch_name}.Perm.Value']) for perm_patch_name in perm_patch_names}
print('Perm')
print(Perm)

poro_patch_names = PFNamelistDict['Geom.Porosity.GeomNames']
tcl_poro_keys = [ ]
Poro = {poro_patch_name:float(PFNamelistDict[f'Geom.{poro_patch_name}.Porosity.Value']) for poro_patch_name in poro_patch_names}
print('Poro')
print(Poro)

vanGA_patch_names = PFNamelistDict['Phase.RelPerm.GeomNames']
tcl_vanGA_keys = [ ]
vanG_RelPerm_a = {vanGA_patch_name:float(PFNamelistDict[f'Geom.{vanGA_patch_name}.RelPerm.Alpha']) for vanGA_patch_name in vanGA_patch_names}
print('vanG_RelPerm_a')
print(vanG_RelPerm_a)

vanGn_patch_names = PFNamelistDict['Phase.RelPerm.GeomNames']
tcl_vanGn_keys = [ ]
vanG_RelPerm_n = {vanGn_patch_name:float(PFNamelistDict[f'Geom.{vanGn_patch_name}.RelPerm.N']) for vanGn_patch_name in vanGn_patch_names}
print('vanG_RelPerm_n')
print(vanG_RelPerm_n)

vanGA_patch_names = PFNamelistDict['Phase.Saturation.GeomNames']
tcl_vanGA_keys = [ ]
vanG_Saturation_a = {vanGA_patch_name:float(PFNamelistDict[f'Geom.{vanGA_patch_name}.Saturation.Alpha']) for vanGA_patch_name in vanGA_patch_names}
print('vanG_Saturation_a')
print(vanG_Saturation_a)

vanGn_patch_names = PFNamelistDict['Phase.Saturation.GeomNames']
tcl_vanGn_keys = [ ]
vanG_Saturation_n = {vanGn_patch_name:float(PFNamelistDict[f'Geom.{vanGn_patch_name}.Saturation.N']) for vanGn_patch_name in vanGn_patch_names}
print('vanG_Saturation_n')
print(vanG_Saturation_n)
