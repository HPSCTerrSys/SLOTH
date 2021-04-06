import numpy as np
import heat as ht
import matplotlib as mpl
import netCDF4 as nc
import cftime
import sys
import glob
import copy

src_path='../src'
sys.path.append(src_path)
import SanityCheck
import ParFlow_IO as pio
import VanG
import pars_ParFlowTCL as ppfl

def mappIndicator(indiRoot, alphaD=2.0, nD=3.0, SresD=0.1):
    Indi3D = pio.read_pfb(f'{indiRoot}/DE-0055_INDICATOR_regridded_rescaled_SoilGrids250-v2017_BGR3_allv.pfb')
    shape3D= Indi3D.shape
    print(f'shape3D: {shape3D}')
    alpha  = np.full(shape3D,alphaD)
    nvg    = np.full(shape3D,nD)
    sres   = np.full(shape3D,SresD)

    IndicatorInput = [     1,     2,     3,     4,     5,     6,     7,     8,     9,    10,    11,    12,    13,    14,    15,    16,    17,    18,    19,   20,   21]
    GeomInput      = ['TC01','TC02','TC03','TC04','TC05','TC06','TC07','TC08','TC09','TC10','TC11','TC12','BGR1','BGR2','BGR3','BGR4','BGR5','BGR6','Allv','Lake','Sea']
    Alpha          = [ 3.523709, 3.475362, 2.666859, 1.111732, 0.505825, 0.657658, 2.108628, 1.581248,  0.83946,  3.34195,  1.62181, 1.496236, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
    Nvg            = [ 3.176874, 2.745822, 2.448772, 2.472313, 2.663413, 2.678804, 2.330454, 2.415794, 2.520548, 2.207814, 2.321296, 2.253141, 4.,   4.,  4.,  4.,  4.,  4.,  4.,  4.,  4.]
    # Nvg            = [ 4.176874, 2.745822, 2.448772, 2.472313, 2.663413, 2.678804, 2.330454, 2.415794, 2.520548, 2.207814, 2.321296, 2.253141, 4.,   4.,  4.,  4.,  4.,  4.,  4.,  4.,  4.]
    Sres           = [ 0.141333, 0.125641, 0.100775, 0.152882, 0.148064, 0.102249, 0.164063, 0.178733, 0.186722, 0.303896, 0.230769, 0.213508, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]

    for k in range(len(IndicatorInput)):
        alpha    = np.where(Indi3D == IndicatorInput[k], Alpha[k], alpha)
        nvg      = np.where(Indi3D == IndicatorInput[k], Nvg[k],   nvg)
        sres     = np.where(Indi3D == IndicatorInput[k], Sres[k],  sres)

    return alpha, nvg, sres

def mappIndicator(ParFlowNamelist, IndicatorFile):
    ###############################################################################
    ### parse ParFlow-Namelist to get indicators values
    ###############################################################################
    PFNamelistDict = ppfl.pars_ParFlowNamelist(f'{ParFlowNamelist}')
    
    nz = int(PFNamelistDict['ComputationalGrid.NZ'])
    tcl_dz_keys = [f'Cell.{i}.dzScale.Value' for i in range(nz)]
    dz_mult = [float(PFNamelistDict[tcl_dz_key]) for tcl_dz_key in tcl_dz_keys]
    
    GeomInputs = PFNamelistDict['GeomInput.indinput.GeomNames']
    
    IndicatorInput = {GeomInput:int(PFNamelistDict[f'GeomInput.{GeomInput}.Value']) for GeomInput in GeomInputs}
    # print(f'IndicatorInput: {IndicatorInput}')
    vanGA_GeomNames = PFNamelistDict['Phase.Saturation.GeomNames']
    # print(f'vanGA_GeomNames: {vanGA_GeomNames}')
    tcl_vanGA_keys = [ ]
    vanG_a = {vanGA_GeomName:float(PFNamelistDict[f'Geom.{vanGA_GeomName}.Saturation.Alpha']) for vanGA_GeomName in vanGA_GeomNames}
    # print(f'vanG_a: {vanG_a}')
    vanG_n = {vanGA_GeomName:float(PFNamelistDict[f'Geom.{vanGA_GeomName}.Saturation.N']) for vanGA_GeomName in vanGA_GeomNames}
    vanG_sres = {vanGA_GeomName:float(PFNamelistDict[f'Geom.{vanGA_GeomName}.Saturation.SRes']) for vanGA_GeomName in vanGA_GeomNames}

    ###############################################################################
    ### prepare arrays to return
    ###############################################################################
    Indi3D = pio.read_pfb(f'{IndicatorFile}')
    shape3D= Indi3D.shape
    print(f'shape3D: {shape3D}')
    alpha  = np.full(shape3D,vanG_a['domain'])
    nvg    = np.full(shape3D,vanG_n['domain'])
    sres   = np.full(shape3D,vanG_sres['domain'])

    ###############################################################################
    ### mapp indicators to VanGenuchten parameter according to ParFlow-Namelist
    ###############################################################################
    for GeomName in vanGA_GeomNames:
        if GeomName == 'domain': 
            continue
        alpha    = np.where(Indi3D == IndicatorInput[GeomName], vanG_a[GeomName], alpha)
        nvg      = np.where(Indi3D == IndicatorInput[GeomName], vanG_n[GeomName],   nvg)
        sres     = np.where(Indi3D == IndicatorInput[GeomName], vanG_sres[GeomName],  sres)

    return alpha, nvg, sres

simresDir = f'/p/scratch/cesmtst/tsmpforecast/DE05Clima_DE05_FZJ-IBG3-mgrowa_clima_FZJ-IBG3-ParFlow/simres'
ParFlowNamelist = '/p/scratch/cjibg35/tsmpforecast/DE05Clima_DE05_FZJ-IBG3-mgrowa_clima_FZJ-IBG3-ParFlow/ctrl/ParFlowCLM_DE05.tcl'
IndicatorFile   = '/p/scratch/cjibg35/tsmpforecast/DE05Clima_DE05_FZJ-IBG3-mgrowa_clima_FZJ-IBG3-ParFlow/run/DE-0055_INDICATOR_regridded_rescaled_SoilGrids250-v2017_BGR3_allv.pfb'

# get a sorted list of all dirs within rootDir
subDirs = []
for year in range(1965,1971):
    subDirs += sorted(glob.glob(f'{simresDir}/{year}_*/'))
#subDirBlacklist = ['1960','1961', '1962'. '1963', '1964']
#subDirs = [item for item in subDirs if subDirBlacklist not in item]
#1961_01/modeloutput

lvl = -1
varName = 'pressure'

for subDir in subDirs:
    ###############################################################################
    ### Read in individual files and merge (choose lvl) to one ndarray
    ###############################################################################
    print(f'subDir: {subDir}')
    fileNames = sorted(glob.glob(f'{subDir}/modeloutput/ParFlowCLM_DE05.out.?????.nc'))
    VarRaw = None
    for fileName in fileNames:
        print(f'handling: {fileName}')
        tmpMask = pio.read_pfb(f'{subDir}/modeloutput/ParFlowCLM_DE05.out.mask.pfb')[lvl,...]
        tmpMask = np.where(tmpMask>0, 1, 0)
        print(f'tmpMask.shape: {tmpMask.shape}')
        with nc.Dataset(f'{fileName}', 'r') as tmp_nc:
            # ParFlow netCDF output is in (t, z, y, x= always!
            masked_tmpVar = tmp_nc.variables[varName][:,lvl,...]
            tmpVar = np.ma.getdata(masked_tmpVar)

            if VarRaw is None:
                VarRaw = tmpVar
                Var_mask = tmpMask[np.newaxis]
            else:
                VarRaw = np.append(VarRaw, tmpVar, axis=0)
                Var_mask = np.append(Var_mask, tmpMask[np.newaxis], axis=0)


    ###############################################################################
    ### Do some calculations e.g. press--> satur
    ###############################################################################
    # alpha, nvg, sres = mappIndicator(f'/p/scratch/cjibg35/tsmpforecast/DE05Clima_DE05_FZJ-IBG3-mgrowa_clima_FZJ-IBG3-ParFlow/run')
    alpha, nvg, sres = mappIndicator(ParFlowNamelist=ParFlowNamelist, IndicatorFile=IndicatorFile)
    # alpha, nvg, sres = mappIndicator(f'{rootDir}/run')
    satur = VanG.vanGenuchten(refP=VarRaw, sSat=1.0, sRes=sres[lvl], nVanG=nvg[lvl], aVanG=alpha[lvl])

    ###############################################################################
    ### Start SanityCheck 3D
    ###############################################################################
    # define some title for plot, which can be passed via function arguments
    # see function definition for full potential
    varPlotName = 'satur'
    Var         = satur
    tmp_title_str = [
        f'Sanity-Check for {subDir.split("/")[-2]}',
        f'Var: {varPlotName}',
        ]
    fig_title    = '\n'.join(tmp_title_str)
    figname      = f'{simresDir}/{subDir.split("/")[-2]}_satur_SanityCheck.png'
    minax_title  = f'{varPlotName} min'
    maxax_title  = f'{varPlotName} max'
    kinax_title  = f'{varPlotName} mean'
    hisax_title  = f'{varPlotName} mean - distribution'

    # use 'plot_SanityCheck_3D' from script ../src/SanityCheck.py imported above
    SanityCheck.plot_SanityCheck_3D(data=Var, 
            # below is optional
            data_mask=Var_mask, kind='mean', figname=figname,
            lowerP=2, upperP=98, interactive=False,
            # below is even more optional (**kwargs)
            fig_title=fig_title, minax_title=minax_title, maxax_title=maxax_title, 
            kinax_title=kinax_title, hisax_title=hisax_title)
