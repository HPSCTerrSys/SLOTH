import numpy as np
import argparse

def vanGenuchten(refP, sSat, sRes, nVanG, aVanG):
    """  Calculates the degree of saturation as a function of the pressure head.

    The degree of saturation is calculated as a function of the pressure head
    according to M. Th. van Genuchten.
    Name: A Closed‐form Equation for Predicting the Hydraulic Conductivity of 
            Unsaturated Soils
    DOI: https://doi.org/10.2136/sssaj1980.03615995004400050002x

    Parameters:
    refP:       Pressure head [L]
    sSat:       Relative saturated water content [-]
    sRes:       Relative residual saturation [-]
    nVanG:      Non linearity coefficient of the soil [-]
    aVanG:      Air entry values of the soil [L^-1]

    Returns:
    vanG:

    """
    mVanG = 1 - (1 / nVanG)

    vanG = ( (sSat - sRes) / ( 1 + (aVanG * np.absolute(refP))**(nVanG) )**mVanG ) + sRes

    vanG = np.where(refP>=0., 1., vanG) # avoid unsaturated values where water is ponding

    return vanG

def invers_vanGnuchten(vanG, sSat, sRes, nVanG, aVanG):
    mVanG = 1 / (1 - (1/nVanG))

    refP = ( ((sSat - sRes) / (vanG - sRes))**mVanG - 1 )**(1./nVanG) / aVanG

    return refP

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Tell me what this script can do!.')
    parser.add_argument('--refX', type=float,
                    help='refP or vanG')
    parser.add_argument('--sRes', type=float, default=0,
                    help='sRes')
    parser.add_argument('--sSat', type=float, default=1,
                    help='sSat')
    parser.add_argument('--aVanG', type=float,
                    help='aVanG')
    parser.add_argument('--nVanG', type=float,
                    help='nVanG')
    parser.add_argument('--invers', type=int, default=0,
                    help='invers-function')
    args = parser.parse_args()


    refX = args.refX
    sRes = args.sRes
    sSat = args.sSat
    aVanG = args.aVanG
    nVanG = args.nVanG
    invers = args.invers

    if not invers:
        vanG = vanGenuchten(refX, sSat, sRes, nVanG, aVanG)
        print(f'vanG: {vanG}')
    if invers:
        refP = invers_vanGnuchten(refX, sSat, sRes, nVanG, aVanG)
        print(f'refP: {refP}')
