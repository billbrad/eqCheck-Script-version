#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 10:05:44 2020

@author: williamstone
Analysis of a Membrane Equilibration
Calculates and plots RMSD of Equilibration
Calculates and plots energy potential of energy minimization
Plots changes temp for NVT,
and pressure for NPT

If specified as a CHARMM-GUI multistep Equilibration, will automatically
concatenate trajectory and energy files, if their prefix has not been changed
from the original CHARMM-GUI format.

"""

import sys
import argparse
import pexpect
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

###################### Reading Command Line Arguments #########################

parser = argparse.ArgumentParser(
    description='Generate Plots to visualize Equilibration data',
    prog='eqCheck')

#Required Args
parser.add_argument('-e', '--input1', type=str, required=False,
                    metavar='', help='prefix of em files')
parser.add_argument('-t', '--input2', type=str, required=False,
                    metavar='', help='prefix of NVT eq file')
parser.add_argument('-p', '--input3', type=str, required=False,
                    metavar='', help='prefix of NPT eq file')

# Bonus Args
group = parser.add_mutually_exclusive_group()
group.add_argument('-q', '--quiet', action='store_true',
                    help='Only display the exit message')
group.add_argument('-v', '--verbose', action='store_false',
                    help='Watch the words fly across the screen (default)')

# CHARMM-GUI or no
group2 = parser.add_mutually_exclusive_group()
group2.add_argument('-C', '--CHARMMGUI', action='store_true',
                    help='eq files have original CHARMM-GUI prefix' +
                    ' (em files can still have custom prefix)')
group2.add_argument('-N', '--notCHARMM', action='store_false',
                    help='files have custom prefix (default)')

args = parser.parse_args()

############################### FUNCTIONS ######################################
def call_gmx(inFile, outFile, keyword):
    # Open the dialog and write its output to a temp file
    pexpect.run('touch templog.txt')
    with open('templog.txt', 'w') as fout:
        gmx_check = pexpect.spawn(f'gmx energy -f {inFile}.edr -o {outFile}.xvg',
                            encoding='utf-8')
        gmx_check.logfile = fout
        gmx_check.expect('\r\n')
        gmx_check.send(f' 0\r')
        gmx_check.expect(pexpect.EOF)
    # Parse temp file for desired energy term
    num = None
    with open('templog.txt', 'r') as fout:
        for line in fout:
            split_line = line.split()
            if split_line != []:
                for word in split_line:
                    if word == keyword:
                        pos = split_line.index(word) - 1
                        num = split_line[pos]
    #rm from the temp file
    pexpect.run('rm templog.txt')

    if num != None:
        gmx = pexpect.spawn(f'gmx energy -f {inFile}.edr -o {outFile}.xvg',
                            encoding='utf-8')
        if args.verbose:
            gmx.logfile = sys.stdout
        gmx.expect('\r\n')
        gmx.send(f'{num} 0\r')
        gmx.expect(pexpect.EOF)
    else:
        print('ERROR: Could not parse energy terms')
        sys.exit()


def glue(bin, ext, output, cnt, last):
    cmd = f'gmx {bin} -f'
    while cnt <= last:
        cmd += f' step6.{cnt}_equilibration{ext}'
        cnt += 1
    cmd += f'  -o {output}{ext}'
    glue = pexpect.spawn(cmd)
    if args.verbose:
        glue.logfile = sys.stdout.buffer
    glue.expect(pexpect.EOF)


def read_xvg(filename, col1, col2):
    output = pd.read_csv(filename,
                 sep='\s+',  # columns are separated by spaces
                 skiprows=14,  # GROMACS xvgs normally have ~12 lines intro
                 comment='@',  # skip xvg plot stuff
                 header=None,  # don't get column names from file
                 names=[col1,  # x-axis column title
                        col2])  # assuming the file has 2 columns, else add more
    return output


def calc_RMSD(ref, sel, out):
    gmx_rms = pexpect.spawn(f'gmx rms -s {ref} -f {sel} -o {out}',
                                encoding = 'utf-8')
    gmx_rms.expect(' ')
    gmx_rms.sendline(f'7\r')
    gmx_rms.expect(' ')
    gmx_rms.sendline(f'7\r')
    if args.verbose:
        gmx_rms.logfile = sys.stdout
    gmx_rms.expect(pexpect.EOF)

    rmsd_data = read_xvg(out, 'Time (ps)', 'RMSD (nm)')
    return rmsd_data

########################### Sorting Arguments ##################################
if args.CHARMMGUI == True:
    if args.input1 == None and args.input2 == None:
        emPrefix = 'step6.0_minimization'
        NVTprefix = 'NVT_eq'
        NPTprefix = 'NPT_eq'
    elif args.input1 != None and args.input2 == None and args.input3 == None:
        emPrefix = args.input1
        NVTprefix = 'NVT_eq'
        NPTprefix = 'NPT_eq'
    else:
        print('ERROR: Please do not specify eq input prefix if using -C option')
        sys.exit()
    # concatenate NVT trajectory files
    glue('trjcat', '.trr', 'NVT_eq', 1, 2)
    # concatenate NPT trajectory files
    glue('trjcat', '.trr', 'NPT_eq', 3, 6)
    # concatenate NVT energy files (6.3 currently missing)
    glue('eneconv', '.edr', 'NVT_eq', 1, 2)
    # concatenate NVT energy files
    glue('eneconv', '.edr', 'NPT_eq', 3, 6)
else:
    emPrefix = args.input1
    NVTprefix = args.input2
    NPTprefix = args.input3

########################## ENERGY POTENTIAL ####################################
if (args.CHARMMGUI == False and args.input1 != None) or args.CHARMMGUI ==True:
    call_gmx(emPrefix, 'potential', 'Potential')

    em = pd.DataFrame(read_xvg('potential.xvg', 'Steps', 'Potential'),
                  columns=['Steps', 'Potential'])
    em = em[em['Potential'] < 0] # filter out positive numbers

    ax = em.plot(x='Steps', y='Potential', kind='line')
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

    plt.title('Potential During Energy Minimization')

    plt.savefig(f'potential.png')
else:
    print('ERROR: Energy minimization .edr file must be specified or -C flag must be used.')

############################## RMSD ############################################
if (args.CHARMMGUI == False and args.input1 != None) or args.CHARMMGUI ==True:

    nvt_data = calc_RMSD('step6.1_equilibration.tpr', 'NVT_eq.trr', 'RMSD_NVT.xvg')
    npt_data = calc_RMSD('step6.3_equilibration.tpr', 'NPT_eq.trr', 'RMSD_NPT.xvg')

    bx = nvt_data.plot(x='Time (ps)', y='RMSD (nm)', kind='line')
    bx.set_ylabel('RMSD (nm')
    plt.title(f'RMSD of NVT Relative to First Frame')
    plt.savefig(f'RMSD_NVT.png')

    cx = npt_data.plot(x='Time (ps)', y='RMSD (nm)', kind='line')
    cx.set_ylabel('RMSD (nm')
    plt.title(f'RMSD of NPT Relative to First Frame')
    plt.savefig(f'RMSD_NPT.png')

############################# BOX SIZE ########################################
# for step in steps:
#     call_gmx('energy', f'step6.{step}_equilibration', '.edr', f'box_{step}', 12)
#     if step == 1:
#         box = pd.DataFrame(mda.auxiliary.XVG.XVGReader(f'box_{step}.xvg')._auxdata_values,
#                    columns = ['Frame', 'x', 'y', 'z'])
#     else:
#         box_add = pd.DataFrame(mda.auxiliary.XVG.XVGReader(f'box_{step}.xvg')._auxdata_values,
#                         columns = ['Frame', 'x', 'y', 'z'])
#         box = pd.concat([box, box_add], axis=1, sort=False)
#
# print(box)
# #    cb = box.plot(x = 'Frame', y = ['x', 'y', 'z'])

############################# TEMPERATURE ########################################

    call_gmx(NVTprefix, 'temperature', 'Temperature')


    temp = pd.DataFrame(read_xvg('temperature.xvg', 'Frame', 'Kelvin'),
                   columns = ['Frame', 'Kelvin'])

    dx = temp.plot('Frame', 'Kelvin', kind='line')

    plt.title('Temperature of NVT')
    plt.ylabel('Temperature (K)')
    plt.xlabel('Frame')
    plt.savefig('temperature.png')

############################# PRESSURE ########################################

    call_gmx(NPTprefix, 'pressure', 'Pressure')

    bar = pd.DataFrame(read_xvg('pressure.xvg', 'Frame', 'Bar'),
                            columns = ['Frame', 'Bar'])

    ex = bar.plot(x='Frame', y='Bar', kind='line')

    plt.title('Pressure of NPT')
    plt.ylabel('Pressure')
    plt.xlabel('Frame')

    plt.savefig('pressure.png')
else:
    print('ERROR: Equilibration files must be specified or -C flag must be used.')

############################# OPEN ########################################
stuff = pexpect.run('open potential.png RMSD_NVT.png temperature.png RMSD_NPT.png pressure.png ')
