import sys
import os.path
nuc = raw_input( 'Enter nucleus name: (He4,O16,Ca40) ' ) 

if (nuc == 'He4'): 
    nprot = 2
    nneut = 2
elif (nuc == 'O16'): 
    nprot = 8
    nneut = 8
elif (nuc == 'Ca40'): 
    nprot = 20
    nneut = 20
else:
    print 'Invalid Entry'
    sys.exit()
    
hwstr = raw_input( 'Enter hw values, seperated by commas: ( 20,22,24 ) ') 
hwlist = hwstr.strip().split(',') 

Rstr = raw_input('Enter eMax values, seperated by commas: (3,5,7,9) ' ) 
Rlist = Rstr.strip().split(',') 

lam = raw_input( 'Momentum cutoff in inverse fermi: ') 

hamtype = raw_input('For CM hamiltonian type: "1" for harmonic trap: "2", full: "3": ') 
hf = raw_input( 'For HF type: "HF". Otherwise type: "HO": ') 
mag = raw_input( 'For magnus type: "mag". Otherwise type: "trad": ' ) 
fq = open('run_all.bat','w')

fq.write( '#!/bin/bash \n\n') 

if hf == 'HF':
    HFint = '1'
else:
    HFint = '2'
    
if mag == 'mag':
    magint = '1'
else:
    magint = '2'
    
mem = ['500mb','1gb','2gb','3gb','4gb','5gb','6gb','7gb','8gb'] 
wtime = [ '00:20:00','00:40:00','01:00:00','02:00:00','05:00:00',  
'10:00:00','20:00:00','30:00:00','35:00:00']

for R in Rlist:
    for hw in hwlist:
        
        hwx = int(hw) 
        Rx = int(R) 
        
        memreq = mem[Rx - 3] 
        timreq = wtime[Rx-3]
        
        jobname = nuc+'_'+mag+'_vsrg'+lam+'_emax'+R+'_hw'+hw 
        spfile = 'nl'+R+'.sps'
        TBMEfile = 'vsrg'+lam+'_n3lo500_w_coulomb_emax'+R+'_hw'+hw+'.int' 
        initfile = nuc+'_'+mag+'_vsrg'+lam+'_emax'+R+'_hw'+hw+'.ini'

        # write pbs file ===========================        
        fx = open('pbs_'+jobname,'w') 
        
        fx.write('#!/bin/sh \n\n')
        fx.write('#PBS -l walltime='+timreq+'\n')
        fx.write('#PBS -l nodes=1:ppn=1\n')
        fx.write('#PBS -l mem='+memreq+'\n') 
        fx.write('#PBS -j oe\n')
        fx.write('#PBS -N '+jobname+'\n') 
        fx.write('#PBS -M parzuchowski@frib.msu.edu\n')
        fx.write('#PBS -m a\n\n')
        fx.write('cd $HOME/nuclear_IMSRG/src/im-SRG\n\n')
        fx.write('export OMP_NUM_THREADS=1\n\n')
        fx.write('time ./rum_IMSRG '+initfile+'\n')
        
        fx.close()
        
        fq.write('qsub pbs_'+jobname+'\n') 
        
        # write ini file ==========================
        
        fx = open('inifiles/'+initfile,'w') 
        
        fx.write('##########################\n')
        fx.write('#### IMSRG INPUT FILE ####\n')
        fx.write('#### KEEP THIS FORMAT ####\n')
        fx.write('##########################\n')
        fx.write('##########################\n')
        fx.write('# ENTER OUTPUT FILE PREFIX\n')
        fx.write(jobname +'\n') 
        fx.write('# ENTER INTERACTION FILE NAME\n')
        fx.write(TBMEfile + '\n') 
        fx.write('# ENTER SINGLE PARTICLE INPUT FILE NAME\n')
        fx.write(spfile + '\n')
        fx.write('# ENTER HAMILTONIAN TYPE\n')
        fx.write('# 1: T-V - Tcm -Vcm  2: harmonic trap T+U+V  3. T+V\n')
        fx.write(hamtype+'\n') 
        fx.write('# ENTER HO SPACING hw\n')
        fx.write(hw +'\n')
        fx.write('# ENTER NUMBER OF PROTONS\n')
        fx.write(str(nprot)+'\n')
        fx.write('# ENTER NUMBER OF NEUTRONS\n')
        fx.write(str(nneut)+'\n')
        fx.write('# ENTER 1 for HF basis\n')
        fx.write('# OR 2 for HO basis\n')
        fx.write(HFint+'\n')
        fx.write('# ENTER 1 for magnus method\n')
        fx.write('# or 2 for traditional ode\n') 
        fx.write(magint+'\n') 
        fx.write('########################################################\n')
        fx.write('# NOTES \n')
        fx.write('#\n')
        fx.write("# 1. THIS FILE'S NAME SHOULD BE SUPPLIED AS THE ONLY\n") 
        fx.write('# COMMAND ARGUMENT WITH THE EXECUTABLE "./run_IMSRG"\n')
        fx.write('#\n')
        fx.write('# 2. THIS FILE IS READ BY SUBROUTINE "read_main_input_file"\n')
        fx.write('# which is found in "basic_IMSRG.f90" \n')
        fx.write('# \n')
        fx.write('# 3. FILENAMES ARE READ TO COMMON BLOCK "files" Sorry. \n')
        fx.write('########################################################\n')

        fx.close()
        
        if not os.path.isfile('sp_inputs/'+spfile):
            print 'WARNING!!'
            print 'file sp_inputs/'+spfile+' not present!\n'
            
        if not os.path.isfile('TBME_input/'+TBMEfile):
            print 'WARNING!!'
            print 'file TBME_input/'+TBMEfile+' not present!\n'

    
