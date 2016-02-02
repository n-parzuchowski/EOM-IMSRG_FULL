import sys
import os.path
nuc = raw_input( 'Enter nucleus name: (He4,O16,Ca40,etc...) ' ) 

if (nuc == 'He4'): 
    nprot = 2
    nneut = 2
elif (nuc == 'O16'): 
    nprot = 8
    nneut = 8
elif (nuc == 'O22'): 
    nprot = 8
    nneut = 14
elif (nuc == 'O24'): 
    nprot = 8
    nneut = 16
elif (nuc == 'Ca40'): 
    nprot = 20
    nneut = 20
elif (nuc == 'Ca48'): 
    nprot = 20
    nneut = 28
elif ( nuc == 'C14'):
    nprot = 6
    nneut = 8
else:
    print 'Invalid Entry'
    sys.exit()
    
hwstr = raw_input( 'Enter hw values, seperated by commas: ( 20,22,24 ) ') 
hwlist = hwstr.strip().split(',') 

Rstr = raw_input('Enter eMax values, seperated by commas: (3,5,7,9) ' ) 
Rlist = Rstr.strip().split(',') 

lam = raw_input( 'Momentum cutoff in inverse fermi: ') 
lamlist = lam.strip().split(',')
me2jlam = 'x'
if (lam=='me2j'):
    print 'USING ME2J MATRIX ELEMENTS NOW'
    me2jlam = raw_input( 'Enter four digit srg cutoff 1/lamda^4: (0625) ')
    lamlist = me2jlam.strip().split(',')


threebodyfile = raw_input('enter a threebody file, or type "none": ')

if threebodyfile=="none": 
    E3Max = '0'
else:
    E3Max = raw_input('enter E3Max cutoff: ')

 
hamtype = raw_input('For intrinsic hamiltonian type: "1" for harmonic trap: "2", full: "3": ') 
hf = raw_input( 'For HF type: "HF". Otherwise type: "HO": ') 
mag = raw_input( 'For magnus type: "mag". Otherwise type: "trad" or "disc": ' ) 
tda = raw_input( 'select calculation: GS,EOM,TDA: ' )

leveltag = 'x'

if tda.lower() == 'tda':
    tdaint = '2'
#    Jtarg = raw_input( 'Input target total J: ')
#    Ptarg = raw_input( 'For even parity enter "0", for odd enter "1": ')
    
    Jlist = raw_input( 'enter list of states: "2+,3-,1-": ')
    Jpis = Jlist.strip().split(',')
    
    valshell = raw_input( 'Enter highest major shell in valence space: ')
    lawstring = raw_input('Lawson beta factors: ')
    lawbetas = lawstring.strip().split(',') 
    com_freq = raw_input('Center of Mass frequency: ') 
    
elif tda.lower() == 'eom':
    tdaint = '1'
    Jlist = raw_input( 'enter list of states: "2+,3-,1-": ')
    Jpis = Jlist.strip().split(',')
    valshell = '1s0d' #default
    lawstring = raw_input('Lawson beta factors: ')
    lawbetas = lawstring.strip().split(',') 
    com_freq = raw_input('Center of Mass frequency: ')     
else:
    tdaint = '0'
    Jpis = ['']
    valshell = '1s0d' #default
    lawbetas = ['0.0']
    com_freq = '0.0' 
    
hw_com = '0.0'

contin = raw_input('For other observable options, type "more", else press enter: ') 

transtype = '0'
transval = '0'

CMint = '0'
RRMSint = '0'
if contin.lower() ==  'more':
    print 'To calculate CM expectation value, and factorization frequency'
    CMint = raw_input('enter 1, otherwise enter 0: ') 
    if CMint != '1':
        print 'To calculate Point Nucleon RMS radius,' 
        RRMSint = raw_input('enter 1, otherwise enter 0: ') 
        
        if RRMSint != '1': 
            print 'To calculate transition matrix elements'
            transint = raw_input('enter 1, otherwise enter 0: ')
            if transint == '1':
                transtype=raw_input('enter transition type: (E,M,GT,F): ')
                transval=raw_input('enter transition rank: (0,1,2...): ')
                
                if transtype != 'E':
                    print 'unavailable...' 
                    transtype = '0'
                    transval = '0'

elif contin.lower() == 'com':
    print 'Including COM factorization calculation' 
    CMint = '1' 
elif contin.lower() == 'rsq':
    print 'Including RMS radius calculation' 
    RRMSint = '1'




    
fq = open('run_all.bat','w')

fq.write( '#!/bin/bash \n\n') 

if hf == 'HF':
    HFint = '1'
else:
    HFint = '2'
    
if mag == 'mag':
    magint = '1'
    quads = raw_input('Do you want to correct for quadrupoles? (Y/N) ')
    if (quads.lower() == 'y'):
        trips = raw_input('Do you want to correct for triples? (Y/N) ')
        if (trips.lower() == 'y'):
            magint='5'
            mag = mag+'TQ'
        else:
            magint='4'
            mag = mag+'Q'
    else:
        trips = raw_input('Do you want to correct for triples? (Y/N) ')
        if (trips.lower() == 'y'):
            magint='6'
            mag = mag+'T'

elif mag == 'disc':
    magint = '3'
else:
    magint = '2'
    
mem = ['500mb','1gb','2gb','3gb','4gb','5gb','6gb','9gb','18gb','20gb','25gb'] 
wtime = [ '00:20:00','00:40:00','01:00:00','02:00:00','03:00:00', \
'04:00:00','04:00:00','04:00:00','04:00:00','04:00:00','04:00:00']
ompnum = ['8','8','8','8','8','8','8','4','2','1','1']

for R in Rlist:
    for hw in hwlist:
        for lawbeta in lawbetas:
            for lam in lamlist:
                for Jpi in Jpis: 
                    
                    if Jpi == '':
                        Jtarg = '0'
                        Ptarg = '0' 
                        par = '+'
                        leveltag = 'x'
                    else:
                        Jtarg = Jpi[0]
                        if  Jpi[1]=='+':
                            Ptarg = '0'
                        else:
                            Ptarg = '1'
                        leveltag = Jpi
                        if (tda.lower() == 'tda'):
                            leveltag+='_TDA'
                            
                    Jtarg = str(2*int(Jtarg))
                    
                    hwx = int(hw) 
                    Rx = int(R) 

                    if magint=='5':
                        if Rx > 5: 
                            resubmit=True
                        else:
                            resubmit=False
                    else:
                        if Rx > 8: 
                            resubmit=True
                        else:
                            resubmit=False

                    memreq = mem[Rx - 3] 
                    timreq = wtime[Rx-3]

                    if (leveltag=='x'):
                        if (me2jlam != 'x'):
                            TBMEfile = 'chi2b_srg'+lam+'_eMax'+(2-len(R))*'0'+R+'_hwHO0'+hw+'.me2j.gz'
                            spfile = 'hk'+R+'.sps'
                            jobname = nuc+'_'+mag+'_chi2b_srg'+lam+'_eMax'+R+'_hw'+hw 
                            initfile = nuc+'_'+mag+'_chi2b_srg'+lam+'_eMax'+R+'_hw'+hw+'.ini'
                            prefix = jobname = nuc+'_'+mag+'_chi2b_srg'+lam+'_eMax'+R+'_hw'+hw 
                        else:    
                            TBMEfile = 'vsrg'+lam+'_n3lo500_w_coulomb_emax'+R+'_hw'+hw+'.int.gz' 
                            spfile = 'nl'+R+'.sps'
                            jobname = nuc+'_'+mag+'_vsrg'+lam+'_emax'+R+'_hw'+hw 
                            initfile = nuc+'_'+mag+'_vsrg'+lam+'_emax'+R+'_hw'+hw+'.ini'
                            prefix = jobname = nuc+'_'+mag+'_srg'+lam+'_eMax'+R+'_hw'+hw 
                    else:
                        if (me2jlam != 'x'):
                            TBMEfile = 'chi2b_srg'+lam+'_eMax'+(2-len(R))*'0'+R+'_hwHO0'+hw+'.me2j.gz'
                            spfile = 'hk'+R+'.sps'
                            jobname = nuc+'_'+mag+'_chi2b_srg'+lam+'_eMax'+R+'_hw'+hw+'_'+leveltag+'_law'+lawbeta
                            initfile = nuc+'_'+mag+'_chi2b_srg'+lam+'_eMax'+R+'_hw'+hw+'_'+leveltag+'_law'+lawbeta+'.ini'
                            prefix = jobname = nuc+'_'+mag+'_chi_2b_srg'+lam+'_eMax'+R+'_hw'+hw+'_'+leveltag 
                        else:    
                            TBMEfile = 'vsrg'+lam+'_n3lo500_w_coulomb_emax'+R+'_hw'+hw+'.int.gz' 
                            spfile = 'nl'+R+'.sps'
                            jobname = nuc+'_'+mag+'_vsrg'+lam+'_emax'+R+'_hw'+hw+'_'+leveltag+'_law'+lawbeta
                            initfile = nuc+'_'+mag+'_vsrg'+lam+'_emax'+R+'_hw'+hw+'_'+leveltag+'_law'+lawbeta+'.ini'
                            prefix = jobname = nuc+'_'+mag+'_srg'+lam+'_eMax'+R+'_hw'+hw +'_'+leveltag
                    # write pbs file ===========================        
                    fx = open('pbs_'+jobname,'w') 

                    fx.write('#!/bin/sh \n\n')
                    fx.write('#PBS -l walltime='+timreq+'\n')
                    fx.write('#PBS -l nodes=1:ppn='+ompnum[Rx-3]+'\n')
                    fx.write('#PBS -l mem='+memreq+'\n') 
                    fx.write('#PBS -j oe\n')
                    fx.write('#PBS -N '+jobname+'\n') 
                    fx.write('#PBS -M parzuchowski@frib.msu.edu\n')
                    fx.write('#PBS -m a\n\n')
                    fx.write('cd $HOME/nuclear_IMSRG/src/im-SRG\n\n')
                    fx.write('export OMP_NUM_THREADS='+ompnum[Rx-3]+'\n\n')
                    if resubmit: 
                        fx.write('export BLCR_CHECKFILE="/mnt/scratch/parzuch6/blcr_checkpoints/'+jobname+'.blcr"\n')
                        fx.write('export BLCR_OUTFILE="../../output/'+jobname+'.out"\n')
                        fx.write('export BLCR_WAIT_TIME=12600\n')
                        fx.write('export BLCR_PBS_FILE="pbs_'+jobname+'"\n')
                        fx.write('export PATH_TO_PBS_FILE="../../"\n\n')
                        fx.write('sh resubmit_qsub.sh ./run_IMSRG '+initfile+'\n')
                        print 'WARNING: CHECK "blcr_checkpoints" for .blcr file...'
                        print 'IF AN OLD ONE IS PRESENT, IT WILL BE RE-USED'
                        print
                    else:
                        fx.write('./run_IMSRG '+initfile+'\n')

                    fx.write('qstat -f ${PBS_JOBID}\nexit 0\n')
                    fx.close()

                    fq.write('qsub pbs_'+jobname+'\n') 

                    # write ini file ==========================


                    if (lam == 'me2j'):
                        if not os.path.isfile('sp_inputs/'+spfile):
                            print 'WARNING!!'
                            print 'file sp_inputs/'+spfile+' not present!\n'

                        if not os.path.isfile('/mnt/research/imsrg/nsuite/me/'+TBMEfile):
                            print 'WARNING!!'
                            print 'file TBME_input/'+TBMEfile+' not present!\n'


                    else:
                        if not os.path.isfile('sp_inputs/'+spfile):
                            print 'WARNING!!'
                            print 'file sp_inputs/'+spfile+' not present!\n'

                        if not os.path.isfile('TBME_input/'+TBMEfile):
                            if not os.path.isfile('TBME_input/'+TBMEfile[:-3]+'bin'):
                                print 'WARNING!!'
                                print 'file TBME_input/'+TBMEfile+' not present!\n'
                            else:
                                TBMEfile = TBMEfile[:-3]+'bin'


                    if (magint=='4'):
                        magout = "1,'y','n'"
                    elif (magint=='5'):
                        magout = "1,'y','y'"
                    elif (magint=='6'):
                        magout = "1,'y','y'"
                    else:
                        magout = magint+",'n','n'"

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
                    fx.write('# ENTER 3B INTERACTION FILE NAME\n')
                    fx.write(threebodyfile + '\n') 
                    fx.write('# ENTER E3Max (enter "0" for no three body force)\n')
                    fx.write(E3Max + '\n') 
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
                    fx.write(magout+'\n') 
                    fx.write('# 0: gs only, 1: EOM, 2: TDA\n')
                    fx.write('# ENTER 0 for ground state only\n') 
                    fx.write(tdaint+'\n')
                    fx.write('# ENTER total 2J of target states\n') 
                    fx.write(Jtarg + '\n')
                    fx.write('# ENTER total PARITY (0-even,1-odd) of target states\n') 
                    fx.write(Ptarg+'\n')
                    fx.write('# ENTER HIGHEST MAJOR VALENCE SHELL\n')
                    fx.write(valshell+'\n') 
                    fx.write('# ENTER 1 TO CALCULATE Hcm, 0 otherwise\n')
                    fx.write(CMint+'\n')
                    fx.write('# ENTER 1 TO CALCULATE Rrms, 0 otherwise\n') 
                    fx.write(RRMSint+'\n') 
                    fx.write('# Transition\n') 
                    fx.write("'"+transtype+"',"+transval+"\n") 
                    fx.write('# Lawson beta value\n') 
                    fx.write(lawbeta+','+com_freq+'\n') 
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
        
        
        
        
fq.close()

os.system("chmod 0755 run_all.bat")
