from subprocess import call

srg=raw_input('srg cutoff: ')
highest_emax = int(raw_input('highest emax: '))
lowest_emax = int(raw_input('lowest emax: '))
emaxspace = int(raw_input('spacing of emax: ')) 
highest_hw = int(raw_input('highest hw: ')) 
lowest_hw = int(raw_input('lowest hw: '))
hwspace = int(raw_input('spacing of hw: '))

for hw in range(lowest_hw,highest_hw+1,hwspace):
    for emax in range(lowest_emax+emaxspace,highest_emax+1,emaxspace)[::-1]:

        infile = 'vsrg'+srg+'_n3lo500_w_coulomb_emax'+str(emax)
        infile += '_hw'+str(hw)+'.int' 
        
        outfile = 'vsrg'+srg+'_n3lo500_w_coulomb_emax'+str(emax-emaxspace)
        outfile += '_hw'+str(hw)+'.int' 

        num = str(emax-emaxspace)
        print srg,infile,outfile
        call(["./truncate",infile,outfile,num])

        
