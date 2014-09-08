fl = open('fort.38','r')
fx = open('fort.75','r') 
fy = open('../sp_inputs/nl4.sps') 

jx = [] 

for line in fy:
    a = line.strip().split()
    jx.append(int(a[3]))

fy.close()

g_kosh = [] 
g_nat = []
for i in range(30):
    g_kosh.append([])
    g_nat.append([])
    for j in range(30):
        g_kosh[i].append([])
        g_nat[i].append([])
        for k in range(30):
            g_kosh[i][j].append([])
            g_nat[i][j].append([]) 
            for l in range(30):
                g_kosh[i][j][k].append(-1000.)
                g_nat[i][j][k].append(-1000.)
          
for line in fl:
    a = line.strip().split()
    i = int(a[0])
    j = int(a[1])
    k = int(a[2])
    l = int(a[3]) 
    g = float(a[4])
  #  print i,j,k,l
    g_kosh[i-1][j-1][k-1][l-1] = g
    
fl.close()

for line in fx:
    a = line.strip().split()
    i = int(a[0])
    j = int(a[1])
    k = int(a[2])
    l = int(a[3]) 
    g = float(a[4])
    
    g_nat[i-1][j-1][k-1][l-1] = g

fx.close()

for i in range(1,31):
    for j in range(1,31):
        for k in range(1,31): 
            for l in range(1,31):
              #  print i,j,k,l
                g1 = g_kosh[i-1][j-1][k-1][l-1]
                if (g1 < -999.9):
                    continue
                
                g2 = g_nat[i-1][j-1][k-1][l-1]
                g3 = g2
                if (g2 < -999.9): 
                    g2 = g_nat[i-1][j-1][l-1][k-1]
                    g3 = -g2* (-1.0)**((jx[l-1]+jx[k-1])/2)
                    if (g2 < -999.9):
                        g2 = g_nat[j-1][i-1][k-1][l-1]
                        g3 = -g2* (-1.0)**((jx[j-1]+jx[i-1])/2)
                        if (g2 < -999.9):
                            g2 = g_nat[j-1][i-1][l-1][k-1]
                            g2 = g2* (-1.0)**((jx[l-1]+jx[k-1])/2)
                            g3 = g2* (-1.0)**((jx[i-1]+jx[j-1])/2)
                            if (g2 < -999.9):
                                print 'fuck' 
                                
                if (abs(g1 - g3) > 1e-4):
                    print i,j,k,l,g1,g3
                
