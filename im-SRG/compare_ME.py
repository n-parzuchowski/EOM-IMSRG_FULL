fl = open('fort.45','r')   # kosh
fx = open('fort.75','r')   # nath
#fy = open('../sp_inputs/nl4.sps') 

jx = [] 

#for line in fy:
 #   a = line.strip().split()
  #  jx.append(int(a[3]))

#fy.close()

g_kosh = [0.]*810000
g_nat =  [0.]*810000

          
for line in fl:
    a = line.strip().split()
    i = min(int(a[0]),int(a[1])) - 1
    j = max(int(a[0]),int(a[1])) - 1
    k = min(int(a[2]),int(a[3])) - 1
    l = max(int(a[2]),int(a[3])) - 1
    g = float(a[4])
    
    g_kosh[l + 30 * k + 900 * j + 27000*i] = g 
    
fl.close()
#print sum(g_kosh)
for line in fx:
    a = line.strip().split()
    i = min(int(a[0]),int(a[1])) - 1
    j = max(int(a[0]),int(a[1])) - 1
    k = min(int(a[2]),int(a[3])) - 1
    l = max(int(a[2]),int(a[3])) - 1
    g = float(a[4])

    g_nat[l + 30 * k + 900 * j + 27000*i]= g

fx.close()
#print sum(g_kosh)
#print sum(g_nat)
#ass = raw_input('ass')
for i in range(30):
    for j in range(30):
        for k in range(30): 
            for l in range(30):
              #  print i,j,k,l
                g1 = max(g_kosh[l + 30 * k + 900 * j + 27000*i], \
                         g_kosh[j + 30 * i + 900 * l + 27000*k])
                g2 = max(g_nat[l + 30 * k + 900 * j + 27000*i], \
                    g_nat[j + 30 * i + 900 * l + 27000*k])
                #print g2
                if (abs(g1 - g2) > 1e-2):
                    print i+1,j+1,k+1,l+1,g1,g2
                    qdubs=raw_input('')
                
