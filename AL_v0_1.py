'''
author's:
@L-daria
@Desert_Fox_Fenek
@Michello077
'''

from itertools import zip_longest
from multiprocessing import process
import os
from os import cpu_count, remove
import pandas as pd
import threading
from scipy.stats import chi2_contingency
from scipy.stats import chisquare
import numpy as np
import multiprocessing
from time import perf_counter


#import danych
def import_data(name_s,name_h):

    global data
    data = pd.DataFrame()
    ch_size = 750000
    batch_num = 1

    #dzielimy nasz duzy plik na mniejsze fragmenty
    for chunk_s,chunk_h in zip(pd.read_csv(name_s,chunksize = ch_size,low_memory=False),pd.read_csv(name_h,chunksize = ch_size,low_memory=False)):
        chunk_s.to_csv('chunk_s'+str(batch_num)+'.csv', index = False)
        chunk_h.to_csv('chunk_h'+str(batch_num)+'.csv', index = False)

        df_s = pd.read_csv('chunk_s'+str(batch_num)+'.csv', low_memory=False)
        df_h = pd.read_csv('chunk_h'+str(batch_num)+'.csv', low_memory=False)

        df_h.rename(columns = {'X238':'X182'}, inplace = True)
        df_s.drop('Chr1',axis = 1, inplace = True)

        df_h = pd.merge(df_h, df_s, how='inner', on = 'X182')
        data = pd.concat([data,df_h])

        del df_h
        del df_s

        remove('chunk_s'+str(batch_num)+'.csv')
        remove('chunk_h'+str(batch_num)+'.csv')

        batch_num += 1

def clear():
    _ = os.call('clear' if os.name =='posix' else 'cls')

def chunk_al(name):
    wynik = pd.DataFrame()

    df_p = pd.read_csv(name, low_memory=False)

    df_p = df_p.groupby(by = 'Chr1')
    lista = ['0/0','0/1','1/1']
    for l in df_p:
        y = pd.DataFrame(l[1])
        for row in range(0,len(y)): 
            chore = pd.DataFrame(y.iloc[row].iloc[4:19].value_counts()).reset_index()
            zdrowe = pd.DataFrame(y.iloc[row].iloc[21:].value_counts()).reset_index()
            zdrowe.columns = ['genotyp', 'ilosc']
            chore.columns = ['genotyp', 'ilosc'] 
        
            if len(zdrowe) != 3:
                for i in lista:
                    n = len(zdrowe.loc[zdrowe['genotyp'] == i])
                    if n == 0:
                        data = [[i,0]]
                        df = pd.DataFrame(data,columns=['genotyp','ilosc'])
                        zdrowe = pd.concat([zdrowe,df])
                    
            if len(chore) != 3:
                for i in lista:
                    n = len(chore.loc[chore['genotyp'] == i])
                    if n == 0:
                        data = [[i,0]]
                        df = pd.DataFrame(data,columns=['genotyp','ilosc'])
                        chore = pd.concat([chore,df])

            zdrowe = zdrowe.sort_values(by = 'genotyp')
            chore = chore.sort_values(by = 'genotyp')
            zdrowe = zdrowe['ilosc'].reset_index(drop = True)
            chore = chore['ilosc'].reset_index(drop = True)

            test1 = [[zdrowe[0] + zdrowe[1], zdrowe[2]], [chore[0]+ chore[1], chore[2]]]
            test2 = [[zdrowe[0] + zdrowe[2], zdrowe[1]], [chore[0]+ chore[2], chore[1]]]
            test3 = [[zdrowe[1] + zdrowe[2], zdrowe[0]], [chore[1]+ chore[2], chore[0]]]
        
            if (test1[0][1] == 0 and test1[1][1] == 0) or (test1[1][0] == 0 and test1[0][0] == 0):
                p1 = 1
                pass
            else:
                p1 = chi2_contingency(test1, correction=False)[1]
            
            if (test2[0][1] == 0 and test2[1][1] == 0) or (test2[1][0] == 0 and test2[0][0] == 0):
                pass
                p2 = 1
            else:
                p2 = chi2_contingency(test2, correction=False)[1]            
            
            if (test3[0][1] == 0 and test3[1][1] == 0) or (test3[1][0] == 0 and test3[0][0] == 0):
                pass
                p3 = 1
            else:
                p3 = chi2_contingency(test3, correction=False)[1]            
               

            if p1 < 0.05 or p2 < 0.05 or p3 < 0.05:
                wynik = pd.concat([wynik,pd.DataFrame([l[0]],[1])])
    return wynik.value_counts()

if __name__ == "__main__":
    wyniki = []
    ch_size = 500000
    batch_num = 1

    print('Rozdrobnienie pliku')

    for chunk in pd.read_csv('Obrobione_dane.csv', sep='\t',chunksize = ch_size,low_memory=False):
        chunk.to_csv('TempData'+str(batch_num)+'.csv', index = False)

        batch_num += 1

    in_data = []
    for i in range(1,batch_num):
        in_data.append('TempData'+str(i)+'.csv')

    pool = multiprocessing.Pool()
    pool = multiprocessing.Pool(processes = cpu_count() - 2)
    
    print('Rozpoczecie obliczen')
    start = perf_counter()
    out_data = pool.map(chunk_al,in_data)

    pool.close()
    stop = perf_counter()
    #for num in range(batch_num):
    #    name = 'TempData'+str(num + 1)+'.csv'

    #    wyniki.append(chunk_al(name))

    print('W czasie: ',stop-start)
    print('Wynik:')
    print(out_data)

    f = open('Wynik.txt', 'w')
    f.writelines(str(out_data))