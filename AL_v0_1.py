'''
author's:
@L-daria
@Desert_Fox_Fenek
@Michello077
'''

from itertools import zip_longest
import os
from os import remove
import pandas as pd
import threading
from scipy.stats import chi2_contingency
from scipy.stats import chisquare
import numpy as np

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

if __name__ == "__main__":

    print("Import danych...")

    import_data('NADIR_sick_genotypes.csv','NADIR_healthy_genotypes.csv')

    print("Zaimportowano")

    print(data.head())

    #do poprawy czasowo
    for item in data:
        data = data.replace({item:{'2/2': np.NAN, '0/2' : np.NAN,'1/2' : np.NAN}}).dropna()

    print("Zmieniono dane")
    
    data = data.groupby(by = 'Chr1')

    print("Generowaie wynikow...")

    wynik = pd.DataFrame()
    lista = ['0/0','0/1','1/1']
    for l in data:
        print(l)
        y = pd.DataFrame(l[1])
        for row in range(0,len(y)): 
            chore = pd.DataFrame(y.iloc[row].iloc[3:19].value_counts()).reset_index()
            zdrowe = pd.DataFrame(y.iloc[row].iloc[19:].value_counts()).reset_index()
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
            
            print(zdrowe) # to zle
            print(chore) #liczy dobrze

            zdrowe = zdrowe.sort_values(by = 'genotyp')
            chore = chore.sort_values(by = 'genotyp')
            zdrowe = zdrowe['ilosc'].reset_index(drop = True)
            chore = chore['ilosc'].reset_index(drop = True)

            test1 = [[zdrowe[0] + zdrowe[1], zdrowe[2]], [chore[0]+ chore[1], chore[2]]]
            test2 = [[zdrowe[0] + zdrowe[2], zdrowe[1]], [chore[0]+ chore[2], chore[1]]]
            test3 = [[zdrowe[1] + zdrowe[2], zdrowe[0]], [chore[1]+ chore[2], chore[0]]]
        
            if (test1[0][1] == 0 and test1[1][1] == 0) or (test1[1][0] == 0 and test1[0][0] == 0):
                pass
            else:
                p1 = chi2_contingency(test1)[1]
            
            if (test2[0][1] == 0 and test2[1][1] == 0) or (test2[1][0] == 0 and test2[0][0] == 0):
                pass
            else:
                p2 = chi2_contingency(test2)[1]            
            
            if (test3[0][1] == 0 and test3[1][1] == 0) or (test3[1][0] == 0 and test3[0][0] == 0):
                pass
            else:
                p3 = chi2_contingency(test3)[1]            
            

            if p1 < 0.05 or p2 < 0.05 or p3 < 0.05:
            
                wynik = pd.concat([wynik,pd.DataFrame([l[0]],[1])])
    print(wynik.value_counts())