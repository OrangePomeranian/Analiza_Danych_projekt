'''
author's:
@L-daria
@Desert_Fox_Fenek
@Michello077
'''


import os
from os import remove
import pandas as pd
import threading
from scipy.stats import chi2_contingency
from scipy.stats import chisquare
import numpy as np

def import_data_s(name):
    global data_s
    
    data_s = pd.DataFrame()
    ch_size = 500000
    batch_num = 1

    #dzielimy nasz duzy plik na mniejsze fragmenty
    for chunk in pd.read_csv(name,chunksize = ch_size, low_memory=False):
        chunk.to_csv('chunk_s'+str(batch_num)+'.csv', index = False)

        df = pd.read_csv('chunk_s'+str(batch_num)+'.csv', low_memory=False)
        data_s = pd.concat([data_s,df])

        remove('chunk_s'+str(batch_num)+'.csv')

        batch_num += 1


#import danych ze zdrowymi osobnikami
def import_data_h(name):

    global data_h
    data_h = pd.DataFrame()
    ch_size = 500000
    batch_num = 1

    #dzielimy nasz duzy plik na mniejsze fragmenty
    for chunk in pd.read_csv(name,chunksize = ch_size,low_memory=False):
        chunk.to_csv('chunk_h'+str(batch_num)+'.csv', index = False)

        df = pd.read_csv('chunk_h'+str(batch_num)+'.csv', low_memory=False)
        data_h = pd.concat([data_h,df])

        remove('chunk_h'+str(batch_num)+'.csv')

        batch_num += 1

    data_h.rename(columns = {'X238':'X182'})
    data_h.rename(columns = {'Chr1':'ChrCH'})

#to do 'NADIR_sick_genotypes.csv'
def merg_import(name): 
    df = pd.DataFrame()
    ch_size = 500000
    batch_num = 1

    for chunk in pd.read_csv(name,chunksize = ch_size, low_memory = False):
        chunk.to_csv('chunk' + str(batch_num)+'.csv', index = False)

        df = pd.read_csv('chunk' + str(batch_num)+'.csv',low_memory=False)
        df.dropna('Chr1')

        data_h = pd.merge(data_h, df, how='inner', on=['X182'])

def clear():
    _ = os.call('clear' if os.name =='posix' else 'cls')

if __name__ == "__main__":

    print("Import danych...")

    load_data_th1 = threading.Thread(target=import_data_h, args=('NADIR_healthy_genotypes.csv',))
    load_data_th2 = threading.Thread(target=import_data_s, args=('NADIR_sick_genotypes.csv',))

    load_data_th1.start()
    load_data_th2.start()

    load_data_th1.join()
    load_data_th2.join()

    print("Zaimportowano")

    print(data_s.head())
    print(data_h.head())

    input()

    header_zdrowy = ['X1.1ZD', 'X1.1.1ZD', 'X1.1.2ZD', 'X0.1ZD', 'X0.1.1ZD', 'X1.1.3ZD', 'X0.1.2ZD', 'X1.1.4ZD', 'X1.1.5ZD', 'X1.1.6ZD', 'X1.1.7ZD', 'X1.1.8ZD', 'X1.1.9ZD', 'X1.1.10ZD', 'X1.1.11ZD', 'X1.1.12ZD']
    header_chory = ['X1.1CH', 'X1.1.1CH', 'X1.1.2CH', 'X0.1CH', 'X0.1.1CH', 'X1.1.3CH', 'X0.1.2CH', 'X1.1.4CH', 'X1.1.5CH', 'X1.1.6CH', 'X1.1.7CH', 'X1.1.8CH', 'X1.1.9CH', 'X1.1.10CH', 'X1.1.11CH', 'X1.1.12CH']
    header_razem = header_zdrowy + header_chory

    data_h.rename(columns = {'X238':'X182'})
    data_h.rename(columns = {'Chr1':'ChrCH'})
    

    #zbior_testowy = pd.merge(healthy, sick, how='inner', on=['X182']).drop(columns = ['ChrCH']).rename(columns = {'ChrZD':'Chr'})


    for item in header_razem:
        data_s = data_s.replace({item:{'2/2': np.NAN, '0/2' : np.NAN,'1/2' : np.NAN}}).dropna()
    
    data_s = data_s.groupby(by = 'Chr')

    wynik = pd.DataFrame()
    lista = ['0/0','0/1','1/1']
    for l in data_s:
        y = pd.DataFrame(l[1])
        for row in range(0,len(y)): 
            chore = pd.DataFrame(y.iloc[row].iloc[2:18].value_counts()).reset_index()
            zdrowe = pd.DataFrame(y.iloc[row].iloc[18:].value_counts()).reset_index()
            zdrowe.columns = ['genotyp', 'ilosc']
            chore.columns = ['genotyp', 'ilosc'] 
        
            if len(zdrowe) != 3:
                for i in lista:
                    n = len(zdrowe.loc[zdrowe['genotyp'] == i])
                    if n == 0:
                        data_s = [[i,0]]
                        df = pd.DataFrame(data_s,columns=['genotyp','ilosc'])
                        zdrowe = pd.concat([zdrowe,df])      
                    
            if len(chore) != 3:
                for i in lista:
                    n = len(chore.loc[chore['genotyp'] == i])
                    if n == 0:
                        data_s = [[i,0]]
                        df = pd.DataFrame(data_s,columns=['genotyp','ilosc'])
                        chore = pd.concat([chore,df])
        
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