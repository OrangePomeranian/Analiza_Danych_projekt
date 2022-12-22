#Projekt - analiza danych

#porownanie ilosci SNP na chromosomach,
#w ktorych obserwujemy roznice czestosci u osob zdrowych i osob chorych


#wczytanie danych
NADIR_healthy_genotypes <- fread("~/Documents/Studia/Bioinformatyka/Semestr V/Analiza danych/Ćwiczenia/Projekt/NADIR_healthy_genotypes.vcf", sep="")
NADIR_healthy_genotypes
NADIR_sick_genotypes <- fread("~/Documents/Studia/Bioinformatyka/Semestr V/Analiza danych/Ćwiczenia/Projekt/NADIR_sick_genotypes.vcf", sep="")
NADIR_sick_genotypes


#załadowanie paczek
#install.packages("dplyr")
library(dplyr)
#install.packages('parallel')
library(parallel)
#install.packages('microbenchmark')
library(microbenchmark)

#==============================================================================

#analiza


#okreslenie badanych chromosomów
nazwy <- unique(NADIR_healthy_genotypes$Chr1)
chromosomy <- grep("^[Chr].*", nazwy, value=TRUE)
chromosomy


#zmiana nazw kolumn
names(NADIR_healthy_genotypes)[1] <- "Chr"
names(NADIR_healthy_genotypes)[2] <- "Genotype"
names(NADIR_sick_genotypes)[1] <- "Chr"
names(NADIR_sick_genotypes)[2] <- "Genotype"


#podzial tabeli na chromosomy i znalezienie wspolnych SNP
nr_chromosomu = 12
for (element in chromosomy) {
  chromosom_zdrowi <- NADIR_healthy_genotypes %>% filter(Chr == element)
  chromosom_chorzy <- NADIR_sick_genotypes %>% filter(Chr == element)
  wspolne_SNP <- inner_join(x = chromosom_zdrowi, y = chromosom_chorzy[2:18], by = "Genotype")
  assign(paste0("SNP_chromosom_", nr_chromosomu), wspolne_SNP %>% filter(wspolne_SNP[3:34] != "0/2", wspolne_SNP[3:34] != "1/2", wspolne_SNP[3:34] != "2/2", wspolne_SNP[3:34] != "0/3", wspolne_SNP[3:34] != "1/3", wspolne_SNP[3:34] != "2/3", wspolne_SNP[3:34] != "3/3"))
  nr_chromosomu <- nr_chromosomu + 1
}


#zliczanie poszczegolnych genotypow w wierszach
#wykonywanie testow statystycznych
#zliczanie wszystkich istotnych mutacji

zbior_testowy <- SNP_chromosom_12[369910:369915,]
zbior_testowy

mutacje <- 0
start_time <- Sys.time()
for (wiersz in 1:nrow(zbior_testowy)) {
  zdrowi_00 <- rowSums(zbior_testowy[wiersz, 3:18] == "0/0")
  zdrowi_01 <- rowSums(zbior_testowy[wiersz, 3:18] == "0/1")
  zdrowi_11 <- rowSums(zbior_testowy[wiersz, 3:18] == "1/1")
  chorzy_00 <- rowSums(zbior_testowy[wiersz, 19:34] == "0/0")
  chorzy_01 <- rowSums(zbior_testowy[wiersz, 19:34] == "0/1")
  chorzy_11 <- rowSums(zbior_testowy[wiersz, 19:34] == "1/1")
  
  tabela_0001 <- rbind(c(zdrowi_00 + zdrowi_01, zdrowi_11), c(chorzy_00 + chorzy_01, chorzy_11))
  test_0001 <- chisq.test(tabela_0001, correct = FALSE)
  
  tabela_0011 <- rbind(c(zdrowi_00 + zdrowi_11, zdrowi_01), c(chorzy_00 + chorzy_11, chorzy_01))
  test_0011 <- chisq.test(tabela_0011, correct = FALSE)
  
  tabela_0111 <- rbind(c(zdrowi_01 + zdrowi_11, zdrowi_00), c(chorzy_01 + chorzy_11, chorzy_00))
  test_0111 <- chisq.test(tabela_0111, correct = FALSE)
  
  wartosci_p <- c(test_0001$p.value, test_0011$p.value,test_0111$p.value)
  wartosci_p <- ifelse(is.nan(wartosci_p),1,wartosci_p)
  #print(wartosci_p)
  
  if (wartosci_p[1] < 0.05 | wartosci_p[2] < 0.05 | wartosci_p[3] < 0.05) {
    mutacje <- mutacje + 1
  }
  #print(mutacje)
}
end_time <- Sys.time()
mutacje
czas = end_time - start_time
czas


rownoleglenie <- function(zbior_testowy) {
  mutacje <- 0
  for (wiersz in 1:nrow(zbior_testowy)) {
    zdrowi_00 <- rowSums(zbior_testowy[wiersz, 3:18] == "0/0")
    zdrowi_01 <- rowSums(zbior_testowy[wiersz, 3:18] == "0/1")
    zdrowi_11 <- rowSums(zbior_testowy[wiersz, 3:18] == "1/1")
    chorzy_00 <- rowSums(zbior_testowy[wiersz, 19:34] == "0/0")
    chorzy_01 <- rowSums(zbior_testowy[wiersz, 19:34] == "0/1")
    chorzy_11 <- rowSums(zbior_testowy[wiersz, 19:34] == "1/1")
    
    tabela_0001 <- rbind(c(zdrowi_00 + zdrowi_01, zdrowi_11), c(chorzy_00 + chorzy_01, chorzy_11))
    test_0001 <- chisq.test(tabela_0001, correct = FALSE)
    
    tabela_0011 <- rbind(c(zdrowi_00 + zdrowi_11, zdrowi_01), c(chorzy_00 + chorzy_11, chorzy_01))
    test_0011 <- chisq.test(tabela_0011, correct = FALSE)
    
    tabela_0111 <- rbind(c(zdrowi_01 + zdrowi_11, zdrowi_00), c(chorzy_01 + chorzy_11, chorzy_00))
    test_0111 <- chisq.test(tabela_0111, correct = FALSE)
    
    wartosci_p <- c(test_0001$p.value, test_0011$p.value,test_0111$p.value)
    wartosci_p <- ifelse(is.nan(wartosci_p),1,wartosci_p)
    
    if (wartosci_p[1] < 0.05 | wartosci_p[2] < 0.05 | wartosci_p[3] < 0.05) {
      mutacje <- mutacje + 1
    }
  }
  print(mutacje)
}

#==============================================================================

#porownanie czasow dwoch funkcji
obliczenie_czasu_1 <- function(dane) {
  mutacje <- 0
  for (wiersz in 1:nrow(dane)) {
    zdrowi_00 <- rowSums(wiersz, dane[3:18] == "0/0")
    zdrowi_01 <- rowSums(wiersz, dane[3:18] == "0/1")
    zdrowi_11 <- rowSums(wiersz, dane[3:18] == "1/1")
    chorzy_00 <- rowSums(wiersz, dane[19:34] == "0/0")
    chorzy_01 <- rowSums(wiersz, dane[19:34] == "0/1")
    chorzy_11 <- rowSums(wiersz, dane[19:34] == "1/1")
    
    tabela_0001 <- rbind(c(zdrowi_00 + zdrowi_01, zdrowi_11), c(chorzy_00 + chorzy_01, chorzy_11))
    test_0001 <- chisq.test(tabela_0001, correct = FALSE)
    
    tabela_0011 <- rbind(c(zdrowi_00 + zdrowi_11, zdrowi_01), c(chorzy_00 + chorzy_11, chorzy_01))
    test_0011 <- chisq.test(tabela_0011, correct = FALSE)
    
    tabela_0111 <- rbind(c(zdrowi_01 + zdrowi_11, zdrowi_00), c(chorzy_01 + chorzy_11, chorzy_00))
    test_0111 <- chisq.test(tabela_0111, correct = FALSE)
    
    wartosci_p <- c(test_0001$p.value, test_0011$p.value,test_0111$p.value)
    wartosci_p <- ifelse(is.nan(wartosci_p),1,wartosci_p)
    
    if (wartosci_p[1] < 0.05 | wartosci_p[2] < 0.05 | wartosci_p[3] < 0.05) {
      mutacje <- mutacje + 1
    }
  }
  print(mutacje)
}

obliczenie_czasu_2 <- function(dane) {
  mutacje <- 0
  for (wiersz in 1:nrow(dane)) {
    zdrowi_00 <- rowSums(wiersz, dane[3:18] == "0/0")
    zdrowi_01 <- rowSums(wiersz, dane[3:18] == "0/1")
    zdrowi_11 <- rowSums(wiersz, dane[3:18] == "1/1")
    chorzy_00 <- rowSums(wiersz, dane[19:34] == "0/0")
    chorzy_01 <- rowSums(wiersz, dane[19:34] == "0/1")
    chorzy_11 <- rowSums(wiersz, dane[19:34] == "1/1")
    
    tabela_0001 <- rbind(c(zdrowi_00 + zdrowi_01, zdrowi_11), c(chorzy_00 + chorzy_01, chorzy_11))
    test_0001 <- chisq.test(tabela_0001, correct = FALSE)
    
    if (test_0001$p.value > 0.05 | is.na(test_0001$p.value)) {
      tabela_0011 <- rbind(c(zdrowi_00 + zdrowi_11, zdrowi_01), c(chorzy_00 + chorzy_11, chorzy_01))
      test_0011 <- chisq.test(tabela_0011, correct = FALSE)
    } else {
      mutacje <- mutacje + 1
    }
    
    if (test_0001$p.value > 0.05 | is.na(test_0001$p.value)) {
      tabela_0011 <- rbind(c(zdrowi_00 + zdrowi_11, zdrowi_01), c(chorzy_00 + chorzy_11, chorzy_01))
      test_0011 <- chisq.test(tabela_0011, correct = FALSE)
    } else {
      mutacje <- mutacje + 1
    }
    
    if (test_0011$p.value > 0.05 | is.na(test_0011$p.value)) {
      tabela_0111 <- rbind(c(zdrowi_01 + zdrowi_11, zdrowi_00), c(chorzy_01 + chorzy_11, chorzy_00))
      test_0111 <- chisq.test(tabela_0111, correct = FALSE)
    } else {
      mutacje <- mutacje + 1
    }
    
    if (test_0011$p.value > 0.05 | is.na(test_0011$p.value)) {
      NULL
    } else {
      mutacje <- mutacje + 1
    }
  }
  print(mutacje)
}

dane_testowe_1 <- SNP_chromosom_1[1:10,]
dane_testowe_2 <- SNP_chromosom_1[1:100,]
dane_testowe_3 <- SNP_chromosom_1[1:1000,]
dane_testowe_4 <- SNP_chromosom_1[1:10000,]

#microbenchmark(obliczenie_czasu_1(dane_testowe_1),
#               obliczenie_czasu_1(dane_testowe_2),
#               obliczenie_czasu_1(dane_testowe_3),
#               obliczenie_czasu_1(dane_testowe_4),
#               times = 10)

microbenchmark(obliczenie_czasu_1(dane_testowe_3),
               obliczenie_czasu_2(dane_testowe_3),
               times = 10)

#==============================================================================

#rownoleglizacja danych
rownoleglenie <- function(dane) {
  mutacje <- 0
  for (wiersz in 1:nrow(dane)) {
    zdrowi_00 <- rowSums(wiersz, dane[3:18] == "0/0")
    zdrowi_01 <- rowSums(wiersz, dane[3:18] == "0/1")
    zdrowi_11 <- rowSums(wiersz, dane[3:18] == "1/1")
    chorzy_00 <- rowSums(wiersz, dane[19:34] == "0/0")
    chorzy_01 <- rowSums(wiersz, dane[19:34] == "0/1")
    chorzy_11 <- rowSums(wiersz, dane[19:34] == "1/1")
    
    tabela_0001 <- rbind(c(zdrowi_00 + zdrowi_01, zdrowi_11), c(chorzy_00 + chorzy_01, chorzy_11))
    test_0001 <- chisq.test(tabela_0001, correct = FALSE)
    
    tabela_0011 <- rbind(c(zdrowi_00 + zdrowi_11, zdrowi_01), c(chorzy_00 + chorzy_11, chorzy_01))
    test_0011 <- chisq.test(tabela_0011, correct = FALSE)
    
    tabela_0111 <- rbind(c(zdrowi_01 + zdrowi_11, zdrowi_00), c(chorzy_01 + chorzy_11, chorzy_00))
    test_0111 <- chisq.test(tabela_0111, correct = FALSE)
    
    wartosci_p <- c(test_0001$p.value, test_0011$p.value,test_0111$p.value)
    wartosci_p <- ifelse(is.nan(wartosci_p),1,wartosci_p)
    
    if (wartosci_p[1] < 0.05 || wartosci_p[2] < 0.05 || wartosci_p[3] < 0.05) {
      mutacje <- mutacje + 1
    }
  }
  print(mutacje)
}

start_1 = Sys.time()
parSapply(cl, 1:00, function(n) rownoleglenie(dane_testowe_4))
koniec_1 = Sys.time()
stopCluster(cl)
czas_1 = koniec1 - start_1
czas_1

start_2 = Sys.time()
wynik <- rownoleglenie(dane_testowe_4)
koniec_2 = Sys.time()
czas_2 = koniec_2 - start_2
czas_2

#==============================================================================

#praca na wiekszej ilości klastrow

dane_testowe = SNP_chromosom_1[1,]
dlugosc = nrow(dane_testowe)

zdrowi_00 <- rowSums(dane_testowe[3:18] == "0/0")
zdrowi_01 <- rowSums(dane_testowe[3:18] == "0/1")
zdrowi_11 <- rowSums(dane_testowe[3:18] == "1/1")
chorzy_00 <- rowSums(dane_testowe[19:34] == "0/0")
chorzy_01 <- rowSums(dane_testowe[19:34] == "0/1")
chorzy_11 <- rowSums(dane_testowe[19:34] == "1/1")

rownoleglenie <- function(z00, z01, z11, c00, c01, c11) {
  mutacje <- 0
  for (element in seq(dlugosc)) {
    tabela_0001 <- rbind(c(z00[element] + z01[element], z11[element]), c(c00[element] + c01[element], c11[element]))
    test_0001 <- chisq.test(tabela_0001, correct = TRUE)
    
    tabela_0011 <- rbind(c(z00[element] + z11[element], z01[element]), c(c00[element] + c11[element], c01[element]))
    test_0011 <- chisq.test(tabela_0011, correct = TRUE)
    
    tabela_0111 <- rbind(c(z01[element] + z11[element], z00[element]), c(c01[element] + c11[element], c00[element]))
    test_0111 <- chisq.test(tabela_0111, correct = TRUE)
    
    wartosci_p <- c(test_0001$p.value, test_0011$p.value,test_0111$p.value)
    wartosci_p <- ifelse(is.nan(wartosci_p),1,wartosci_p)

    if (wartosci_p[1] < 0.05 | wartosci_p[2] < 0.05 | wartosci_p[3] < 0.05) {
      mutacje <- mutacje + 1
    }
  }
  print(mutacje)
}

start_1 <- Sys.time()
test_1 <- rownoleglenie(z00=zdrowi_00, z01=zdrowi_01, z11=zdrowi_11, c00=chorzy_00, c01=chorzy_01, c11=chorzy_11)
end_1 <- Sys.time()
czas_1 = end_1 - start_1
czas_1

rdzenie <- detectCores()
klaster <- makeCluster(6)
clusterExport(klaster, c('zdrowi_00', 'zdrowi_00', 'zdrowi_00', 'chorzy_00', 'chorzy_01', 'chorzy_11', 'rownoleglenie'))
start_2 <- Sys.time()
clusterApply(klaster, x = c(zdrowi_00, zdrowi_01, zdrowi_11, chorzy_00, chorzy_01, chorzy_11), fun = rownoleglenie(z00=zdrowi_00, z01=zdrowi_01, z11=zdrowi_11, c00=chorzy_00, c01=chorzy_01, c11=chorzy_11))
end_2 <- Sys.time()
stopCluster()
czas_2 = end_2 - start_2
czas_2

#==============================================================================

#wlasciwy kod
names(NADIR_healthy_genotypes)[1] <- "Chr"
names(NADIR_healthy_genotypes)[2] <- "Genotype"
names(NADIR_sick_genotypes)[1] <- "Chr"
names(NADIR_sick_genotypes)[2] <- "Genotype"

nazwy <- unique(NADIR_healthy_genotypes$Chr)
chromosomy <- grep("^[Chr].*", nazwy, value=TRUE)

SNP <- function (dane, tabela_1, tabela_2) {
  start <- Sys.time()
  SNP_o_znaczeniu_biologicznym <- numeric(length(dane))
  stosunek_SNP <- numeric(length(dane))
  nr_chromosomu <- 1
  for (chromosom in chromosomy) {
    print(paste0("Zaczynamy analize chromosomu: ", nr_chromosomu))
    chromosom_zdrowi <- NADIR_healthy_genotypes %>% filter(Chr == chromosom)
    chromosom_chorzy <- NADIR_sick_genotypes %>% filter(Chr == chromosom)
    wspolne_SNP <- inner_join(x = chromosom_zdrowi, y = chromosom_chorzy[2:18], by = "Genotype")
    ilosc_SNP <- nrow(wspolne_SNP)
    SNP_chromosom <- wspolne_SNP %>% filter(wspolne_SNP[3:34] != "0/2", wspolne_SNP[3:34] != "1/2", wspolne_SNP[3:34] != "2/2", wspolne_SNP[3:34] != "0/3", wspolne_SNP[3:34] != "1/3", wspolne_SNP[3:34] != "2/3", wspolne_SNP[3:34] != "3/3", wspolne_SNP[3:34] != "./.")
    
    print(paste0("Wykonano tabele dla chromosomu: ", nr_chromosomu))
  
    dlugosc <- nrow(SNP_chromosom)
 
    zdrowi_00 <- rowSums(SNP_chromosom[3:18] == "0/0")
    zdrowi_01 <- rowSums(SNP_chromosom[3:18] == "0/1")
    zdrowi_11 <- rowSums(SNP_chromosom[3:18] == "1/1")
    chorzy_00 <- rowSums(SNP_chromosom[19:34] == "0/0")
    chorzy_01 <- rowSums(SNP_chromosom[19:34] == "0/1")
    chorzy_11 <- rowSums(SNP_chromosom[19:34] == "1/1")
    
    print("Rozpoczynamy poszukiwanie istotnych SNP")
    czas_1 <- Sys.time()
  
    mutacje <- 0
    for (element in seq(dlugosc)) {
      tabela_0001 <- rbind(c(zdrowi_00[element] + zdrowi_01[element], zdrowi_11[element]), c(chorzy_00[element] + chorzy_01[element], chorzy_11[element]))
      test_0001 <- chisq.test(tabela_0001, correct = TRUE)
    
      tabela_0011 <- rbind(c(zdrowi_00[element] + zdrowi_11[element], zdrowi_01[element]), c(chorzy_00[element] + chorzy_11[element], chorzy_01[element]))
      test_0011 <- chisq.test(tabela_0011, correct = TRUE)
    
      tabela_0111 <- rbind(c(zdrowi_01[element] + zdrowi_11[element], zdrowi_00[element]), c(chorzy_01[element] + chorzy_11[element], chorzy_00[element]))
      test_0111 <- chisq.test(tabela_0111, correct = TRUE)
    
      wartosci_p <- c(test_0001$p.value, test_0011$p.value,test_0111$p.value)
      wartosci_p <- ifelse(is.nan(wartosci_p),1,wartosci_p)
    
      if (wartosci_p[1] < 0.05 | wartosci_p[2] < 0.05 | wartosci_p[3] < 0.05) {
        mutacje <- mutacje + 1
      }
    }
    czas_2 <- Sys.time()
    roznica_czasu <- czas_2 - czas_1
    SNP_o_znaczeniu_biologicznym[nr_chromosomu] <- mutacje
    stosunek_SNP[nr_chromosomu] <- mutacje/ilosc_SNP
    print(paste0("Ilosc mutacji na chromosomie ", nr_chromosomu, " wynosi: ", mutacje))
    print(paste0("Czas poszukiwan istotnych SNP: ", roznica_czasu))
    nr_chromosomu <- nr_chromosomu + 1
  }
  koniec <- Sys.time()
  czas <- koniec - start
  print(paste0("Czas wykonania funkcji: ", czas))
  print(SNP_o_znaczeniu_biologicznym)
}

system.time(SNP(chromosomy, NADIR_healthy_genotypes, NADIR_sick_genotypes))

rdzenie <- detectCores()
klaster <- makeCluster(4)
clusterExport(klaster, c('chromosomy', 'NADIR_healthy_genotypes', 'NADIR_sick_genotypes'))
start_2 <- Sys.time()
clusterApply(klaster, x = c(chromosomy), fun = SNP(chromosomy, NADIR_healthy_genotypes, NADIR_sick_genotypes))
end_2 <- Sys.time()
stopCluster()
czas_2 <- end_2 - start_2
czas_2

#==============================================================================

#test czasu

#==============================================================================

#dodatkowe funkcje

licz <- function (dane, tabela_1, tabela_2) {
  nr_chromosomu <- 1
  for (chromosom in chromosomy) {
    print(paste0("Zaczynamy analize chromosomu: ", nr_chromosomu))
    chromosom_zdrowi <- NADIR_healthy_genotypes %>% filter(Chr == chromosom)
    chromosom_chorzy <- NADIR_sick_genotypes %>% filter(Chr == chromosom)
    print(nrow(chromosom_zdrowi))
    print(nrow(chromosom_chorzy))
    nr_chromosomu <- nr_chromosomu + 1
  }
}

licz(chromosomy, NADIR_healthy_genotypes, NADIR_sick_genotypes)


liczba_SNP <- function (dane, tabela_1, tabela_2) {
  nr_chromosomu <- 1
  for (chromosom in chromosomy) {
    print(paste0("Zaczynamy analize chromosomu: ", nr_chromosomu))
    chromosom_zdrowi <- NADIR_healthy_genotypes %>% filter(Chr == chromosom)
    chromosom_chorzy <- NADIR_sick_genotypes %>% filter(Chr == chromosom)
    wspolne_SNP <- inner_join(x = chromosom_zdrowi, y = chromosom_chorzy[2:18], by = "Genotype")
    print(nrow(wspolne_SNP))
    nr_chromosomu <- nr_chromosomu + 1
  }
}
liczba_SNP(chromosomy, NADIR_healthy_genotypes, NADIR_sick_genotypes)
