#funkcja zapisujaca dane wyjsciowe do pliku .txt
def zapis_do_pliku(output, zawartosc):
    plik=open(output+'.txt', 'w')
    plik.write(str(zawartosc))
    plik.close()

#random pozwala na generowanie losowych wynikow
import random

base=["A","C","G","T"]

#generowanie losowej sekwencji 10 nukleotydow
def losowa_sekwencja(lista):
    seq=""
    for i in range(0,10):
        i=random.choice(lista) #tworzy losowy ciag 10-literowy z liter zawartych w liscie wejsciowej
        seq+=i
    return(seq)

nic_1=(losowa_sekwencja(base))

#slownik na potrzeby utworzenia nici komplementarnej
rep={'A':'T', 'C':'G', 'G':'C', 'T':'A'}

#replikacja nici komplementarnej(bez odwracania)
def replikacja(strand, slownik):
    seq_2=""
    for litera in strand:
        seq_2+=(slownik[litera]) #tworzy ciag 10-literowy z wartosci odpowiadajacych podanym kluczom wedlug slownika
    return(seq_2)

nic_2=(replikacja(nic_1, rep))

#zamiana DNA na mRNA
def transkrypcja(strand):
    seq_3=""
    for litera in strand:
        if litera=="T":
            seq_3+="U" #zamiana kazdego wystapienia "T" na "U" w celu utworzenia mRNA
        else:
            seq_3+=litera #A, C i G pozostaja bez zmian
    return(seq_3)
    
nic_3=(transkrypcja(nic_2))

#slownik kodonow niezbedny do translacji
Aminokwasy={
'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys',
'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys',
'UUA': 'Leu', 'UCA': 'Ser', 'UAA': 'STOP!', 'UGA': 'STOP!',
'UUG': 'Leu', 'UCG': 'Ser', 'UAG': 'STOP!', 'UGG': 'Trp',

'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg',
'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg',
'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg',
'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',

'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser',
'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser',
'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg',
'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',

'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly',
'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly',
'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',
'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'}

#translacja mRNA do oligopeptydow
def translacja(strand, slownik, start):
    peptyd=""
    for i in range(0,(len(strand)-1)):
        kodon=strand[start:start+3] #tzw. ramka odczytu, czyli zczytywanie 3 kolejnych liter
        if (len(kodon))==3 and kodon in slownik: #jezeli dlugosc kodonu=3 i kodon ma odpowiednik w slowniku
            peptyd+=slownik[kodon] #mozna dodac aminokwas do "lancucha"
            start=start+3 #przesuniecie ramki odczytu
        else:
            break #jezeli warunki if nie sa spelnione, nalezy zakonczyc odczyt
    return peptyd

pep_1=(translacja(nic_3, Aminokwasy, 0)) #start ramki odczytu od pierwszej pozycji
pep_2=(translacja(nic_3, Aminokwasy, 1)) #start ramki odczytu od drugiej pozycji
pep_3=(translacja(nic_3, Aminokwasy, 2)) #start ramki odczytu od trzeciej pozycji

peptydy=(pep_1+"\n"+pep_2+"\n"+pep_3)

#wywolania wszystkich funkcji i zapis do plikow
zapis_do_pliku("wygenerowana sekwencja", nic_1)
zapis_do_pliku("nic komplementarna", nic_2)
zapis_do_pliku("mRNA",nic_3)
zapis_do_pliku("translacja", peptydy)


#testy do wszystkich funkcji
#kazdy test w momencie kompilacji wyrzuca blad
#prawidlowe wartosci/odpowiedzi znajduja sie w komentarzach po prawej stronie

#assert losowa_sekwencja(base)==" ", 'error: lista jest pusta' #losowa lista wygeneruje zawsze inny wynik
#assert replikacja("AACC", rep)=="TTAG", 'error: zla zamiana liter' #"TTGG"
#assert transkrypcja("ATATCCT")=="ATAUCCU", 'error:nieprawidlowa zamiana T na U' #"AUAUCCU"
#assert translacja("GGGCCCAAA", Aminokwasy, 0)=="GlyProLeu", 'error: blad w odczycie mRNA' #"GlyProLys"