
# uvoz modulov
import csv, random, sys, os       

# Uporaba v terminalu: python split_dataset.py <input_csv> <train_csv> <test_csv>
if len(sys.argv) != 4:      # 4 argumenti klicanja
    script = os.path.basename(sys.argv[0])     # pot do skripte
    print(f"Usage: python {script} <input_csv> <train_csv> <test_csv>")  
    sys.exit(1)  
    
# definiramo input podatke, podatke za ucenje (80%) in testne podatke(20%)
infile, trainfile, testfile = sys.argv[1], sys.argv[2], sys.argv[3]  

# odprtje datoteke shranjeno v infile in branje po vrsticah
# list naredi Python seznam
data = list(csv.reader(open(infile))) 
header, rows = data[0], data[1:] 
random.seed(42)            
random.shuffle(rows)       
cut = int(0.8 * len(rows))      
# vse vrstice od 0 do indeksa cut-1 bodo ucni set podatkov
# vse vrstice od indeksa cut, cut+1 do konca t.j. len(rows)-1 bodo testni set
train, test = rows[:cut], rows[cut:]      

# locen zapis podatkov v train.csv in test.csv
for name, subset in [(trainfile, train), (testfile, test)]:
    with open(name, 'w', newline='') as f: 
        writer = csv.writer(f)            
        writer.writerow(header)           # posebej zapises glavo
        writer.writerows(subset)          # zapis ostalih vrstic train ali test podatkov
    print(f"Wrote {len(subset)} rows to {name}") # izpis kolk vrstic v kteri file
