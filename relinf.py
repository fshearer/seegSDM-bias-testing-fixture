import csv
diseases = ("cchf", "chik", "deng", "hat", "melio", "nwcl", "owcl", "scrub")
variants = ("1A", "1B", "1C", "2C", "2D", "3B", "4A", "5A", "6A")

string = ""
count = 0
skip = True
for disease in diseases:
  for variant in variants:
    with open('all/results/' + disease + "_" + variant + '/relative_influence.csv', 'r') as csvfile:
        lines = csv.reader(csvfile, delimiter=',', quotechar='"')
        for row in lines:
            if skip:
                skip = False
                continue
            string += ((row[0] +",%.2f,") % float(row[1]))
        print((disease + "_" + variant  + ",") + string)
        
        string = ""
        count =0
        skip = True
