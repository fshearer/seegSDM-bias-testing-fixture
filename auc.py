import csv
import os.path
import sys

diseases = ("cchf", "chik", "deng", "hat", "melio", "nwcl", "owcl", "scrub");
variant = ("1A", "1B", "1C", "2C", "2D", "3B", "4A", "5A", "6A")

sum = 0.0
count = 0
skip = True
for d in diseases:
    print (d)
    sys.stdout.write("[")
    for v in variant:
        with open("C:/Temp/r/" + d + "_" + v + '/statistics.csv', 'r') as csvfile:
            lines = csv.reader(csvfile, delimiter=',', quotechar='"')
            for row in lines:
                if skip:
                    skip = False
                    continue
                sum += float(row[0])
                count += 1
            sys.stdout.write("\"%.4f\", " % (sum/count))
            sum = 0.0
            count =0
            skip = True
    sys.stdout.write("],\n")