from BioFile import *

acc = "NC_011916.1"

record = SeqIO.read(fetch_file(acc, "test@example.com", None, "gbwithparts"), "gb")

for feature in record.features:
    name = ""
    if feature.type == 'CDS' and 'gene' in feature.qualifiers:
        name = feature.qualifiers['gene'][0] 
        if len(name) < 4:
            name = name + " " * (4 - len(name))
        print(f"{name}\t {int(feature.location.start)} : {int(feature.location.end)}")
    if feature.type == 'rep_origin':
        name = "oriC"
        print(f"{name}\t {int(feature.location.start)} : {int(feature.location.end)}")