import coreapi
import json

remap = {}

client = coreapi.Client()
codec = coreapi.codecs.CoreJSONCodec()

url = "http://remap.univ-amu.fr:80/api/v1/list/biotypes/taxid=9606"
biotypes = json.loads(codec.encode(client.get(url)))

for bt in biotypes["targets"]:
    url = "http://remap.univ-amu.fr:80/api/v1/datasets/findByBiotype/biotype=%s&taxid=9606" % bt["biotype_name"]
    datasets = json.loads(codec.encode(client.get(url)))
    for ds in datasets[bt["biotype_name"]]["datasets"]:
        remap.setdefault((bt["biotype_name"], ds["biotype_modification"]), 0)
        remap[(bt["biotype_name"], ds["biotype_modification"])] += 1

for biotype_name, biotype_modification in sorted(remap):
    print("{}\t{}\t{}".format(biotype_name, biotype_modification, remap[(biotype_name, biotype_modification)]))