import random, os, argparse

def construct_set(input_seq, clust_data, query, output_dir):
    input = open(input_seq)
    clust = open(clust_data)
    query_seq = ""
    query_num = 0
    sequences = dict()

    header = ""
    counter = 0
    seq = ""
    store_query = False
    for line in input.readlines():
        line = line.strip()
        if line[0] == '>':
            if seq != "":
                sequences[counter] = [header, seq]
                counter += 1
                if store_query:
                    store_query = False
                    query_seq = seq
                seq = ""
                if header == ">" + query:
                    store_query = True
                    query_num = counter - 1
            header = line
        else:
            seq += line
    sequences[counter] = [header, seq]
    counter += 1

    clusters = dict()
    for line in clust.readlines():
        line = line.strip()
        parts = line.split(",")
        if parts[1] in clusters:
            temp = clusters[parts[1]]
            temp.append(parts[0])
            clusters[parts[1]] = temp
        else:
            clusters[parts[1]] = [parts[0]]

    used = dict()
    biggest = 0
    for clust in clusters:
        used[clust] = list()
        if len(clusters[clust]) > biggest:
            biggest = len(clusters[clust])

    for i in range(biggest):
        current = list()
        for clust in clusters:
            if len(used[clust]) == len(clusters[clust]):
                used[clust] = list()
            while True:
                act = random.randint(0, len(clusters[clust]) - 1)
                if not act in used[clust]:
                    used[clust].append(act)
                    current.append(clusters[clust][act])
                    break

        found = False
        output = open(output_dir + "/set_" + str(i), "w")
        for cur in current:
            try:
                output.write(sequences[int(cur)][0] + "\n" + sequences[int(cur)][1] + "\n")
            except:
                continue
            if query in sequences[int(cur)][0]:
                found = True
        if not found:
             output.write(">" + query + "\n" + query_seq)

parser = argparse.ArgumentParser()
parser.add_argument("--input_seq", dest="input_seq", type=str, help="set of sequences in fasta format")
parser.add_argument("--clusters", dest="clusters", type=str, help="output of SigClust program")
parser.add_argument("--output_dir", dest="output", type=str, help="output folder")
parser.add_argument("--query", dest="query", type=str, help="name of query sequence")
args = parser.parse_args()

if args.input_seq != None and args.clusters != None and args.output != None and args.query != None:
    if args.output[-1] == "/":
        args.output = args.output[0:-1]
    construct_set(args.input_seq, args.clusters, args.query, args.output)
else:
    print("Missing arguments.\nUsage: create_sets.py --input_seq [fasta sequences] --clusters [SigCLust output] --query [query name] --output_dir [output dir]")

