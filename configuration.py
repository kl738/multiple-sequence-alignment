import json

'''
Config class for heuristics
'''
class Config:
    def __init__(self, scores, d, e, seq_weighting, init_gap_penalties, pos_gap_penalties, weight_m, divergent_seq):
        self.scores = scores
        self.d = d
        self.e = e
        self.seq_weighting = seq_weighting 
        self.init_gap_penalties = init_gap_penalties
        self.pos_gap_penalties = pos_gap_penalties
        self.weight_m = weight_m
        self.divergent_seq = divergent_seq

def parseConfigFile(config_file):
    with open(config_file) as f:
        lines = [line.strip() for line in f.readlines()]
        scores_path = lines[0]
        scores = {}
        matrices = ['BLOSUM30','BLOSUM45','BLOSUM62','BLOSUM80']
        for m in matrices:
            with open(scores_path+m+'.json') as f:
                s = json.loads(f.readlines()[0])
                scores[m] = s
        d = int(lines[1])
        e = int(lines[2])
        seq_weighting = int(lines[3])
        init_gap_penalties = int(lines[4])
        pos_gap_penalties = int(lines[5])
        weight_m = int(lines[6])
        divergent_seq = int(lines[7])
        config = Config(scores,d,e,seq_weighting,init_gap_penalties,pos_gap_penalties,weight_m,divergent_seq)
    return config