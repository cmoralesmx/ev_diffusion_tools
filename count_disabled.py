import argparse
from collections import Counter

def count_states(xml_file):
    td = '<disabledReason>'
    tl = '<last_ev_id>'
    ntd = len(td)
    ntl = len(tl)
    cnt = Counter()
    total_evs = 0
    
    with open(xml_file, 'r') as xin:
        for line in xin:
            if len(line) > ntd and line[: ntd] == td:
                # check the value
                cnt[line[ntd : ntd + 1]] += 1
            elif len(line) > ntl and line[: ntl] == tl:
                #print(line)
                #print(line[ntl: -(ntl + 2)])
                #break
                total_evs += int(line[ntl: -(ntl + 2)])
    print('Counts per state:')
    t = 0
    for k, v in cnt.items():
        print(k, v)
        t += v
    print('Total EVs in current state', t)
    print('Total EVs at any point during simulation', total_evs)
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('xmlFile', help='XML state file to read')
    args = p.parse_args()

    if args.xmlFile:
        print('Counting states from', args.xmlFile)
        count_states(args.xmlFile)