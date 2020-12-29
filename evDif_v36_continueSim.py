import argparse
import re
from os import listdir, rename
from os.path import isfile, join

constants = ['const_Boltzmann',
                'const_Boltzmann_x_Temp_K',
                'const_boltzmann_x_temp',
                'const_mass_per_volume_unit',
                'const_mass_p_vol_u_x_4div3_pi',
                'minMaxCoords']
evi = re.compile(r'closest_ev_initial')
tiis = re.compile(r'time_in_initial_state')
otini = re.compile(r'origTimeInInitial')
vums = re.compile(r'velocity_ums')
ba = re.compile(r'<becameAt>\d*,')

def rename_variables(line, of):
    nts = ['defaultAt','disabledAt', 'disabledReason']

    if evi.search(line):
        of.write(re.sub('initial', 'biogenesis', line, 2))
    elif tiis.search(line):
        of.write(re.sub('time_in_initial_state', 'time_in_biogenesis_state', line, 2))
    elif otini.search(line) or vums.search(line):
        # skip this line
        pass
    elif ba.match(line):
        # produces 3 lines
        for t,v in zip(nts, line[10:-12].split(', ')):
            of.write(f'<{t}>{v}</{t}>\n')
    else:
        of.write(line)

def remove_constants(line, of):
    # replace constants
    tag = line[line.index('<') + 1 : line.index('>')]
    if tag in constants:
        print('Removed:',line.strip())
    else:
        of.write(line)

def prepare_xml(xml_file):
    env = re.compile(r'\s*</environment>')
    
    reading_environment = True
    with open(xml_file, 'r') as inp:
        n = xml_file.rindex('/')+1
        new_file = xml_file[:n] + 'new_0_' + xml_file[n:]
        with open(new_file, 'w') as of:
            for line in inp:
                if reading_environment:
                    # check we are still reading environment entries
                    if env.match(line):
                        print('found closing ENVIRONMENT tag')
                        reading_environment = False
                    else:
                        remove_constants(line, of)
                else:
                    # rename variables
                    rename_variables(line, of)
        print('Done. New 0.xml file:', new_file)

def rename_files(path, rlog):
    print('Browsing files in', path)
    fs = [f for f in listdir(path)]
    
    iterations = []
    for f in fs:
        if f[:6] == 'new_0_':
            print('0.xml:', f)
            offset = int(f[6 : -4])
            print('offset:', offset)
            rlog.write(f'0.xml: {f}. File names will be offset by: {offset}\n')
        elif f[-4:] == '.xml':
            iterations.append(int(f[:-4]))
    for i in reversed(sorted(iterations)):
        f = f'{i}.xml'
        nf = f'{offset + i}.xml'
        print(f'{f} -> {nf}')
        rlog.write(f'iteration: {f} -> {nf}\n')
        rename(join(path, f), join(path, nf))
                    
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('path', help='XML file with the last state of the system to process')
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument('--prepare', action='store_true', help='Triggers the preparation process for a new 0.xml file from a previous state file')
    group.add_argument('--rename', action='store_true', help='Fixes the naming of the XML (state) files produced by a continued simulation')

    group.add_argument('--download', action='store_true', help='Gets the files from the remote server')
    group.add_argument('--upload', action='store_true', help='Gets the prepared files to the remote server')
    p.add_argument('--server')
    p.add_argument('--username')
    p.add_argument('--password')
    p.add_argument('--remote_path')
    p.add_argument('--repeats')

    args = p.parse_args()

    if args.prepare:
        print('Will produce a new 0.xml file from', args.path)
        prepare_xml(args.path)
    elif args.rename:
        with open(join(args.path, 'rename.log'), 'w') as rlog:
            #rename_files(args.path, rlog)
            rename_files(args.path, rlog)
            rlog.write(f'DONE\n')
    elif args.download:
        # path provides the directory where the EV files should be stored
        # needs: remote path
        pass
    elif args.upload:
        # path provides the directory where the prepared EV files are stored
        pass
    else:
        print('ERROR')