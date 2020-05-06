

sections = ['isthmus', 'ampulla']

classes2 = {'ordering':[('narrow_end','#30a2da'), ('wide_end','#fc4f30'),
            ('narrow_lumen','#e5ae38'),('wide_lumen','green')]}
#classes2['utj'] = {'narrow_end':[1,2,3, 7,8, 13,21], 'wide_end':[4,5,6, 9,14,19], 'narrow_lumen':[15, 16, 17,18], 'wide_lumen':[10,11,12, 20]}
classes2['isthmus'] = {'narrow_end':[1,2,3,4,5,6], 
                        'wide_end':[7,8, 9],
                        'narrow_lumen':[10, 11, 12, 13, 14, 15],
                        'wide_lumen':[16,17,18,19,20,21]}
#classes2['ia-junction'] = {'narrow_end':[1,2,3,4,5,6,7,8], 'wide_end':[9,10,11,12,13,14,15,16], 'narrow_lumen':[17,18,19,20,21,22,23,24], 'wide_lumen':[25,26,27,28,29,30,31,32]}
classes2['ampulla'] = {'narrow_end':[1,2,3,4,5,6],
                        'wide_end':[7,8,9,10, 11, 12],
                        'narrow_lumen':[13,14,15,16,17,18], 
                        'wide_lumen':[19, 20, 21, 22, 23, 24]}



def compare_single(section1, path1, iter1, rep1):
    print('Comparing single section:')
    print(f'   section: {section1}')
    print(f' iteration: {iter1}')
    print(f'replicates: {rep1}')
    print(f'      path: {path1}')

def compare_double(section1, path1, iter1, rep1, section2, path2, iter2, rep2):
    print( 'Comparing->    [1]   |   [2]')
    print(f'   section:  {section1} | {section2}')
    print(f' iteration:    {int(iter1)/1e6:>.1f}M | {int(iter2)/1e6:>.1f}M')
    print(f'replicates:        {rep1} | {rep2}')
    print(f'path [1]: {path1}')
    print(f'path [2]: {path2}')
    


if __name__ == '__main__':
    import argparse

    """
    comparisons possible:
    1 section - 1) within classes, 2) across classes
    2 sections - 1) within classes, 2) across classes
    """

    parser = argparse.ArgumentParser(add_help=True)
    s1 = parser.add_argument_group('section1', 'Settings for the first section to compare')
    s1.add_argument('section', choices=sections, help='The 1st section to analyse.')
    s1.add_argument('path', metavar='base_path', help='Directory where the target files are located (1st section to compare)')
    s1.add_argument('--iteration1', nargs='?', type=int, help='Iteration number to read - 1st section. Default 14.4k')
    s1.add_argument('--replicates1', nargs='?', type=int, help='The number of replicates to read - 1st section. Default 3')
    
    s2 = parser.add_argument_group('section2', 
        """Settings for a second section to compare.
        Section2 is optional, however, it's required arguments are:
        --other_section, --other_roi_class, --other_path
        """)
    s2.add_argument('--section2', choices=sections, help='The 2nd section to analyse.')
    s2.add_argument('--path2', help='Directory where the target files are located (2nd section to compare)') 
    s2.add_argument('--iteration2', nargs='?', type=int, help='Iteration number to read - 2nd section. Default 14.4k')
    s2.add_argument('--replicates2', nargs='?', type=int, help='The number of replicates to read - 2nd section. Default 3')
    args = parser.parse_args()


    iter1 = args.iteration1 if args.iteration1 else 14400000
    rep1 = args.replicates1 if args.replicates1 else 3
    
    if args.section2 or args.path2:
        
        iter2 = args.iteration2 if args.iteration2 else 14400000
        rep2 = args.replicates2 if args.replicates2 else 3
        compare_double(args.section, args.path, iter1, rep1,
                args.section2, args.path2, iter2, rep2)
    else:
        compare_single(args.section, args.path, iter1, rep1)