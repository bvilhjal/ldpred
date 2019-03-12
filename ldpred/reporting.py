"""
Bunch of functions for beautifying output.
"""

def print_summary(summary_dict, title, log_file=None):
    print('')
    title_len= len(title)
    num_dashes = 80-title_len-2
    left_dashes = '-'*int(num_dashes/2)
    right_dashes = '-'*int(num_dashes/2+num_dashes%2)
    
    print( left_dashes  + ' '+title+' '+right_dashes)
    for k in sorted(summary_dict.keys()):
        sd = summary_dict[k]
        if sd['name']=='dash':
            print ('-' * 80)
        else:
            val_str= str(sd['value'])
            if len(val_str)>29:
                print ('{:<50}'.format(sd['name']))
                print ('{:>80}'.format(val_str))
            else:
                print ('{:<50}{:>30}'.format(sd['name'],str(sd['value'])))
    print ('-' * 80)
    print('')
        