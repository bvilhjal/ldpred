"""
Bunch of functions for beautifying output.
"""
window_length = 80
name_length = 50

def print_summary(summary_dict, title, log_file=None):
    print('')
    title_len= len(title)
    num_dashes = window_length-title_len-2
    left_dashes = '='*int(num_dashes/2)
    right_dashes = '='*int(num_dashes/2+num_dashes%2)
    val_length=window_length-name_length
    
    print( left_dashes  + ' '+title+' '+right_dashes)
    for k in sorted(summary_dict.keys()):
        sd = summary_dict[k]
        if sd['name']=='dash':
            if 'value' not in sd:
                print ('-' * window_length)
            else:
                val_str= str(sd['value'])
                str_len= len(val_str)
                num_dashes = window_length-str_len-2
                left_dashes = '-'*int(num_dashes/2)
                right_dashes = '-'*int(num_dashes/2+num_dashes%2)
                print( left_dashes  + ' '+val_str+' '+right_dashes)
        else:
            val_str= str(sd['value'])
            if len(val_str)>val_length-1:
                print (('{:<%d}'%name_length).format(sd['name']))
                print (('{:>%d}'%window_length).format(val_str))
            else:
                print (('{:<%d}{:>%d}'%(name_length,val_length)).format(sd['name'],str(sd['value'])))
    print ('=' * window_length)
    print('')
    
