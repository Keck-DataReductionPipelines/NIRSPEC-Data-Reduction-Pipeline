import argparse

def log2csv(content):
    
    d = None    # single data dictionary
    ds = []     # array of data dictionaries
    
    for line in content:
        
        if line.find('starting reduction of') >= 0:
            # start of a new record
            
            if d is not None:
                if dict_is_complete(d):
                    ds.append(d)
                    
            d = get_empty_dict()
        
        if d is not None:            
            read_line(line, d)
        
    if dict_is_complete(d):
        print_record(d)
        
    ds.sort()
    for d in ds:
        print_record(d)
        
    
    return

def read_line(line, d):
    
    if line.find('starting reduction of') >= 0:
        d['fn'] = line.split()[-1:][0]
    elif line.find('date of observation =') >= 0:
        d['date'] = line[line.find('=') + 2:]
    elif line.find('target name =') >= 0:
        d['target'] = line[line.find('=') + 2:]
    elif line.find('filter name =') >= 0:
        d['filter'] = line.split()[-1:][0]
    elif line.find('slit name =') >= 0:
        d['slit'] = line.split()[-1:][0]
#     elif line.find('n orders on the detector =') >= 0:
#         d['n_orders'] = line.split()[-1:][0]
#     elif line.find('n orders reduced =') >= 0:
        d['n_reduced'] = line.split()[-1:][0]
    elif line.find('mean signal-to-noise ratio = ') >= 0:
        d['snr_mean'] = line.split()[-1:][0]
    elif line.find('minimum signal-to-noise ratio = ') >= 0:
        d['snr_min'] = line.split()[-1:][0]
    elif line.find('mean spatial peak width = ') >= 0:
        if line.find('unknown') >= 0:
            d['w_mean'] = ''
        else:
            d['w_mean'] = line.split()[-2:-1][0]
    elif line.find('maximum spatial peak width = ') >= 0:
        d['w_max'] = line.split()[-2:-1][0]
    elif line.find('n sky lines identified = ') >= 0:
        d['n_lines_found'] = line.split()[-1:][0]
    elif line.find('n lines used in wavelength fit = ') >= 0:
        d['n_lines_used'] = line.split()[-1:][0]
    elif line.find('rms wavelength fit residual = ') >= 0:
        d['r_rms'] = line.split()[-1:][0]
    elif line.find('integration time = ') >= 0:
        d['itime'] = line.split()[-2:-1][0]
#         d['itime'] = line.split()[-1:][0]
    elif line.find('n coadds = ') >= 0:
        d['n_coadds'] = line.split()[-1:][0]
    return
        
def get_empty_dict():
    d = dict()
    d['fn'] = None
    d['date'] = None
    d['target'] = None
    d['filter'] = None
    d['slit'] = None
#     d['n_orders'] = None
#     d['n_reduced'] = None
    d['snr_mean'] = None
    d['snr_min'] = None
    d['w_mean'] = None
    d['w_max'] = None
    d['n_lines_found'] = None
    d['n_lines_used'] = None
    d['r_rms'] = None
    d['itime'] = None
    d['n_coadds'] = None
    return d

def dict_is_complete(d):
    for k, v in d.items():
        if v is None:
            if k.lower() == 'n_lines_used':
                continue
            if k.lower() == 'r_rms':
                if d['n_lines_used'] is None:
                    continue
            else:
                return False
    return True

i = 0
def print_record(d):
    global i
    if i == 0:
#        print('index, fn, date, target, filter, slit, itime, coadds, n orders, n reduced, '),
        print(('index, fn, date, target, filter, slit, itime, coadds, '), end=' ')
        print(('snr mean, snr min, width mean, width max, n lines found, n lines used, '), end=' ')
        print('r rms')
    i += 1
    
    for k, v in d.items():
        if v is None:
            d[k] = ''
            
    print(('{}, '.format(i)), end=' ')
    print(('{}, '.format(d['fn'].strip())), end=' ')
    print(('{}, '.format(d['date'].strip())), end=' ')
    print(('{}, '.format(d['target'].strip())), end=' ')
    print(('{}, '.format(d['filter'])), end=' ')
    print(('{}, '.format(d['slit'])), end=' ')
    print(('{}, '.format(d['itime'])), end=' ')
    print(('{}, '.format(d['n_coadds'])), end=' ')
#     print('{}, '.format(d['n_orders'])),
#     print('{}, '.format(d['n_reduced'])),
    print(('{:.2f}, '.format(float(d['snr_mean']))), end=' ')
    print(('{}, '.format(d['snr_min'])), end=' ')
    print(('{}, '.format(d['w_mean'])), end=' ')
    print(('{}, '.format(d['w_max'])), end=' ')
    print(('{}, '.format(d['n_lines_found'])), end=' ')
    print(('{}, '.format(d['n_lines_used'])), end=' ')
    print(('{}, '.format(d['r_rms'])), end=' ')
    print('')

def main():
    parser = argparse.ArgumentParser(description="NIRSPEC DRP - log2csv")
    parser.add_argument('log_fn', help='log file name')
    args = parser.parse_args()
    
    try:
        f = open(args.log_fn)
    except IOError as e:
        print(('can\'t open {}: {}'.format(args.log_fn, e)))
        exit(1)
    
    log2csv(f.readlines())

if __name__ == "__main__":
    """
    line id test
    """
    main()   