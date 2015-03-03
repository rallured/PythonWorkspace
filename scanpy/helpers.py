#Define a few functions for plotting histogram of automatic scanning runs.

from scipy import *
from numpy import *
from matplotlib.pyplot import *
import curvefitting as c
import sys


def plot_slice(info, x = None, y = None):
    if x != None:
        d = info[where(info['x'] == x)]
        d.sort(order = 'y')
        plot(d['y'], d['centroid'])
        return
    if y != None:
        d = info[where(info['y'] == y)]
        d.sort(order = 'x')
        plot(d['x'], d['centroid'])
        return

def plot_hist(data, info, x, y, x_offset = 10, y_offset = 5):
    x = x - x_offset
    y = y - y_offset
    print('Plotting at positions (' + str(x) + ', ' + str(y) + ')')
    print('Since an offset of y = ' + str(y_offset) + ' and x = ' +
    str(x_offset) +' is applied.')
    temp = info[where(info['x'] == x)]
    c_hist = data[where(temp['y'] == y)]
    c_hist[0]
    len(c_hist[0])
    absc = arange(0, 256)
    print absc
    print c_hist[0]
    print len(absc), len(c_hist[0])
    plot(absc, c_hist[0])

def manual_fitting(fname, filename_pos = True, stepsize = 1):
    """Routine for setting up a whole bunch of fitting. Basically, this program displays a graph, then asks the user for input on help to generate a proper fit. These inputs include:
    Left bound region of interest:
    Right bound region of interest:
    nterms for fit: 3 or 4
    
    """
    ion()
    bounds = []
    absc = arange(256)
    final_list = []
    data_list = []
    for x in open(fname):
        while True:
            hold(False)
            data = hist(loadtxt(x.split()[0]), bins = 256, range = (0, 255))[0]
            #These next two lines assume a filename of the form:
            #'specx_[-+][0-9][0-9]y_[-+][0-9][0-9][ab].txt\n'
            if filename_pos ==True:
                curr_absc = x.split()[0].split('_')[1][:-1]
                curr_ord = x.split()[0].split('_')[2][:-5]
            else:
                curr_absc = 0
                curr_ord = 0
            
            #Display the current spectrum name.
            print('Currently on spectrum: ' + x)

            #Give the user the option to mark this image bad and move on.
            if raw_input('Do you want to mark this image bad and move on? ') in ('Y', 'y'):
                #Append to the list in the order of
                #Specname, x, y, Cent, left, right, nterms, 
                #cent value of -1 agreed upon with Ryan Allured.
                final_list.append((x.split()[0], curr_absc, curr_ord, -1, 0, 0, 0))
                data_list.append(data)
                break

            #Take in region of interest and fit nterms data.
            try:
                left = int(raw_input('Enter the left bound of the desired ROI for ' + x + ':\n'))
                right = int(raw_input('Enter the right bound of the desired ROI for ' + x + ':\n'))
                fit_terms = int(raw_input('Enter the number of terms for the gaussian fit. 3 for just the peak, 4 to add a constant offset:'))
            #Catch any typos so that the whole data run isn't thrown away.
            except ValueError:
                print('One of the inputs was not correct! Start over.\n')
                #Restart the innermost while loop.
                continue

            #Run the fit, after slicing the input data array according to the
            # ROI.
            rois = (left, right)
            d_slice = data[left:right]
            a_slice = absc[left:right]
            print('Running fit with region of interest at (' + str(left) + 
            ', ' + str(right) + ' and nterms = ' + str(fit_terms))
            try:
                c_fit = c.gaussguess(a_slice, d_slice, nterms = fit_terms)
            except ValueError:
                print('Something went wrong with fitting! Restarting the loop.')
                print('Try adjusting the ROI.')
                print('If it keeps on not working, mark it bad.')
                #TODO: Fix this issue in gaussguess.
                continue
            plot(absc, data)
            hold(True)

            #Attempt to plot the gaussian. Since sometimes the fitting #routines return arrays rather than single numbers, 
            #I've put in a branch to give the user control over which fit #result is used.
            try:
                plot(a_slice, c.gaussian(c_fit, a_slice))
            except ValueError:
                print('Something went wrong with plotting the fit. Printing'
                + ' found fit parameters...\n')
                print('Fit parameters array is ' + str(c_fit) + '\n')
                if raw_input('Restart the loop? ') in ('y','Y'): 
                    continue
                else:
                    try:
                        index = int(raw_input('Enter the desired index to work with:'))
                    except ValueError:
                        print('Wrong input! Repeating the loop...')
                        continue
                    #This loop replaces all the elements of c_fit that are 
                    #arrays with single elements, according to the user input 
                    #index #.
                    for i in range(len(c_fit)):
                        if len(c_fit[i]) >= index + 1:
                            c_fit[i] = array([c_fit[i][index]])
                try:
                    plot(a_slice, c.gaussian(c_fit, a_slice))
                except ValueError:
                    print('Still didn\'t work! Continuing...')
                    continue

            print('Fitted centroid is ' + str(c_fit[1][0]) + '\n')
            next = raw_input('Enter N, n if you are happy and want to go on to the next spectrum, otherwise press enter:')
            if next == 'N' or next == 'n':
                #Append to the list in the order of
                #Cent, left, right, nterms, 
                final_list.append((x.split()[0], curr_absc, curr_ord, c_fit[1][0], rois[0], rois[1], fit_terms))
                data_list.append(data)
                break

    
    #Stick everything into arrays with column names.
    info_arr = array(final_list, dtype = [('spec name', 'a25'), ('x', '<i4'), ('y', '<i4'), ('centroid', '<i4'), ('roi left', '<i4'), ('roi right', '<i4'), ('nterms', '<i4')])    
    data_arr = array(data_list)
    
    xs = info_arr['x']
    ys = info_arr['y']
    
    mat2d = zeros(((max(xs) - min(xs)) / stepsize + 1, (max(ys) - min(ys)) / stepsize + 1))
    
    #NB: This array is output with respect to the Python row, column standard.
    #Make sure to transpose it if you wish to use it in IDL for making a 
    #contour plot or surface.
    for temp in zip(xs, ys):
        t = info_arr[where(xs == temp[0])]
        val = t[where(ys[where(xs == temp[0])] == temp[1])]
        mat2d[temp[0] / stepsize, temp[1] / stepsize] = val['centroid'][0]
        
    save('info_array.npy', info_arr)
    save('data_array.npy', data_arr)
    savetxt('2dmatrix.txt', mat2d)
        
    
    return (info_arr, data_arr, mat2d)

#

def gen_2d_matrix(info, x_interval = 1, y_interval = 1):
    """Quick function to generate a 2d matrix from an information array as given by manual_fitting().
    
    Inputs:
        info -      Info array from manual_fitting.
    Optional Inputs:
        x_interval -    Corresponds to x_jog_size on the automatic acquisition             
                            routine.
        y_interval -    Corresponds to y_jog_size on the automatic             
                            acquisition.
        
    """
    xs = info['x']
    xs -= min(xs)
    ys = info['y']
    ys -= min(ys)
    
    bounds = ((max(xs) - min(xs) + x_interval) / x_interval, (max(ys) - min(ys) + y_interval) /
                y_interval)
    m = zeros(bounds)
    
    c = 0
    
    for temp in zip(xs, ys):
        a = info[where(xs == temp[0])]
        val = a[where(ys[where(xs == temp[0])] == temp[1])]
        m[temp[0] / x_interval, temp[1]/ y_interval] = val['centroid'][0]
        print val
        c += 1
        print c
    return m

    
 
def gen_gain_curve(fname, volt_info = True):
    """
    Quick routine to perform fitting on the spectra acquired from a typical gain curve run.
    
    Inputs:
        fname -         List of filenames of spectra.
    
    Optional Inputs:
        volt_info -     Boolean. If True, filenames conform to the format hv_%i.txt, where %i is an integer denoting the number of volts. EG, a spectrum at 2.52 kv is hv_2520.txt
    
    """
    info, data = manual_fitting(fname, filename_pos = False)
    volts = []
    print info
    print data
    for x in info['spec name']:
        if volt_info == True:
            volts.append(int(x.split('_')[1].split('.')[0]))
    return array(zip(volts, info['centroid']))
 
 
 
 
 
 
 

        
"""The following functions were part of a suite that I threw together one afternoon in order to generate the PNG, along with some colormap stuff that I was messing around with. I wouldn't pay too much attention to it at this point.
"""


def hist_datas(fname, man_rois = False):
    #Reads the spectra from the file and loads them into a hash map.
    hists = []
    names = []
    for x in open(fname):
        if man_rois == True:
            bounds = x[1]
            hists.append(histogram(loadtxt(x.split()[0])[bounds[0]:bounds[1]], bins = 256, range = (0, 255))[0])
            names.append(x.replace('spec', '').replace('.txt','').split()[0])
        else:
            hists.append(hist(loadtxt(x.split()[0]), bins = 256, range = (0, 255))[0])
            names.append(x.replace('spec', '').replace('.txt','').split()[0])
    x = arange(256)
    hist_d = dict(zip(names, hists))
    return hist_d
    
def manual_filtering(fname):
    #This function exists primarily to facilitate the fitting of data using
    #user defined input in order to cut down on the bad data mucking up a
    #particular fit.
    return None
        
def cut_bgs(hist_dict):
    #Helper function to cut background from spectrae.
    #Really, this just zeroes the first two channels.
    for x in hist_dict:
        hist_dict[x][0] = 0
        hist_dict[x][1] = 0
        hist_dict[x][2] = 0
        hist_dict[x][3] = 0
        hist_dict[x][4] = 0

    

def run_fits(hist_dict, fit_terms = 4):
    fits = []
    names = []
    xs = arange(256)
    for h in hist_dict:
        print h
        try:
            fits.append(c.gaussguess(xs, hist_dict[h], nterms = fit_terms))
            names.append(h)
        except:
            print sys.exc_info()[0]
            #None means the fit failed to run.
            fits.append(None)
            names.append(h)
    fits_d = dict(zip(names, fits))
    return(fits_d)
        
def vis_fits(fits, hists, nterms = 4, prefix = '', xmin = 0, ymin = 0):
    """Writes a series of images to the cwd that is a plot of the data with a fit overlaid.
    """
    absc = arange(256)
    output_vals = []
    prefix = prefix + 'vis_' + str(nterms) + '_'
    print 'test'
    for x in fits:
        if fits[x] == None:
            continue
        if nterms == 3:
            pars = (fits[x][0], fits[x][1][0], fits[x][2][0])
        elif nterms == 4:
            pars = (fits[x][0], fits[x][1][0], fits[x][2][0], fits[x][3])
        ord = c.gaussian(pars, absc)
        output_vals.append((x, absc, ord))
        plot(absc, hists[x],hold = False)
        plot(absc, ord)
        draw()
        savefig(prefix + adjust_vals(str(x), x_min = xmin, y_min = ymin) + '.png', format='png')
    return output_vals
    
def adjust_vals(name_str, x_min = 0, y_min = 0):
    x = int(name_str.split('y')[0].split('_')[1])
    y = int(name_str.split('y')[1].split('_')[1][:-1])
    x += x_min
    y += y_min
    return 'x_ ' + str(x) + 'y_ ' + str(y) + name_str[-1]
    
    
        
def plot_with_weights(hist_dict, fit_dict):
    #Generates a plot of data points that are colored according to a colormap
    #from 0x00ff00 to 0xff, eg. black to light blue along the blue spectrum.
    colormax = float(int('0x00ff00', 16)) #Green
    colormin = float(int('0x0000ff', 16)) #Blue
    plot_infos = []
    cents = []
    for x in fit_dict:
        if fit_dict[x] == None:
            continue
        cents.append(fit_dict[x][1][0])
    w_max = float(max(cents))
    w_min = float(min(cents))
    w_range = float(w_max - w_min)
    for x in fit_dict:
        if fit_dict[x] == None:
            continue
        absc = int(x.split('_')[1][:-1])
        ord = int(x.split('_')[2][:-1])
        cent = fit_dict[x][1][0]
        color = (cent - w_min) / w_range * (colormax - colormin) + colormin
        plot_infos.append((x,(absc, ord, cent, hex(int(color)))))
    p_inf = dict(plot_infos)
    alt = ([], [])
    for x in p_inf:
        abs = p_inf[x][0]
        ord = p_inf[x][1]
        c = '#' + p_inf[x][3].replace('0x', '').zfill(6)
        #plot([abs], [ord], color = c, marker = 'o')
        alt[0].append(abs)
        alt[1].append(p_inf[x][2])
    plot(alt[0], alt[1], marker = 'o', linewidth = 0)
    return p_inf
    
def gen_2dcentmap(fit_dict):
    xs = []
    ys = []
    cents = []
    for el in fit_dict:
        s = el.split('_')
        curr_x = int(s[1][:-1])
        curr_y = int(s[2][:-1])
        if fit_dict[el] == None:
            cent = None
        else:
            cent = fit_dict[el][1]
        xs.append(curr_x)
        ys.append(curr_y)
        cents.append(cent)
    output = zeros((max(xs) - min(xs) + 1, max(ys) - min(ys) + 1))
    for x in zip(xs, ys, cents):
        if x[2] == None:
            output[x[0], x[1]] = None
        else:
            output[x[0], x[1]] = x[2][0]
    return output        
