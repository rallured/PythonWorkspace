import matplotlib.pyplot as plt

def plot_labels(labels, xdata, ydata):
    """"Add fancy labels at positions xdata,ydata."""
    for label, x, y in zip(labels, xdata, ydata):
        plt.annotate(
            label, 
            xy = (x, y), xytext = (-20, 20),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
            
if __name__=="__main__":
    #this is what I want to pass as label, it must be a 3 element list with [string, xpos, ypos]
    xdata=(1.,2.,3.)
    ydata=(1.,4.,9.)
    labels=('uno','due','tre')
    
    plt.plot(xdata,ydata,'o')
    plt.plot(xdata,ydata)
    plot_labels(labels, xdata, ydata)
    plt.show()