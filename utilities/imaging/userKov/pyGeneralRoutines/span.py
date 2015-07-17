import numpy as np

def span(array,size=False):
	"""return min and max of an array as a tuple.
	If size is set to True, return the size of the interval instead."""
	min=np.nanmin(array)
	max=np.nanmax(array)
	if size:
		return max-min
	else:
		return (min,max)