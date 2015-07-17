import os

def fn_add_subfix(fileName,subfix="",newext=None):
    """Add a subfix to a filename. A new extension can be defined (dot must be included)."""
    #2014/07/28 exclude dot from extension (it must be explicitily included in the extension string).
    #   This is useful to remove extension by calling fn_add_subfix(fileName,"","") and also more
    #   consistent with os.path.splitext that returns the extension string with dot included.
    import os
    fn,ext=os.path.splitext(fileName)
    if newext!=None:
        ext=newext
    dir,base=os.path.split(fn)
    return os.path.join(dir,base+subfix+ext)