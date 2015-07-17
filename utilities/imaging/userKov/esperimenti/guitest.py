#!/usr/bin/env python

import wx

ID_ABOUT=101
ID_EXIT=110

class MiaApp(wxApp):
    def OnInit(self):
        frame = wxFrame(None, -1, "Ciao mondo")
        # viene mostrata la finestra
        frame.Show(1)
        # imposta la finestra principale
        self.SetTopWindow(frame)
        return 1
    # crea un'istanza della classe MiaApp

#----------------------------------------------------------------------

def runTest(frame, nb, log):
    win = TestPanel(nb, log)
    return win

#----------------------------------------------------------------------

if __name__ == '__main__':
    import sys,os
    import run
    app = MiaApp(0)
    app.MainLoop()
    cc = wx.combo.ComboCtrl(style=style, size=(250,-1))
    app.control=cc
    # Create a Popup
    popup = ListCtrlComboPopup()
    # Associate them with each other.  This also triggers the
    # creation of the ListCtrl.
    cc.SetPopupControl(popup)
    # Add some items to the listctrl.
    for x in range(75):
         popup.AddItem("Item-%02d" % x)
    MainWindow.control=cc

