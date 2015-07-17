import wx

class MyFrame(wx.Frame):
	def __init__(self, parent, ID, title, pos=wx.DefaultPosition,size=wx.DefaultSize, style=wx.DEFAULT_FRAME_STYLE):

		#inizializza il frame chiamando il metodo init
		wx.Frame.__init__(self, parent, ID, title, pos, size, style)

		#piazza un panel e un button
		panel = wx.Panel(self, -1)
		button = wx.Button(panel, 1003, "Plot!!")
		button.SetPosition((15, 15))
		pl = wx.lib.plot.TestFrame(panel, -1, "PlotCanvas Demo")

		#collega gli eventi alle routine da chiamare
		self.Bind(wx.EVT_BUTTON, self.plotta, button)
		self.Bind(wx.EVT_CLOSE, self.OnCloseWindow)
	
	def _draw1Objects():
		# 100 points sin function, plotted as green circles
		data1 = 2.*_Numeric.pi*_Numeric.arange(200)/200.
		data1.shape = (100, 2)
		data1[:,1] = _Numeric.sin(data1[:,0])
		markers1 = PolyMarker(data1, legend='Green Markers', colour='green', marker='circle',size=1)
	
		# 50 points cos function, plotted as red line
		data1 = 2.*_Numeric.pi*_Numeric.arange(100)/100.
		data1.shape = (50,2)
		data1[:,1] = _Numeric.cos(data1[:,0])
		lines = PolyLine(data1, legend= 'Red Line', colour='red')
	
		# A few more points...
		pi = _Numeric.pi
		markers2 = PolyMarker([(0., 0.), (pi/4., 1.), (pi/2, 0.),
							  (3.*pi/4., -1)], legend='Cross Legend', colour='blue',
							  marker='cross')
		
		return PlotGraphics([markers1, lines, markers2],"Graph Title", "X Axis", "Y Axis")

	#routine da chiamare
	def plotta(self, event):
		pass

	def OnCloseWindow(self, event):
		self.Destroy()

if __name__=="__main__":
    win = MyFrame(None, -1, "Plottiamo!!", size=(550, 400),style = wx.DEFAULT_FRAME_STYLE)
    win.Show(True)