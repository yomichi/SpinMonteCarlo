import sys
import os
import wx
from wx import glcanvas
from OpenGL.GL import *
import matplotlib.colors as colors

import numpy as np
import numpy.random as random

import loop

class MyFrame(wx.Frame):
  ID_TIMER = 1

  def __init__(self, parent, id, size, title):
    wx.Frame.__init__(self, None, -1, title=title, size=size)
    self.menubar=self.MenuItems() # create instance of wx.MenuBar
    self.SetMenuBar(self.menubar) # attach Menubar to this frame
    self.Bind(wx.EVT_MENU,self.OnMenu) # Bind - the method of wx.Frame class
    self.Bind(wx.EVT_PAINT, self.OnPaint)
    self.Bind(wx.EVT_TIMER, self.OnTimer)

    self.speed = 100
    self.isAuto = False
    self.timer = wx.Timer(self, MyFrame.ID_TIMER)
    self.timer.Stop()

    self.L = 16
    self.T = 0.5
    self.loop = loop.Looper(L=self.L, T=self.T, periodic_boundary = False)
    self.initUI()

  def initUI(self):
    self.canvas = MyCanvas(self)
    update_button = wx.Button(self, label='update')
    update_button.Bind(wx.EVT_BUTTON, self.OnMCUpdate)

    self.T_txt = wx.StaticText(self, label="T = {}".format(self.T))
    Tslider = wx.Slider(self, value = 50, style=wx.SL_HORIZONTAL)
    Tslider.Bind(wx.EVT_SLIDER, self.OnTempSlider)
    Tslider.SetMin(1)
    Tslider.SetMax(100)

    L_txt = wx.StaticText(self, label = "L = ")

    Ls = map(str, xrange(4,129,2))
    LBox = wx.ComboBox(self, value = "16", choices=Ls, style=wx.CB_READONLY)
    LBox.Bind(wx.EVT_COMBOBOX, self.OnLBox)

    autoBox = wx.CheckBox(self, label = 'auto update')
    autoBox.SetValue(False)
    autoBox.Bind(wx.EVT_CHECKBOX, self.OnAutoBox)

    vbox = wx.BoxSizer(wx.VERTICAL)
    hbox = wx.BoxSizer(wx.HORIZONTAL)

    vbox.Add(self.canvas, 1, flag = wx.EXPAND | wx.ALIGN_CENTER_HORIZONTAL | wx.TOP, border=10)
    vbox.Add(hbox)

    hbox.Add(update_button, 1, flag = wx.EXPAND | wx.RIGHT, border = 5)
    hbox.Add(autoBox, 1, flag = wx.EXPAND)
    hbox.Add(self.T_txt, 1, flag = wx.EXPAND | wx.RIGHT | wx.ALIGN_BOTTOM)
    hbox.Add(Tslider, 1 )
    hbox.Add(L_txt, 1, flag = wx.ALIGN_CENTER_VERTICAL)
    hbox.Add(LBox, 1, flag = wx.EXPAND | wx.LEFT)
    self.SetSizer(vbox)

  def OnSize(self, event):
    self.canvas.ReSize()

  def OnPaint(self, event):
    self.Paint()

  def Paint(self):
    self.canvas.Paint()

  def OnTempSlider(self, event):
    obj = event.GetEventObject()
    val = obj.GetValue()
    self.T = val/100.0
    self.T_txt.SetLabel("T = {}".format(self.T))
    self.loop = loop.Looper(L=self.L, T=self.T, periodic_boundary = False)

  def OnLBox(self, event):
    obj = event.GetEventObject()
    val = obj.GetValue()
    self.L = int(val)
    self.loop = loop.Looper(L=self.L, T=self.T, periodic_boundary = False)

  def OnMCUpdate(self, event):
    self.MCUpdate()

  def MCUpdate(self):
    self.loop.update()
    self.Paint()

  def OnAutoBox(self, event):
    obj = event.GetEventObject()
    val = obj.GetValue()
    if val:
      self.start()
    else:
      self.pause()

  def MenuItems(self):
    """ Create ManuBar/Menu items """
    menubar=wx.MenuBar() # manubar instance
    # File menu
    submenu=wx.Menu()  # menu instance
    submenu.Append(-1,'Quit','Close this panel')
    menubar.Append(submenu,'File')    

    return menubar

  def OnMenu(self,event):
    """ Menu event handler """
    menuid=event.GetId()  # get id number of the clicked menuitem
    item=self.menubar.GetLabel(menuid) # get the label of menuid
    # File menu
    if item == "Quit":
      self.Destroy()   # Destroy the frame to quit the program

  def start(self):
    if not self.isAuto:
      self.isAuto = True
      self.timer.Start(self.speed)
  def pause(self):
    if self.isAuto:
      self.isAuto = False
      self.timer.Stop()
  def OnTimer(self, event):
    if event.GetId() == MyFrame.ID_TIMER:
      self.MCUpdate()
    else:
      event.Skip()


class MyCanvas(glcanvas.GLCanvas):
  def __init__(self, parent):
    super(MyCanvas, self).__init__(parent, 1, attribList=[])
    self.parent = parent
    self.context = glcanvas.GLContext(self)
    self.initialized = False
    self.Bind(wx.EVT_SIZE, self.OnSize)
    self.Bind(wx.EVT_PAINT, self.OnPaint)

  def InitGL(self):
    glClearColor(1.0, 1.0, 1.0, 1.0)

  def OnSize(self, event):
    self.ReSize()

  def ReSize(self):
    w, h = self.GetClientSize()
    glViewport(0, 0, w, h)
    
  def OnPaint(self, event):
    self.Paint()

  def Paint(self):
    self.SetCurrent(self.context)
    cc = colors.ColorConverter()
    if not self.initialized:
      self.InitGL()
      self.initialized = True
    self.ReSize()
    L = self.parent.L
    ts = np.zeros(L)
    xs = np.linspace(-0.9, 0.9, L)
    ss = self.parent.loop.spins[:]
    glClear(GL_COLOR_BUFFER_BIT)
    glBegin(GL_LINES)
    glColor3d(0.0, 0.0, 0.0)
    for op in self.parent.loop.operators:
      ls = self.parent.loop.leftsite(op.bond)
      rs = self.parent.loop.rightsite(op.bond)
      if ss[ls] == 1:
        glColor3d(1.0, 0.0, 0.0)
      else:
        glColor3d(0.0, 0.0, 1.0)
      glVertex2d(xs[ls], 1.8*ts[ls]-0.9)
      ts[ls] = op.time
      glVertex2d(xs[ls], 1.8*ts[ls]-0.9)

      if ss[rs] == 1:
        glColor3d(1.0, 0.0, 0.0)
      else:
        glColor3d(0.0, 0.0, 1.0)
      glVertex2d(xs[rs], 1.8*ts[rs]-0.9)
      ts[rs] = op.time
      glVertex2d(xs[rs], 1.8*ts[rs]-0.9)

      glColor3d(0.0, 0.0, 0.0)
      glVertex2d(xs[ls], 1.8*ts[ls]-0.9)
      glVertex2d(xs[rs], 1.8*ts[rs]-0.9)

      if not op.isdiagonal :
        ss[ls] *= -1
        ss[rs] *= -1

    for x in range(L):
      if ss[x] == 1:
        glColor3d(1.0, 0.0, 0.0)
      else:
        glColor3d(0.0, 0.0, 1.0)
      glVertex2d(xs[x], 1.8*ts[x]-0.9)
      glVertex2d(xs[x], 0.9)

    glEnd()
    glFlush()

if __name__ == '__main__':
    app = wx.App()
    frame = MyFrame(parent = None, id = -1, title = sys.argv[0], size=(480,480))
    frame.Show()
    app.MainLoop()
    app.Destroy()
