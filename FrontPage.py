'''
Created on 27 Apr 2014

@author: NUT67271
'''
import wx,time,os

class Fpanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        #panel = wx.Panel(self, -1)
        self.title_label = wx.StaticText(self, -1, "Title" ) 
        purpose_label = wx.StaticText(self, -1, "Purpose of assessment" )         

        self.author_signature = wx.TextCtrl(self, -1, "Not signed",size=(400,25), style =wx.TE_READONLY)
        self.checker_signature = wx.TextCtrl(self, -1, "Not signed",size=(400,25), style =wx.TE_READONLY)
        
          
        self.author_sgn_btn = wx.Button(self, -1, ' Sign as author ')
        self.author_clr_btn = wx.Button(self, -1, ' Clear ')
        self.chker_sgn_btn = wx.Button(self, -1, ' Sign as checker ')
        self.chker_clr_btn = wx.Button(self, -1, ' Clear ')
        
        self.title = wx.TextCtrl(self, -1, "Catchment @ location", size=(500,50),style=wx.TE_MULTILINE)
        self.purpose = wx.TextCtrl(self, -1, "State purpose of assessment", size=(500,50),style=wx.TE_MULTILINE)
        self.author_notes = wx.TextCtrl(self, -1, "Author's notes", size=(500, 150), style=wx.TE_MULTILINE)
        self.checkers_notes = wx.TextCtrl(self, -1, "Checkers's notes", size=(500, 150), style=wx.TE_MULTILINE)
        
        #  Assign actions to buttons
        self.author_sgn_btn.Bind(wx.EVT_BUTTON, self.authorSign)
        self.author_clr_btn.Bind(wx.EVT_BUTTON, self.authorClear)
        self.chker_sgn_btn.Bind(wx.EVT_BUTTON, self.chkerSign)
        self.chker_clr_btn.Bind(wx.EVT_BUTTON, self.chkerClear)
        
        self.title.Bind(wx.EVT_BUTTON, self.updateTitle)
        
        
        # use gridbagsizer for layout of widgets
        sizer = wx.GridBagSizer(vgap=10, hgap=10)
        sizer.Add(self.title_label, pos=(0,0))
        sizer.Add(purpose_label, pos=(2,0))

        
        sizer.Add(self.author_signature,pos=(5,0),span=(1,2))
        sizer.Add(self.checker_signature,pos=(8,0),span=(1,2))

        
        sizer.Add(self.author_sgn_btn,pos=(5,3), span=(1,1))
        sizer.Add(self.author_clr_btn,pos=(5,4), span=(1,1))
        sizer.Add(self.chker_sgn_btn,pos=(8,3), span=(1,1))
        sizer.Add(self.chker_clr_btn,pos=(8,4), span=(1,1))
        
        sizer.Add(self.title, pos=(1,0), span=(1,5))
        sizer.Add(self.purpose, pos=(3,0), span=(1,5))
        sizer.Add(self.author_notes, pos=(4, 0), span=(1,5))
        sizer.Add(self.checkers_notes, pos=(7, 0), span=(1,5))
        
        
        # use boxsizer to add border around sizer
        border = wx.BoxSizer()
        border.Add(sizer, 0, wx.ALL, 20)
        self.SetSizerAndFit(border)
        self.Fit()
    
    def generateSignature(self):
      username=os.environ['USERNAME']
      sign_time = time.asctime(time.localtime())
      qmed = 0.0
      signature = "Adopted QMED="+str(qmed)+" "+str(username)+" "+str(sign_time)
      return signature
    
    
    def updateTitle(self,event):
      self.store['title']=self.title.GetValue()
    
    def authorSign(self,event):
      signature=self.generateSignature()
      self.author_signature.SetLabel(str(signature))
      self.Refresh()
      self.Update()
    
    def authorClear(self,event):
      signature="Not signed"
      self.author_signature.SetLabel(str(signature))
      self.Refresh()
      self.Update()
    
    def chkerSign(self,evet):
      signature=self.generateSignature()
      self.checker_signature.SetLabel(str(signature))
      self.Refresh()
      self.Update()
    
    def chkerClear(self,event):
      signature="Not signed"
      self.checker_signature.SetLabel(str(signature))
      self.Refresh()
      self.Update()
