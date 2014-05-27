'''
Created on 27 Apr 2014

@author: NUT67271
'''
import wx
import sys
from wx.lib.mixins.listctrl import CheckListCtrlMixin, ListCtrlAutoWidthMixin, ColumnSorterMixin
import feh_statistical,AMAX

class CheckListCtrl(wx.ListCtrl, CheckListCtrlMixin, ListCtrlAutoWidthMixin):
    def __init__(self, parent):
        wx.ListCtrl.__init__(self, parent, -1,size=(750,200), style=wx.LC_REPORT | wx.SUNKEN_BORDER)
        CheckListCtrlMixin.__init__(self)
        ListCtrlAutoWidthMixin.__init__(self)
        #ColumnSorterMixin.__init__(self,8)

class Fpanel(wx.Panel):
    def __init__(self, parent,p):
      wx.Panel.__init__(self, parent)
      self.p=p
      
      self.data_series = None
      self.amax_data_series = None
      
      self.calc_cds2008_btn = wx.Button(self, -1, ' QMED CDS 2008 ')
      self.calc_cds1999_btn = wx.Button(self, -1, ' QMED CDS 1999 ')
      self.calc_areaOnly_btn = wx.Button(self, -1, ' QMED by area ')
      self.calc_obs_amax_btn = wx.Button(self, -1, ' QMED AMAX OBS ')
      self.calc_obs_pot_btn = wx.Button(self, -1, ' QMED POT OBS ')  
      self.calc_geom_btn = wx.Button(self, -1, ' Channel Width ')
      
      self.user_qmed_label = wx.StaticText(self, -1, "User QMED")
      self.selected_unadj_qmed_label = wx.StaticText(self, -1, "Selected unadjusted QMED")
      self.station_search_distance_label = wx.StaticText(self, -1, "Station search distance")
      self.selected_unadj_qmed = wx.TextCtrl(self, -1, "-", style =wx.TE_READONLY)
      
      self.qmed_cds2008 = wx.TextCtrl(self, -1, "-", style =wx.TE_READONLY)
      self.qmed_cds1999 = wx.TextCtrl(self, -1, "-", style =wx.TE_READONLY)
      self.qmed_areaOnly = wx.TextCtrl(self, -1, "-", style =wx.TE_READONLY)
      self.qmed_obs_amax = wx.TextCtrl(self, -1, "-", style =wx.TE_READONLY)
      self.qmed_obs_pot = wx.TextCtrl(self, -1, "-", style =wx.TE_READONLY)
      self.qmed_geom = wx.TextCtrl(self, -1, "-", style =wx.TE_READONLY) 
      self.qmed_user = wx.TextCtrl(self, -1, "-") 
      
      self.rb1 = wx.RadioButton(self, -1, '', style=wx.RB_GROUP)
      self.rb2 = wx.RadioButton(self, -1, '')
      self.rb3 = wx.RadioButton(self, -1, '')
      self.rb4 = wx.RadioButton(self, -1, '')
      self.rb5 = wx.RadioButton(self, -1, '')
      self.rb6 = wx.RadioButton(self, -1, '')
      self.rb7 = wx.RadioButton(self, -1, '')
      
      self.Bind(wx.EVT_RADIOBUTTON, self.SetVal, id=self.rb1.GetId())
      self.Bind(wx.EVT_RADIOBUTTON, self.SetVal, id=self.rb2.GetId())
      self.Bind(wx.EVT_RADIOBUTTON, self.SetVal, id=self.rb3.GetId())
      self.Bind(wx.EVT_RADIOBUTTON, self.SetVal, id=self.rb4.GetId())
      self.Bind(wx.EVT_RADIOBUTTON, self.SetVal, id=self.rb5.GetId())
      self.Bind(wx.EVT_RADIOBUTTON, self.SetVal, id=self.rb6.GetId())
      self.Bind(wx.EVT_RADIOBUTTON, self.SetVal, id=self.rb7.GetId())

      self.station_search_distance = wx.TextCtrl(self, -1, "25.0")
      self.refresh_stations_btn = wx.Button(self, -1, ' Refresh stations ')
      self.add_station_btn = wx.Button(self, -1, ' Add station ')  
      self.remove_station_btn = wx.Button(self, -1, ' Remove station ')
      self.clear_user_stations_btn = wx.Button(self, -1, ' Clear add/removes ')
      
      self.refresh_stations_btn.Bind(wx.EVT_BUTTON, self.OnRefreshStations)
      self.add_station_btn.Bind(wx.EVT_BUTTON, self.OnAddStation)
      self.remove_station_btn.Bind(wx.EVT_BUTTON, self.OnRemoveStation)
      self.clear_user_stations_btn.Bind(wx.EVT_BUTTON, self.OnClearUser)
      
      self.distance_decay_update= wx.RadioButton(self, -1, "Decaying with distance", style=wx.RB_GROUP)
      self.direct_transfer_update  = wx.RadioButton(self, -1, 'Direct transfer')
      self.dont_update  = wx.RadioButton(self, -1, "Don't update")
      
      self.Bind(wx.EVT_RADIOBUTTON, self.SetUpdate, id=self.dont_update.GetId())
      self.Bind(wx.EVT_RADIOBUTTON, self.SetUpdate, id=self.distance_decay_update.GetId())
      self.Bind(wx.EVT_RADIOBUTTON, self.SetUpdate, id=self.direct_transfer_update.GetId())
      
      self.locally_adjusted_qmed_label = wx.StaticText(self, -1, "Locally adjusted QMED")
      self.locally_adjusted_qmed = wx.TextCtrl(self, -1, "-", style =wx.TE_READONLY)
      
      sizer = wx.GridBagSizer(vgap=5, hgap=10)
      sizer.Add(self.calc_cds2008_btn, pos=(0,0))
      sizer.Add(self.calc_cds1999_btn, pos=(1,0))
      sizer.Add(self.calc_areaOnly_btn, pos=(2,0))
      sizer.Add(self.calc_obs_amax_btn, pos=(3,0))
      sizer.Add(self.calc_obs_pot_btn, pos=(4,0))
      sizer.Add(self.calc_geom_btn, pos=(5,0))
      
      sizer.Add(self.user_qmed_label,pos=(6,0))
      sizer.Add(self.selected_unadj_qmed_label,pos=(7,0))
      
      sizer.Add(self.qmed_cds2008, pos=(0,1))
      sizer.Add(self.qmed_cds1999, pos=(1,1))
      sizer.Add(self.qmed_areaOnly, pos=(2,1))
      sizer.Add(self.qmed_obs_amax, pos=(3,1))
      sizer.Add(self.qmed_obs_pot, pos=(4,1))
      sizer.Add(self.qmed_geom, pos=(5,1))
      sizer.Add(self.qmed_user, pos=(6,1))
      sizer.Add(self.selected_unadj_qmed, pos=(7,1))
      
      sizer.Add(self.rb1, pos=(0,2))
      sizer.Add(self.rb2, pos=(1,2))
      sizer.Add(self.rb3, pos=(2,2))
      sizer.Add(self.rb4, pos=(3,2))
      sizer.Add(self.rb5, pos=(4,2))
      sizer.Add(self.rb6, pos=(5,2))
      sizer.Add(self.rb7, pos=(6,2))
      
      self.table = wx.Panel(self, -1)
      self.list = CheckListCtrl(self.table)
      self.list.InsertColumn(0, 'STATION')
      self.list.InsertColumn(1, 'DISTANCE')
      self.list.InsertColumn(2, 'AREA')
      self.list.InsertColumn(3, 'SAAR')
      self.list.InsertColumn(5, 'BFIHOST')
      self.list.InsertColumn(5, 'FARL')
      self.list.InsertColumn(6, 'QMED OBS')
      self.list.InsertColumn(7, 'QMED ERROR')
      self.list.InsertColumn(8, 'ASG')
      
      self.findStations()
      self.refreshStations()
      
      sizer.Add(self.table, pos=(9,0),span=(1,6),flag=wx.EXPAND)
      sizer.Add(self.station_search_distance_label, pos=(10,0),span=(1,1))
      sizer.Add(self.station_search_distance, pos=(10,1),span=(1,1))
      sizer.Add(self.refresh_stations_btn, pos=(10,2),span=(1,1))
      sizer.Add(self.add_station_btn, pos=(10,3),span=(1,1))
      sizer.Add(self.remove_station_btn, pos=(10,4),span=(1,1))
      sizer.Add(self.clear_user_stations_btn, pos=(10,5),span=(1,1))      
      sizer.Add(self.dont_update, pos=(11,0),span=(1,2))
      sizer.Add(self.distance_decay_update, pos=(12,0),span=(1,2))
      sizer.Add(self.direct_transfer_update, pos=(13,0),span=(1,2))
      
      sizer.Add(self.locally_adjusted_qmed_label, pos=(14,0),span=(1,1))
      sizer.Add(self.locally_adjusted_qmed, pos=(14,1),span=(1,1))
      
      #  Assign actions to buttons
      self.calc_cds2008_btn.Bind(wx.EVT_BUTTON, self.calcCds2008)
      self.calc_cds1999_btn.Bind(wx.EVT_BUTTON, self.calcCds1999)
      self.calc_areaOnly_btn.Bind(wx.EVT_BUTTON, self.calcAreaOnly)
      self.calc_obs_amax_btn.Bind(wx.EVT_BUTTON, self.calcAMAX)
      self.calc_obs_pot_btn.Bind(wx.EVT_BUTTON, self.calcPOT)
      self.calc_geom_btn.Bind(wx.EVT_BUTTON, self.calcChannel)
      
      border = wx.BoxSizer()
      border.Add(sizer, 0, wx.ALL, 20)
      self.SetSizerAndFit(border)
      self.Fit()
      
    def findStations(self):
          self.qmed_databse = "C:\\Users\\nut67271\\workspace\\RandomNumberGenerator\\qmed_db.csv"
          f = open(self.qmed_databse,'r')
          lines = f.readlines()
          
          headers = lines[0].split(',')
          station_cds_index=list()
          i=0
          
          for header in headers:
            label = header.replace(' ','').lower()
            station_cds_index.append(label)
            i=i+1
          
          self.stations = list()
          
          
          for line in lines[1:]:
            entries = line.split(',')
            
            station_cds=dict()
            i=0
            for field in station_cds_index:
              try:
                station_cds[field]=float(entries[i])
              except:
                station_cds[field]=entries[i]
              i=i+1
            
            
            station_centroid_x = float(entries[6])
            station_centroid_y = float(entries[7])    
            self.location_centroid_x =float(self.p.centroid_x.GetValue())
            self.location_centroid_y = float(self.p.centroid_y.GetValue())
            
            station_cds['distance'] = 0.001*(((station_cds['centroid_x']-self.location_centroid_x)**2.0+(station_cds['centroid_y']-self.location_centroid_y))**2.0)**0.5        
            station_cds['user_added'] = None
            
            self.stations.append(station_cds)
       
    def OnClearUser(self,event):
        i=0
        for station_cds in self.stations:
          self.stations[i]['user_added'] = None
          i=i+1
        self.list.DeleteAllItems()     
        self.refreshStations()
        self.Refresh()
        self.Update()  
              
    def OnRefreshStations(self,event): 
      self.list.DeleteAllItems()     
      self.refreshStations()
      self.Refresh()
      self.Update()
      
    def refreshStations(self):
      self.search_distance =float(self.station_search_distance.GetValue())
      for station_cds in self.stations:
            station_cds['distance'] = 0.001*(((station_cds['centroid_x']-self.location_centroid_x)**2.0+(station_cds['centroid_y']-self.location_centroid_y)**2.0))**0.5
            if station_cds['distance'] > self.search_distance and station_cds['user_added'] != True:
              pass
            elif station_cds['user_added'] == False:
              pass
            else:
              if bool(self.rb1.GetValue()) == True:
                qmed_cds = feh_statistical.qmed_cds2008(station_cds['dtm_area'],station_cds['saar'],station_cds['farl'],station_cds['bfihost'])
              if bool(self.rb2.GetValue()) == True:
                qmed_cds = feh_statistical.qmed_cds1999(station_cds['dtm_area'],station_cds['saar'],station_cds['farl'],station_cds['sprhost'],station_cds['bfihost'])
              if bool(self.rb3.GetValue()) == True:
                qmed_cds = feh_statistical.area_based_qmed(station_cds['dtm_area'])
              if bool(self.rb4.GetValue()) == True:
                qmed_cds = station_cds['qmed_obs']
              if bool(self.rb5.GetValue()) == True:
                qmed_cds = station_cds['qmed_obs']
              if bool(self.rb6.GetValue()) == True:
                qmed_cds = station_cds['qmed_obs']
              if bool(self.rb7.GetValue()) == True:
                qmed_cds = station_cds['qmed_obs']
              
              try:
                station_cds['qmed_error'] = station_cds['qmed_obs']/qmed_cds
              except:
                qmed_error = 0.0
              station_cds['asg'] = feh_statistical.calc_asg(station_cds['distance'])
  
              index = self.list.InsertStringItem(sys.maxint, str(int(station_cds['station'])))
              self.list.SetStringItem(index, 1, str(station_cds['distance']))
              self.list.SetStringItem(index, 2, str(station_cds['dtm_area']))
              self.list.SetStringItem(index, 3, str(station_cds['saar']))
              self.list.SetStringItem(index, 4, str(station_cds['bfihost']))
              self.list.SetStringItem(index, 5, str(station_cds['farl']))
              self.list.SetStringItem(index, 6, str(station_cds['qmed_obs']))
              self.list.SetStringItem(index, 7, str(station_cds['qmed_error']))
              self.list.SetStringItem(index, 8, str(station_cds['asg']))

    def OnAddStation(self,event):
        stations_txt = list()
        stations_ids =list()
        for station_cds in self.stations:
          stations_txt.append(str(station_cds['station']))
          stations_ids.append(station_cds['station'])
        dlg = wx.SingleChoiceDialog(self, 'Add station', 'Which station do you want to add?', stations_txt, wx.CHOICEDLG_STYLE)
        if dlg.ShowModal() == wx.ID_OK:
          i = 0
          for station_cds in self.stations:
            if dlg.GetStringSelection() == str(station_cds['station']):
              self.stations[i]['user_added'] = True
            i=i+1
          dlg.Destroy()
          self.list.DeleteAllItems()     
          self.refreshStations()
          self.Refresh()
          self.Update()


    
    def OnRemoveStation(self,event):
        stations_txt = list()
        stations_ids =list()
        stations_count = self.list.GetItemCount()
        for i in range(stations_count):
          stations_txt.append(str(self.list.GetItem(i,0).GetText()))
          stations_ids.append(self.list.GetItem(i,0).GetText())
        
        dlg = wx.SingleChoiceDialog(self, 'Remove station', 'Which station do you want to remove?', stations_txt, wx.CHOICEDLG_STYLE)
        if dlg.ShowModal() == wx.ID_OK:
          i = 0
          for station_cds in self.stations:
            if dlg.GetStringSelection() == str(int(station_cds['station'])):
              self.stations[i]['user_added'] = False
            i=i+1
          dlg.Destroy()
          self.list.DeleteAllItems()     
          self.refreshStations()
          self.Refresh()
          self.Update()
    
    def SetUpdate(self,event):
      qmed_cds = float(self.selected_unadj_qmed.GetValue())
      f=list()
      stations = self.list.GetItemCount()
      
      if bool(self.distance_decay_update.GetValue()) == True: 
        for i in range(stations):
          if bool(self.list.IsChecked(i)) == True:
            d = float(self.list.GetItem(i,1).GetText())
            qmed_error = float(self.list.GetItem(i,7).GetText())
            asg = feh_statistical.calc_asg(d)
            f.append(feh_statistical.distance_decaying_weighting_factor(qmed_error, asg))
        if len(f) == 0:
          updated_qmed = qmed_cds
        else:
          updated_qmed = qmed_cds * sum(f)/len(f)
        self.locally_adjusted_qmed.SetLabel(str(updated_qmed))
        
      if bool(self.direct_transfer_update.GetValue()) == True:
        for i in range(stations):
          if bool(self.list.IsChecked(i)) == True:
            qmed_error = float(self.list.GetItem(i,7).GetText())
            f.append(qmed_error)
        if len(f) == 0:
          updated_qmed = qmed_cds
        else:
          updated_qmed = qmed_cds * sum(f)/len(f)
        self.locally_adjusted_qmed.SetLabel(str(updated_qmed))
      
      if bool(self.dont_update.GetValue()) == True:
        self.locally_adjusted_qmed.SetLabel(str(qmed_cds))
      self.Refresh()
      self.Update()
    
    def SetVal(self, event):
      if bool(self.rb1.GetValue()) == True:
        self.selected_unadj_qmed.SetLabel(str(self.qmed_cds2008.GetValue()))
      if bool(self.rb2.GetValue()) == True:
        self.selected_unadj_qmed.SetLabel(str(self.qmed_cds1999.GetValue()))
      if bool(self.rb3.GetValue()) == True:
        self.selected_unadj_qmed.SetLabel(str(self.qmed_areaOnly.GetValue()))
      if bool(self.rb4.GetValue()) == True:
        self.selected_unadj_qmed.SetLabel(str(self.qmed_amax.GetValue()))
      if bool(self.rb5.GetValue()) == True:
        self.selected_unadj_qmed.SetLabel(str(self.qmed_pot.GetValue()))
      if bool(self.rb6.GetValue()) == True:
        self.selected_unadj_qmed.SetLabel(str(self.qmed_geom.GetValue()))
      if bool(self.rb7.GetValue()) == True:
        self.selected_unadj_qmed.SetLabel(str(self.qmed_user.GetValue()))
      
      self.list.DeleteAllItems() 
      self.refreshStations()
      
      self.Refresh()
      self.Update()
    
    def calcCds2008(self,event):
      try:
        saar = float(self.p.saar.GetValue())
        carea = float(self.p.carea.GetValue())
        bfihost = float(self.p.bfihost.GetValue())
        farl = float(self.p.farl.GetValue())
        qmed = feh_statistical.qmed_cds2008(carea,saar,farl,bfihost)
        if bool(self.rb1.GetValue()) == True:
          self.selected_unadj_qmed.SetLabel(str(qmed))
      except:
        qmed="invalid cds"
      self.qmed_cds2008.SetLabel(str(qmed))
      self.Refresh()
      self.Update()
      
    def calcCds1999(self,event):
      import math
      try:
        bfihost = float(self.p.bfihost.GetValue())
        sprhost = float(self.p.sprhost.GetValue())
        farl = float(self.p.farl.GetValue())
        saar = float(self.p.saar.GetValue())
        carea = float(self.p.carea.GetValue())
        
        qmed = feh_statistical.qmed_cds1999(carea,saar,farl,sprhost,bfihost)
        
        if bool(self.rb2.GetValue()) == True:
          self.selected_unadj_qmed.SetLabel(str(qmed))
      except:
        qmed="invalid cds"
      self.qmed_cds1999.SetLabel(str(qmed))
      self.Refresh()
      self.Update()    
      
    def calcAreaOnly(self,event):
      try:
        carea = float(self.p.carea.GetValue())        
        
        qmed = feh_statistical.area_based_qmed(carea)
        
        if bool(self.rb3.GetValue()) == True:
          self.selected_unadj_qmed.SetLabel(str(qmed))
      except:
        qmed="invalid cds"
      self.qmed_areaOnly.SetLabel(str(qmed))
      self.Refresh()
      self.Update()   

    def calcAMAX(self,event):
      '''
      To be pop up window for user to load/see the amax series, possibly with graph showing trend in data
      '''
      AMAX.AmaxFrame(self).Show()
      #self.qmed_obs_amax.SetLabel("Not implemented")
      self.Refresh()
      self.Update()   
    
    def calcPOT(self,event):
      '''
      To be pop up window for user to load/see flow series and set a POT threshold, possibly with graph showing trend in data
      Possibly a means of assisting the user in selecting a POT threshold and a means of detecting clustering
      '''
      self.qmed_obs_amax.SetLabel("Not implemented")
      self.Refresh()
      self.Update()   
    
    
    def calcChannel(self,event):
      try:
        width = float(self.p.chnl_width.GetValue())
        qmed=feh_statistical.channel_width_qmed(width)

        if bool(self.rb6.GetValue()) == True:
          self.selected_unadj_qmed.SetLabel(str(qmed))
      except:
        qmed="invalid cds"
      self.qmed_geom.SetLabel(str(qmed))
      self.Refresh()
      self.Update()  
      