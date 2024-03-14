#!/usr/bin/env Python3  
import os
import shutil

import PySimpleGUI as sg  
import xml.dom.minidom as xmlobj

import subprocess
import sys

sg.ChangeLookAndFeel('GreenTan')  

# ------ coppy files to desti ----- #
def coppyfiles(source,target):

    if len(source) != 0:
        #assert not os.path.isabs(source)
        #target = os.path.join(target, os.path.dirname(source))

	    # create the folders if not already exists
        #os.makedirs(target)

	    # adding exception handling
        try:
           shutil.copy(source, target)
        except IOError as e:
           print("Unable to copy file. %s" % e)
           return str(e)
        except:
           print("Unexpected error:", sys.exc_info())
           return sys.exc_info()

        return ' '+source
    else:
    	return ''

# ------ load params from xml ------#
def readXmlGiveTags(xmlfilename):
	doc = xmlobj.parse(xmlfilename)

	
	proteinname = doc.getElementsByTagName('protein_name')
	#print(proteinname[0].firstChild.nodeValue)
	seqfile = doc.getElementsByTagName('seq_file')
	#print(seqfile[0].firstChild.nodeValue)
	upls = doc.getElementsByTagName('upl_file')
	#print(upls[0].firstChild.nodeValue)
	hfile = doc.getElementsByTagName('hbond_file')
	#print(hfile[0].firstChild.nodeValue)	
	angfile = doc.getElementsByTagName('ang_file')
	#print(angfile[0].firstChild.nodeValue)
	proteinpath = doc.getElementsByTagName('protein_path')
	#print(proteinpath[0].firstChild.nodeValue)

	homit = doc.getElementsByTagName('hydrogen_omission')
	#print(homit[0].firstChild.nodeValue)
	augbounds = doc.getElementsByTagName('aug_bounds')
	#print(augbounds[0].firstChild.nodeValue)	
	etalo = doc.getElementsByTagName('eta_lo')
	#print(etalo[0].firstChild.nodeValue)	
	includeneighbour = doc.getElementsByTagName('include_neighbour')
	#print(includeneighbour[0].firstChild.nodeValue)
	grpexpand = doc.getElementsByTagName('grp_expand')
	#print(grpexpand[0].firstChild.nodeValue)

	proteinname_     = proteinname[0].firstChild.nodeValue 
	seqfile_          = seqfile[0].firstChild.nodeValue 
	upls_             = upls[0].firstChild.nodeValue 
	print(upls[0].firstChild)
	print(hfile[0].firstChild)
	if hfile[0].firstChild:
		hfile_            = hfile[0].firstChild.nodeValue 
	else:
		hfile_ = ''
	if angfile[0].firstChild:
		angfile_          = angfile[0].firstChild.nodeValue 
	else:
		angfile_ = ''
	proteinpath_     = proteinpath[0].firstChild.nodeValue 
	homit_          = homit[0].firstChild.nodeValue 
	if augbounds[0].firstChild:
		augbounds_        = augbounds[0].firstChild.nodeValue 
	else:
		augbounds=''
	etalo_            = etalo[0].firstChild.nodeValue 
	if includeneighbour[0].firstChild:
		includeneighbour_ = includeneighbour[0].firstChild.nodeValue 
	else:
		includeneighbour_ = ''
	if grpexpand[0].firstChild:
		grpexpand_        = grpexpand[0].firstChild.nodeValue
	else :
		grpexpand_ = ''

	#return proteinname[0].firstChild.nodeValue, seqfile[0].firstChild.nodeValue, upls[0].firstChild.nodeValue, hfile[0].firstChild.nodeValue, angfile[0].firstChild.nodeValue, proteinpath[0].firstChild.nodeValue, homit[0].firstChild.nodeValue, augbounds[0].firstChild.nodeValue, etalo[0].firstChild.nodeValue, includeneighbour[0].firstChild.nodeValue, grpexpand[0].firstChild.nodeValue
	return proteinname_ ,seqfile_,upls_,hfile_,angfile_,proteinpath_,homit_,augbounds_,etalo_,includeneighbour_,grpexpand_

def popUpForLoad():	
	text = sg.PopupGetFile('Please enter a file name')	
	proteinname, seqfile, upls, hfile, angfile, proteinpath, homit, augbounds, etalo, includeneighbour, grpexpand = readXmlGiveTags(text)
	upl1=''
	upl2=''
	upl3=''
	upl4=''
	upllist = upls.split(',')

	if len(upllist)==1:
		upl1 = upllist[0]
	if len(upllist)==2:
		upl1 = upllist[0]
		upl2 = upllist[1]
	if len(upllist)==3:
		upl1 = upllist[0]
		upl2 = upllist[1]
		upl3 = upllist[2]		
	if len(upllist)==4:	
		upl1 = upllist[0]
		upl2 = upllist[1]
		upl3 = upllist[2]
		upl4 = upllist[3]	
	
	return proteinname, seqfile, upl1, upl2, upl3, upl4, hfile, angfile, proteinpath, homit, augbounds, etalo, includeneighbour, grpexpand



# ------ writeinto xml file -------
def writeIntoXml(xmlfilename, proteinname, seqfile, upl1, upl2, upl3, upl4, hfile, angfile, proteinpath, homit, augbounds, etalo, includeneighbour, grpexpand):
	upl1_ = upl1.split('/')[-1]
	upls=upl1_
	if len(upl2) != 0:
		upl2_ = upl2.split('/')[-1]
		upls = upls+','+upl2_
	if len(upl3) != 0:
		upl3_ = upl3.split('/')[-1]
		upls = upls+','+upl3_
	if len(upl4) != 0:
		upl4_ = upl4.split('/')[-1]
		upls = upls+','+upl4_

	seqfile_ = seqfile.split('/')[-1]
	hfile_ = hfile.split('/')[-1]
	acofile_= acofile.split('/')[-1]

	f = open(xmlfilename,"w")
	f.write('<set_up>\n')
	f.write('	<inputs>\n')
	f.write('		<protein_name>%s</protein_name>\n'%(proteinname))
	f.write('		<seq_file>%s</seq_file>\n'%(seqfile_))
	f.write('		<upl_file>%s</upl_file>\n'%(upls))
	f.write('		<hbond_file>%s</hbond_file>\n'%(hfile_))
	f.write('		<ang_file>%s</ang_file>\n'%(acofile_))
	f.write('		<protein_path>%s</protein_path>\n'%(proteinpath))
	f.write('	</inputs>\n')
	f.write('	<params>\n')
	f.write('		<hydrogen_omission>%s</hydrogen_omission>\n'%(homit))
	f.write('		<aug_bounds>%s</aug_bounds>\n'%(augbounds))
	f.write('		<break_graph>\n')
	f.write('			<eta_lo>%s</eta_lo>\n'%(etalo))
	f.write('			<eta_hi>1</eta_hi>\n')
	f.write('		</break_graph>\n')
	f.write('		<include_neighbour>%s</include_neighbour>\n'%(includeneighbour))
	f.write('		<multi_expand>\n')
	f.write('			<k_2>3</k_2>\n')
	f.write('			<size_cutoff>0.1</size_cutoff>\n')
	f.write('			<grp_min>250</grp_min>\n')
	f.write('			<incr_min>2</incr_min>\n')
	f.write('		</multi_expand>\n')
	f.write('		<grp_expand>%s</grp_expand>\n'%(grpexpand))
	f.write('	</params>\n')
	f.write('</set_up>\n')


def runCommand(cmd, timeout=None, window=None):
    """ run shell command
    @param cmd: command to execute
    @param timeout: timeout for command execution
    @param window: the PySimpleGUI window that the output is going to (needed to do refresh on)
    @return: (return code from command, command output)
    """

    print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = ''
    for line in p.stdout:
        line = line.decode(errors='replace' if (sys.version_info) < (3, 5) else 'backslashreplace').rstrip()
        output += line
        print(line)
        window.Refresh() if window else None        # yes, a 1-line if, so shoot me

    retval = p.wait(timeout)
    return (retval, output)


# ------ Menu Definition ------ #  
menu_def = [['File', ['Open', 'Save', 'Exit', 'Properties']],  
['Edit', ['Paste', ['Special', 'Normal', ], 'Undo'], ],  
['Help', 'About...'], ]  

# ------ Column Definition ------ #  
column1 = [[sg.Text('Column 1', background_color='#F7F3EC', justification='center', size=(10, 1))],  
   [sg.Spin(values=('Spin Box 1', '2', '3'), initial_value='Spin Box 1')],  
   [sg.Spin(values=('Spin Box 1', '2', '3'), initial_value='Spin Box 2')],  
   [sg.Spin(values=('Spin Box 1', '2', '3'), initial_value='Spin Box 3')]]  

layout = [
#[sg.Menu(menu_def, tearoff=True)],
[sg.Frame(layout=[
[sg.Text('Protein Name', size=(15, 1)), sg.InputText(key='proteinname')],
[sg.Text('Sequence File:', size=(20, 1)), sg.Input(key='seqfile'), sg.FileBrowse()],
[sg.Frame(layout=[
                 [sg.Text('upl file-1', size=(8, 1)), sg.Input(key='upl1'), sg.FileBrowse()],
                 [sg.Text('upl file-2', size=(8, 1)), sg.Input(key='upl2'), sg.FileBrowse()],
                 [sg.Text('upl file-3', size=(8, 1)), sg.Input(key='upl3'), sg.FileBrowse()],
                 [sg.Text('upl file-4', size=(8, 1)), sg.Input(key='upl4'), sg.FileBrowse()],
                 ],title='Upper Bound File(s):')],
[sg.Text('H-Bond File:', size=(20, 1)), sg.Input(key='hfile'), sg.FileBrowse()],
[sg.Text('Angular bound File:', size=(20, 1)), sg.Input(key='acofile'), sg.FileBrowse()],
[sg.Text('Protein Path', size=(15, 1)), sg.InputText(key='proteinpath')]], title='Inputs',title_color='red', relief=sg.RELIEF_SUNKEN,tooltip='set lo only')],
[sg.Frame(layout=[   
[sg.Checkbox('H-omission', default=True, size=(10,1),key='homit'),
sg.Frame(layout=[
                  [sg.Text('Eta lo', size=(5, 1)), sg.Slider(range=(0,1), resolution=.05, size=(20,15), orientation='horizontal', key='etalo'),sg.Text('Eta hi : 1', size=(15, 1))],
                  ], title='Break Graph:', title_color='red', relief=sg.RELIEF_SUNKEN,tooltip='set lo only')],
[sg.Text('Augment bounds', size=(15, 1)), sg.InputText('',size=(5, 1),key='augbounds'),
sg.Text('Include Neibhors:', size=(15, 1)), sg.InputText('',size=(5, 1),key='includeneigh'),
sg.Frame(layout=[
                  [sg.Text('k_2: 3', size=(15, 1))],
                  [sg.Text('Size cutoff: 0.1', size=(15, 1))],
                  [sg.Text('Group min: 250', size=(15, 1))],
                  [sg.Text('Increament min: 2', size=(15, 1))],
                  ], title='Multi-Expand', relief=sg.RELIEF_SUNKEN,tooltip='set lo only'), 
sg.Text('Group Expand:', size=(15, 1)), sg.InputText('',size=(5, 1),key='grpexp')],
], title='Parameters',title_color='red', relief=sg.RELIEF_SUNKEN,tooltip='set lo only')],
[sg.Button("Load", button_color=("white","green"), size=(6, 1)),sg.Submit(tooltip='Click to submit this window'), sg.Cancel(key='_BUTTON_KEY3_')]                  
]



#window = sg.Window('Our algo', layout, default_element_size=(40, 1), grab_anywhere=False)  
window = sg.Window('Window Title', layout, location=(0,0), keep_on_top=True)  

while True:
    event, values = window.Read()

    if event in (None, 'Exit'):
        break

    if event == 'Load':
        proteinname, seqfile, upl1, upl2, upl3, upl4, hfile, angfile, proteinpath, homit, augbounds, etalo, includeneighbour, grpexpand=popUpForLoad()        
        window.FindElement('proteinname').Update(proteinname)
        window.FindElement('seqfile').Update(proteinpath+'/'+seqfile)
        window.FindElement('upl1').Update(proteinpath+'/'+upl1)
        if len(upl2)!=0:
           window.FindElement('upl2').Update(proteinpath+'/'+upl2)
        if len(upl3)!=0:
           window.FindElement('upl3').Update(proteinpath+'/'+upl3)
        if len(upl4)!=0:
           window.FindElement('upl4').Update(proteinpath+'/'+upl4)
        if len(hfile)!=0:
           window.FindElement('hfile').Update(proteinpath+'/'+hfile)
        if len(angfile)!=0:
           window.FindElement('acofile').Update(proteinpath+'/'+angfile)
        window.FindElement('proteinpath').Update(proteinpath)
        if int(homit) == 1:
        	window.FindElement('homit').Update(True)
        else:
            window.FindElement('homit').Update(False)
        window.FindElement('augbounds').Update(augbounds)
        window.FindElement('includeneigh').Update(includeneighbour)
        window.FindElement('grpexp').Update(grpexpand)
        window.FindElement('etalo').Update(etalo)

    if event == 'Submit':        
        xmlfilename = values.get('proteinname')+'.xml'

        proteinname = values.get('proteinname')
        seqfile    = values.get('seqfile')
        upl1 = values.get('upl1')
        upl2 = values.get('upl2')
        upl3 = values.get('upl3')
        upl4 = values.get('upl4')
        hfile        = values.get('hfile') if values.get('hfile') else ''
        acofile      = values.get('acofile') if values.get('acofile') else ''
        proteinpath  = values.get('proteinpath')
        if values.get('homit') == True:
             homit = 1
        else:
        	homit = 0        
        augbounds    = values.get('augbounds')
        etalo        = values.get('etalo')
        includeneigh = values.get('includeneigh')
        grpexp       = values.get('grpexp')


        if os.path.isdir(proteinpath):
        	print('ok')
        else:
        	os.makedirs(proteinpath)


        seqmsg=coppyfiles(seqfile,proteinpath)
        upl1msg=coppyfiles(upl1,proteinpath)
        upl2msg=coppyfiles(upl2,proteinpath)
        upl3msg=coppyfiles(upl3,proteinpath)
        upl4msg=coppyfiles(upl4,proteinpath)
        hfilemsg=coppyfiles(hfile,proteinpath)
        acofilemsg=coppyfiles(acofile,proteinpath)

        runcmd='./genRunScripts_bckend.sh '+proteinname+' > ../log/'+proteinname+'_log.txt'
        print('Run: %s'%(runcmd))

        ## write the run script into the proteinpath
        f = open(proteinpath+'/'+proteinname+'.sh',"w")
        f.write('#!/bin/bash')
        f.write('\n%s'%(runcmd))
        f.close()

        writeIntoXml(proteinpath+'/'+xmlfilename, proteinname, seqfile, upl1, upl2, upl3, upl4, hfile, acofile, proteinpath, homit, augbounds, etalo, includeneigh, grpexp)

        sg.Popup('Title',
      	'Coppied files:',
        'msg'+seqmsg+upl1msg+upl2msg+upl3msg+upl4msg+hfilemsg+acofilemsg,
       'Run: '+runcmd
      	)

        layout2 = [[sg.Output(size=(60,15))],
                     [sg.Button('Run'), sg.Button('Exit')] ]

        window2 = sg.Window( runcmd, layout2)
        while True:             # Event Loop
             event, values = window2.Read()
             # print(event, values)
             if event in (None, 'Exit'):
                 break
             elif event == 'Run':
             	#####p = subprocess.Popen(runcmd,shell=True)

             	runCommand(cmd='tail -1 ../log/'+proteinname+'_log.txt', window=window2)
             	#runCommand(cmd=runcmd+" > /dev/null &")
             	#runCommand(cmd='tail -f ../log/'+proteinname+'_log.txt', window=window2)


             	 #runCommand(cmd='echo '+runcmd, window=window2)
                 #runCommand(cmd=values['_IN_'], window=window2)
        window2.Close()



#sg.Popup('Title',  
# 'The results of the window.',  
# 'The button clicked was "{}"'.format(event),  
# 'The values are', values.get('etalo'))

#print(values)
#print(values.get('etalo'))


