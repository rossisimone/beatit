
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

#alpha=[0.05, 0.1, 0.15, 0.20, 0.25, 0.3, 0.35, 0.4, 0.45]
alpha=[0.1, 0.2, 0.3, 0.4]
alpha=[0.25]
tau=[0, 0.1, 0.2, 0.5]
tau=[0.5]
#el=[12, 25, 50, 100, 200, 400]
el=[50, 100, 200, 400, 800, 1600]
#el=[50, 100, 200, 400, 800]
base_folder='/not_backed_up/srossi/dev/BeatIt/build-opt/examples/example_monowave_1D/'
#test_name='mic'
test_name='monowave'
for a in alpha:
    for t in tau:
        for m in el:
            base_folder='/not_backed_up/srossi/dev/BeatIt/build-opt/examples/example_monowave_1D/'
	    base_folder+=test_name
            base_folder+='_TAU'
            base_folder+=str(t)
            base_folder+='_ALPHA'
            base_folder+=str(a)
            base_folder+='_M'
            base_folder+=str(m)
            base_folder+='/'
            file=base_folder+'parameters.exo'
            print base_folder

            save_file=base_folder
            save_file+='tau'
            save_file+=str(t)
            save_file+='_alpha'
            save_file+=str(a)
            save_file+='_M'
            save_file+=str(m)
            save_file+='_t.csv'
            # create a new 'ExodusIIReader'
            parametersexo = ExodusIIReader(FileName=[file])
            parametersexo.ElementVariables = []
            parametersexo.PointVariables = []
            parametersexo.SideSetArrayStatus = []

# get active view
#renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [666, 700]

# destroy renderView1
#Delete(renderView1)
#del renderView1

# get layout
#layout1 = GetLayoutByName("Layout #1")

# close an empty frame
#layout1.Collapse(1)

            # find view
            spreadSheetView1 = FindViewOrCreate('SpreadSheetView1', viewtype='SpreadSheetView')
            # uncomment following to set a specific view size
            spreadSheetView1.ViewSize = [400, 400]

            # set active view
            SetActiveView(spreadSheetView1)

            #Properties modified on parametersexo
            parametersexo.GenerateObjectIdCellArray = 0
            parametersexo.GenerateGlobalElementIdArray = 0
            parametersexo.GenerateGlobalNodeIdArray = 0
            parametersexo.PointVariables = ['activation_times']
            parametersexo.ElementBlocks = ['Unnamed block ID: 0 Type: EDGE2']

            # show data in view
            parametersexoDisplay = Show(parametersexo, spreadSheetView1)
            # trace defaults for the display properties.
            parametersexoDisplay.CompositeDataSetIndex = [2]

            #create a new 'Calculator'
            calculator1 = Calculator(Input=parametersexo)
            # Properties modified on calculator1
            calculator1.ResultArrayName = 'X'
            calculator1.Function = 'coordsX'

            # show data in view
            calculator1Display = Show(calculator1, spreadSheetView1)
            # trace defaults for the display properties.
            calculator1Display.CompositeDataSetIndex = [2]

# hide data in view
#Hide(parametersexo, spreadSheetView1)

            # save data
            SaveData(save_file, proxy=calculator1)

            # destroy calculator1
            Delete(calculator1)
            del calculator1

            # destroy parametersexo
            Delete(parametersexo)
            del parametersexo
#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
