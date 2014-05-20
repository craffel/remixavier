# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

'''
GUI App for running remixavier
'''

# <codecell>

import sys, os
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import numpy as np
import librosa
import estimate

# <codecell>

class Remixavier( QThread ) :
    '''
    Class which runs the subtraction algorithm, in a thread.
    '''
    percentDoneSignal = pyqtSignal( int )
    doneSignal = pyqtSignal( bool )
    statusSignal = pyqtSignal( str )

    def __init__( self ):
        '''
        Get ready to run - initialize thread and GUI elements
        '''
        super(Remixavier, self).__init__()
        self.mix_file = None
        self.source_file = None
        self.wiener_threshold = None
    
    def loadNewValues( self, mix_file, source_file, wiener_threshold ):
        '''
        Load in a new parameter setting from the GUI
        
        Input:
            mix_file - path to mixure .wav file
            source_file - path to source .wav file
            wiener_threshold - threshold for wiener filtering, in dB
        '''
        self.mix_file = mix_file
        self.source_file = source_file
        self.wiener_threshold = wiener_threshold
    
    def run( self ):
        # Initialize signals
        self.percentDoneSignal.emit( 0 )
        percent_scale = 1000.0/5
        self.doneSignal.emit( 0 )
        self.statusSignal.emit( "" )
        # Load in audio data
        self.statusSignal.emit( "Loading {}".format( os.path.split(self.mix_file)[1] ) )
        mix, self.fs = librosa.load( self.mix_file, sr=None )
        self.percentDoneSignal.emit( 1*percent_scale )
        self.statusSignal.emit( "Loading {}".format( os.path.split(self.source_file)[1] ) )
        source, self.fs = librosa.load( self.source_file, sr=self.fs )
        self.percentDoneSignal.emit( 2*percent_scale )
        # Fix any gross timing offset
        self.statusSignal.emit( "Aligning..." )
        mix, source = estimate.align( mix, source, self.fs )
        self.percentDoneSignal.emit( 3*percent_scale )
        self.statusSignal.emit( "Subtracting..." )
        source = estimate.reverse_channel(mix, source)
        mix, source = estimate.pad(mix, source)
        self.percentDoneSignal.emit( 4*percent_scale )
        self.statusSignal.emit( "Enhancing..." )
        self.subtracted = estimate.wiener_enhance( mix - source, source, self.wiener_threshold )
        self.percentDoneSignal.emit( 5*percent_scale )
        self.doneSignal.emit( 1 )

# <codecell>

# The GUI app
class AppForm(QMainWindow):
    
    # Initialize
    def __init__(self, parent=None):
        # Initialize window
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('Remixavier')
                
        # Create GUI elements
        self.createMenu()
        self.createMainFrame()
        self.createStatusBar()
    
        # Instance of remixaiver
        self.remixavierInstance = Remixavier()
        self.remixavierInstance.percentDoneSignal.connect( self.updateProgressBar )
        self.remixavierInstance.doneSignal.connect( self.remixavierFinished )
        self.remixavierInstance.statusSignal.connect( self.updateStatusLabel )
    
    def updateProgressBar( self, percent ):
        self.progressBar.setValue( percent )

    def remixavierFinished( self, finished ):
        if finished:
            librosa.output.write_wav( self.output_file, self.remixavierInstance.subtracted, self.remixavierInstance.fs )
            self.statusBar().showMessage( 'Done.', 0 )
            self.progressBar.setValue( 1000 )
            self.startStopButton.setEnabled( 1 )
    
    def updateStatusLabel( self, text ):
        self.statusBar().showMessage( text, 0 )
    
    # Convenience function for creating sliders
    def createSlider( self, name, minimum, maximum, default, scale ):
        slider = QSlider( Qt.Horizontal, self )
        slider.setFocusPolicy( Qt.NoFocus )
        slider.setMinimum( minimum*scale )
        slider.setMaximum( maximum*scale )
        slider.setValue( default*scale )
        slider.connect( slider, SIGNAL('valueChanged(int)'), self.updateLabels )
        return QLabel( name ), slider

    # Create main GUI window
    def createMainFrame( self ):
        # Initialize main widget
        self.mainFrame = QWidget()
        # Make widgets less squished
        self.mainFrame.setMinimumWidth( 400 )
        
        # Parameter sliders
        self.wiener_thresholdLabel, self.wiener_thresholdSlider = self.createSlider( "Wiener Threshold", -10, 10, 0, 1 )
        
        # Update the slider labels with their default values
        self.updateLabels()
        
        # Progress bar for displaying progress
        self.progressBar = QProgressBar()
        self.progressBar.setMinimum(1)
        self.progressBar.setMaximum(1000)
        
        # Open button for starting analysis
        self.startStopButton = QPushButton( "&Start" )
        self.connect(self.startStopButton, SIGNAL( 'clicked()'), self.startStopButtonClicked )
        
        # VBox for snippet length controls
        parametersVBox = QVBoxLayout()
        parametersVBox.addWidget( self.wiener_thresholdLabel )
        parametersVBox.addWidget( self.wiener_thresholdSlider )
        
        # Box for status bar and open button
        statusHBox = QHBoxLayout()
        statusHBox.addWidget( self.progressBar )
        statusHBox.addWidget( self.startStopButton )
        
        # Add all to main layout
        mainLayout = QVBoxLayout()
        mainLayout.addLayout( parametersVBox )
        mainLayout.addLayout( statusHBox )
        
        # Add the main layout to the frame
        self.mainFrame.setLayout( mainLayout )
        # Main widget is the frame
        self.setCentralWidget( self.mainFrame )
        self.setFixedSize( self.mainFrame.minimumSize() )
    
    # Updates all slider labels with their values
    def updateLabels( self ):
        # Convert lengths to float (in seconds)
        self.wiener_thresholdLabel.setText( "Wiener Threshold: {}".format( self.wiener_thresholdSlider.value() ) )

    # Create menus
    def createMenu( self ):
        # Open item
        openFile = QAction('Start Analysis', self)
        # Shortcut
        openFile.setShortcut('Ctrl+O')
        # Connect open action to show dialog
        self.connect(openFile, SIGNAL('triggered()'), self.showDialog)
        # Exit (like quit)
        exitAction = QAction('Exit - like quitting', self)
        # Shortcut
        exitAction.setShortcut('Ctrl+E')
        self.connect(exitAction, SIGNAL('triggered()'), self.exit)
        
        # Create menubar
        menubar = self.menuBar()
        # Add "file" menu
        fileMenu = menubar.addMenu('&File')
        # Add open and exit
        fileMenu.addAction(openFile)
        fileMenu.addAction(exitAction)
    
    def exit(self):
        if self.remixavierInstance.isRunning():
            self.remixavierInstance.quit()
        sys.exit()
    
    # Callback for when the start/stop button (aka "open" or "cancel") is clicked
    def startStopButtonClicked( self ):
        if self.remixavierInstance.isRunning():
            # I don't know how to stop the thread... this doesn't work.
            self.remixavierInstance.quit()
            self.startStopButton.setEnabled( 1 )
        else:
            self.showDialog()
    
    
    # Show a dialog box for a set of audio files and start running remixavier
    def showDialog(self):
        # Get mixture file name from dialog box
        mix_file = str(QFileDialog.getOpenFileName(self, "Select a mixture file", ".", "Audio Files (*.mp3 *.wav)"))
        # If the user didn't hit cancel
        if mix_file is not '':
            # ... source file
            source_file = str(QFileDialog.getOpenFileName(self, "Select a source file", ".", "Audio Files (*.mp3 *.wav)"))
            if source_file is not '':
                # Where to save the resulting file
                self.output_file = str(QFileDialog.getSaveFileName(self, "Save file", ".", "Audio Files (*.mp3 *.wav)"))
                if self.output_file is not '':
                    self.progressBar.reset()
                    self.remixavierInstance.loadNewValues( mix_file, source_file, self.wiener_thresholdSlider.value() )
                    self.remixavierInstance.start()
                    self.startStopButton.setEnabled( 0 )
    
    def createStatusBar(self):
        self.statusText = QLabel('Click "Start" to begin.')
        self.statusBar().addWidget(self.statusText, 1)


def main():
    app = QApplication(sys.argv)
    form = AppForm()
    form.show()
    form.raise_()
    app.exec_()

if __name__ == "__main__":
    main()

