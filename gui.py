from datetime import datetime
import threading
import time
import os
import glob
import fnmatch
import subprocess as sp
from PIL import Image as imopen
from PIL import ImageTk


try:
    from Tkinter import *
    from tkFileDialog import askopenfilename
    from tkFileDialog import askdirectory
except ImportError:
    from tkinter import *
    from tkinter.filedialog import askopenfilename
    from tkinter.filedialog import askdirectory


class nsdrpGui():
    def __init__(self):
        self.gui = Tk()
        self.gui.title("NSDRP GUI")

        try:
            import ktl
            self.inputDir = ktl.cache('nirspec', 'outdir').read()
            self.batchInVar = StringVar(value=self.inputDir)
        except BaseException:
            self.inputDir = os.getcwd()
            self.batchInVar = StringVar(value=self.inputDir)

        self.dataProducts = []
        self.checkedOptions = []
        self.modes = [1]

        self.mainFrame = Frame(self.gui, width=500, height=500)
        self.radioButtonFrame = Frame(self.gui)
        # , width=600, height=500)
        self.plotFrame = Frame(self.gui, bg='lavender')

        self.radioButtonFrame.grid(row=0, column=0, sticky='NW')
        self.mainFrame.grid(row=1, column=0, sticky='N')
        self.plotFrame.grid(row=1, column=1)

        self.missingFile = ''
        self.objectVar = StringVar()
        self.flatVar = StringVar()
        self.aVar = StringVar()
        self.bVar = StringVar()
        self.batchVar = StringVar()
        self.pairFlatVar = StringVar()
        self.modeVar = IntVar()
        self.var1 = StringVar()
        self.var2 = StringVar()
        self.productVar = StringVar()

        if datetime.now().day != datetime.utcnow().day:
            self.utDate = datetime.utcnow().strftime('%Y%m%d')
        else:
            self.utDate = datetime.utcnow().strftime('%Y%m%d')

        self.defaultOutDir = '/h/disks/scratch/nsdrp'
        if os.path.isdir(self.defaultOutDir):
            if os.path.isdir(self.defaultOutDir + '/' + self.utDate):
                self.outDirVar = StringVar(
                    value=self.defaultOutDir + '/' + self.utDate)
                self.batchOutVar = self.outDirVar
            else:
                self.outDirVar = StringVar(value=os.makedirs(
                    self.defaultOutDir + '/' + self.utDate))
                self.batchOutVar = self.outDirVar
        else:
            self.outDirVar = StringVar(value=os.getcwd())
            self.batchOutVar = self.outDirVar

        self.previewsListBox = Listbox(self.mainFrame, selectmode='SINGLE')

        self.make_labels()
        self.make_entries()
        self.make_dropboxes()
        self.make_buttons()

        self.gui.mainloop()

    def make_labels(self):
        self.statusLabel = Label(self.mainFrame, text='')
        self.fileLabel = Label(self.mainFrame, text='')
        self.labelsFont = ('times', 14, 'bold')

        self.objectFileLabel = Label(
            self.mainFrame,
            text='Object filename:',
            font=self.labelsFont)
        self.objectFileLabel.grid(row=1, column=0, padx=5, sticky='W')

        self.flatLabel = Label(
            self.mainFrame,
            text='Flat filename:',
            font=self.labelsFont)
        self.flatLabel.grid(row=3, column=0, padx=5, sticky='W')

        self.aFrameFileLabel = Label(
            self.mainFrame,
            text='A Frame Filename:',
            font=self.labelsFont)
        self.bFrameFileLabel = Label(
            self.mainFrame,
            text='B Frame Filename:',
            font=self.labelsFont)
        self.pairFlatLabel = Label(
            self.mainFrame,
            text='Flat Filename',
            font=self.labelsFont)

        self.batchInputLabel = Label(
            self.mainFrame,
            text='Input Directory:',
            font=self.labelsFont)
        self.batchOutputLabel = Label(
            self.mainFrame,
            text='Output Directory:',
            font=self.labelsFont)

        self.cmdOptionsLabel = Label(
            self.mainFrame,
            text='Options:',
            font=self.labelsFont)

        self.outDirLabel = Label(
            self.mainFrame,
            text='Output Directory:',
            font=self.labelsFont)
        self.outDirLabel.grid(row=5, column=0, sticky='W')

        self.cmdOptionsLabel.grid(row=7, column=0, sticky='W')

    def make_entries(self):
        self.outDirEntry = Entry(
            self.mainFrame,
            textvariable=self.outDirVar,
            width=len(
                self.outDirVar.get()))
        self.outDirEntry.grid(row=6, column=0, sticky='W')

        self.batchInputEntry = Entry(
            self.mainFrame,
            textvariable=self.batchInVar,
            width=len(
                self.batchInVar.get()))
        self.batchOutputEntry = Entry(
            self.mainFrame,
            textvariable=self.batchOutVar,
            width=len(
                self.batchOutVar.get()))

    def make_dropboxes(self):
        self.files = glob.glob(self.inputDir + '/*.fits*')
        if self.files:
            self.fitsFiles = []
            for file in self.files:
                self.fitsFiles.append(file.split('/')[-1])
            self.width = max(len(string) for string in self.fitsFiles)
        else:
            self.missingFile = 'No FITS files found!'
            self.objectVar = StringVar(value=self.missingFile)
            self.flatVar = StringVar(value=self.missingFile)
            self.aVar = StringVar(value=self.missingFile)
            self.bVar = StringVar(value=self.missingFile)
            self.pairFlatVar = StringVar(value=self.missingFile)
            self.width = len(self.objectVar.get())
            self.fitsFiles = ['']

        self.objectFileDropBox = OptionMenu(
            self.mainFrame, self.objectVar, *self.fitsFiles)
        self.objectFileDropBox.configure(width=self.width)
        self.objectFileDropBox.grid(row=2, column=0, sticky='W')

        self.aFrameDropBox = OptionMenu(
            self.mainFrame, self.aVar, *self.fitsFiles)
        self.aFrameDropBox.configure(width=self.width)

        self.bFrameDropBox = OptionMenu(
            self.mainFrame, self.bVar, *self.fitsFiles)
        self.bFrameDropBox.configure(width=self.width)

        self.flatFileDropBox = OptionMenu(
            self.mainFrame, self.flatVar, *self.fitsFiles)
        self.flatFileDropBox.configure(width=self.width)
        self.flatFileDropBox.grid(row=4, column=0, sticky='W')

        self.pairFlatDropBox = OptionMenu(
            self.mainFrame, self.pairFlatVar, *self.fitsFiles)
        self.pairFlatDropBox.configure(width=self.width)

#       self.gui.after(5000, self.make_dropboxes)

    def update_outdir(self):
        self.defaultDir = self.outDirEntry.get()
        self.outDirName = askdirectory()

        if self.outDirName:
            self.outDirEntry.delete(0, END)
            self.outDirEntry.insert(0, self.outDirName)
            self.outDirEntry.configure(width=len(self.outDirName))
        else:
            self.outDirEntry.delete(0, END)
            self.outDirEntry.insert(0, self.defaultDir)
            self.outDirEntry.configure(width=len(self.defaultDir))

        if self.modeVar.get() == 1:
            self.outDirEntry.grid(row=6, column=0, sticky='W')
        else:
            self.outDirEntry.grid(row=8, column=0, sticky='W')

    def update_batch_indir(self):
        self.defaultInDir = self.batchInputEntry.get()
        self.inDirName = askdirectory()

        if self.inDirName:
            self.batchInputEntry.delete(0, END)
            self.batchInputEntry.insert(0, self.inDirName)
            self.batchInputEntry.configure(width=len(self.inDirName))
        else:
            self.batchInputEntry.delete(0, END)
            self.batchInputEntry.insert(0, self.defaultInDir)
            self.batchInputEntry.configure(width=len(self.defaultInDir))

        self.batchInputEntry.grid(row=2, column=0, sticky='W')

    def update_batch_outdir(self):
        self.defaultOutDir = self.batchOutputEntry.get()
        self.outDirName = askdirectory()

        if self.outDirName:
            self.batchOutputEntry.delete(0, END)
            self.batchOutputEntry.insert(0, self.outDirName)
            self.batchOutputEntry.configure(width=len(self.outDirName))
        else:
            self.batchOutputEntry.delete(0, END)
            self.batchOutputEntry.insert(0, self.defaultInDir)
            self.batchOutputEntry.configure(width=len(self.defaultOutDir))

    def make_buttons(self):
        self.objectButton = Button(
            self.mainFrame,
            text='Select Object',
            command=self.select_object)
        self.objectButton.grid(row=2, column=1, sticky='W')
        self.objectFileDropBox.grid(row=2, column=0, sticky='W')

        self.flatButton = Button(
            self.mainFrame,
            text='Select Flat',
            command=self.select_flat)
        self.flatButton.grid(row=4, column=1, sticky='W')
        self.flatFileDropBox.grid(row=4, column=0, sticky='W')

        self.aFrameButton = Button(
            self.mainFrame,
            text='Select A Frame',
            command=self.select_a)
        self.bFrameButton = Button(
            self.mainFrame,
            text='Select B Frame',
            command=self.select_b)

        self.pairFlatButton = Button(
            self.mainFrame,
            text='Select Flat',
            command=self.select_pair_flat)

        self.outputDirButton = Button(
            self.mainFrame,
            text='Select',
            command=self.update_outdir)
        self.outputDirButton.grid(row=6, column=1, sticky='W')

        self.batchInputButton = Button(
            self.mainFrame,
            text='Select',
            command=self.update_batch_indir)
        self.batchOutputButton = Button(
            self.mainFrame,
            text='Select',
            command=self.update_batch_outdir)

        self.cosmicOption = Checkbutton(
            self.mainFrame,
            text='Disable Cosmic',
            onvalue='-no_cosmic',
            offvalue='',
            variable=self.var1)
        self.cosmicOption.grid(row=8, column=0, sticky='W')
        self.jpegOption = Checkbutton(
            self.mainFrame,
            text='JPG',
            onvalue='-jpg',
            offvalue='',
            variable=self.var2)
        self.jpegOption.grid(row=9, column=0, sticky='W')

        self.runButton = Button(
            self.mainFrame,
            text='Run NSDRP',
            command=self.go)
        self.runButton.configure(width=10)
        self.runButton.grid(row=10, column=0, sticky='W')

        self.plotButton = Button(
            self.mainFrame,
            text='Show Plot',
            command=self.show_plot)

        self.normalModeButton = Radiobutton(
            self.radioButtonFrame,
            text='Single Frame',
            variable=self.modeVar,
            value=1,
            command=self.get_mode).grid(
            row=0,
            column=0,
            sticky='W')
        self.pairModeButton = Radiobutton(
            self.radioButtonFrame,
            text='Pair Subtraction',
            variable=self.modeVar,
            value=2,
            command=self.get_mode).grid(
            row=0,
            column=1,
            sticky='W')
        self.modeVar.set(1)
        self.batchButton = Radiobutton(
            self.radioButtonFrame,
            text='Batch',
            variable=self.modeVar,
            value=3,
            command=self.get_mode).grid(
            row=0,
            column=2,
            sticky='W')

    def get_mode(self):
        # try:
            # self.productsDropBox.grid_forget()
            # self.previewsListBox.grid_forget()
            # self.productsLabel.grid_forget()
            # self.plotButton.grid_forget()
        # except:
            # pass

        self.modes.append(self.modeVar.get())
        self.statusLabel.grid_forget()
        if self.modeVar.get() == 1:
            # Remove pair
            if self.modes[-2] == 2:
                self.aFrameFileLabel.grid_forget()
                self.aFrameDropBox.grid_forget()
                self.aFrameButton.grid_forget()

                self.bFrameFileLabel.grid_forget()
                self.bFrameDropBox.grid_forget()
                self.bFrameButton.grid_forget()

                self.pairFlatLabel.grid_forget()
                self.pairFlatDropBox.grid_forget()
                self.pairFlatButton.grid_forget()

                self.outDirLabel.grid_forget()
                self.outDirEntry.grid_forget()
                self.outputDirButton.grid_forget()
            # Remove batch
            else:
                self.batchInputLabel.grid_forget()
                self.batchInputEntry.grid_forget()
                self.batchInputButton.grid_forget()

                self.batchOutputLabel.grid_forget()
                self.batchOutputEntry.grid_forget()
                self.batchOutputButton.grid_forget()

            self.objectFileLabel.grid(row=1, column=0, sticky='W')
            self.objectFileDropBox.grid(row=2, column=0, sticky='W')
            self.objectButton.grid(row=2, column=1, sticky='W')

            self.flatLabel.grid(row=3, column=0, stick='W')
            self.flatFileDropBox.grid(row=4, column=0, sticky='W')
            self.flatButton.grid(row=4, column=1, sticky='W')

            self.outDirLabel.grid(row=5, column=0, sticky='W')
            self.outDirEntry.grid(row=6, column=0, sticky='W')
            self.outputDirButton.grid(row=6, column=1, sticky='W')

#                self.cmdOptionsLabel.grid(row=7, column=0, stick='W')
#                self.cosmicOption.grid(row=8, column=0, sticky='W')
#                self.jpegOption.grid(row=9, column=0, sticky='W')
#                self.runButton.grid(row=10, column=0, sticky='W')

        elif self.modeVar.get() == 2:
            if self.modes[-2] == 1:
                # Remove  normal mode
                self.objectFileLabel.grid_forget()
                self.objectButton.grid_forget()
                self.objectFileDropBox.grid_forget()

                self.flatLabel.grid_forget()
                self.flatButton.grid_forget()
                self.flatFileDropBox.grid_forget()

                self.outDirLabel.grid_forget()
                self.outDirEntry.grid_forget()
                self.outputDirButton.grid_forget()

            else:
                # Remove batch mode
                self.batchInputLabel.grid_forget()
                self.batchInputEntry.grid_forget()
                self.batchInputButton.grid_forget()

                self.batchOutputLabel.grid_forget()
                self.batchOutputEntry.grid_forget()
                self.batchOutputButton.grid_forget()

                self.outDirLabel.grid_forget()
                self.outDirEntry.grid_forget()
                self.outputDirButton.grid_forget()

            # Make pair
            self.aFrameFileLabel.grid(row=1, column=0, sticky='W')
            self.aFrameDropBox.grid(row=2, column=0, sticky='W')
            self.aFrameButton.grid(row=2, column=1, sticky='W')

            self.bFrameFileLabel.grid(row=3, column=0, sticky='W')
            self.bFrameDropBox.grid(row=4, column=0, sticky='W')
            self.bFrameButton.grid(row=4, column=1, sticky='W')

            self.pairFlatLabel.grid(row=5, column=0, sticky='W')
            # self.pairFlatDropBox = OptionMenu(self.mainFrame, self.pairFlatVar, *self.fitsFiles)
            # self.pairFlatDropBox.configure(width=self.width)
            self.pairFlatDropBox.grid(row=6, column=0, stick='W')
            self.pairFlatButton.grid(row=6, column=1, sticky='W')

            self.outDirLabel.grid(row=7, column=0, sticky='W')
            self.outputDirButton.grid(row=8, column=1, sticky='W')
            self.outDirEntry.grid(row=8, column=0, sticky='W')

            self.cmdOptionsLabel.grid(row=9, column=0, sticky='W')
            self.cosmicOption.grid(row=10, column=0, sticky='W')
            self.jpegOption.grid(row=11, column=0, sticky='W')
            self.runButton.grid(row=12, column=0, sticky='W')

        else:
            # Remove object
            if self.modes[-2] == 1:
                self.objectFileLabel.grid_forget()
                self.objectFileDropBox.grid_forget()
                self.objectButton.grid_forget()

                self.flatLabel.grid_forget()
                self.flatFileDropBox.grid_forget()
                self.flatButton.grid_forget()

                self.outDirLabel.grid_forget()
                self.outDirEntry.grid_forget()
                self.outputDirButton.grid_forget()
            # Remove pair
            else:
                self.aFrameFileLabel.grid_forget()
                self.aFrameDropBox.grid_forget()
                self.aFrameButton.grid_forget()

                self.bFrameFileLabel.grid_forget()
                self.bFrameDropBox.grid_forget()
                self.bFrameButton.grid_forget()

                self.pairFlatLabel.grid_forget()
                self.pairFlatDropBox.grid_forget()
                self.pairFlatButton.grid_forget()

                self.outDirLabel.grid_forget()
                self.outputDirButton.grid_forget()
                self.outDirEntry.grid_forget()
            # Make batch
            self.batchInputLabel.grid(row=1, column=0, sticky='W')
            self.batchInputEntry.grid(row=2, column=0, sticky='W')
            self.batchInputButton.grid(row=2, column=1, sticky='W')

            self.batchOutputLabel.grid(row=3, column=0, sticky='W')
            self.batchOutputEntry.grid(row=4, column=0, sticky='W')
            self.batchOutputButton.grid(row=4, column=1, sticky='W')

    def select_object(self):
        self.selectedObject = askopenfilename(initialdir=self.inputDir)
        self.files.append(self.selectedObject)
        if isinstance(self.selectedObject, tuple):
            self.selectedObject = ''
        else:
            self.selectedObject = self.selectedObject.split('/')[-1]
            self.objectVar.set(self.selectedObject)

    def select_flat(self):
        self.selectedFlat = askopenfilename(initialdir=self.inputDir)
        self.files.append(self.selectedFlat)
        if isinstance(self.selectedFlat, tuple):
            self.selectedFlat = ''
        else:
            self.selectedFlat = self.selectedFlat.split('/')[-1]
            self.flatVar.set(self.selectedFlat)

    def select_a(self):
        self.selectedA = askopenfilename(initialdir=self.inputDir)
        self.files.append(self.selectedA)
        if isinstance(self.selectedA, tuple):
            self.selectedA = ''
        else:
            self.selectedA = self.selectedA.split('/')[-1]
            self.aVar.set(self.selectedA)

    def select_b(self):
        self.selectedB = askopenfilename(initialdir=self.inputDir)
        self.files.append(self.selectedB)
        if isinstance(self.selectedB, tuple):
            self.selectedB = ''
        else:
            self.selectedB = self.selectedB.split('/')[-1]
            self.bVar.set(self.selectedB)

    def select_pair_flat(self):
        self.selectedPairFlat = askopenfilename(initialdir=self.inputDir)
        self.files.append(self.selectedPairFlat)
        if isinstance(self.selectedPairFlat, tuple):
            self.selectedPairFlat = ''
        else:
            self.selectedPairFlat = self.selectedPairFlat.split('/')[-1]
            self.pairFlatVar.set(self.selectedPairFlat)

    def show_plot(self):
        try:
            self.fileLabel.grid_forget()
            self.img = imopen.open(self.productsDir +
                                   '/previews/' +
                                   self.productVar.get() +
                                   '/' +
                                   self.previews[self.previewsListBox.curselection()[0]])
            self.photo = ImageTk.PhotoImage(self.img)
            self.photoLabel = Label(self.plotFrame, image=self.photo)
            self.photoLabel.grid(row=0, column=0)
        except IndexError:
            self.fileLabel.configure(text='Select a data product!', fg='red')
            self.fileLabel.grid(row=16, column=0, sticky='E')

    def func(self, value):
        self.previews = (
            os.listdir(
                self.productsDir +
                '/previews/' +
                self.productVar.get()))
        self.previewFileLength = max(len(s) for s in self.previews)
        self.previewsListBox = Listbox(self.mainFrame, selectmode='SINGLE')
        self.previewsListBox.configure(width=self.previewFileLength)
        for item in self.previews:
            self.previewsListBox.insert(END, item)
            self.previewsListBox.grid(row=15, column=0, sticky='W')

        self.plotButton.grid(row=16, column=0, sticky='W')

    def products(self):
        self.t.join()

        self.runButton.configure(state=NORMAL, relief=RAISED)
        self.dataProducts = []

        if os.path.isdir(self.productsDir):
            for root, dirs, files in os.walk(self.productsDir + '/previews'):
                for name in dirs:
                    self.dataProducts.append(name)

        if self.err:
            self.err = str(self.err.splitlines()[-1])
            self.err = self.err.split(':')[1:]
            self.statusLabel.configure(
                text='Critical Error:' +
                ' '.join(
                    self.err),
                fg='red')
            self.statusLabel.grid(row=self.statusRow, column=1, sticky='W')
        elif not self.err:
            # self.dataProducts:
            self.productsLabel = Label(
                self.mainFrame,
                text='Data Products:',
                font=self.labelsFont)
            self.productsLabel.grid(row=13, column=0, sticky='W')
            self.productsDropBox = OptionMenu(
                self.mainFrame,
                self.productVar,
                *self.dataProducts,
                command=self.func)
            self.productsDropBox.configure(width=10)
            self.productsDropBox.grid(row=14, column=0, sticky='W')
        else:
            self.statusLabel.configure(
                text='No data products. Check the log.', fg='black')
            self.statusLabel.grid(row=self.statusRow, column=1, sticky='W')

    def go(self):
        self.statusLabel.grid_forget()

        if self.modeVar.get() == 1:
            self.statusRow = 12
        else:
            self.statusRow = 10

        if self.modeVar.get() == 1 or self.modeVar.get() == 2:
            self.productsDir = self.outDirEntry.get()
        else:
            self.productsDir = self.batchOutVar.get()

        if not os.path.exists(self.productsDir):
            os.makedirs(self.productsDir)

        if os.path.exists(self.productsDir):
            if self.modeVar.get() == 2:
                self.pairFiles = [
                    self.aVar.get(),
                    self.bVar.get(),
                    self.pairFlatVar.get()]
                if len(self.aVar.get()) == 0 or len(self.bVar.get()) == 0 or len(self.pairFlatVar.get()) == 0 or self.aVar.get(
                ) == self.missingFile or self.bVar.get() == self.missingFile or self.pairFlatVar.get() == self.missingFile:
                    self.statusLabel.configure(
                        text='Missing a filename!', fg='red')
                    self.statusLabel.grid(row=12, column=0, sticky='E')
                elif len(self.pairFiles) != len(set(self.pairFiles)):
                    self.statusLabel.configure(text='Same filename!', fg='red')
                    self.statusLabel.grid(row=12, column=0, sticky='E')
                else:
                    self.statusLabel.configure(text='Reducing...')
                    self.statusLabel.grid(row=12, column=0, sticky='E')

                    self.runButton.configure(state=DISABLED, relief=SUNKEN)

                    self.objectPattern = '*' + self.aVar.get()
                    self.bFramePattern = '*' + self.bVar.get()
                    self.pairFlatPattern = '*' + self.pairFlatVar.get()

                    self.pairObject = fnmatch.filter(
                        self.files, self.objectPattern)[0]
                    self.bFrame = fnmatch.filter(
                        self.files, self.bFramePattern)[0]
                    self.pairFlat = fnmatch.filter(
                        self.files, self.pairFlatPattern)[0]

                    self.pairCmd = 'python nsdrp.py -out_dir ' + self.productsDir + ' ' + self.var1.get() + \
                        ' ' + self.var2.get() + ' ' + self.pairFlat + ' ' + \
                        self.pairObject + ' -b ' + self.bFrame

                    self.execute = sp.Popen(
                        [self.pairCmd], stdout=sp.PIPE, stderr=sp.PIPE, shell=True)

                    def run():
                        self.out, self.err = self.execute.communicate()

                    self.t = threading.Thread(target=run)
                    self.t.start()

                    self.products()
            elif self.modeVar.get() == 1:
                if len(self.objectVar.get()) == 0 or len(self.flatVar.get()) == 0 or self.objectVar.get(
                ) == self.missingFile or self.flatVar.get() == self.missingFile:
                    self.statusLabel.configure(
                        text='Missing a filename!', fg='red')
                    self.statusLabel.grid(row=10, column=0, sticky='E')
                elif self.objectVar.get() == self.flatVar.get():
                    self.statusLabel.configure(text='Same filename!', fg='red')
                    self.statusLabel.grid(row=10, column=0, sticky='E')
                else:
                    # self.statusLabel.configure(text='Reducing...')
                    # self.statusLabel.grid(row=8, column=0, sticky='E')

                    self.runButton.configure(state=DISABLED, relief=SUNKEN)

                    self.objectPattern = '*' + self.objectVar.get()
                    self.flatPattern = '*' + self.flatVar.get()

                    self.object = fnmatch.filter(
                        self.files, self.objectPattern)[0]
                    self.flat = fnmatch.filter(self.files, self.flatPattern)[0]

                    self.cmd = 'python nsdrp.py -out_dir ' + self.productsDir + ' ' + \
                        self.var1.get() + ' ' + self.var2.get() + ' ' + self.object + ' ' + self.flat

                    self.execute = sp.Popen(
                        [self.cmd], stdout=sp.PIPE, stderr=sp.PIPE, shell=True)

                    def run():
                        self.out, self.err = self.execute.communicate()

                    self.t = threading.Thread(target=run)
                    self.t.start()

                    self.products()
            else:
                if len(self.batchInVar.get()) == 0 or len(
                        self.batchOutVar.get()) == 0:
                    self.statusLabel.configure(
                        text='Missing a directory', fg='red')
                    self.statusLabel.grid(row=10, column=0, sticky='E')
                # elif not os.path.isdir(self.batchInVar.get()) or not
                # os.path.isdir(self.batchOutVar.get()):
                elif not os.path.isdir(self.batchInVar.get()):
                    self.statusLabel.configure(
                        text='Input directory does not exists', fg='red')
                    self.statusLabel.grid(row=10, column=0, sticky='E')
                else:
                    self.runButton.configure(state=DISABLED, relief=SUNKEN)

                    self.batchCmd = 'python nsdrp.py ' + self.var1.get() + ' ' + self.var2.get() + \
                        ' ' + self.batchInVar.get() + ' ' + self.batchOutVar.get()

                    self.execute = sp.Popen(
                        [self.batchCmd], stdout=sp.PIPE, stderr=sp.PIPE, shell=True)

                    def run():
                        self.out, self.err = self.execute.communicate()

                    self.t = threading.Thread(target=run)
                    self.t.start()

                    self.products()
    #    else:
        #    pass
        #    self.statusLabel.configure(text='Output directory does not exist', fg='red')
        #    self.statusLabel.grid(row=self.statusRow, column=1, sticky='E')


if __name__ == "__main__":
    my_gui = nsdrpGui()
